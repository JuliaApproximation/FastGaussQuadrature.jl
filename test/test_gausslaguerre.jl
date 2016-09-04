# Test gausslaguerre.jl
import FastGaussQuadrature.getUQ
import FastGaussQuadrature.asyAiry
import FastGaussQuadrature.asyAirygen
import FastGaussQuadrature.asyBessel
import FastGaussQuadrature.asyBesselgen
import FastGaussQuadrature.asyBulk
import FastGaussQuadrature.asyBulkgen
import FastGaussQuadrature.polyAsyRH
import FastGaussQuadrature.polyAsyRHgen

# Test integration
tol = 3.e-13
tolGen = 2.e-9
# Possibly set srand(...) to get the same results for integration of a random polynomial.
ns = [42; 251; 5000];
for wei = 1:3
    if wei == 1; alpha = 0.0;qm = 1.0; m = 1;
    elseif wei == 2; alpha = pi+0.1; qm = 1.0; m = 1;
    elseif wei == 3; alpha = 0.67; qm = 0.34; m = 4;
    end

    for ni = 1:length(ns)
        n = ns[ni]
	# Make sure the integrand does not overflow at x[n]
        d = min(2*n-1, ceil(Int64, log(realmax(Float64)-4*n)/n ) )
        ra = rand(d+1)
        function fct(x)
            if length(x) == 1
                return sum(ra.*(x.^(0:d)))
            end
            result = NaN*x
            for i = 1:length(x)
                result[i] = sum(ra.*(x[i].^(0:d)))
            end
            return result
        end
        fcti(x) = fct(x).*exp(-qm*x.^m).*x^(alpha)
        # Test integration of a random polynomial which should be integrated exactly:
        exa = quadgk(fcti, 0, Inf, reltol = tol)[1]
        for methi = 1:ifelse(ni == 1, 1, 5)
            if (ni == 1)
                # Golub-Welsch is the accurate algorithm for low n but does not allow m or qm != 1.
                if wei == 3; continue; end
                println("    α = $alpha, n = $ni, Golub-Welsch")
                x, w = gausslaguerre( n, alpha )
            elseif methi == 1
                println("    α = $alpha, n = $ni, default")
                x, w = gausslaguerre( n, alpha, "default", qm, m)
            elseif methi == 2
                if wei == 3; continue; end # General weight
                println("    α = $alpha, n = $ni, RH")
                x, w = gausslaguerre( n, alpha, "RH", qm, m)
            elseif methi == 3
                if wei == 3; continue; end # General weight
                println("    α = $alpha, n = $ni, RHW")
                x, w = gausslaguerre( n, alpha, "RHW", qm, m)
            elseif methi == 4
                println("    α = $alpha, n = $ni, gen")
                x, w = gausslaguerre( n, alpha, "gen", qm, m)
            elseif methi == 5
                println("    α = $alpha, n = $ni, genW")
                x, w = gausslaguerre( n, alpha, "genW", qm, m)
            end
            # Multiply the tolerance by sqrt(n) because the error on each node/meaningful weight should be about 1e-14
            @test abs(dot(x,w) - gamma((alpha+2)/m)*qm^(-(alpha+2)/m)/m )/gamma((alpha+2)/m)/qm^(-(alpha+2)/m)*m <= sqrt(n)*ifelse(wei == 3, tolGen, tol)
            @test abs(dot(x.^2,w) - gamma((alpha+3)/m)*qm^(-(alpha+3)/m)/m )/gamma((alpha+3)/m)/qm^(-(alpha+3)/m)*m <= sqrt(n)*ifelse(wei == 3, tolGen, tol)
            if (n == 42) && (wei == 1)
                @test abs(x[37] - 98.388267163326702) < tol
                @test abs(w[7] - 0.055372813167092) < tol
            elseif (n == 251) && (wei == 1)
                @test abs(x[37] - 13.309000189442097) < tol
                @test abs(w[3] - 0.050091759039996) < tol
            end
            @test abs(dot(fct(x), w) - exa)/abs(exa) < sqrt(n)*ifelse(wei == 3, tolGen, tol)
	    end
    end
end

# Test whether higher-order terms correspond between "RH(W)", "gen(W)" and exact results computed in high precision.
n = NaN
tolHot = eps(Float64)*10
np = 200
tolEx = np*4e-15
T = ceil(Int64, 34/log(np) )
@test T == 7
for alpha = [0.0; 4.15]
    println("    α = $alpha")
    UQ = getUQ(alpha, 1.0, 1, T+1)

    z = 0.95 # > 3.7/4
    mnxi = 2*np*( sqrt(z)*sqrt(1 - z) - acos(sqrt(z) ) )
    fn = (mnxi*3im/2)^(2/3)
    aRHa = asyAiry(np, 4*np*z, alpha, T)
    @test abs(asyAirygen(np, z, alpha, T, 1.0, 1, UQ, fn, true) -aRHa)./abs(aRHa) < tolHot
    # Not using Q has an explicit substraction of ceil(3*T/2)-order poles, so loosen the tolerance.
    @test abs(asyAirygen(np, z, alpha, T, 1.0, 1, UQ, fn, false, mnxi*1im/np) -aRHa)./abs(aRHa) < 300*tolHot

    z = 0.01 # < sqrt(np)/4/np
    mnxi = 2*np*( sqrt(z)*sqrt(1 - z) - acos(sqrt(z) ) )
    aRHe = asyBessel(np, 4*np*z, alpha, T)
    @test abs(asyBesselgen(np, z, alpha, T, 1.0, 1, UQ, mnxi + pi*np, true) -aRHe)./abs(aRHe) < 20*tolHot
    @test abs(asyBesselgen(np, z, alpha, T, 1.0, 1, UQ, mnxi + pi*np, false) -aRHe)./abs(aRHe) < 20*tolHot

    z = 0.35
    mnxi = 2*np*( sqrt(z)*sqrt(1 - z) - acos(sqrt(z) ) )
    aRHu = asyBulk(np, 4*np*z, alpha, T)
    @test abs(asyBulkgen(np, z, alpha, T, 1.0, 1, UQ, mnxi) -aRHu)./abs(aRHu) < tolHot

    a = alpha
    factorp = (1/3840*a^10 - 5/2304*a^9 + 11/2304*a^8 + 7/1920*a^7 - 229/11520*a^6 + 107/34560*a^5 + 2653/103680*a^4 - 989/155520*a^3 -         3481/311040*a^2 + 139/103680*a + 9871/6531840)/np^5
    factorp = factorp + (1/384*a^8 - 1/96*a^7 + 1/576*a^6 + 43/1440*a^5 - 5/384*a^4 - 23/864*a^3 + 163/25920*a^2 + 31/6480*a - 139/155520)/np^4
    factorp = factorp + (1/48*a^6 - 1/48*a^5 - 1/24*a^4 + 5/144*a^3 + 1/36*a^2 - 1/144*a - 31/6480)/np^3
    factorp = factorp + (1/8*a^4 + 1/12*a^3 - 1/24*a^2 + 1/72)/np^2
    factorp = factorp + (1/2*a^2 + 1/2*a + 1/6)/np^1 +1
    # Add an asymptotic approximation of gamma(alpha+n+1)/gamma(n+1)
    factorp = (4*np)^(-alpha/2-1/2)*sqrt(1/2/pi)/sqrt(factorp)*sqrt(np^alpha*(1+alpha*(alpha+1)/2/np + alpha*(3*alpha^3+2*alpha^2-3*alpha-2)/24/np^2 + alpha^2*(alpha+1)^2*(alpha^2-3*alpha+2)/48/np^3 + alpha*(15*alpha^7-60*alpha^6-10*alpha^5+192*alpha^4 -25*alpha^3-180*alpha^2+20*alpha+48)/5760/np^4 + alpha^2*(alpha+1)^2*(3*alpha^6-31*alpha^5+109*alpha^4-125*alpha^3-88*alpha^2+276*alpha-144)/11520/np^5) )
    if alpha == 0.0
        @test norm([+6.975493821675553e-05  ,  -8.294589753829749e-05/1im]-vec(UQ[1,:,6,1,2]))/norm(UQ[1,:,6,1,2]) < tolHot
        @test norm([-6.975493821675553e-05  ,  +1.401453805748335e-05/1im]-vec(UQ[1,:,6,1,1]))/norm(UQ[1,:,6,1,1]) < tolHot*10
        # We choose to test the constants most prone to errors.
        @test norm([-0.01196102063075393    ,  +0.1904949571852236/1im]-vec(UQ[1,:,6,3,4]))/norm(UQ[1,:,6,1,4]) < tolHot*1e3
        @test norm([+1.309628097160176e-05  ,  -0.001326131973043531/1im]-vec(UQ[1,:,6,3,3]))/norm(UQ[1,:,6,1,3]) < tolHot*10

        # Compare with results obtained from a high precision calculation of the Laguerre polynomials with standard normalisation.
        @test abs(-4.9327289759964085492326419516911524182405e163/exp(0.95*4*np/2)/factorp -aRHa)/abs(aRHa) < tolEx
        @test abs(-1.6224971613979232593218524827791333158448e59/exp(0.35*4*np/2)/factorp -aRHu)/abs(aRHu) < tolEx
        @test abs(-3.9144031220161646274486036731469804751328/exp(0.01*4*np/2)/factorp -aRHe)/abs(aRHe) < tolEx
    else
        @test abs(5/49152*alpha^8 - 35/49152*alpha^7 + 67/49152*alpha^6 + 11/36864*alpha^5 - 1529/589824*alpha^4 + 1891/2359296*alpha^3 + 26827/26542080*alpha^2 - 109/524288*alpha - 190867/1698693120  -UQ[1,1,4,1,2])/abs(UQ[1,1,4,1,2]) < 6*tolHot
        @test abs(-1/24576*alpha^8 - 17/24576*alpha^7 - 107/147456*alpha^6 + 1751/983040*alpha^5 + 1535/1179648*alpha^4 - 853/884736*alpha^3 - 7583/26542080*alpha^2 + 14131/141557760*alpha - 12829/1274019840 -1im*4^alpha*UQ[1,2,4,2,1])/abs(4^alpha*UQ[1,2,4,2,1]) < tolHot
        @test abs(-23/259187003712000*alpha^21 + 353/76796149248000*alpha^20 - 333791/4561691265331200*alpha^19 + 87631/600222534912000*alpha^18 + 24928777/4001483566080000*alpha^17 - 205079/4803701760000*alpha^16 - 121008773/882680198400000*alpha^15 + 221720011/117690693120000*alpha^14 - 11808705347/16991593819200000*alpha^13 - 119431297679/3485455142400000*alpha^12 + 100119618541/1742727571200000*alpha^11 + 3526062948193/12357522777600000*alpha^10 - 83939498901133/123575227776000000*alpha^9 - 48402313811119/41191742592000000*alpha^8 + 4645283107358809/1050389436096000000*alpha^7 - 91277362389443/40399593696000000*alpha^6 - 49328830551721603/7622617782780000000*alpha^5 + 23985609069512791/1960101715572000000*alpha^4 + 2660362470248081/3587114250720000000*alpha^3 - 779839198817230403/60980942262240000000*alpha^2 + 361055235303455299/3347022171893400000000*alpha + 91136762727842675269/36817243890827400000000 - 1im*4^alpha*UQ[1,2,4,5,3])/abs(4^alpha*UQ[1,2,4,5,3]) < 3e-11 # First term = 8.9e-14 times alpha^21=9.5e12, and some integers in the formula cannot be represented in Int64, so lower accuracy can be expected
        @test abs(1/362880*alpha^13 - 1/453600*alpha^12 - 71/435456*alpha^11 + 1361/7257600*alpha^10 + 625/217728*alpha^9 - 21569/7257600*alpha^8 - 1337/64800*alpha^7 + 291101/21772800*alpha^6 + 9329/170100*alpha^5 + 3224237/59875200*alpha^4 - 272927/3628800*alpha^3 - 6823639/20528640*alpha^2 + 138793/2721600*alpha + 40018577/179625600 - UQ[1,1,3,4,4])/abs(UQ[1,1,3,4,4]) < 40*tolHot

        @test abs(-3.7536142288312445691754139973055317263278e162/exp(0.95*4*np/2)/factorp -aRHa)/abs(aRHa) < tolEx
        @test abs(1.3002334410655323780033551365563985116651e59/exp(0.35*4*np/2)/factorp -aRHu)/abs(aRHu) < tolEx
        @test abs(-2395.51952258921326919097391744358588135/exp(0.01*4*np/2)/factorp -aRHe)/abs(aRHe) < tolEx
    end

    for n = ceil(Int64, linspace(111, 5000, 5) )
        for z = linspace(1e-4, 0.99, 10)
            aRH = polyAsyRH(n, 4*n*z, alpha, T)
            @test abs(polyAsyRHgen(n, 4*n*z, alpha, T, 1.0, 1, UQ) -aRH)/abs(aRH) < 2*tolEx
        end
    end
end
