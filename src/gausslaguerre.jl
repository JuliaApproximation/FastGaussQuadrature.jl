# (x,w) = gausslaguerre(n) returns n Gauss-Laguerre nodes and weights.
# (x,w) = gausslaguerre(n,alpha) allows generalized Gauss-Laguerre quadrature.
# (x,w) = gausslaguerre(n,alpha,method) allows the user to select which method to use.
# (x,w) = gausslaguerre(n,alpha,method,qm,m) allows the generalised weight function x^alpha*exp(-qm*x^m) where qm and m should be one when method = GW, rec(W) or exp(W)
# METHOD = "GW" will use the traditional Golub-Welsch eigenvalue method, which is best for when N is small.
# METHOD = "gen" can generate an arbitrary number of terms of the asymptotic expansion of Laguerre-type polynomials, orthogonal with respect to x^alpha*exp(-qm*x^m). "genW" does the same, but stops as the weights underflow: it is O(sqrt(n)) when m is one. gausslaguerre(round(Int64, (n/17)^2), alpha, "genW") returns about n nodes and weights above realmin(Float64) for large n.
# METHOD = "exp" will use explicit expansions and METHOD = "expW" will only compute the weights above realmin.
# METHOD = "rec" uses forward recurrence and METHOD = "recW" stops when the weights underflow.
# METHOD = "default" uses "gen" when m or qm are not one, explicit formulae when n <= 2, "rec" when 2 < n < 128 and else "exp".

function gausslaguerre( n::Int64, alpha::Float64=0.0, method::AbstractString="default",  qm::Float64=1.0, m::Int64=1 )
    
    if ( imag(alpha) != 0 ) || ( alpha < -1 )
        error(string("alpha = ", alpha, " is not allowed.") )
    elseif ( (m != 1) || (qm != 1) ) && ( (method != "gen") && (method != "genW") && (method != "default") )
        error(string("Method ", method, " is not implemented for generalised weights.") )
    end

    if ( (method == "default") && ( (m != 1) || (qm != 1.0) ) ) || (method == "gen") || (method == "genW")
        (x,w) = asyRHgen(n, method == "genW", alpha, m, qm)
        w = w/sum(w)*gamma((alpha+1)/m)*qm^(-(alpha+1)/m)/m # We left out a constant factor while computing w in case we use finite differences
        (x,w)
    elseif (method == "default") && (n == 0)
        Float64[], Float64[]
    elseif (method == "default") && (n == 1)
        [1.0+alpha], [1.0]
    elseif (method == "default") && (n == 2)
        [alpha+2.-sqrt(alpha+2.),alpha+2.+sqrt(alpha+2.)], [((alpha-sqrt(alpha+2)+2)*gamma(alpha+2))/(2*(alpha+2)*(sqrt(alpha+2)-1)^2),((alpha+sqrt(alpha+2)+2)*gamma(alpha+2))/(2*(alpha+2)*(sqrt(alpha+2)+1)^2)]
    elseif method == "GW"
        laguerreGW( n, alpha )         # Use Golub-Welsch
    elseif ( (method == "default") && (n < 128 ) ) || (method == "rec") || (method == "recW")
        laguerreRec( n, method == "recW", alpha)   # Use forward recurrence
    elseif ( (method == "default") && (n >= 128) ) || (method == "exp") || (method == "expW")
        laguerreExp( n, method == "expW", alpha)   # Use explicit expansions and possibly only compute the representable weights
    else
        error(string("Wrong method string, got ", method) )
    end
end

function laguerreGW( n::Int64, alpha::Float64 )
# Calculate Gauss-Laguerre nodes and weights based on Golub-Welsch

    alph = 2*(1:n)-1+alpha           # 3-term recurrence coeffs
    beta = sqrt.( (1:n-1).*(alpha .+ (1:n-1) ) )
    T = SymTridiagonal(alph, beta)  # Jacobi matrix
    x, V = eig( T )                  # eigenvalue decomposition
    w = gamma(alpha+1)*V[1,:].^2     # Quadrature weights
    x, vec(w)
end



########################## Routines for the forward recurrence ##########################

function laguerreRec( n::Int64, compRepr::Bool, alpha::Float64)

    print("myRec")
    if compRepr
        # Get a heuristic for the indices where the weights are about above realmin.
        mn = min(ceil(Int64, 17*sqrt(n)), n)
    else
        mn = n
    end
    itric = min(mn, 7);
    bes = besselroots(alpha, itric).^2/(4*n + 2*alpha+2) # [DLMF 18.16.10] says this is lower than the zero, so we do not risk skipping a zero (if none of the results coincide)

    factorw = (n^2 +alpha*n)^(-1/2) # Ratio of leading order coefficients
    w = zeros(mn)
    x = [ bes; zeros(mn - itric)]
    noUnderflow = true
    for k = 1:mn
        if ( x[k] == 0 ) # Use sextic extrapolation for the initial guesses.
            x[k] = 7*x[k-1] -21*x[k-2] +35*x[k-3] -35*x[k-4] +21*x[k-5] -7*x[k-6] +x[k-7]
        end
        step = x[k]
        l = 0 # Newton-Raphson iteration number
        ov = realmax(Float64) # Previous/old value
        ox = x[k] # Old x
        while ( ( abs(step) > eps(Float64)*40*x[k] ) && ( l < 20) )
            l = l + 1
            # pe = Float64(lagpnRec(big(n), big(alpha), big(x[k]) ) ) #Calculations in BigFloat
            pe = lagpnRec(n, alpha, x[k] )
            if (abs(pe) >= abs(ov)*(1-35*eps(Float64)) )
                # The function values do not decrease enough any more due to roundoff errors.
                x[k] = ox # Set to the previous value and quit.
                break
            end
            # poly' = (p*exp(-Q/2) )' = exp(-Q/2)*(p' -p/2) with orthonormal p.
            step = pe/lagpnRecDer(n, alpha, x[k])
            ox = x[k]
            x[k] = x[k] -step
            ov = pe
    end
    if ( x[k] < 0 ) || ( x[k] > 4*n + 2*alpha + 2 ) || ( l == 20 ) || ( ( k != 1 ) && ( x[k - 1] >= x[k] ) )
        error("Newton method may not have converged.")
    end
    if noUnderflow
        w[k] = factorw/lagpnRec(n -1, alpha, x[k] )/lagpnRecDer(n, alpha, x[k] )
    end
    if noUnderflow && ( w[k] == 0 ) && ( k > 1 ) && ( w[k-1] > 0 ) # We could stop now as the weights underflow.
        if compRepr
		x = x[1:k-1]
		w = w[1:k-1]
		return (x,w)
	    else
                noUnderflow = false
	    end
        end
    end
    x, w
end

# Orthonormal associated Laguerre polynomial with positive leading coefficient, allows BigFloats
function lagpnRec(n,alpha,x)
    pnprev = 0*alpha;
    pn= 1/sqrt(gamma(alpha+1) );
    for k=1:n
        pnold = pn; 
        pn = (x -2*k -alpha+1)/sqrt(k*(alpha+k))*pn-sqrt((k-1+alpha)*(k-1)/k/(k+alpha))*pnprev;
        pnprev = pnold; 
    end; 
    return pn; 
end;

# Derivative of the orthonormal associated Laguerre polynomial, allows BigFloats
function lagpnRecDer(n,alpha,x)
    pnprev = 0*alpha;
    pn= 1/sqrt(gamma(alpha+1) );
    pndprev = 0*alpha;
    pnd = 0*alpha;
    for k=1:n
        pnold = pn;
        pn = (x -2*k -alpha+1)/sqrt(k*(alpha+k))*pnold-sqrt((k-1+alpha)*(k-1)/k/(k+alpha))*pnprev;
        pnprev = pnold;
        pndold = pnd;
        pnd = (pnold+(x-2*k-alpha+1)*pndold)/sqrt(k*(alpha+k)) -sqrt((k-1+alpha)*(k-1)/k/(alpha+k))*pndprev;
        pndprev = pndold;
    end;
    return pnd;
end;



############## Routines for the "gen(W)" algorithm for computing an arbitrary number of terms with general w(x) = x^alpha*exp(-qm*x^m) ##########################

function asyRHgen(n, compRepr, alpha, m, qm)

    T = ceil(Int64, 34/log(n) ) # Heuristic for number of terms, should be scaled by the logarithm of eps(Float64) over the machine precision.
    UQ0 = getUQ(alpha, qm ,m, T)
    if compRepr
        mn = min(n,ceil(Int64, exp(exp(1/m)*1.05)*n^(1-1/2/m) ))
    else
        mn = n
    end
    itric = max(ceil(Int64, 3.6*n^0.188), 7)
    # Heuristics to switch between Bessel, extrapolation and Airy initial guesses.
    igatt = ceil(Int64, mn + 1.31*n^0.4 - n)

    A = zeros(m+1)
    for k =0:m
        A[k+1] = prod((2*(1:k)-1)/2./(1:k))
    end
    softEdge = (n*2/m/qm/A[m+1] )^(1/m)
    # Use finite differences for derivative of polynomial when not x^alpha*exp(-x) and use other initial approximations
    useFinDiff = (m != 1) || (qm != 1.0)
    bes = besselroots(alpha, itric).^2 # [Tricomi 1947 pg. 296]
    w = zeros(mn)
    if useFinDiff
        x = [bes*(2*m-1)^2/16/m^2/n^2*softEdge ; zeros(mn-itric) ]
    else
        ak = [-13.69148903521072; -12.828776752865757; -11.93601556323626;    -11.00852430373326; -10.04017434155809; -9.02265085340981; -7.944133587120853;    -6.786708090071759; -5.520559828095551; -4.08794944413097; -2.338107410459767]
        t = 3*pi/2*( (igatt:-1:12)-0.25) # [DLMF (9.9.6)]
        ak = [-t.^(2/3).*(1 + 5/48./t.^2 - 5/36./t.^4 + 77125/82944./t.^6     -10856875/6967296./t.^8); ak[max(1,12-igatt):11] ]
        nu = 4*n+2*alpha+2 # [Gatteshi 2002 (4.9)]
        air = (nu+ak*(4*nu)^(1/3)+ ak.^2*(nu/16)^(-1/3)/5 + (11/35-alpha^2-12/175*ak.^3)/nu + (16/1575*ak+92/7875*ak.^4)*2^(2/3)*nu^(-5/3) -(15152/3031875*ak.^5+1088/121275*ak.^2)*2^(1/3)*nu^(-7/3))
        x = [ bes/(4*n + 2*alpha+2).*(1 + (bes + 2*(alpha^2 - 1) )/(4*n + 2*alpha+2)^2/3 ) ; zeros(mn - itric -max(igatt,0) ) ; air]
    end

    if !useFinDiff
        UQ1 = getUQ(alpha+1, qm, m, T)
        factor0 = 1-4im*4^alpha*sum((UQ0[1,2,1:(T-1),1, 2] + UQ0[1,2,1:(T-1),1])./n.^reshape(1:(T-1), (1,1,T-1)) )
        factor1 = 1-4im*4^(alpha+1)*sum((UQ1[1,2,1:(T-1),1, 2] + UQ1[1,2,1:(T-1),1])./n.^reshape(1:(T-1), (1,1,T-1)) )
        factorx = real(sqrt(factor0/factor1 )/2/(1 - 1/n)^(1+alpha/2))
        factorw = real( -(1 - 1/(n + 1) )^(n + 1+ alpha/2)*(1 - 1/n)^(1 + alpha/2)*exp(1 + 2*log(2) )*4^(1+alpha)*pi*n^alpha*sqrt(factor0*factor1)*(1 + 1/n)^(alpha/2) )
    end

    if ( alpha^2/n > 1 )
        warn("A large alpha may lead to inaccurate results because the weight is low and R(z) is not close to identity.")
    end
    noUnderflow = true
    for k = 1:mn
        if ( x[k] == 0 ) && useFinDiff # Use linear extrapolation for the initial guesses for robustness in generalised weights.
            x[k] = 2*x[k-1] -x[k-2]
        elseif ( x[k] == 0 ) # Use sextic extrapolation for the initial guesses.
            x[k] = 7*x[k-1] -21*x[k-2] +35*x[k-3] -35*x[k-4] +21*x[k-5] -7*x[k-6] +x[k-7]
        end
        step = x[k]
        l = 0 # Newton-Raphson iteration number
        ov = realmax(Float64) # Previous/old value
        ox = x[k] # Old x
        # Accuracy of the expansions up to machine precision would lower this bound.
        while ( ( abs(step) > eps(Float64)*40*x[k] ) && ( l < 20) )
            l = l + 1
            pe = polyAsyRHgen(n, x[k], alpha, T, qm, m, UQ0)
            if (abs(pe) >= abs(ov)*(1-35*eps(Float64)) )
                # The function values do not decrease enough any more due to roundoff errors.
                x[k] = ox # Set to the previous value and quit.
                break
            end
            if useFinDiff
                hh = max(sqrt(eps(Float64))*x[k], sqrt(eps(Float64)) )
                step = pe*hh/(polyAsyRHgen(n, x[k]+hh, alpha, T, qm, m, UQ0) - pe)
            else
                # poly' = (p*exp(-Q/2) )' = exp(-Q/2)*(p' -p/2) with orthonormal p.
                step = pe/(polyAsyRHgen(n-1, x[k], alpha+1, T, qm, m, UQ1)*factorx - pe/2)
            end
            ox = x[k]
            x[k] = x[k] -step
            ov = pe
        end
        if ( x[k] < 0 ) || ( l == 20 ) || ( ( k != 1 ) && ( x[k - 1] >= x[k] ) ) || isnan(x[k])
            print(x[k], "=x[k], k=", k, ", l=", l, ", x[k-1]=", x[k-1], ", x[k-2]=", x[k-1], ", step=", step, ", ox=", ox, ", ov=", ov, ".\n") # Print some debugging information and throw an error.
            error("Newton method may not have converged.")
        elseif ( x[k] > softEdge)
            warn("Node is outside the support of the measure: inaccuracy is expected.");
        end
        if noUnderflow&& useFinDiff
            hh = max(sqrt(eps(Float64))*x[k], sqrt(eps(Float64)) )
            w[k] = hh/(polyAsyRHgen(n,x[k]+hh,alpha,T,qm,m,UQ0) -polyAsyRHgen(n,x[k],alpha,T,qm,m,UQ0))/polyAsyRHgen(n-1,x[k],alpha,T,qm,m,UQ0)/exp(qm*x[k]^m) # This leaves out a constant factor, given by a ratio of leading order coefficients and normalising constants
        elseif noUnderflow
            w[k] = factorw/polyAsyRHgen(n-1, x[k], alpha+1, T, qm, m, UQ1)/polyAsyRHgen(n+1, x[k], alpha, T,qm,m,UQ0)/exp( x[k] )
        end
        if noUnderflow && ( w[k] == 0 ) && ( k > 1 ) && ( w[k-1] > 0 ) # We could stop now as the weights underflow.
            if compRepr
	        x = x[1:k-1]
	        w = w[1:k-1]
	        return (x,w)
	    else
                noUnderflow = false
	    end
        end
    end
    x, w
end

# Compute the expansion of the orthonormal polynomial without e^(qm*x^m/2) nor a constant factor based on some heuristics.
function polyAsyRHgen(np, y, alpha, T::Int64, qm, m::Int64, UQ)
    if (qm == 1) && (m == 1)
        z = y/4/np
        mnxi = 2*np*( sqrt(z).*sqrt(1 - z) - acos(sqrt(z) ) ) # = -n*xin/i
    else
        A = zeros(m+1)
        for k =0:m
            A[k+1] = prod((2*(1:k)-1)/2./(1:k))
        end
        z = y/(np*2/m/qm/A[m+1] )^(1/m)
        # Also correct but much slower: Hn = 4*m/(2*m-1)*double(hypergeom([1, 1-m], 3/2-m, z))/m
        Hn = 2/A[m+1]*sum(z.^(0:m-1).*A[m-(0:m-1)])/m
        mnxi = np*(sqrt(z+0im).*sqrt(1-z+0im).*Hn/2 -2*acos(sqrt(z+0im)))
    end
    # We could avoid these tests by splitting the loop k=1:mn into three parts with heuristics for the bounding indices.
    if y < sqrt(np)
        # The fixed delta in the Riemann-Hilbert Problem would mean this bound has to be proportional to n, but x(1:k) are O(1/n) so choose the bound in between them to make more use of the (cheap) expansion in the bulk.
        return asyBesselgen(np, z, alpha, T, qm, m, UQ, mnxi + pi*np, true)
    elseif y > 3.7*np
        # Use the expansion in terms of the (expensive) Airy function, although the corresponding weights will start underflowing for n >= 186 for standard associated Laguerre polynomials.
        return asyAirygen(np, z, alpha, T, qm, m, UQ, (mnxi*3im/2)^(2/3), true)
    end
    asyBulkgen(np, z, alpha, T, qm, m, UQ, mnxi)
end

function asyBulkgen(np, z, alpha, T::Int64, qm, m::Int64, UQ, mnxi)
    if T == 1
        return real( 2/(z+0im)^(1/4 + alpha/2)/(1 - z+0im)^(1/4)*cos(acos(2*z - 1+0im)*(1/2 + alpha/2) - mnxi - pi/4) )
    end
    R = [1 0]
    for k = 1:T-1
        for i = 1:ceil(Int64, 3*T/2)
            R = R + reshape((UQ[1,:,k,i,1]/(z-1)^i + UQ[1,:,k,i,2]/z^i)/np^k,(1,2))
        end
    end
    p = real( 2/(z+0im)^(1/4 + alpha/2)*(cos(acos(2*z-1+0im)*(1/2+alpha/2) - mnxi-pi/4)*R[1]-cos(acos(2*z-1+0im)*(-1/2+alpha/2)-mnxi-pi/4)*R[2]*1im*4^alpha)/(1 - z+0im)^(1/4) )
end

function asyBesselgen(np, z, alpha, T::Int64, qm, m::Int64, UQ, npb, useQ::Bool)
    if T == 1
        return real( sqrt(2*pi)*(-1)^np*sqrt(npb)/(z+0im)^(1/4+alpha/2)/(1 - z+0im)^(1/4)*(sin( (alpha + 1)/2*acos(2*z - 1+0im) - pi*alpha/2)*besselj(alpha,npb) + cos( (alpha + 1)/2*acos(2*z - 1+0im) - pi*alpha/2)*(besselj(alpha-1,npb) - alpha/(npb)*besselj(alpha, npb) ) ) )
    end
    R = (1+0im)*[1 0]
    # Use the series expansion of R because it is faster and we use asyBessel only very close to zero to have less calls to besselj.
    if useQ
        for k = 1:T-1
            for i = 1:min(size(UQ,4),9-k)
                R = R + reshape(UQ[1, :, k, i, 4]*z^(i-1)/np^k,(1,2))
            end
        end
    else
        d = z-1+0im
        phi = 2*z-1+2*sqrt(z)*sqrt(d)
        Rko = (1+0im)zeros(T-1,2)
        sL = (1+0im)zeros(2,2,T-1)
	for m = 1:T-1
            for i = 1:ceil(Int64,3*m/2)
                Rko[m,:] += UQ[1, :, m, i, 1]/d^i+UQ[1, :, m, i, 2]/z^i
            end
            sL[:,:,m] = brac(m-1,alpha)/2^(1+2*m)/(npb/2im/np)^m*( [2^(-alpha)  0 ; 0   2^(alpha)]*[sqrt(phi)    1im/sqrt(phi)  ;  -1im/sqrt(phi)   sqrt(phi)]/2/z^(1/2)/d^(1/2)*[(-phi)^(alpha/2)  0 ; 0   (-phi)^(-alpha/2) ]*[((-1)^m)/m*(alpha^2+m/2-1/4)     (m-1/2)*1im  ;  -((-1)^m)*(m-1/2)*1im    (alpha^2+m/2-1/4)/m]*[(-phi)^(-alpha/2)  0 ; 0   (-phi)^(alpha/2) ]*[sqrt(phi)     -1im/sqrt(phi) ;  1im/sqrt(phi)    sqrt(phi)]*[2^(alpha)   0 ; 0 2^(-alpha)] -mod(m+1,2)*(4*alpha^2+2*m-1)/m*eye(2,2) )
        end
        for k = 1:T-1
            R += reshape((Rko[k,:] - sL[1,:,k])/np^k,(1,2))
            for m = 1:k-1
                R -= reshape(reshape(Rko[k-m,:],(1,2))*sL[:,:,m]/np^k,(1,2))
            end
        end
    end
    p = real( sqrt(2*pi)*(-1)^np*sqrt(npb)/(z+0im)^(1/4+alpha/2)/(1 - z+0im)^(1/4)*( (sin( (alpha + 1)/2*acos(2*z - 1+0im) - pi*alpha/2)*R[1] -sin( (alpha - 1)/2*acos(2*z - 1+0im) - pi*alpha/2)*R[2]*1im*4^alpha)*besselj(alpha, npb) + (cos( (alpha + 1)/2*acos(2*z - 1+0im)- pi*alpha/2)*R[1] - cos( (alpha - 1)/2*acos(2*z - 1+0im) - pi*alpha/2)*R[2]*1im*4^alpha)*(besselj(alpha-1, npb) - alpha/npb*besselj(alpha, npb) ) ) )
end

function asyAirygen(np, z, alpha, T::Int64, qm, m::Int64, UQ, fn, useQ::Bool, xin=NaN+NaN*1im)
    d = z - 1.0 +0im
    if T == 1
        return real( 4*sqrt(pi)/(z+0im)^(1/4+alpha/2)/d^(1/4)*(cos( (alpha + 1)/2*acos(2*z - 1+0im) )*fn^(1/4)*airyai(fn) -1im*sin( (alpha + 1)/2*acos(2*z - 1+0im) )*ifelse(angle(z-1) <= 0, -one(z), one(z) )*fn^(-1/4)*airyaiprime(fn) ) )
    end
    R = (1+0im)*[1 0]
    if useQ
        for k = 1:T-1
            for i = 1:min(size(UQ,4),9-k)
                R = R + reshape(UQ[1, :, k, i, 3]*d^(i-1)/np^k,(1,2))
            end
        end
    else
        phi = 2*z-1+2*sqrt(z)*sqrt(d)
        Rko = (1+0im)zeros(T-1,2)
        sR = (1+0im)zeros(2,2,T-1)
	for m = 1:T-1
            for i = 1:ceil(Int64,3*m/2)
                Rko[m,:] += UQ[1, :, m, i, 1]/d^i+UQ[1, :, m, i, 2]/z^i
            end
            sR[:,:,m] = nuk(m)/xin^m*( [2^(-alpha)  0 ; 0   2^(alpha)]*[sqrt(phi)    1im/sqrt(phi)  ;  -1im/sqrt(phi)   sqrt(phi)]/8/z^(1/2)/d^(1/2)*[phi^(alpha/2)  0 ; 0   phi^(-alpha/2) ]*[(-1.0)^m  -m*6im ; 6im*m*(-1)^m    1.0]*[phi^(-alpha/2)  0 ; 0   phi^(alpha/2) ]*[sqrt(phi)     -1im/sqrt(phi) ;  1im/sqrt(phi)    sqrt(phi)]*[2^(alpha)   0 ; 0 2^(-alpha)] - mod(m+1,2)*eye(2,2) )
        end
        for k = 1:T-1
            R += reshape((Rko[k,:] -sR[1,:,k])/np^k,(1,2))
            for m = 1:k-1
                R -= reshape(reshape(Rko[k-m,:],(1,2))*sR[:,:,m]/np^k,(1,2))
            end
        end
    end
    p = real( 4*sqrt(pi)/(z+0im)^(1/4+alpha/2)/d^(1/4)*( (R[1]*cos( (alpha + 1)/2*acos(2*z - 1+0im) ) -cos( (alpha - 1)/2*acos(2*z - 1+0im) )*R[2]*1im*4^alpha)*fn^(1/4)*airyai(fn) + 1im*(-sin( (alpha + 1)/2*acos(2*z - 1+0im) )*R[1] +sin( (alpha - 1)/2*acos(2*z - 1+0im) )*R[2]*1im*4^alpha)*ifelse(angle(z-1) <= 0, -one(z), one(z) )*fn^(-1/4)*airyaiprime(fn) ) )
end

# Additional short functions
function poch(x,n) # pochhammer does not seem to exist yet in Julia
    p = prod(x+(0:(n-1)) )
end
function binom(x,n) # binomial only works for integer x
    b = 1.0
    for i = 1:n
        b *= (x-(n-i))/i
    end
    b
end
function nuk(n)
    nu = -gamma(3*n-1/2)*2^n/27^n/2/n/sqrt(pi)/gamma(n*2)
end
function brac(n,alpha)
    b = prod(4*alpha^2-(2*(1:n)-1).^2 )/(2^(2*n)*gamma(1.0+n))
end

# Compute the W or V-matrices to construct the asymptotic expansion of R.
# Input
#   alpha, qm, m - Factors in the weight function w(x) = x^alpha*exp(-qm*x^m)
#   maxOrder     - The maximum order of the error
#   r            - 1 when computing Wright, -1 when computing Wleft
#   isW          - Whether we compute W(iso V)-matrices
# Output
#   WV           - Coefficient matrices for (z + 1/2 \pm 1/2)^m of Delta_k(z) or s_k(z)
function getV(alpha,qm,m::Int64,maxOrder::Int64,r)
    mo = ceil(Int64, 3*maxOrder/2) + 4
    ns = 0:mo
    f = NaN*zeros(mo+1) # Coefficients in the expansion of \bar{phi}_n(z) or \xi_n(z)
    g = NaN*zeros(maxOrder-1,mo+1)

    A = zeros(m+1)
    for k =0:m
        A[k+1] = prod((2*(1:k)-1.0)/2./(1:k))
    end
    if (r == 1) # Right disk: near z=1
        f = NaN*zeros(mo+2)
        ns = [ns; ns[mo+1]+1] # Extend by one because f(1) = 0 while not for left
        u = zeros(ns[mo+2]+1,ns[mo+2]+2)
        v = zeros(ns[mo+2]+1,ns[mo+2]+2)
        u[1,1] = 1.0
        v[1,1] = 1.0
        for n = [ns; ns[mo+2]+1]
            u[2,n+1] = binom(1/2,n+1)
            v[2,n+1] = binom(1/2,n+2)
        end
        for kt = 2:ns[mo+2]
            for n = ns
                u[kt+1,n+1] = sum(u[kt,(0:n)+1].*u[2,n-(0:n)+1])
                v[kt+1,n+1] = sum(v[kt,(0:n)+1].*v[2,n-(0:n)+1])
            end
        end
        q = zeros(ns[mo+2]+1)
        rr = zeros(ns[mo+2]+1) # Coeffs in the expansion of sqrt(2-2*sqrt(1-w))
        for kt = ns
            for l = 0:kt
                q[kt+1] = q[kt+1] + poch(1/2,kt-l)*u[kt-l+1,l+1]/(-2)^(kt-l)/gamma(1.0+kt-l)/(1+2*(kt-l) )
                rr[kt+1] = rr[kt+1] + binom(1/2,kt-l)*v[kt-l+1,l+1]*2^(kt-l)
            end
        end
        if (m == 1)
            for n = ns
                f[n+1,1] = -2*binom(1/2,n)
                for l = 0:n
                    f[n+1,1] = f[n+1,1] + 2*q[l+1]*rr[n-l+1]
                end
            end
        else
            for j = ns
                f[j+1,1] = 0
                for i=0:min(j,m-1)
                    f[j+1,1] = f[j+1,1] + binom(1/2,j-i)*(-1)^(m-i-1)*gamma(-1/2-i)/gamma(1/2-m)/gamma(m-i)
                end
                f[j+1,1] = -f[j+1,1]/m/A[m+1]
                for l = 0:j
                    f[j+1,1] = f[j+1,1] + 2*q[l+1]*rr[j-l+1]
                end
            end
        end
        if(abs(f[1]) > 10*eps(Float64) )
            error("xi_n should be O( (z-1)^(3/2) ): Expected f[1] to be zero")
        end
        ns = ns[1:mo+1] # Reset ns to its value before computing f's
        g[1,1,1] = -1/f[2,1]
        for n = 1:mo
            g[1,n+1,1] = -sum(reshape(g[1,1:n,1],(n,1)).*reshape(f[(n+2):-1:3,1],(n,1) ) )/f[2,1]
        end
    else # Left disk: near z=0
        if (m == 1)
            for n = ns
                f[n+1,1] = -(binom(1/2,n)*(-1)^n + poch(1/2,n)./(1+2*n)./gamma(1.0+n))
            end
        else
            for n = ns
                f[n+1,1] = 0.0
                for k = 0:min(m-1,n)
                    f[n+1,1] = f[n+1,1] + binom(1/2,n-k)*(-1)^(n-k)*A[m-k]
                end
                f[n+1,1] = -f[n+1,1]/2/m/A[m+1]-poch(1/2,n)./(1+2*n)./gamma(1.0+n)
            end
        end
        g[1,1,1] = 1/f[1,1]
        for n = 1:mo
            g[1,n+1,1] = -sum(reshape(g[1,1:n,1],(n,1)).*reshape(f[(n+1):-1:2,1],(n,1)) )/f[1,1]
        end
    end
    rho = (1+1im)*zeros(2*mo+3,mo+1)
    for n = ns
        rho[2,n+1] = poch(1/2,n)/gamma(1.0+n)/(1+2*n)*(-r)^n
    end
    rho[1,1] = 1
    for i = 2:(maxOrder-1)
        for n = ns
            g[i,n+1] = sum(g[i-1,1:(n+1) ].*g[1,(n+1):-1:1] )
        end
    end
    for i = 2:(mo*2+2)
        for n = ns
            rho[i+1,n+1] = sum(rho[i,1:(n+1) ].*rho[2,(n+1):-1:1] )
        end
    end
    OmOdd = (1+1im)*zeros(mo+1); OmEven = (1+1im)*zeros(mo+1)
    XiOdd = (1+1im)*zeros(mo+1); XiEven = (1+1im)*zeros(mo+1)
    ThOdd = (1+1im)*zeros(mo+1); ThEven = (1+1im)*zeros(mo+1)
    OmO = (1+1im)*zeros(mo+1); OmE = (1+1im)*zeros(mo+1)
    XiO = (1+1im)*zeros(mo+1); XiE = (1+1im)*zeros(mo+1)
    ThO = (1+1im)*zeros(mo+1); ThE = (1+1im)*zeros(mo+1)
    for n = ns
        js = 0:n
        for j = js
            OmOdd[n+1] = OmOdd[n+1] + (-1)^j/gamma(1.0+2.0*j)*(-2*alpha/sqrt(-r+0.0im))^(2*j)*rho[2*j+1,n-j+1]
            XiOdd[n+1] = XiOdd[n+1] + (-1)^j/gamma(1.0+2.0*j)*(-2*(alpha+1)/sqrt(-r+0.0im))^(2*j)*rho[2*j+1,n-j+1]
            ThOdd[n+1] = ThOdd[n+1] + (-1)^j/gamma(1.0+2.0*j)*(-2*(alpha-1)/sqrt(-r+0.0im))^(2*j)*rho[2*j+1,n-j+1]
            OmEven[n+1] = OmEven[n+1] + (-1)^j/gamma(1.0+2*j+1.0)*(-2*alpha/sqrt(-r+0.0im))^(2*j+1)*rho[2*j+2,n-j+1]
            XiEven[n+1] = XiEven[n+1] + (-1)^j/gamma(1.0+2*j+1.0)*(-2*(alpha+1)/sqrt(-r+0.0im))^(2*j+1)*rho[2*j+2,n-j+1]
            ThEven[n+1] = ThEven[n+1] + (-1)^j/gamma(1.0+2*j+1.0)*(-2*(alpha-1)/sqrt(-r+0.0im))^(2*j+1)*rho[2*j+2,n-j+1]
        end
        for j = js
            OmO[n+1] = OmO[n+1] + binom(-1/2,j)*(r)^j*OmOdd[n-j+1]
            XiO[n+1] = XiO[n+1] + binom(-1/2,j)*(r)^j*XiOdd[n-j+1]
            ThO[n+1] = ThO[n+1] + binom(-1/2,j)*(r)^j*ThOdd[n-j+1]
            OmE[n+1] = OmE[n+1] + binom(-1/2,j)*(r)^j*OmEven[n-j+1]
            XiE[n+1] = XiE[n+1] + binom(-1/2,j)*(r)^j*XiEven[n-j+1]
            ThE[n+1] = ThE[n+1] + binom(-1/2,j)*(r)^j*ThEven[n-j+1]
        end
    end
    Ts = (1+1im)*zeros(2,2,mo+1) # = G_{k,n}^{odd/even} depending on k, overwritten on each new k
    WV = (1+1im)*zeros(2,2,maxOrder-1,mo+1)
    for k = 1:(maxOrder-1)
        Ts[:,:,:] = 0
        if r == 1
            if mod(k,2) == 1
                for n = 0:mo
                    Ts[:,:,n+1] = nuk(k)*[-2*(2*binom(-1/2,n-1)*(n>0)+binom(-1/2,n))     2im*4^(-alpha)*binom(-1/2,n)    ;    2im*4^(alpha)*binom(-1/2,n)    (2*(2*binom(-1/2,n-1)*(n>0) +binom(-1/2,n)))] -6*k*nuk(k)*[-2*OmO[n+1]   4^(-alpha)*2im*XiO[n+1]  ;   4^(alpha)*2im*ThO[n+1]    2*OmO[n+1]]
                    WV[:,:,k,n+1] = sum(repeat(reshape(g[k,1:(n+1) ], (1,1,n+1) ), outer=[2,2,1]).*Ts[:,:,(n+1):-1:1],3)/8
                end
            else
                for n = 0:mo
                     Ts[:,:,n+1] = nuk(k)*4*(n==0)*eye(2) +6*k*nuk(k)*[-2im*OmE[n+1]    -2*4^(-alpha)*XiE[n+1]  ;   -2*4^alpha*ThE[n+1]   2im*OmE[n+1]]
                     WV[:,:,k,n+1] = sum(repeat(reshape(g[k,1:(n+1) ], (1,1,n+1) ), outer=[2,2,1]).*Ts[:,:,(n+1):-1:1],3)/8
                end
            end
        else
            if mod(k,2) == 1
                for n = 0:mo
                    Ts[:,:,n+1] = -(alpha^2+k/2-1/4)/k*[-(-1)^n*(2*binom(-1/2,n-1)*(n>0)+binom(-1/2,n))*2    -1im*4^(-alpha)*2*(-1)^n*binom(-1/2,n)  ;  -1im*4^(alpha)*2*(-1)^n*binom(-1/2,n)     ( (-1)^n*(2*binom(-1/2,n-1)*(n>0) +binom(-1/2,n))*2)] - (k-1/2)*[2*OmO[n+1]   4^(-alpha)*2im*XiO[n+1]  ;   4^(alpha)*2im*ThO[n+1]   -2*OmO[n+1]] # binom(-1/2,-1) should be zero
                    WV[:,:,k,n+1] = -(-1)^(ceil(Int64, k/2)+1)*(1im*sqrt(2))^k*(-2+0im)^(-k/2)/4^(k+1)*brac(k-1,alpha)*sum(repeat(reshape(g[k,1:(n+1) ], (1,1,n+1) ), outer=[2,2,1]).*Ts[:,:,(n+1):-1:1],3)
                end
            else
                for n = 0:mo
                    Ts[:,:,n+1] = (alpha^2+k/2-1/4)/k*4*(n==0)*eye(2)  -2*(k-1/2)*[ OmE[n+1]   4^(-alpha)*1im*XiE[n+1]  ;   4^alpha*1im*ThE[n+1]   -OmE[n+1] ]
                    WV[:,:,k,n+1] = -(-1)^(ceil(Int64, k/2)+1)*(1im*sqrt(2))^k*(-2)^(-k/2)/4^(k+1)*brac(k-1,alpha)*sum(repeat(reshape(g[k,1:(n+1) ], (1,1,n+1) ), outer=[2,2,1]).*Ts[:,:,(n+1):-1:1],3)
                end
            end
        end
    end
    WV
end

# Get the U-matrices to construct the asymptotic expansion of R using the procedure with the convolutions with a specified method.
# Input
#   alpha, qm, m - Parts of the weight function
#   maxOrder     - The maximal order of the error
# Output
#   UQ           - Coefficient matrices of R_k(z) for (z-1)^(-m) [Uright], or z^(-m) [Uleft] of R_k^{right}(z) for (z-1)^n [Qright] and of R_k^{left}(z) for z^n [Qleft]
function getUQ(alpha, qm, m::Int64, maxOrder::Int64)
    Vr = getV(alpha, qm, m, maxOrder, 1)
    Vl = getV(alpha, qm, m, maxOrder, -1)
    UQ = (1+1im)*zeros(2,2,maxOrder-1,ceil(Int64,3*maxOrder/2)+2, 4)
    for kt = 0:(maxOrder-2)
        # Uright(:,:,(maxOrder-1)+1,:) will not be used later on because first term in expansions is without U's
        for mt = 0:(ceil(Int64,3*(kt+1)/2)-1)
            UQ[:,:,kt+1,mt+1,1] = Vr[:,:,kt+1,ceil(Int64,3*(kt+1)/2)-mt]
            for j = 0:(kt-1)
                for l = 0:(ceil(Int64,3*(j+1)/2)-mt-1)
                    UQ[:,:,kt+1,mt+1,1] = UQ[:,:,kt+1,mt+1,1] + UQ[:,:,kt-j,l+1,3]*Vr[:,:,j+1,ceil(Int64,3*(j+1)/2)-l-mt]
                end
            end
        end
        for mt = 0:(ceil(Int64,(kt+1)/2)-1)
            UQ[:,:,kt+1,mt+1,2] = Vl[:,:,kt+1,ceil(Int64,(kt+1)/2)-mt]
            for j= 0:(kt-1)
                for l = 0:(ceil(Int64,(j+1)/2)-mt-1)
                    UQ[:,:,kt+1,mt+1,2] = UQ[:,:,kt+1,mt+1,2] + UQ[:,:,kt-j,l+1,4]*Vl[:,:,j+1,ceil(Int64,(j+1)/2)-l-mt]
                end
            end
        end
        for n = 0:(ceil(Int64,3*(maxOrder-kt+1)/2)-1)
            UQ[:,:,kt+1,n+1,3] = -Vr[:,:,kt+1,ceil(Int64,3*(kt+1)/2)+1+n]
            UQ[:,:,kt+1,n+1,4] = -Vl[:,:,kt+1,ceil(Int64,(kt+1)/2)+1+n]
            for i = 0:(ceil(Int64,(kt+1)/2)-1)
                UQ[:,:,kt+1,n+1,3] = UQ[:,:,kt+1,n+1,3] + binom(-i-1,n)*UQ[:,:,kt+1,i+1,2]
            end
            for i = 0:(ceil(Int64,3*(kt+1)/2)-1)
                UQ[:,:,kt+1,n+1,4] = UQ[:,:,kt+1,n+1,4] + binom(-i-1,n)*(-1.0)^(-i-1-n)*UQ[:,:,kt+1,i+1,1]
            end
            for j = 0:(kt-1)
                for l = 0:(ceil(Int64, (j+1)/2)+n)
                    UQ[:,:,kt+1,n+1,4] = UQ[:,:,kt+1,n+1,4] -UQ[:,:,kt-j,l+1,4]*Vl[:,:,j+1,n-l+1+ceil(Int64,(j+1)/2) ]
                end
                for l = 0:(ceil(Int64, 3*(j+1)/2)+n)
                    UQ[:,:,kt+1,n+1,3] = UQ[:,:,kt+1,n+1,3] -UQ[:,:,kt-j,l+1,3]*Vr[:,:,j+1,n-l+1+ceil(Int64, 3*(j+1)/2) ]
                end
            end
        end
    end
    UQ
end



########################## Routines for the explicit expansions ##########################

function laguerreExp( n::Int64, compRepr::Bool, alpha::Float64)

    if compRepr
        # Get a heuristic for the indices where the weights are about above realmin.
        mn = min(ceil(Int64, 17*sqrt(n)), n)
    else
        mn = n
    end
    if ( alpha^2/n > 1 )
        warn("A large alpha may lead to inaccurate results because the weight is low and R(z) is not close to identity.")
    end
    # This is a heuristic for the number of terms in the expansions that follow.
    T = ceil(Int64, 34/log(n) )
    d = 1/(4*n+2*alpha+2);

    itric = max(ceil(Int64, 3.6*n^0.188), 7)
    # Heuristics to switch between Bessel, extrapolation and Airy initial guesses.
    igatt = ceil(Int64, mn + 1.31*n^0.4 - n)

    pt = (4*n -4*((itric+1):(mn -max(igatt,0))) +3)*d;
    t = pi^2/16*(pt -1).^2;
    # This t is not a very good initial approximation of the inverse function of f(y) = (4*n -4*k +3)*d*pi +2*sqrt(y)*sqrt(1-y) -acos(2*y-1);
    for it = 1:6
	t = t - (pt*pi +2*sqrt(t-t.^2) -acos(2*t-1) ).*sqrt(t./(1-t))/2;
    end
    # A second option is to use the Roots package and t = fzero(f, 0, 1); for each k
    # A third option is to improve the initial approximations for k starting from ichebSer = Int64(round(n/10-0.45*alpha+0.3)) with barycentric interpolation on a chebyshev grid as 
#    xk = 0.5+ 0.49*cos((2*(1:16) -1)/2/16*pi)';
#    baryc = [8.392374553903038e+02, -1.573904756067641e+04,  1.217801459223059e+05,  -5.241187241822941e+05,   1.561671400871400e+06, -3.632334022840342e+06,   7.056947602530619e+06,  -1.192343707574069e+07,   1.795854669340786e+07,  -2.446733743951402e+07,   3.036493740995492e+07,  -3.430589960928102e+07,   3.489836980231517e+07,  -3.097800484853177e+07,   2.192179297293409e+07,  -8.079248247153265e+06];
#    for mix = 1:16
#	t[ichebSer+1:end] = t[ichebSer+1:end].*(pt[ichebSer+1:end] -xk[mix]);
#	asum = asum + baryc[mix]./(pt[ichebSer+1:end] -xk[end]);
#    end
#    t[ichebSer+1:end] = t[ichebSer+1:end].*asum;
    # The third option ensures that only 3 Newton iterations are required. However, both are slower in total than the uncommented version, which still takes more than half of the total execution time of laguerreExp. Maybe a planned Clenshaw algorithm can improve this ...


    jak = besselroots(alpha, itric)
    bes = jak*0;
    wbes = 0*bes;

    bulk = t*0;
    wbulk = t*0;

    ak = [-13.69148903521072; -12.828776752865757; -11.93601556323626;    -11.00852430373326; -10.04017434155809; -9.02265085340981; -7.944133587120853;    -6.786708090071759; -5.520559828095551; -4.08794944413097; -2.338107410459767]
    tair = 3*pi/2*( (igatt:-1:12)-0.25) # [DLMF (9.9.6)]
    ak = [-tair.^(2/3).*(1 + 5/48./tair.^2 - 5/36./tair.^4 + 77125/82944./tair.^6     -10856875/6967296./tair.^8); ak[max(1,12-igatt):11] ]
    air = 0*ak;

    if (T >= 7)
        # These higher order terms in the left and bulk region are derived in [Opsomer 2017, in preparation]
	bes = bes + (657*jak.^6 +36*jak.^4*(73*alpha^2-181) +2*jak.^2*(2459*alpha^4 -10750*alpha^2 +14051) + 4*(1493*alpha^6 -9303*alpha^4 +19887*alpha^2 - 12077) )*d^6/2835;
	wbes = wbes + (11944*alpha^6 + 5256*jak.^6 - (5061*alpha^5 + 5085*alpha^4 + 4830*alpha^3 -22724*alpha^2 - 22932*alpha + 39164)*jak.^4 - 74424*alpha^4 + 8*(2459*alpha^4 -10750*alpha^2 + 14051)*jak.^2 + 159096*alpha^2 - 96616)/2835/2*d^6;

        bulk = bulk -d^5/181440*(9216*(21*alpha^6 - 105*alpha^4 + 147*alpha^2 - 31)*t.^10 -69120*(21*alpha^6 - 105*alpha^4 + 147*alpha^2 - 31)*t.^9 + 384*(12285*alpha^6 -61320*alpha^4 + 85785*alpha^2 - 18086)*t.^8 - 64*(136080*alpha^6 - 675675*alpha^4 +943110*alpha^2 - 198743)*t.^7 + 144*(70560*alpha^6 - 345765*alpha^4 +479850*alpha^2 - 101293)*t.^6 + 72576*alpha^6 - (8128512*alpha^6 - 38656800*alpha^4+ 52928064*alpha^2 - 13067711)*t.^5 + 5*(1016064*alpha^6 - 4581360*alpha^4 +6114528*alpha^2 + 113401)*t.^4 - 317520*alpha^4 - 10*(290304*alpha^6 -1245888*alpha^4 + 1620864*alpha^2 - 528065)*t.^3 + 5*(290304*alpha^6 -1234800*alpha^4 + 1598688*alpha^2 - 327031)*t.^2 + 417312*alpha^2 -5*(96768*alpha^6 - 417312*alpha^4 + 544320*alpha^2 - 111509)*t -85616)./(t-1).^8./t.^2;
        wbulk = wbulk + d^6/362880*(9216*(21*alpha^6 - 105*alpha^4 + 147*alpha^2 - 31)*t.^10 -1536*(945*alpha^6 - 4830*alpha^4 + 6825*alpha^2 - 1444)*t.^9 + 384*(11340*alpha^6 -60165*alpha^4 + 86310*alpha^2 - 18289)*t.^8 - 2*(2903040*alpha^6 - 17055360*alpha^4+ 25401600*alpha^2 - 5*alpha - 5252997)*t.^7 - (11753280*alpha^4 - 23506560*alpha^2+ 67*alpha - 13987519)*t.^6 - 290304*alpha^6 + 12*(1016064*alpha^6 -3578400*alpha^4 + 4108608*alpha^2 + 16*alpha + 7100871)*t.^5 - 5*(4064256*alpha^6 -16559424*alpha^4 + 20926080*alpha^2 + 61*alpha - 15239393)*t.^4 + 1270080*alpha^4 +10*(1741824*alpha^6 - 7386624*alpha^4 + 9547776*alpha^2 + 29*alpha - 1560107)*t.^3- 15*(580608*alpha^6 - 2503872*alpha^4 + 3265920*alpha^2 + 11*alpha - 669051)*t.^2- 1669248*alpha^2 + 4*(604800*alpha^6 - 2630880*alpha^4 + 3447360*alpha^2 + 13*alpha- 706850)*t - 7*alpha + 342463)./(t-1).^9./t.^3;
    end
    if (T >= 5)
	bes = bes + (11*jak.^4 +3*jak.^2.*(11*alpha^2-19) +46*alpha^4 -140*alpha^2 +94)*d^4/45;
	wbes = wbes + (46*alpha^4 + 33*jak.^4 +6*jak.^2*(11*alpha^2 -19) -140*alpha^2 +94)/45*d^4;

	air = air -(15152/3031875*ak.^5+1088/121275*ak.^2)*2^(1/3)*d^(7/3) # [Gatteshi 2002 (4.9)], Gives an O(n^{-4}) relative error

        bulk = bulk - d^3/720*(32*(15*alpha^4 - 30*alpha^2 + 7)*t.^6 -144*(15*alpha^4 - 30*alpha^2 + 7)*t.^5 + 16*(225*alpha^4 - 450*alpha^2 +104)*t.^4 - 240*alpha^4 - 480*(5*alpha^4 - 10*alpha^2 + 1)*t.^3 + 480*alpha^2 +45*(16*alpha^4 - 32*alpha^2 + 7)*t + 990*t.^2 - 105)./(t-1).^5./t;
        wbulk = wbulk + d^4/720*(16*(15*alpha^4 - 30*alpha^2 + 7)*t.^6 - 32*(45*alpha^4 - 90*alpha^2 +22)*t.^5 + 48*(75*alpha^4 - 150*alpha^2 + 74)*t.^4 + 240*alpha^4 - 600*(8*alpha^4- 16*alpha^2 - 5)*t.^3 + 45*(80*alpha^4 - 160*alpha^2 + 57)*t.^2 - 480*alpha^2 -90*(16*alpha^4 - 32*alpha^2 + 7)*t + 105)./(t-1).^6./t.^2;
    end
    if (T >= 3)
        # From here, the terms are also in [Tricomi 1947 pg. 296] and 
	bes = bes + (jak.^2 + 2*alpha^2 - 2)*d^2/3;
	wbes = wbes + (alpha^2 + jak.^2 -1)*2/3*d^2;

	air = air +  ak.^2*(d*16)^(1/3)/5 + (11/35-alpha^2-12/175*ak.^3)*d +  (16/1575*ak+92/7875*ak.^4)*2^(2/3)*d^(5/3);
        bulk = bulk - d/12*(4*(3*alpha^2 - 1)*t.^2 +12*alpha^2 - 12*(2*alpha^2 - 1)*t - 3)./(t-1).^2;
        wbulk = wbulk  + d^2/6*(2*t + 3)./(t-1).^3;
    end
    bes = jak.^2*d.*(1 + bes )
    air = 1/d +ak*(d/4)^(-1/3) + air; 
    bulk = bulk + t/d;

    w = [ 4*d*bes.^alpha.*exp(-bes)./(besselj(alpha-1, jak)).^2.*(1+ wbes); bulk.^alpha.*exp(-bulk)*2*pi.*sqrt(t./(1-t)).*(1 +wbulk) ; 2^(2*alpha+1/3)*n^(alpha+1/3)*exp.(-air)./(airyaiprime.(ak)).^2];
    # For the Airy region, only O(n^{-4}) relative error for x and O(n^{-2/3}) for w as the latter are extremely small or even underflow
    x = [ bes; bulk ; air]


    if ( minimum(x) < 0.0 ) || ( maximum(x) > 4*n + 2*alpha + 2 ) ||  ( minimum(diff(x)) <= 0.0 ) || (minimum(w) < 0.0)
	# Temp. disable for testing: error("Wrong node.")
    end

    if compRepr
	k = findlast(w);
	x = x[1:k]
	w = w[1:k]
    end
    x, w
end
