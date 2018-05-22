
"""
Compute the n-point Gaussian quadrature rule for the Freud-type weight function
`x^α exp(-q_m x^m)`.
The computations are based on asymptotic expansions for the corresponding orthogonal
polynomials.
"""
function gaussfreud(n, alpha = 0.0, m = 1, qm = 1.0; reduced = false)
    @assert alpha > -1
    @assert n >= 0
    T = ceil(Int64, 34/log(n) ) # Heuristic for number of terms, should be scaled by the logarithm of eps(Float64) over the machine precision.
    UQ0 = getUQ(alpha, qm, m, T)
    if reduced
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
            warn("Node is outside the support of the measure: inaccuracy is expected.")
        end
        if noUnderflow&& useFinDiff
            hh = max(sqrt(eps(Float64))*x[k], sqrt(eps(Float64)) )
            w[k] = hh/(polyAsyRHgen(n,x[k]+hh,alpha,T,qm,m,UQ0) -polyAsyRHgen(n,x[k],alpha,T,qm,m,UQ0))/polyAsyRHgen(n-1,x[k],alpha,T,qm,m,UQ0)/exp(qm*x[k]^m) # This leaves out a constant factor, given by a ratio of leading order coefficients and normalising constants
        elseif noUnderflow
            w[k] = factorw/polyAsyRHgen(n-1, x[k], alpha+1, T, qm, m, UQ1)/polyAsyRHgen(n+1, x[k], alpha, T,qm,m,UQ0)/exp( x[k] )
        end
        if noUnderflow && ( w[k] == 0 ) && ( k > 1 ) && ( w[k-1] > 0 ) # We could stop now as the weights underflow.
            if reduced
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
