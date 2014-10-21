function gaussjacobi(n::Int64, a, b)
#GAUSS-JACOBI QUADRATURE NODES AND WEIGHTS

    if ( a == 0 && b == 0 )
        x = gausslegendre(n)
    elseif ( a == -0.5 && b == -0.5 )
        x = gausschebyshev(n,1)
    elseif ( a == 0.5 && b == 0.5)
        x = gausschebyshev(n,2)
    elseif ( a == -0.5 && b == 0.5)
        x = gausschebyshev(n,3)
    elseif ( a == 0.5 && b == -0.5 )
        x = gausschebyshev(n,4)
    elseif ( n == 0 )
        x = (Float64[],Float64[])
    elseif ( n == 1 )
        x = ([(b-a)/(a+b+2)], [2^(a+b+1)*beta(a+1, b+1)])
    elseif ( n <= 100 )
        x = JacobiRec(n, a, b)
    elseif ( n > 100 )
        x = JacobiAsy(n, a, b)
    else 
        error("1st argument must be a positive integer.")
    end
    return x
end


function JacobiRec(n::Int64, a::Float64, b::Float64)
#Compute nodes and weights using recurrrence relation.
    x1 = HalfRec(n, a, b, 1)
    x2 = HalfRec(n, b, a, 0) 
    x = vcat( -flipud(x2[1]) , x1[1] )
    ders = vcat( flipud(x2[2]) , x1[2] )
    w = 1./((1-x.^2).*ders.^2)
    w = 2^(a+b+1)*gamma(2+a) * gamma(2+b) / (gamma(2+a+b)*(a+1)*(b+1)) * w ./ sum(w)
    return x, w
end

function HalfRec(n::Int64, a::Float64, b::Float64, flag)
#HALFREC   Jacobi polynomial recurrence relation.
    # Asymptotic formula - only valid for positive x.
    r = (flag==1) ? [ceil(n/2):-1:1] : [floor(n/2):-1:1]
    C = (2*r+a-.5)*pi/(2*n+a+b+1)
    x = cos( C + 1/(2*n+a+b+1)^2 * ((.25-a^2)*cot(.5*C) - (.25-b^2)*tan(.5*C)) )
    dx = 1.0 
    counter = 0
    # Loop until convergence:
    while ( norm(dx,Inf) > sqrt(eps(Float64))/1000 && counter < 10 )
        counter = counter + 1
        P = innerJacobiRec(n, x, a, b)
        dx = -P[1]./P[2] 
        x += dx
    end
    # Once more for derivatives:
    P = innerJacobiRec(n, x, a, b)     
    return x, P[2]
end

function innerJacobiRec(n::Int64, x::Array{Float64}, a::Float64, b::Float64)
# EVALUATE JACOBI POLYNOMIALS AND ITS DERIVATIVE USING THREE-TERM RECURRENCE.
    #P, Pm1, PP, PPm1 = [.5*(a-b+(a+b+2)*x)], ones(x), .5*(a+b+2)*ones(x), zeros(x)
    N = length(x)
    P = Array(Float64,N); Pm1 = Array(Float64,N); PP = Array(Float64,N); PPm1 = Array(Float64,N); 
    for j = 1:N
        P[j] = .5*(a-b+(a+b+2)*x[j]); Pm1[j] = 1.0; PP[j] = .5*(a+b+2); PPm1[j] = 0.0
        for k = 1:n-1
            A = 2*(k + 1)*(k + a+b + 1)*(2*k + a+b)
            B = (2*k + a+b + 1)*(a^2 - b^2)
            C = prod(2*k + a+b + [0:2])
            D = 2*(k + a)*(k + b)*(2*k + a+b + 2)
            Pm1[j], P[j] = P[j], ( (B+C*x[j])*P[j] - D*Pm1[j] ) / A
            PPm1[j], PP[j] = PP[j], ( (B+C*x[j])*PP[j] + C*Pm1[j] - D*PPm1[j] ) / A
        end
    end
    return P, PP
end

function weightsConstant( n::Int64, a::Float64, b::Float64 )
    # Compute the constant for weights:
    M = min(20, n-1)
    C = 1.0 
    p::Float64 = -a*b/n
    for m = 1:M
        C += p
        p *= -(m+a)*(m+b)/(m+1)/(n-m)
        if ( abs(p / C) < eps(Float64)/100 ) break end
    end
    return 2^(a+b+1)*C
end

function JacobiAsy(n::Int64, a::Float64, b::Float64)
#ASY   Compute nodes and weights using asymptotic formulae.

    # Determine switch between interior and boundary regions:
    nbdy = 10
    bdyidx1 = n-(nbdy-1):n
    bdyidx2 = nbdy:-1:1

    # Interior formula:
    x = asy1(n, a, b, nbdy) 
    w = x[2]; x = x[1]

    # Boundary formula (right):
    xbdy = boundary(n, a, b, nbdy);  
    x[bdyidx1] = xbdy[1]
    w[bdyidx1] = xbdy[2]
    
    # Boundary formula (left):
    if ( !(a == b) )
        xbdy = boundary(n, b, a, nbdy)  
    end
    x[bdyidx2] = -xbdy[1] 
    w[bdyidx2] = xbdy[2]
    
    w *= weightsConstant(n, a, b)
    return x[:], w[:]
end

function asy1(n::Int64, a::Float64, b::Float64, nbdy)
# Algorithm for computing nodes and weights in the interior.

    # Approximate roots via asymptotic formula: (Gatteschi and Pittaluga, 1985)
    K = (2*[n:-1:1]+a-.5)*pi/(2*n+a+b+1)
    tt = K + 1/(2*n+a+b+1)^2*((.25-a^2)*cot(.5*K)-(.25-b^2)*tan(.5*K))

    # First half (x > 0):
    t = tt[tt .<= pi/2]'
    mint = t[end-nbdy+1]
    idx = 1:max(findfirst(float64(t .< mint))-1, 1)

    dt = 1.0; counter = 0
    # Newton iteration
    while ( norm(dt,Inf) > sqrt(eps(Float64))/100 && counter < 10 )
        vals = feval_asy1(n, a, b, t, idx, 0)         # Evaluate
        dt = vals[1]./vals[2]                                 # Newton update
        t += dt                                     # Next iterate
        counter += 1
        dt = dt[idx]
    end
    vals = feval_asy1(n, a, b, t, idx, 0)      # Once more for luck
    t += vals[1]./vals[2]                                 # Newton update.

    # Store:
    x = cos(t)
    w = 1./vals[2].^2

    # Second half (x < 0):
    tmp = a; a = b; b = tmp
    t = pi - tt[1:(n-length(x))]'
    mint = t[nbdy]
    idx = max(findfirst(float64(t .> mint)), 1):length(t)

    dt = 1.0; counter = 0;
    # Newton iteration
    while ( norm(dt,Inf) > sqrt(eps(Float64))/100 && counter < 10 )
        vals = feval_asy1(n, a, b, t, idx, 0)  # Evaluate.
        dt = vals[1]./vals[2]                                # Newton update.
        t += dt                                     # Next iterate.
        counter += 1
        dt = dt[idx]
    end
    vals = feval_asy1(n, a, b, t, idx, 0)     # Once more for luck.
    t += vals[1]./vals[2]                                 # Newton update.

    # Store:
    x = [-cos(t) x];
    w = [1./vals[2].^2 w];
    return x, w
end

function feval_asy1(n::Int64, a::Float64, b::Float64, t, idx, flag)
# Evaluate the interior asymptotic formula at x = cos(t).
    
    # Number of terms in the expansion:
    M = 20

    # Some often used vectors/matrices:
    onesT = ones(1,length(t)); onesM = ones(M); MM = ([0:M-1]')';

    # The sine and cosine terms:
    alpha = (.5*(2*n+a+b+1+MM))*onesT .* (onesM*t) - .5*(a+.5)*pi
    cosA = cos(alpha); sinA = sin(alpha)

    if ( flag == 1 )
        k = ( idx[1] == 1 ) ? [length(t):-1:1]' : [1:length(t)]'
        ta = float64(float32(t)) 
        tb = t - ta; hi = n*ta; lo = n*tb + (a+b+1)*.5*t
        pia = float64(float32(pi))
        pib = -8.742278000372485e-08; #pib = pi - pia
        dh = (hi - (k-.25)*pia) + lo - .5*a*pia - (k-.25+.5*a)*pib
        tmp = 0; sgn = 1; fact = 1; DH = dh; dh2 = dh.*dh
        for j = 0:20
            dc = sgn*DH/fact
            tmp = tmp + dc
            sgn = -sgn
            fact = fact*(2*j+3)*(2*j+2)
            DH = DH.*dh2
            if ( norm(dc,Inf) < eps(Float64)/2000 )
                break
            end
        end
        tmp = -tmp; tmp[1] = -tmp[1]
        tmp = sign(cosA[1,2]*tmp[2])*tmp
        cosA[1,:] = tmp
    end

    sinT = onesM*sin(t)
    cosT = onesM*cos(t)
    cosA2 = cosA.*cosT + sinA.*sinT
    sinA2 = sinA.*cosT - cosA.*sinT

    one = ones(t)
    sinT = vcat( one , cumprod(onesM[2:end]*(.5*csc(.5*t))))
    cosT = .5*sec(.5*t)

    j = transpose([0:M-2])
    vec = (.5+a+j).*(.5-a+j)./(j+1)./(2*n+a+b+j+2)
    P1 = [1 cumprod(vec,2)]
    P1[3:4:end] = -P1[3:4:end]
    P1[4:4:end] = -P1[4:4:end]
    P2 = eye(M)
    for l = 0:M-1
        j = transpose([0:(M-l-2)])
        vec = (.5+b+j).*(.5-b+j)./(j+1)./(2*n+a+b+j+l+2)
        P2[l+1+(1:length(j)),l+1] = cumprod(vec,2)
    end
    PHI = repmat(P1,M,1).*P2

    j = transpose([0:M-2])
    vec = (.5+a+j).*(.5-a+j)./(j+1)./(2*(n-1)+a+b+j+2)
    P1 = [1 cumprod(vec,2)]
    P1[3:4:end] = -P1[3:4:end]
    P1[4:4:end] = -P1[4:4:end]
    P2 = eye(M)
    for l = 0:M-1
        j = transpose([0:(M-l-2)])
        vec = (.5+b+j).*(.5-b+j)./(j+1)./(2*(n-1)+a+b+j+l+2)
        P2[l+1+(1:length(j)),l+1] = cumprod(vec,2)
    end
    PHI2 = repmat(P1,M,1).*P2

    S = 0; S2 = 0;
    SC = sinT
    for m = 0:M-1
        l = [0:2:m]
        phi = PHI[m+1,l+1]
        dS1 = phi*SC[l+1,:].*cosA[m+1,:]
        phi2 = PHI2[m+1,l+1]
        dS12 = phi2*SC[l+1,:].*cosA2[m+1,:]
        l = [1:2:m]
        phi = PHI[m+1,l+1]
        dS2 = phi*SC[l+1,:].*sinA[m+1,:]
        phi2 = PHI2[m+1,l+1]
        dS22 = phi2*SC[l+1,:].*sinA2[m+1,:]
        if m > 10 && norm(dS1[idx] + dS2[idx], Inf) < eps(Float64)/100 break end
        S = S + dS1 + dS2
        S2 = S2 + dS12 + dS22
        SC[1:m+1,:] = broadcast(.*,SC[1:m+1,:],cosT)
    end

    # Constant out the front:
    dsa = .5*(a^2)/n; dsb = .5*(b^2)/n; dsab = .25*(a+b)^2/n
    ds = dsa + dsb - dsab; s = ds; j = 1
    dsold = ds # to fix a = -b bug.
    while ( (abs(ds/s) + dsold) > eps(Float64)/10 )
        dsold = abs(ds/s)
        j += 1
        tmp = -(j-1)/(j+1)/n
        dsa = tmp*dsa*a
        dsb = tmp*dsb*b
        dsab = .5*tmp*dsab*(a+b)
        ds = dsa + dsb - dsab
        s = s + ds
    end
    p2 = exp(s)*sqrt(2*pi)*sqrt((n+a)*(n+b)/(2*n+a+b))/(2*n+a+b+1)
    g = [1, 1/12, 1/288, -139/51840, -571/2488320, 163879/209018880, 
         5246819/75246796800, -534703531/902961561600,
         -4483131259/86684309913600, 432261921612371/514904800886784000]
    f(g,z) = sum(g.*[1,cumprod(ones(9)./z)])
    C = p2*(f(g,n+a)*f(g,n+b)/f(g,2*n+a+b))*2/pi
    C2 = C*(a+b+2*n).*(a+b+1+2*n)./(4*(a+n).*(b+n))

    vals = C*S
    S2 = C2*S2

    # Use relation for derivative:
    ders = (n*(a-b-(2*n+a+b)*cos(t)).*vals + 2*(n+a)*(n+b)*S2)/(2*n+a+b)./sin(t)
    denom = 1./real(sin(t/2).^(a+.5).*cos(t/2).^(b+.5))
    vals = vals.*denom
    ders = ders.*denom
     
    return (vals, ders)
end

function boundary(n::Int64, a::Float64, b::Float64, npts)
# Algorithm for computing nodes and weights near the boundary.

    # Use Newton iterations to find the first few Bessel roots:
    smallK = min(30, npts);
    jk = besselRoots(a, min(npts, smallK))
    # Use asy formula for larger ones (See NIST 10.21.19, Olver 1974 p247)
    if ( npts > smallK )
        mu = 4*a^2;
        a8 = 8*([length(jk)+1:npts]'+.5*a-.25)*pi
        jk2 = .125*a8-(mu-1)./a8 - 4*(mu-1)*(7*mu-31)/3./a8.^3 - 
          32*(mu-1)*(83*mu.^2-983*mu+3779)/15./a8.^5 -
          64*(mu-1)*(6949*mu^3-153855*mu^2+1585743*mu-6277237)/105./a8.^7
        jk = [jk ; jk2]
    end
    jk = real(jk[1:npts])

    # Approximate roots via asymptotic formula: (see Olver 1974)
    phik = jk/(n + .5*(a + b + 1))
    x = cos( phik + ((a^2-.25)*(1-phik.*cot(phik))./(8*phik) - .25*(a^2-b^2)*tan(.5*phik))/(n + .5*(a + b + 1))^2 )

    dx = 1.0; counter = 0;
    # Newton iteration:
    while ( norm(dx,Inf) > sqrt(eps(Float64))/200 && counter < 10)
        vals = innerJacobiRec(n, x, a, b)   # Evaluate via asymptotic formula.
        dx = -vals[1]./vals[2]                   # Newton update.
        x += dx                             # Next iterate.
        counter += 1
    end
    vals = innerJacobiRec(n, x, a, b);     # Evaluate via asymptotic formula.
    dx = -vals[1]./vals[2]                        # Newton update
    x += dx    
    
    # flip:
    x = x[npts:-1:1]
    ders = vals[2]; vals = vals[1]
    ders = ders[npts:-1:1]

    # Revert to x-space:     
    w = 1./((1-x.^2).*ders.^2)
    return x, w
end

function besselRoots(nu::Float64, m::Int64)
# BESSELROOTS(NU, M) returns the first M roots of besselj(nu, x).
# Find an approximation:
    jk = Array(Float64, m) 
    m1 = 3
    if ( nu == 0 )
        xs = 2.404825557695773
    elseif ( nu > 0 )
        xs = nu + 1.8557*nu^(1/3)
    else
        nu1 = nu + 1;
        # See Piessens 1984:
        xs = 2*sqrt(nu+1)*(1 + nu1/4 - 7*nu1^2/96 + 49*nu1^3/1152 - 8363*nu1/276480);
        m1 = min(max(2*ceil(abs(log10(nu1))), 3), m);
    end

    jk[1] = besselNewton(nu, xs);
    if ( m > 1 )
        jk[2] = besselNewton(nu, jk[1]+.9*pi)
    end
    if ( m > 2)
        for k = 3:m1
            jk[k] = besselNewton(nu, jk[k-1]+.99*pi)
        end
        for k = m1+1:m
            jk[k] = besselNewton(nu, jk[k-1]+pi);
        end
    end
    return jk
end

function besselNewton(nu::Float64, x::Float64)
# BESSELNEWTON(NU, JK)   Find roots of Bessel function using Newton iteration.
    dx = 1.0
    counter = 0
    while ( (dx > sqrt(eps(Float64)/1000)) && counter < 20 )
        u = besselj(nu, x)
        du = besselj(nu-1, x) - nu./x*u
        dx = u./du
        x -= dx
        counter += 1
    end  
    return x
end
