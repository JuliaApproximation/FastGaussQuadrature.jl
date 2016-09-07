gaussjacobi(n::Number, a::Number, b::Number) =
    gaussjacobi(Int(n),Float64(a),Float64(b))


function gaussjacobi(n::Int, a::Float64, b::Float64)
    #GAUSS-JACOBI QUADRATURE NODES AND WEIGHTS
    if n < 0
        error("gaussjacobi($n,$a,$b) not defined: n must be positive.")
    elseif a == 0. && b == 0.
        gausslegendre(n)
    elseif a == -0.5 && b == -0.5
        gausschebyshev(n, 1)
    elseif a == 0.5 && b == 0.5
        gausschebyshev(n, 2)
    elseif a == -0.5 && b == 0.5
        gausschebyshev(n, 3)
    elseif a == 0.5 && b == -0.5
        gausschebyshev(n, 4)
    elseif n == 0
        Float64[], Float64[]
    elseif n == 1
        [(b - a) / (a + b + 2)], [2^(a + b + 1) * beta(a + 1, b + 1)]
    elseif min(a,b) <= -1.
        error("The Jacobi parameters correspond to a nonintegrable weight function")
    elseif n <= 100 && max(a,b) < 5.
        JacobiRec(n, a, b)
    elseif n > 100 && max(a,b) < 5.
        JacobiAsy(n, a, b)
    elseif n <= 4000 && max(a,b)>=5.
        JacobiGW(n, a, b)
    else
        error("gaussjacobi($n,$a,$b) is not implemented: n must be ≤ 4000 for max(a,b)≥5.")
    end
end

function JacobiRec(n::Int, a::Float64, b::Float64)
    #Compute nodes and weights using recurrrence relation.
    x11, x12 = HalfRec(n, a, b, 1)
    x21, x22 = HalfRec(n, b, a, 0)

    x = Array(Float64, n)
    w = Array(Float64, n)
    m1 = length(x11)
    m2 = length(x21)
    sum_w = 0.0
    @inbounds for i in 1:m2
        idx = m2 + 1 - i
        xi = -x21[i]
        der = x22[i]
        wi = 1 / ((1 - xi^2) * der^2)
        w[idx] = wi
        x[idx] = xi
        sum_w += wi
    end
    @inbounds for i in 1:m1
        idx = m2 + i
        xi = x11[i]
        der = x12[i]
        wi = 1 / ((1 - xi^2) * der^2)
        w[idx] = wi
        x[idx] = xi
        sum_w += wi
    end
    c = (2^(a + b + 1) * gamma(2 + a) *
         gamma(2 + b) / (gamma(2 + a + b) * (a + 1) * (b + 1)))
    scale!(w, c / sum_w)
    x, w
end

function HalfRec(n::Int, a::Float64, b::Float64, flag)
    # HALFREC  Jacobi polynomial recurrence relation.
    # Asymptotic formula - only valid for positive x.
    r = (flag == 1) ? (ceil(n / 2):-1:1) : (floor(n / 2):-1:1)
    m = length(r)
    c1 = 1 / (2 * n + a + b + 1)
    a1 = .25 - a^2
    b1 = .25 - b^2
    c1² = c1^2
    x = Array(Float64, m)
    @inbounds for i in 1:m
        C = muladd(2.0, r[i], a - .5) * (π * c1)
        C_2 = 0.5 * C
        x[i] = cos(muladd(c1², muladd(-b1, tan(C_2), a1 * cot(C_2)), C))
    end

    P1 = Array(Float64, m)
    P2 = Array(Float64, m)
    # Loop until convergence:
    for _ in 1:10
        innerJacobiRec!(n, x, a, b, P1, P2)
        dx2 = 0.0
        @inbounds for i in 1:m
            dx = P1[i] / P2[i]
            _dx2 = abs2(dx)
            dx2 = ifelse(_dx2 > dx2, _dx2, dx2)
            x[i] = x[i] - dx
        end
        dx2 > eps(Float64) / 1e6 || break
    end
    # Once more for derivatives:
    innerJacobiRec!(n, x, a, b, P1, P2)
    x, P2
end

function innerJacobiRec!(n, x, a, b, P, PP)
    # EVALUATE JACOBI POLYNOMIALS AND ITS DERIVATIVE USING THREE-TERM RECURRENCE.
    N = length(x)
    @inbounds for j = 1:N
        xj = x[j]

        Pj = .5 * (a - b + (a + b + 2) * xj)
        Pm1 = 1.0
        PPj = .5 * (a + b + 2)
        PPm1 = 0.0
        for k = 1:n-1
            k0 = muladd(2.0, k, a + b)
            k1 = k0 + 1
            k2 = k0 + 2
            A = 2 * (k + 1) * (k + (a + b + 1)) * k0
            B = k1 * (a^2 - b^2)
            C = k0 * k1 * k2
            D = 2 * (k + a) * (k + b) * k2
            c1 = muladd(C, xj, B)
            Pm1, Pj = Pj, muladd(-D, Pm1, c1 * Pj) / A
            PPm1, PPj = PPj, muladd(c1, PPj, muladd(-D, PPm1, C * Pm1)) / A
        end
        P[j] = Pj
        PP[j] = PPj
    end
    nothing
end

function innerJacobiRec(n, x, a, b)
    # EVALUATE JACOBI POLYNOMIALS AND ITS DERIVATIVE USING THREE-TERM RECURRENCE.
    N = length(x)
    P = Array(Float64, N)
    PP = Array(Float64, N)
    innerJacobiRec!(n, x, a, b, P, PP)
    P, PP
end

function weightsConstant(n, a, b)
    # Compute the constant for weights:
    M = min(20, n - 1)
    C = 1.0
    p = -a * b / n
    for m = 1:M
        C += p
        p *= -(m + a) * (m + b) / (m + 1) / (n - m)
        abs(p / C) < eps(Float64) / 100 && break
    end
    2^(a + b + 1) * C
end

function JacobiAsy(n, a, b)
    # ASY  Compute nodes and weights using asymptotic formulae.

    # Determine switch between interior and boundary regions:
    nbdy = 10
    bdyidx1 = n - (nbdy - 1):n
    bdyidx2 = nbdy:-1:1

    # Interior formula:
    x, w = asy1(n, a, b, nbdy)

    # Boundary formula (right):
    xbdy = boundary(n, a, b, nbdy)
    x[bdyidx1], w[bdyidx1] = xbdy

    # Boundary formula (left):
    if a != b
        xbdy = boundary(n, b, a, nbdy)
    end
    x[bdyidx2] = -xbdy[1]
    w[bdyidx2] = xbdy[2]

    scale!(w, weightsConstant(n, a, b))
    x, w
end

function asy1(n::Int, a::Float64, b::Float64, nbdy)
    # Algorithm for computing nodes and weights in the interior.

    # Approximate roots via asymptotic formula: (Gatteschi and Pittaluga, 1985)
    K = (2*(n:-1:1)+a-.5)*pi/(2*n+a+b+1)
    tt = K + 1/(2*n+a+b+1)^2*((.25-a^2)*cot(.5*K)-(.25-b^2)*tan(.5*K))

    # First half (x > 0):
    t = tt[tt .<= pi/2]'
    mint = t[end-nbdy+1]
    idx = 1:max(findfirst(t .< mint)-1, 1)

    dt = 1.0
    counter = 0
    # Newton iteration
    while ( norm(dt,Inf) > sqrt(eps(Float64))/100 && counter < 10 )
        vals = feval_asy1(n, a, b, t, idx)         # Evaluate
        dt = vals[1]./vals[2]                                 # Newton update
        t += dt                                     # Next iterate
        counter += 1
        dt = dt[idx]
    end
    vals = feval_asy1(n, a, b, t, idx)      # Once more for luck
    t += vals[1]./vals[2]                                 # Newton update.

    # Store:
    x = cos(t)
    w = 1./vals[2].^2

    # Second half (x < 0):
    tmp = a; a = b; b = tmp
    t = pi - tt[1:(n-length(x))]'
    mint = t[nbdy]
    idx = max(findfirst(t .> mint), 1):length(t)

    dt = 1.0; counter = 0;
    # Newton iteration
    while ( norm(dt,Inf) > sqrt(eps(Float64))/100 && counter < 10 )
        vals = feval_asy1(n, a, b, t, idx)  # Evaluate.
        dt = vals[1]./vals[2]                                # Newton update.
        t += dt                                     # Next iterate.
        counter += 1
        dt = dt[idx]
    end
    vals = feval_asy1(n, a, b, t, idx)     # Once more for luck.
    t += vals[1]./vals[2]                                 # Newton update.

    # Store:
    [-cos(vec(t));vec(x)],[1./vec(vals[2]).^2;vec(w)]
end

function feval_asy1(n::Int, a::Float64, b::Float64, t, idx)
    # Evaluate the interior asymptotic formula at x = cos(t).

    # Number of terms in the expansion:
    M = 20

    # Some often used vectors/matrices:
    onesT = ones(1,length(t)); onesM = ones(M); MM = (collect(0:M-1)')';

    # The sine and cosine terms:
    alpha = (.5*(2*n+a+b+1+MM))*onesT .* (onesM*t) - .5*(a+.5)*pi
    cosA = cos(alpha); sinA = sin(alpha)

    sinT = onesM*sin(t)
    cosT = onesM*cos(t)
    cosA2 = cosA.*cosT + sinA.*sinT
    sinA2 = sinA.*cosT - cosA.*sinT

    one = ones(t)
    sinT = vcat( one , cumprod(onesM[2:end]*(.5*csc(.5*t))))
    cosT = .5*sec(.5*t)

    j = transpose(0:M-2)
    vec = (.5+a+j).*(.5-a+j)./(j+1)./(2*n+a+b+j+2)
    P1 = [1 cumprod(vec,2)]
    P1[3:4:end] = -P1[3:4:end]
    P1[4:4:end] = -P1[4:4:end]
    P2 = eye(M)
    for l = 0:M-1
        j = transpose(0:(M-l-2))
        vec = (.5+b+j).*(.5-b+j)./(j+1)./(2*n+a+b+j+l+2)
        P2[l+1+(1:length(j)),l+1] = cumprod(vec,2)
    end
    PHI = repmat(P1,M,1).*P2

    j = transpose(0:M-2)
    vec = (.5+a+j).*(.5-a+j)./(j+1)./(2*(n-1)+a+b+j+2)
    P1 = [1 cumprod(vec,2)]
    P1[3:4:end] = -P1[3:4:end]
    P1[4:4:end] = -P1[4:4:end]
    P2 = eye(M)
    for l = 0:M-1
        j = transpose(0:(M-l-2))
        vec = (.5+b+j).*(.5-b+j)./(j+1)./(2*(n-1)+a+b+j+l+2)
        P2[l+1+(1:length(j)),l+1] = cumprod(vec,2)
    end
    PHI2 = repmat(P1,M,1).*P2

    S = 0; S2 = 0;
    SC = sinT
    for m = 1:M
        l = 1:2:m
        phi = PHI[m:m, l]
        dS1 = phi * SC[l, :] .* cosA[m:m, :]
        phi2 = PHI2[m:m, l]
        dS12 = phi2*SC[l, :] .* cosA2[m:m, :]
        l = 2:2:m
        phi = PHI[m:m, l]
        dS2 = phi * SC[l, :] .* sinA[m:m, :]
        phi2 = PHI2[m:m, l]
        dS22 = phi2 * SC[l, :] .* sinA2[m:m, :]
        if m - 1 > 10 && norm(dS1[idx] + dS2[idx], Inf) < eps(Float64) / 100
            break
        end
        S = S + dS1 + dS2
        S2 = S2 + dS12 + dS22
        SC[1:m,:] = broadcast(.*,SC[1:m,:],cosT)
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
    f(g,z) = sum(g.*[1;cumprod(ones(9)./z)])
    C = p2*(f(g,n+a)*f(g,n+b)/f(g,2*n+a+b))*2/pi
    C2 = C*(a+b+2*n).*(a+b+1+2*n)./(4*(a+n).*(b+n))

    vals = C*S
    S2 = C2*S2

    # Use relation for derivative:
    ders = (n*(a-b-(2*n+a+b)*cos(t)).*vals + 2*(n+a)*(n+b)*S2)/(2*n+a+b)./sin(t)
    t = abs( t )
    denom = 1./real(sin(t/2).^(a+.5).*cos(t/2).^(b+.5))
    vals = vals.*denom
    ders = ders.*denom

    return (vals, ders)
end

function boundary(n::Int, a::Float64, b::Float64, npts)
# Algorithm for computing nodes and weights near the boundary.

    # Use Newton iterations to find the first few Bessel roots:
    smallK = min(30, npts)
    jk = besselroots(a, min(npts, smallK))

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

function JacobiGW( n::Int64, a::Float64, b::Float64 )
    # Golub-Welsh for Gauss--Jacobi quadrature. This is used when max(a,b)>5.
    ab = a + b;
    ii = 2:n-1;
    abi = 2*ii + ab;
    aa = [(b - a)/(2 + ab)
          (b^2 - a^2)./((abi - 2).*abi)
          (b^2 - a^2)./((2*n - 2+ab).*(2*n+ab))]
    bb = [2*sqrt( (1 + a)*(1 + b)/(ab + 3))/(ab + 2) ;
          2*sqrt(ii.*(ii + a).*(ii + b).*(ii + ab)./(abi.^2 - 1))./abi]
    TT = diagm(bb,1) + diagm(aa) + diagm(bb,-1) # Jacobi matrix.
    x, V = eig( TT )                       # Eigenvalue decomposition.
    # Quadrature weights:
    w = V[1,:].^2*( 2^(ab+1)*gamma(a+1)*gamma(b+1)/gamma(2+ab) );
    w = w/sum(w);
    x, vec(w)
end
