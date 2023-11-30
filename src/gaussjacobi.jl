@doc raw"""
    gaussjacobi(n::Integer, α::Real, β::Real) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Jacobi quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Jacobi_quadrature)
for exponents `α` and `β`.

```math
\int_{-1}^{1} f(x) (1-x)^\alpha(1+x)^\beta dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gaussjacobi(3, 1/3, -1/3);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ 268π/729(√3)
true
```
"""
function gaussjacobi(n::Integer, α::Real, β::Real)
    #GAUSS-JACOBI QUADRATURE NODES AND WEIGHTS
    if n < 0
        throw(DomainError(n, "gaussjacobi($n,$α,$β) not defined: n must be non-negative."))
    elseif α == 0. && β == 0.
        gausslegendre(n)
    elseif α == -0.5 && β == -0.5
        gausschebyshev(n, 1)
    elseif α == 0.5 && β == 0.5
        gausschebyshev(n, 2)
    elseif α == -0.5 && β == 0.5
        gausschebyshev(n, 3)
    elseif α == 0.5 && β == -0.5
        gausschebyshev(n, 4)
    elseif n == 0
        Float64[], Float64[]
    elseif n == 1
        [(β - α) / (α + β + 2)], [2^(α + β + 1) * beta(α + 1, β + 1)]
    elseif min(α,β) ≤ -1.
        throw(DomainError((α,β), "The Jacobi parameters correspond to a nonintegrable weight function"))
    elseif n ≤ 100 && max(α,β) < 5.
        jacobi_rec(n, α, β)
    elseif n > 100 && max(α,β) < 5.
        jacobi_asy(n, α, β)
    elseif n ≤ 4000 && max(α,β) ≥ 5.
        jacobi_gw(n, α, β)
    else
        error("gaussjacobi($n,$α,$β) is not implemented: n must be ≤ 4000 for max(α,β) ≥ 5.")
    end
end

# Convenience function: convert any kind of numbers a and b to a joint floating point type
jacobi_rec(n::Integer, α::Real, β::Real) = jacobi_rec(n, promote(float(α), float(β))...)
jacobi_asy(n::Integer, α::Real, β::Real) = jacobi_asy(n, promote(float(α), float(β))...)

function jacobi_rec(n::Integer, α::T, β::T) where {T <: AbstractFloat}
    #Compute nodes and weights using recurrrence relation.
    x11, x12 = HalfRec(n, α, β, 1)
    x21, x22 = HalfRec(n, β, α, 0)

    x = Array{T}(undef,n)
    w = Array{T}(undef,n)
    m1 = length(x11)
    m2 = length(x21)
    sum_w = zero(T)
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
    c = (2^(α+β+1)*gamma(2+α)*gamma(2+β)/(gamma(2+α+β)*(α+1)*(β+1)))
    rmul!(w, c / sum_w)
    return x, w
end

function HalfRec(n::Integer, α::T, β::T, flag) where {T <: AbstractFloat}
    # HALFREC  Jacobi polynomial recurrence relation.
    # Asymptotic formula - only valid for positive x.
    r = (flag == 1) ? (ceil(n / 2):-1:1) : (floor(n / 2):-1:1)
    m = length(r)
    c1 = 1 / (2 * n + α + β + 1)
    a1 = one(T)/4 - α^2
    b1 = one(T)/4 - β^2
    c1² = c1^2
    x = Array{T}(undef,m)
    @inbounds for i in 1:m
        C = muladd(2*one(T), r[i], α - one(T)/2) * (T(π) * c1)
        C_2 = C/2
        x[i] = cos(muladd(c1², muladd(-b1, tan(C_2), a1 * cot(C_2)), C))
    end

    P1 = Array{T}(undef,m)
    P2 = Array{T}(undef,m)
    # Loop until convergence:
    for _ in 1:10
        innerjacobi_rec!(n, x, α, β, P1, P2)
        dx2 = 0.0
        @inbounds for i in 1:m
            dx = P1[i] / P2[i]
            _dx2 = abs2(dx)
            dx2 = ifelse(_dx2 > dx2, _dx2, dx2)
            x[i] = x[i] - dx
        end
        dx2 > eps(T) / 1e6 || break
    end
    # Once more for derivatives:
    innerjacobi_rec!(n, x, α, β, P1, P2)
    return x, P2
end

function innerjacobi_rec!(n, x, α::T, β::T, P, PP) where {T <: AbstractFloat}
    # EVALUATE JACOBI POLYNOMIALS AND ITS DERIVATIVE USING THREE-TERM RECURRENCE.
    N = length(x)
    @inbounds for j = 1:N
        xj = x[j]

        Pj = (α - β + (α + β + 2) * xj)/2
        Pm1 = one(T)
        PPj = (α + β + 2)/2
        PPm1 = zero(T)
        for k in 1:n-1
            k0 = muladd(2*one(T), k, α + β)
            k1 = k0 + 1
            k2 = k0 + 2
            A = 2 * (k + 1) * (k + (α + β + 1)) * k0
            B = k1 * (α^2 - β^2)
            C = k0 * k1 * k2
            D = 2 * (k + α) * (k + β) * k2
            c1 = muladd(C, xj, B)
            Pm1, Pj = Pj, muladd(-D, Pm1, c1 * Pj) / A
            PPm1, PPj = PPj, muladd(c1, PPj, muladd(-D, PPm1, C * Pm1)) / A
        end
        P[j] = Pj
        PP[j] = PPj
    end
    nothing
end

function innerjacobi_rec(n, x, α::T, β::T) where {T <: AbstractFloat}
    # EVALUATE JACOBI POLYNOMIALS AND ITS DERIVATIVE USING THREE-TERM RECURRENCE.
    N = length(x)
    P = Array{T}(undef,N)
    PP = Array{T}(undef,N)
    innerjacobi_rec!(n, x, α, β, P, PP)
    return P, PP
end

function weightsConstant(n, α, β)
    # Compute the constant for weights:
    M = min(20, n - 1)
    C = 1.0
    p = -α * β / n
    for m = 1:M
        C += p
        p *= -(m + α) * (m + β) / (m + 1) / (n - m)
        abs(p / C) < eps(Float64) / 100 && break
    end
    return 2^(α + β + 1) * C
end

function jacobi_asy(n::Integer, α::Float64, β::Float64)
    # ASY Compute nodes and weights using asymptotic formulae.

    # Determine switch between interior and boundary regions:
    nbdy = 10
    bdyidx1 = n - (nbdy - 1):n
    bdyidx2 = nbdy:-1:1

    # Interior formula:
    x, w = asy1(n, α, β, nbdy)

    # Boundary formula (right):
    xbdy = asy2(n, α, β, nbdy)
    x[bdyidx1], w[bdyidx1] = xbdy

    # Boundary formula (left):
    if α ≠ β
        xbdy = asy2(n, β, α, nbdy)
    end
    x[bdyidx2] = -xbdy[1]
    w[bdyidx2] = xbdy[2]

    rmul!(w, weightsConstant(n, α, β))
    return x, w
end

function asy1(n::Integer, α::Float64, β::Float64, nbdy::Integer)
    # Algorithm for computing nodes and weights in the interior.

    # Approximate roots via asymptotic formula: (Gatteschi and Pittaluga, 1985)
    K = π*(2(n:-1:1).+α.-0.5)/(2n+α+β+1)
    tt = K .+ (1/(2n+α+β+1)^2).*((0.25-α^2).*cot.(K/2).-(0.25-β^2).*tan.(K/2))

    # First half (x > 0):
    t = tt[tt .≤ π/2]
    mint = t[end-nbdy+1]
    idx = 1:max(findfirst(t .< mint)-1, 1)

    # Newton iteration
    for _ in 1:10
        vals, ders = feval_asy1(n, α, β, t, idx)  # Evaluate
        dt = vals./ders
        t += dt  # Next iterate
        if norm(dt[idx],Inf) < sqrt(eps(Float64))/100
            break
        end
    end
    vals, ders = feval_asy1(n, α, β, t, idx)  # Once more for luck
    t += vals./ders

    # Store
    x_right = cos.(t)
    w_right = 1 ./ ders.^2

    # Second half (x < 0):
    α, β = β, α
    t = π .- tt[1:(n-length(x_right))]
    mint = t[nbdy]
    idx = max(findfirst(t .> mint), 1):length(t)

    # Newton iteration
    for _ in 1:10
        vals, ders = feval_asy1(n, α, β, t, idx)  # Evaluate.
        dt = vals./ders  # Newton update.
        t += dt
        if norm(dt[idx],Inf) < sqrt(eps(Float64))/100
            break
        end
    end
    vals, ders = feval_asy1(n, α, β, t, idx)  # Once more for luck.
    t += vals./ders  # Newton update.

    # Store
    x_left = cos.(t)
    w_left = 1 ./ ders.^2

    return vcat(-x_left, x_right), vcat(w_left, w_right)
end

"""
Evaluate the interior asymptotic formula at x = cos(t).
Assumption:
* `length(t) == n ÷ 2`
"""
function feval_asy1(n::Integer, α::Float64, β::Float64, t::AbstractVector, idx)
    # Number of terms in the expansion:
    M = 20

    # Number of elements in t:
    N = length(t)

    # The sine and cosine terms:
    A = repeat((2n+α+β).+(1:M),1,N).*repeat(t',M)/2 .- (α+1/2)*π/2  # M × N matrix
    cosA = cos.(A)
    sinA = sin.(A)

    sinT = repeat(sin.(t)',M)
    cosT = repeat(cos.(t)',M)
    cosA2 = cosA.*cosT .+ sinA.*sinT
    sinA2 = sinA.*cosT .- cosA.*sinT

    sinT = hcat(ones(N), cumprod(repeat((csc.(t/2)/2),1,M-1), dims=2))  # M × N matrix
    secT = sec.(t/2)/2

    _vec = [(α+j-1/2)*(-α+j-1/2)/(2n+α+β+j+1)/j for j in 1:M-1]
    P1 = [1;cumprod(_vec)]
    P1[3:4:end] = -P1[3:4:end]
    P1[4:4:end] = -P1[4:4:end]
    P2 = Matrix(1.0I, M, M)
    for l in 1:M
        _vec = [(β+j-1/2)*(-β+j-1/2)/(2n+α+β+j+l)/j for j in 1:M-l-2]
        P2[l,l+1:M-2] = cumprod(_vec)
    end
    PHI = repeat(P1,1,M).*P2

    _vec = [(α+j-1/2)*(-α+j-1/2)/(2n+α+β+j-1)/j for j in 1:M-1]
    P1 = [1;cumprod(_vec)]
    P1[3:4:end] = -P1[3:4:end]
    P1[4:4:end] = -P1[4:4:end]
    P2 = Matrix(1.0I, M, M)
    for l in 1:M
        _vec = [(β+j-1/2)*(-β+j-1/2)/(2n+α+β+j+l-2)/j for j in 1:M-l-2]
        P2[l,l+1:M-2] = cumprod(_vec)
    end
    PHI2 = repeat(P1,1,M).*P2

    S = zeros(N)
    S2 = zeros(N)
    for m in 1:M
        l = 1:2:m
        phi = PHI[l, m]
        dS1 = (sinT[:, l]*phi) .* cosA[m, :]
        phi2 = PHI2[l, m]
        dS12 = (sinT[:, l]*phi2) .* cosA2[m, :]
        l = 2:2:m
        phi = PHI[l, m]
        dS2 = (sinT[:, l]*phi) .* sinA[m, :]
        phi2 = PHI2[l, m]
        dS22 = (sinT[:, l]*phi2) .* sinA2[m, :]
        if m - 1 > 10 && norm(dS1[idx] + dS2[idx], Inf) < eps(Float64) / 100
            break
        end
        S .+= dS1
        S .+= dS2
        S2 .+= dS12
        S2 .+= dS22
        sinT[:,1:m] .*= secT
    end

    # Constant out the front:
    dsa = α^2/2n
    dsb = β^2/2n
    dsab = (α+β)^2/4n
    ds = dsa + dsb - dsab
    s = ds
    i = 1
    dsold = ds # to fix α = -β bug.
    while abs(ds/s)+dsold > eps(Float64)/10
        dsold = abs(ds/s)
        i += 1
        tmp = -(i-1)/(i+1)/n
        dsa = tmp*dsa*α
        dsb = tmp*dsb*β
        dsab = tmp*dsab*(α+β)/2
        ds = dsa + dsb - dsab
        s = s + ds
    end
    p2 = exp(s)*sqrt(2π*(n+α)*(n+β)/(2n+α+β))/(2n+α+β+1)
    # g is a vector of coefficients in
    # ``\Gamma(z) = \frac{z^{z-1/2}}{\exp(z)}\sqrt{2\pi} \left(\sum_{i} B_i z^{-i}\right)``, where B_{i-1} = g[i].
    # https://math.stackexchange.com/questions/1714423/what-is-the-pattern-of-the-stirling-series
    g = [1, 1/12, 1/288, -139/51840, -571/2488320, 163879/209018880,
         5246819/75246796800, -534703531/902961561600,
         -4483131259/86684309913600, 432261921612371/514904800886784000]
    f(z) = dot(g, [1;cumprod(ones(9)./z)])

    # Float constant C, C2
    C = 2*p2*(f(n+α)*f(n+β)/f(2n+α+β))/π
    C2 = C*(α+β+2n)*(α+β+1+2n)/(4*(α+n)*(β+n))

    vals = C*S

    # Use relation for derivative:
    ders = (n*((α-β).-(2n+α+β)*cos.(t)).*vals .+ (2*(n+α)*(n+β)*C2).*S2)./(2n+α+β)./sin.(t)
    denom = 1 ./ (sin.(abs.(t)/2).^(α+0.5).*cos.(t/2).^(β+0.5))
    vals .*= denom
    ders .*= denom

    return vals, ders
end

function asy2(n::Integer, α::Float64, β::Float64, npts::Integer)
    # Algorithm for computing nodes and weights near the boundary.

    # Use Newton iterations to find the first few Bessel roots:
    smallK = min(30, npts)
    jk = approx_besselroots(α, smallK)
    # Use asy formula for larger ones (See NIST 10.21.19, Olver 1974 p247)
    if (npts > smallK)
        μ  = 4*α^2
        a8 = 8*(transpose(length(jk)+1:npts)+.5*α-.25)*pi
        jk2 = .125*a8-(μ-1)./a8 - 4*(μ-1)*(7*μ-31)/3 ./ a8.^3 -
              32*(μ-1)*(83*μ.^2-982*μ+3779)/15 ./ a8.^5 -
              64*(μ-1)*(6949*μ^3-153855*μ^2+1585743*μ-6277237)/105 ./ a8.^7
        jk = [jk; jk2]
    end
    jk = real(jk[1:npts])

    # Approximate roots via asymptotic formula: (see Olver 1974)
    phik = jk/(n + .5*(α + β + 1))
    t = phik .+ ((α^2-0.25).*(1 .-phik.*cot.(phik))./(2*phik) .- 0.25.*(α^2-β^2).*tan.(0.5.*phik))./(n + .5*(α + β + 1))^2

    # Compute higher terms:
    higherterms = asy2_higherterms(α, β, t, n)

    # Newton iteration:
    for _ in 1:10
        vals, ders = feval_asy2(n, α, β, t, higherterms)  # Evaluate via asymptotic formula.
        dt = vals./ders  # Newton update.
        t += dt  # Next iterate.
        if norm(dt,Inf) < sqrt(eps(Float64))/200
            break
        end
    end
    vals, ders = feval_asy2(n, α, β, t, higherterms)  # Evaluate via asymptotic formula.
    dt = vals./ders
    t += dt

    # flip:
    t = reverse(t)
    ders = reverse(ders)

    # Revert to x-space:
    x = cos.(t)
    w = transpose(1 ./ ders.^2)
    v = sin.(t)./ders

    return x, w
end

"""
Evaluate the boundary asymptotic formula at x = cos(t).
Assumption:
* `length(t) == n ÷ 2`
"""
function feval_asy2(n::Integer, α::Float64, β::Float64, t::AbstractVector, higherterms)
    rho  = n + .5*(α + β + 1) 
    rho2 = n + .5*(α + β - 1)
    A = (.25 - α^2)       
    B = (.25 - β^2)
    
    # Evaluate the Bessel functions:
    Ja = besselj.(α, rho*t)
    Jb = besselj.(α + 1, rho*t)
    Jbb = besselj.(α + 1, rho2*t)
    Jab = besselj.(α, rho2*t)

    # Evaluate functions for recursive definition of coefficients:
    gt = A*(cot.(t/2) .- (2 ./ t)) .- B*tan.(t/2)
    gtdx = A*(2 ./ t.^2 .- .5*csc.(t/2).^2) .- .5*B*sec.(t/2).^2
    tB0 = .25*gt
    A10 = α*(A+3*B)/24
    A1 = gtdx/8 .- (1+2*α)/8*gt./t .- gt.^2/32 .- A10
    # Higher terms:
    tB1, A2, tB2, A3 = higherterms
    tB1t = tB1(t) 
    A2t  = A2(t) 

    # VALS:
    vals  = Ja + Jb.*tB0/rho + Ja.*A1/rho^2 + Jb.*tB1t/rho^3 + Ja.*A2t/rho^4
    # DERS:
    vals2 = Jab + Jbb.*tB0/rho2 + Jab.*A1/rho2^2 + Jbb.*tB1t/rho2^3 + Jab.*A2t/rho2^4
    
    # Higher terms (not needed for n > 1000).
    tB2t    = tB2(t) 
    A3t     = A3(t)
    vals  .+= Jb.*tB2t/rho^5 + Ja.*A3t/rho^6
    vals2 .+= Jbb.*tB2t/rho2^5 + Jab.*A3t/rho2^6

    # Constant out the front (Computed accurately!)
    ds = .5*(α^2)/n
    s = ds 
    jj = 1
    while abs(ds/s) > eps(Float64)/10
        jj = jj+1
        ds = -(jj-1)/(jj+1)/n*(ds*α)
        s = s + ds
    end
    p2 = exp(s)*sqrt((n+α)/n)*(n/rho)^α
    g = [1, 1/12, 1/288, -139/51840, -571/2488320, 163879/209018880,
         5246819/75246796800, -534703531/902961561600,
         -4483131259/86684309913600, 432261921612371/514904800886784000]
    f(z) = dot(g, [1;cumprod(ones(9)./z)])
    C = p2*(f(n+α)/f(n))/sqrt(2)

    # Scaling:
    valstmp = C*vals
    denom = sin.(t/2).^(α+.5).*cos.(t/2).^(β+.5)
    vals = sqrt.(t).*valstmp./denom

    # Relation for derivative:
    C2 = C*n/(n+α)*(rho/rho2)^α
    ders = (n*(α-β .- (2n+α+β)*cos.(t)).*valstmp .+ 2*(n+α)*(n+β)*C2*vals2)/(2n+α+β)
    ders = ders.*(sqrt.(t)./(denom.*sin.(t)))

    return vals, ders
end

function asy2_higherterms(α::Float64, β::Float64, theta::AbstractVector, n::Integer)
    # Higher-order terms for boundary asymptotic series.
    # Compute the higher order terms in asy2 boundary formula. See [2]. 
    
    # These constants are more useful than α and β:
    A = (0.25 - α^2)
    B = (0.25 - β^2)
    
    # For now, just work on half of the domain:
    c = max(maximum(theta), .5)
    
    # Scaled 2nd-kind Chebyshev points and barycentric weights:
    t = .5*c*(sin.(pi*(-(N-1):2:(N-1))/(2*(N-1))) .+ 1)
    v = [.5; ones(N-1,1)]
    v[2:2:end] .= -1
    v[end]     *= .5
    
    # The g's:
    g   = A*(cot.(t/2) - 2 ./t) - B*tan.(t/2)
    gp  = A*(2 ./ t.^2 - .5*csc.(t/2).^2) - .5*(.25-β^2)*sec.(t/2).^2
    gpp = A*(-4 ./ t.^3 + .25*sin.(t).*csc.(t/2).^4) - 4*B*sin.(t/2).^4 .* csc.(t).^3
    g[1]   = 0 
    gp[1]  = -A/6-.5*B 
    gpp[1] = 0
    
    # B0:
    B0     = .25*g./t
    B0p    = .25*(gp./t - g./t.^2)
    B0[1]  = .25*(-A/6-.5*B)
    B0p[1] = 0
    
    # A1:
    A10   = α*(A+3*B)/24
    A1    = .125*gp .- (1+2*α)/2*B0 .- g.^2/32 .- A10
    A1p   = .125*gpp .- (1+2*α)/2*B0p .- gp.*g/16
    A1p_t = A1p./t
    A1p_t[1] = -A/720 - A^2/576 - A*B/96 - B^2/64 - B/48 + α*(A/720 + B/48)
    
    # Make f accurately: (Taylor series approx for small t)
    fcos = B./(2*cos.(t/2)).^2
    f = -A*(1/12 .+ t.^2/240+t.^4/6048 + t.^6/172800 + t.^8/5322240 +
        691*t.^10/118879488000 + t.^12/5748019200 +
        3617*t.^14/711374856192000 + 43867*t.^16/300534953951232000)
    idx = t .> .5
    f[idx] = A.*(1 ./ t[idx].^2 - 1 ./ (2*sin.(t[idx]/2)).^2)
    f .-= fcos
    
    # Integrals for B1: (Note that N isn't large, so we don't need to be fancy).
    C = (.5*c)*cumsummatN
    D = (2/c)*diffmatN
    I = (C*A1p_t)
    J = (C*(f.*A1))
    
    # B1:
    tB1    = -.5*A1p - (.5+α)*I + .5*J
    tB1[1] = 0
    B1     = tB1./t
    B1[1]  = A/720 + A^2/576 + A*B/96 + B^2/64 + B/48 +
        α*(A^2/576 + B^2/64 + A*B/96) - α^2*(A/720 + B/48)
    
    # A2:
    K  = C*(f.*tB1)
    A2 = .5*(D*tB1) - (.5+α)*B1 - .5*K
    A2 .-= A2[1]
    
    # A2p:
    A2p = D*A2
    A2p .-= A2p[1]
    A2p_t = A2p./t
    # Extrapolate point at t = 0:
    w = pi/2 .- t[2:end]
    w[2:2:end] = -w[2:2:end]
    w[end]     = .5*w[end]
    A2p_t[1] = sum(w.*A2p_t[2:end])/sum(w)
    
    # B2:
    tB2 = -.5*A2p - (.5+α)*(C*A2p_t) + .5*C*(f.*A2)
    B2  = tB2./t
    # Extrapolate point at t = 0:
    B2[1] = sum(w.*B2[2:end])/sum(w)
    
    # A3:
    K  = C*(f.*tB2)
    A3 = .5*(D*tB2) - (.5+α)*B2 - .5*K
    A3 .-= A3[1]
    
    tB1f(theta) = bary(theta, tB1, t, v)
    A2f(theta)  = bary(theta, A2, t, v)
    tB2f(theta) = bary(theta, tB2, t, v)
    A3f(theta)  = bary(theta, A3, t, v)

    return (tB1f, A2f, tB2f, A3f)
end

function jacobi_jacobimatrix(n, α, β)
    s = α + β
    ii = 2:n-1
    si = 2*ii .+ s
    aa = [(β - α)/(2 + s);
          (β^2 - α^2) ./ ((si .- 2).*si);
          (β^2 - α^2) ./ ((2n - 2+s).*(2n+s))]
    bb = [2*sqrt( (1 + α)*(1 + β)/(s + 3))/(s + 2) ;
          2 .*sqrt.(ii.*(ii .+ α).*(ii .+ β).*(ii .+ s)./(si.^2 .- 1))./si]
    return SymTridiagonal(aa, bb)
end

function jacobimoment(α, β)
    s = α + β
    T = float(typeof(s))
    # Same as 2^(α+β+1) * beta(α+1,β+1)
    return exp((s+1)*log(convert(T,2)) + loggamma(α+1)+loggamma(β+1)-loggamma(2+s))
end

function jacobi_gw(n::Integer, α, β)
    # Golub-Welsh for Gauss--Jacobi quadrature. This is used when max(α,β)>5.
    x, V = eigen(jacobi_jacobimatrix(n, α, β))  # Eigenvalue decomposition.
    # Quadrature weights:
    w = V[1,:].^2 .* jacobimoment(α, β)
    return x, w
end

function bary(x, fvals, xk, vk)
    # barycentric interpolation taken from chebfun/bary.m

    # Initialise return value:
    fx = zeros(length(x))

    # Loop:
    for j in eachindex(x)
        xx = vk ./ (x[j] .- xk)
        fx[j] = dot(xx, fvals) / sum(xx)
    end

    # Try to clean up NaNs:
    for k in findall(isnan.(fx))
        index = findfirst(x[k] .== xk) # Find the corresponding node
        if !isempty(index)
            fx[k] = fvals[index]       # Set entries to the correct value
        end
    end

    return fx
end

# Chebyshev type 2 integration and differentiation matrices from chebfun
const N = 10
const cumsummatN = [
        0 0 0 0 0 0 0 0 0 0;
        0.019080722834519 0.0496969890549313 -0.0150585059796021 0.0126377679164575 -0.0118760811432484 0.0115424841953298 -0.0113725236133433 0.0112812076497144 -0.011235316890839 0.00561063519017238;
        0.000812345683614654 0.14586999854807 0.0976007154946748 -0.0146972757610091 0.00680984376276729 -0.00401953146146086 0.00271970678005437 -0.00205195604894289 0.00172405556686793 -0.000812345683614662;
        0.017554012345679 0.103818185816131 0.249384588781868 0.149559082892416 -0.0321899366961563 0.0210262631520163 -0.0171075837742504 0.0153341224604243 -0.0145160806571407 0.00713734567901234;
        0.00286927716087872 0.136593368810421 0.201074970443365 0.339479954969535 0.164397864607267 -0.0260484364615523 0.0127399306249393 -0.00815620454308202 0.00627037388217603 -0.00286927716087872;
        0.0152149561732244 0.110297082689861 0.233440527881186 0.289200104648429 0.369910942265696 0.179464641196877 -0.0375399196961666 0.0242093528947391 -0.0200259122383839 0.00947640185146695;
        0.00520833333333334 0.131083537229178 0.20995020087768 0.319047619047619 0.322836242652128 0.376052442500301 0.152380952380952 -0.024100265443764 0.0127492707559062 -0.00520833333333333;
        0.0131580246959603 0.114843401005169 0.227336279387047 0.299220328493314 0.347882037265605 0.337052662041377 0.316637311034378 0.12768360784343 -0.0293025419760333 0.011533333328731;
        0.00673504382217329 0.127802773462876 0.21400311568839 0.313312558886712 0.332320021608814 0.355738586947393 0.289302267356911 0.240342829317707 0.0668704675171058 -0.00673504382217329;
        0.0123456790123457 0.116567456572037 0.225284323338104 0.301940035273369 0.343862505804144 0.343862505804144 0.301940035273369 0.225284323338104 0.116567456572037 0.0123456790123457
    ]
const diffmatN = [
        -27.1666666666667 33.1634374775264 -8.54863217041303 4 -2.42027662546121 1.70408819104185 -1.33333333333333 1.13247433143179 -1.03109120412576 0.5;
        -8.29085936938159 4.01654328417507 5.75877048314363 -2.27431608520652 1.30540728933228 -0.898197570222574 0.694592710667722 -0.586256827714545 0.532088886237956 -0.257772801031441;
        2.13715804260326 -5.75877048314363 0.927019729872654 3.75877048314364 -1.68805925749197 1.06417777247591 -0.789861687269397 0.652703644666139 -0.586256827714545 0.283118582857949;
        -1 2.27431608520652 -3.75877048314364 0.333333333333335 3.06417777247591 -1.48445439793712 1 -0.789861687269397 0.694592710667722 -0.333333333333333;
        0.605069156365302 -1.30540728933228 1.68805925749197 -3.06417777247591 0.0895235543024196 2.87938524157182 -1.48445439793712 1.06417777247591 -0.898197570222574 0.426022047760462;
        -0.426022047760462 0.898197570222574 -1.06417777247591 1.48445439793712 -2.87938524157182 -0.0895235543024196 3.06417777247591 -1.68805925749197 1.30540728933228 -0.605069156365302;
        0.333333333333333 -0.694592710667722 0.789861687269397 -1 1.48445439793712 -3.06417777247591 -0.333333333333335 3.75877048314364 -2.27431608520652 1;
        -0.283118582857949 0.586256827714545 -0.652703644666139 0.789861687269397 -1.06417777247591 1.68805925749197 -3.75877048314364 -0.927019729872654 5.75877048314363 -2.13715804260326;
        0.257772801031441 -0.532088886237956 0.586256827714545 -0.694592710667722 0.898197570222574 -1.30540728933228 2.27431608520652 -5.75877048314363 -4.01654328417507 8.29085936938159;
        -0.5 1.03109120412576 -1.13247433143179 1.33333333333333 -1.70408819104185 2.42027662546121 -4 8.54863217041303 -33.1634374775264 27.1666666666667
    ]