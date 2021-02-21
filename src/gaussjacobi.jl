@doc raw"""
    gaussjacobi(n::Integer) -> Tuple{Vector{Float64},Vector{Float64}}

Return nodes and weights of [Gauss-Jacobi quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Jacobi_quadrature).

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

function jacobi_asy(n, α, β)
    # ASY  Compute nodes and weights using asymptotic formulae.

    # Determine switch between interior and boundary regions:
    nbdy = 10
    bdyidx1 = n - (nbdy - 1):n
    bdyidx2 = nbdy:-1:1

    # Interior formula:
    x, w = asy1(n, α, β, nbdy)

    # Boundary formula (right):
    xbdy = boundary(n, α, β, nbdy)
    x[bdyidx1], w[bdyidx1] = xbdy

    # Boundary formula (left):
    if α ≠ β
        xbdy = boundary(n, β, α, nbdy)
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

    # Some often used vectors/matrices:
    onesM = ones(M)

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
    f(g,z) = dot(g, [1;cumprod(ones(9)./z)])

    # Float constant C, C2
    C = 2*p2*(f(g,n+α)*f(g,n+β)/f(g,2n+α+β))/π
    C2 = C*(α+β+2n)*(α+β+1+2n)/(4*(α+n)*(β+n))

    vals = C*S

    # Use relation for derivative:
    ders = (n*((α-β).-(2n+α+β)*cos.(t)).*vals .+ (2*(n+α)*(n+β)*C2).*S2)./(2n+α+β)./sin.(t)
    denom = 1 ./ (sin.(abs.(t)/2).^(α+0.5).*cos.(t/2).^(β+0.5))
    vals .*= denom
    ders .*= denom

    return vals, ders
end

function boundary(n::Integer, α::Float64, β::Float64, npts::Integer)
# Algorithm for computing nodes and weights near the boundary.

    # Use Newton iterations to find the first few Bessel roots:
    smallK = min(30, npts)
    jk = approx_besselroots(α, smallK)

    # Approximate roots via asymptotic formula: (see Olver 1974)
    phik = jk/(n + .5*(α + β + 1))
    x = cos.( phik .+ ((α^2-0.25).*(1 .-phik.*cot.(phik))./(8*phik) .- 0.25.*(α^2-β^2).*tan.(0.5.*phik))./(n + 0.5*(α + β + 1))^2 )

    # Newton iteration:
    for _ in 1:10
        vals, ders = innerjacobi_rec(n, x, α, β)  # Evaluate via asymptotic formula.
        dx = -vals./ders  # Newton update.
        x += dx  # Next iterate.
        if norm(dx,Inf) < sqrt(eps(Float64))/200
            break
        end
    end
    vals, ders = innerjacobi_rec(n, x, α, β)  # Evaluate via asymptotic formula.
    dx = -vals./ders
    x += dx

    # flip:
    x = reverse(x)
    ders = reverse(ders)

    # Revert to x-space:
    w = 1 ./ ((1 .- x.^2) .* ders.^2)
    return x, w
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
