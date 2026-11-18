@doc raw"""
    gaussjacobi(n::Integer, α::Real, β::Real) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Jacobi quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Jacobi_quadrature)
for exponents `α` and `β`.
The parameters must satisfy `n ≥ 0` and `α, β > -1`. The output element type
is the common type of `promote(float(α), float(β))`. For example, pass two
`Float32` parameters for `Float32` outputs, or two `BigFloat` parameters for
arbitrary precision. The iterative method uses the requested precision at 24
significand bits and above. Lower-precision rules compute in `Float32` before
conversion to the output type.

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
    α, β = promote(float(α), float(β))
    T = typeof(α)
    if n < 0
        throw(DomainError(n, "gaussjacobi($n,$α,$β) not defined: n must be non-negative."))
    elseif n == 0
        x, w = T[], T[]
    elseif precision(T) < precision(Float32)
        x, w = gaussjacobi(n, Float32(α), Float32(β))
        x, w = T.(x), T.(w)
    elseif iszero(α) && iszero(β)
        x, w = gausslegendre(T, n)
    elseif abs(α) == abs(β) == 0.5
        if α == -0.5 && β == -0.5
            x, w = gausschebyshevt(T, n)
        elseif α == 0.5 && β == 0.5
            x, w = gausschebyshevu(T, n)
        elseif α == -0.5 && β == 0.5
            x, w = gausschebyshevv(T, n)
        else
            x, w = gausschebyshevw(T, n)
        end
    elseif n == 1
        x, w = [(β - α) / ((α + 1) + (β + 1))], [jacobimoment(α, β)]
    elseif min(α, β) ≤ -1.0
        throw(DomainError((α, β), "The Jacobi parameters correspond to a nonintegrable weight function"))
    elseif T === Float64 && n ≥ 768 && min(α, β) ≥ -0.99 && max(α, β) < 5
        x, w = jacobi_asy(n, α, β)
    elseif n ≤ (T === Float64 ? 24 : 20) && min(α, β) ≥ -0.5 && max(α, β) ≤ 2
        x, w = jacobi_rec(n, α, β)
    else
        x, w = jacobi_gst(n, α, β)
    end
    validate_rule(x, w)
    return x, w
end

# Convenience function: convert any kind of numbers a and b to a joint floating point type
jacobi_rec(n::Integer, α::Real, β::Real) = jacobi_rec(n, promote(float(α), float(β))...)
jacobi_asy(n::Integer, α::Real, β::Real) = jacobi_asy(n, promote(float(α), float(β))...)

function jacobi_rec(n::Integer, α::T, β::T) where {T <: AbstractFloat}
    #Compute nodes and weights using recurrrence relation.
    x11, x12 = HalfRec(n, α, β, 1)
    if α == β
        first = isodd(n) ? 2 : 1
        x21, x22 = @view(x11[first:end]), @view(x12[first:end])
    else
        x21, x22 = HalfRec(n, β, α, 0)
    end

    x = Array{T}(undef, n)
    w = Array{T}(undef, n)
    m2 = length(x21)
    sum_w = zero(T)
    # At a root, the weight is proportional to 1 / ((1-x²)P_n'(x)²).
    # Normalize the completed weights by the exact zeroth moment below.
    @inbounds for i in 1:m2
        idx = m2 + 1 - i
        x[idx] = -x21[i]
        w[idx] = 1 / ((1 - x[idx]^2) * x22[i]^2)
        sum_w += w[idx]
    end
    @inbounds for i in 1:length(x11)
        idx = m2 + i
        x[idx] = x11[i]
        w[idx] = 1 / ((1 - x[idx]^2) * x12[i]^2)
        sum_w += w[idx]
    end
    rmul!(w, (2^(α + β + 1) * gamma(2 + α) * gamma(2 + β) / (gamma(2 + α + β) * (α + 1) * (β + 1))) / sum_w)
    validate_rule(x, w)
    return x, w
end

function HalfRec(n::Integer, α::T, β::T, flag) where {T <: AbstractFloat}
    # HALFREC  Jacobi polynomial recurrence relation.
    # Asymptotic formula - only valid for positive x.
    r = (flag == 1) ? (cld(n, 2):-1:1) : (fld(n, 2):-1:1)
    m = length(r)
    x = Array{T}(undef, m)
    @inbounds for i in 1:m
        # Correct the leading phase by the first asymptotic cotangent/tangent term.
        C = muladd(2 * one(T), r[i], α - one(T) / 2) * (T(π) * (1 / (2 * n + α + β + 1)))
        x[i] = cos(
            muladd(
                (1 / (2 * n + α + β + 1))^2,
                muladd(-(one(T) / 4 - β^2), tan(C / 2), (one(T) / 4 - α^2) * cot(C / 2)),
                C
            )
        )
    end
    if α == β && isodd(n) && flag == 1
        x[1] = zero(T)
    end

    P1 = Array{T}(undef, m)
    P2 = Array{T}(undef, m)
    # Quadratic Newton convergence doubles the correct bits per step.
    for _ in 1:max(10, ndigits(precision(T); base = 2) + 2)
        innerjacobi_rec!(n, x, α, β, P1, P2)
        dx2 = zero(T)
        @inbounds for i in 1:m
            dx = P1[i] / P2[i]
            _dx2 = abs2(dx)
            dx2 = ifelse(_dx2 > dx2, _dx2, dx2)
            x[i] = x[i] - dx
        end
        dx2 > (T === Float64 ? eps(T) / 1.0e6 : eps(T)^2) || break
    end
    # Once more for derivatives:
    innerjacobi_rec!(n, x, α, β, P1, P2)
    return x, P2
end

function innerjacobi_rec!(n, x, α::T, β::T, P, PP) where {T <: AbstractFloat}
    # EVALUATE JACOBI POLYNOMIALS AND ITS DERIVATIVE USING THREE-TERM RECURRENCE.
    @inbounds for j in 1:length(x)
        xj = x[j]

        Pj = (α - β + (α + β + 2) * xj) / 2
        Pm1 = one(T)
        PPj = (α + β + 2) / 2
        PPm1 = zero(T)
        for k in 1:(n - 1)
            k0 = muladd(2 * one(T), k, α + β)
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
    return nothing
end

function weightsConstant(n, α, β)
    # Compute the constant for weights:
    C = 1.0
    p = -α * β / n
    for m in 1:min(20, n - 1)
        C += p
        p *= -(m + α) * (m + β) / (m + 1) / (n - m)
        abs(p / C) < eps(Float64) / 100 && break
    end
    return 2^(α + β + 1) * C
end

function jacobi_asy(n::Integer, α::Float64, β::Float64, transform = Val(:none))
    # ASY  Compute nodes and weights using asymptotic formulae.

    # Determine switch between interior and boundary regions:
    nbdy = 10
    # Interior formula:
    x, w = asy1(n, α, β, nbdy, transform)

    # Boundary formula (right):
    xbdy = asy2(n, α, β, nbdy, transform)
    x[(n - (nbdy - 1)):n], w[(n - (nbdy - 1)):n] = xbdy

    # Boundary formula (left):
    if α ≠ β || transform isa Val{:radau}
        xbdy = asy2(n, β, α, nbdy, transform, true)
    end
    x[nbdy:-1:1] = -xbdy[1]
    w[nbdy:-1:1] = xbdy[2]

    rmul!(w, weightsConstant(n, α, β))
    validate_rule(x, w)
    return x, w
end

# Transform Jacobi weights before discarding the angles. Subtracting a rounded
# endpoint node from ±1 loses relative accuracy in Radau and Lobatto weights.
@inline jacobi_weight_divisor(t, ::Val{:none}, reflected) = one(t)
@inline jacobi_weight_divisor(t, ::Val{:lobatto}, reflected) = sin(t)^2
@inline function jacobi_weight_divisor(t, ::Val{:radau}, reflected)
    return reflected ? 2sin(t / 2)^2 : 2cos(t / 2)^2
end

function asy1(n::Integer, α::Float64, β::Float64, nbdy::Integer, transform = Val(:none))
    # Algorithm for computing nodes and weights in the interior.

    # Approximate roots via asymptotic formula: (Gatteschi and Pittaluga, 1985)
    K = π * (2(n:-1:1) .+ α .- 0.5) / (2n + α + β + 1)
    tt = K .+ (1 / (2n + α + β + 1)^2) .* ((0.25 - α^2) .* cot.(K / 2) .- (0.25 - β^2) .* tan.(K / 2))

    # First half (x > 0):
    t = α == β ? tt[(n ÷ 2 + 1):end] : tt[tt .≤ π / 2]
    boundary = findfirst(t .< t[end - nbdy + 1])
    boundary === nothing && error("Asymptotic Jacobi nodes do not resolve the boundary region.")
    idx = 1:max(boundary - 1, 1)

    # Newton iteration
    for _ in 1:10
        vals, ders = feval_asy1(n, α, β, t, idx)  # Evaluate
        dt = vals ./ ders
        t += dt  # Next iterate
        if maximum(abs, dt[idx]; init = 0.0) < sqrt(eps(Float64)) / 100
            break
        end
    end
    vals, ders = feval_asy1(n, α, β, t, idx)  # Once more for luck
    t += vals ./ ders

    # Store
    x_right = cos.(t)
    w_right = 1 ./ (ders .^ 2 .* jacobi_weight_divisor.(t, Ref(transform), false))

    if α == β && !(transform isa Val{:radau})
        # Reflection is exact for symmetric Jacobi rules, including the centre.
        first = isodd(n) ? 2 : 1
        isodd(n) && (x_right[1] = 0.0)
        return vcat(-reverse(@view(x_right[first:end])), x_right),
            vcat(reverse(@view(w_right[first:end])), w_right)
    end

    # Second half (x < 0):
    α, β = β, α
    t = π .- tt[1:(n - length(x_right))]
    boundary = findfirst(t .> t[nbdy])
    boundary === nothing && error("Asymptotic Jacobi nodes do not resolve the boundary region.")
    idx = max(boundary, 1):length(t)

    # Newton iteration
    for _ in 1:10
        vals, ders = feval_asy1(n, α, β, t, idx)  # Evaluate.
        dt = vals ./ ders  # Newton update.
        t += dt
        if maximum(abs, dt[idx]; init = 0.0) < sqrt(eps(Float64)) / 100
            break
        end
    end
    vals, ders = feval_asy1(n, α, β, t, idx)  # Once more for luck.
    t += vals ./ ders  # Newton update.

    # Store
    x_left = cos.(t)
    w_left = 1 ./ (ders .^ 2 .* jacobi_weight_divisor.(t, Ref(transform), true))

    return vcat(-x_left, x_right), vcat(w_left, w_right)
end

"""
Evaluate the interior asymptotic formula at x = cos(t).
"""
function feval_asy1(n::Integer, α::Float64, β::Float64, t::AbstractVector, idx)
    # Number of terms in the expansion:
    M = 20

    # Number of elements in t:
    N = length(t)

    # The sine and cosine terms:
    A = repeat((2n + α + β) .+ (1:M), 1, N) .* repeat(t', M) / 2 .- (α + 1 / 2) * π / 2  # M × N matrix
    cosA = cos.(A)
    sinA = sin.(A)

    sinT = repeat(sin.(t)', M)
    cosT = repeat(cos.(t)', M)
    cosA2 = cosA .* cosT .+ sinA .* sinT
    sinA2 = sinA .* cosT .- cosA .* sinT

    sinT = hcat(ones(N), cumprod(repeat((csc.(t / 2) / 2), 1, M - 1), dims = 2))  # M × N matrix
    secT = sec.(t / 2) / 2

    _vec = [(α + j - 1 / 2) * (-α + j - 1 / 2) / (2n + α + β + j + 1) / j for j in 1:(M - 1)]
    P1 = [1;cumprod(_vec)]
    P1[3:4:end] = -P1[3:4:end]
    P1[4:4:end] = -P1[4:4:end]
    P2 = Matrix(1.0I, M, M)
    for l in 1:M
        _vec = [(β + j - 1 / 2) * (-β + j - 1 / 2) / (2n + α + β + j + l) / j for j in 1:(M - l - 2)]
        P2[l, (l + 1):(M - 2)] = cumprod(_vec)
    end
    PHI = repeat(P1, 1, M) .* P2

    _vec = [(α + j - 1 / 2) * (-α + j - 1 / 2) / (2n + α + β + j - 1) / j for j in 1:(M - 1)]
    P1 = [1;cumprod(_vec)]
    P1[3:4:end] = -P1[3:4:end]
    P1[4:4:end] = -P1[4:4:end]
    P2 = Matrix(1.0I, M, M)
    for l in 1:M
        _vec = [(β + j - 1 / 2) * (-β + j - 1 / 2) / (2n + α + β + j + l - 2) / j for j in 1:(M - l - 2)]
        P2[l, (l + 1):(M - 2)] = cumprod(_vec)
    end
    PHI2 = repeat(P1, 1, M) .* P2

    S = zeros(N)
    S2 = zeros(N)
    for m in 1:M
        l = 1:2:m
        phi = PHI[l, m]
        dS1 = (sinT[:, l] * phi) .* cosA[m, :]
        phi2 = PHI2[l, m]
        dS12 = (sinT[:, l] * phi2) .* cosA2[m, :]
        l = 2:2:m
        phi = PHI[l, m]
        dS2 = (sinT[:, l] * phi) .* sinA[m, :]
        phi2 = PHI2[l, m]
        dS22 = (sinT[:, l] * phi2) .* sinA2[m, :]
        if m - 1 > 10 && maximum(abs, dS1[idx] + dS2[idx]; init = 0.0) < eps(Float64) / 100
            break
        end
        S .+= dS1
        S .+= dS2
        S2 .+= dS12
        S2 .+= dS22
        sinT[:, 1:m] .*= secT
    end

    # Constant out the front:
    dsa = α^2 / 2n
    dsb = β^2 / 2n
    dsab = (α + β)^2 / 4n
    ds = dsa + dsb - dsab
    s = ds
    i = 1
    dsold = ds # to fix α = -β bug.
    while abs(ds / s) + dsold > eps(Float64) / 10
        dsold = abs(ds / s)
        i += 1
        tmp = -(i - 1) / (i + 1) / n
        dsa = tmp * dsa * α
        dsb = tmp * dsb * β
        dsab = tmp * dsab * (α + β) / 2
        ds = dsa + dsb - dsab
        s = s + ds
    end
    # g is a vector of coefficients in
    # ``\Gamma(z) = \frac{z^{z-1/2}}{\exp(z)}\sqrt{2\pi} \left(\sum_{i} B_i z^{-i}\right)``, where B_{i-1} = g[i].
    # https://math.stackexchange.com/questions/1714423/what-is-the-pattern-of-the-stirling-series
    g = [
        1, 1 / 12, 1 / 288, -139 / 51840, -571 / 2488320, 163879 / 209018880,
        5246819 / 75246796800, -534703531 / 902961561600,
        -4483131259 / 86684309913600, 432261921612371 / 514904800886784000,
    ]
    f(z) = dot(g, [1;cumprod(ones(9) ./ z)])

    # Float constant C, C2
    C = 2 * (exp(s) * sqrt(2π * (n + α) * (n + β) / (2n + α + β)) / (2n + α + β + 1)) *
        (f(n + α) * f(n + β) / f(2n + α + β)) / π
    C2 = C * (α + β + 2n) * (α + β + 1 + 2n) / (4 * (α + n) * (β + n))

    vals = C * S

    # Use relation for derivative:
    ders = (n * ((α - β) .- (2n + α + β) * cos.(t)) .* vals .+ (2 * (n + α) * (n + β) * C2) .* S2) ./ (2n + α + β) ./ sin.(t)
    denom = 1 ./ (sin.(abs.(t) / 2) .^ (α + 0.5) .* cos.(t / 2) .^ (β + 0.5))
    vals .*= denom
    ders .*= denom

    return vals, ders
end

function asy2(n::Integer, α::Float64, β::Float64, npts::Integer, transform = Val(:none), reflected = false)
    # Algorithm for computing nodes and weights near the boundary.

    # Seed the Jacobi angles from Bessel zeros (Olver 1974, NIST 18.16.8).
    phik = approx_besselroots(α, npts) / (n + 0.5 * (α + β + 1))
    t = phik .+ ((α^2 - 0.25) .* (1 .- phik .* cot.(phik)) ./ (2 * phik) .- 0.25 .* (α^2 - β^2) .* tan.(0.5 .* phik)) ./ (n + 0.5 * (α + β + 1))^2

    # Compute higher terms:
    higherterms = asy2_higherterms(α, β, t)

    # Newton iteration:
    for _ in 1:10
        vals, ders = feval_asy2(n, α, β, t, higherterms)  # Evaluate via asymptotic formula.
        dt = vals ./ ders  # Newton update.
        t += dt  # Next iterate.
        if maximum(abs, dt; init = 0.0) < sqrt(eps(Float64)) / 200
            break
        end
    end
    vals, ders = feval_asy2(n, α, β, t, higherterms)  # Evaluate via asymptotic formula.
    t += vals ./ ders

    # flip:
    t = reverse(t)
    ders = reverse(ders)

    # Revert to x-space:
    x = cos.(t)
    w = 1 ./ (ders .^ 2 .* jacobi_weight_divisor.(t, Ref(transform), reflected))

    validate_rule(x, w)
    return x, w
end

"""
Evaluate the boundary asymptotic formula at x = cos(t).
"""
function feval_asy2(n::Integer, α::Float64, β::Float64, t::AbstractVector, higherterms::HT) where {HT <: Tuple{<:Function, <:Function, <:Function, <:Function}}
    rho = n + 0.5 * (α + β + 1)
    rho2 = n + 0.5 * (α + β - 1)
    A = (0.25 - α^2)
    B = (0.25 - β^2)

    # Evaluate the Bessel functions:
    Ja = besselj.(α, rho * t)
    Jb = besselj.(α + 1, rho * t)
    Jbb = besselj.(α + 1, rho2 * t)
    Jab = bessel_taylor(-t, rho * t, α)

    # Evaluate functions for recursive definition of coefficients:
    gt = A * (cot.(t / 2) .- (2 ./ t)) .- B * tan.(t / 2)
    gtdx = A * (2 ./ t .^ 2 .- 0.5 * csc.(t / 2) .^ 2) .- 0.5 * B * sec.(t / 2) .^ 2
    tB0 = 0.25 * gt
    A1 = gtdx / 8 .- (1 + 2 * α) / 8 * gt ./ t .- gt .^ 2 / 32 .- (α * (A + 3 * B) / 24)
    # Higher terms:
    tB1, A2, tB2, A3 = higherterms
    tB1t = tB1(t)
    A2t = A2(t)

    # VALS:
    vals = Ja + Jb .* tB0 / rho + Ja .* A1 / rho^2 + Jb .* tB1t / rho^3 + Ja .* A2t / rho^4
    # DERS:
    vals2 = Jab + Jbb .* tB0 / rho2 + Jab .* A1 / rho2^2 + Jbb .* tB1t / rho2^3 + Jab .* A2t / rho2^4

    # Higher terms (not needed for n > 1000).
    tB2t = tB2(t)
    A3t = A3(t)
    vals .+= Jb .* tB2t / rho^5 + Ja .* A3t / rho^6
    vals2 .+= Jbb .* tB2t / rho2^5 + Jab .* A3t / rho2^6

    # Constant out the front (Computed accurately!)
    ds = 0.5 * (α^2) / n
    s = ds
    jj = 1
    while abs(ds / s) > eps(Float64) / 10
        jj = jj + 1
        ds = -(jj - 1) / (jj + 1) / n * (ds * α)
        s = s + ds
    end
    g = [
        1, 1 / 12, 1 / 288, -139 / 51840, -571 / 2488320, 163879 / 209018880,
        5246819 / 75246796800, -534703531 / 902961561600,
        -4483131259 / 86684309913600, 432261921612371 / 514904800886784000,
    ]
    f(z) = dot(g, [1;cumprod(ones(9) ./ z)])
    C = (exp(s) * sqrt((n + α) / n) * (n / rho)^α) * (f(n + α) / f(n)) / sqrt(2)

    # Scaling:
    valstmp = C * vals
    denom = sin.(t / 2) .^ (α + 0.5) .* cos.(t / 2) .^ (β + 0.5)
    vals = sqrt.(t) .* valstmp ./ denom

    # Relation for derivative:
    C2 = C * n / (n + α) * (rho / rho2)^α
    ders = (n * (α - β .- (2n + α + β) * cos.(t)) .* valstmp .+ 2 * (n + α) * (n + β) * C2 * vals2) / (2n + α + β)
    ders .*= sqrt.(t) ./ (denom .* sin.(t))

    return vals, ders
end

function asy2_higherterms(α::Float64, β::Float64, theta::AbstractVector)
    # Higher-order terms for boundary asymptotic series.
    # Compute the higher order terms in asy2 boundary formula. See [2].

    # These constants are more useful than α and β:
    A = (0.25 - α^2)
    B = (0.25 - β^2)

    # For now, just work on half of the domain:
    c = max(maximum(theta), 0.5)

    # Use 10 Chebyshev modes in order to use precomputed integration and
    # differentiation matrices
    N = 10

    # Scaled 2nd-kind Chebyshev points and barycentric weights:
    t = 0.5 * c * (sin.(pi * (-(N - 1):2:(N - 1)) / (2 * (N - 1))) .+ 1)
    v = [0.5; ones(N - 1)]
    v[2:2:end] .= -1
    v[end] *= 0.5

    # The g's:
    g = A * (cot.(t / 2) - 2 ./ t) - B * tan.(t / 2)
    gp = A * (2 ./ t .^ 2 - 0.5 * csc.(t / 2) .^ 2) - 0.5 * (0.25 - β^2) * sec.(t / 2) .^ 2
    gpp = A * (-4 ./ t .^ 3 + 0.25 * sin.(t) .* csc.(t / 2) .^ 4) - 4 * B * sin.(t / 2) .^ 4 .* csc.(t) .^ 3
    g[1] = 0
    gp[1] = -A / 6 - 0.5 * B
    gpp[1] = 0

    # B0:
    B0 = 0.25 * g ./ t
    B0p = 0.25 * (gp ./ t - g ./ t .^ 2)
    B0[1] = 0.25 * (-A / 6 - 0.5 * B)
    B0p[1] = 0

    # A1:
    A1 = 0.125 * gp .- (1 + 2 * α) / 2 * B0 .- g .^ 2 / 32 .- (α * (A + 3 * B) / 24)
    A1p = 0.125 * gpp .- (1 + 2 * α) / 2 * B0p .- gp .* g / 16
    A1p_t = A1p ./ t
    A1p_t[1] = -A / 720 - A^2 / 576 - A * B / 96 - B^2 / 64 - B / 48 + α * (A / 720 + B / 48)

    # Make f accurately: (Taylor series approx for small t)
    f = -A * (
        1 / 12 .+ t .^ 2 / 240 + t .^ 4 / 6048 + t .^ 6 / 172800 + t .^ 8 / 5322240 +
            691 * t .^ 10 / 118879488000 + t .^ 12 / 5748019200 +
            3617 * t .^ 14 / 711374856192000 + 43867 * t .^ 16 / 300534953951232000
    )
    idx = t .> 0.5
    f[idx] = A .* (1 ./ t[idx] .^ 2 - 1 ./ (2 * sin.(t[idx] / 2)) .^ 2)
    f .-= B ./ (2 * cos.(t / 2)) .^ 2

    # Integrals for B1: (Note that N isn't large, so we don't need to be fancy).
    C = (0.5 * c) * CUMSUMMAT_10
    D = (2 / c) * DIFFMAT_10
    # B1:
    tB1 = -0.5 * A1p - (0.5 + α) * (C * A1p_t) + 0.5 * (C * (f .* A1))
    tB1[1] = 0
    B1 = tB1 ./ t
    B1[1] = A / 720 + A^2 / 576 + A * B / 96 + B^2 / 64 + B / 48 +
        α * (A^2 / 576 + B^2 / 64 + A * B / 96) - α^2 * (A / 720 + B / 48)

    # A2:
    A2 = 0.5 * (D * tB1) - (0.5 + α) * B1 - 0.5 * (C * (f .* tB1))
    A2 .-= A2[1]

    # A2p:
    A2p = D * A2
    A2p .-= A2p[1]
    A2p_t = A2p ./ t
    # Extrapolate point at t = 0:
    w = pi / 2 .- t[2:end]
    w[2:2:end] = -w[2:2:end]
    w[end] = 0.5 * w[end]
    A2p_t[1] = sum(w .* A2p_t[2:end]) / sum(w)

    # B2:
    tB2 = -0.5 * A2p - (0.5 + α) * (C * A2p_t) + 0.5 * C * (f .* A2)
    B2 = tB2 ./ t
    # Extrapolate point at t = 0:
    B2[1] = sum(w .* B2[2:end]) / sum(w)

    # A3:
    A3 = 0.5 * (D * tB2) - (0.5 + α) * B2 - 0.5 * (C * (f .* tB2))
    A3 .-= A3[1]

    tB1f(theta) = bary(theta, tB1, t, v)
    A2f(theta) = bary(theta, A2, t, v)
    tB2f(theta) = bary(theta, tB2, t, v)
    A3f(theta) = bary(theta, A3, t, v)

    return tB1f, A2f, tB2f, A3f
end

function jacobimoment(α, β)
    s = α + β
    T = float(typeof(s))
    a, b = α + 1, β + 1
    # The duplication formula avoids cancellation between terms of size α
    # when the mass of (1 - x^2)^α is only of order 1/sqrt(α).
    α == β && return beta(one(T) / 2, a)
    # Same as 2^(α+β+1) * beta(α+1,β+1)
    return exp((s + 1) * log(convert(T, 2)) + loggamma(a) + loggamma(b) - loggamma(a + b))
end

function bary(x::AbstractVector, fvals::AbstractVector, xk::AbstractVector, vk::AbstractVector)
    # simple barycentric interpolation routine adapted from chebfun/bary.m
    # p(x) = sum(v_k f_k/(x-x_k)) / sum(v_k/(x-x_k)); use f_k at x = x_k.

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
        if index !== nothing
            fx[k] = fvals[index]       # Set entries to the correct value
        end
    end

    return fx
end

function bessel_taylor(t::AbstractVector, z::AbstractVector, α::Float64)
    # Accurate evaluation of Bessel function J_A for asy2. (See [2].)
    # evaluates J_A(Z+T) by a Taylor series expansion about Z.

    npts = length(t)
    kmax = min(ceil(Int64, abs(log(eps(Float64)) / log(maximum(abs, t; init = 0.0)))), 30)
    H = t' .^ (0:kmax)
    # Compute coeffs in Taylor expansions about z (See NIST 10.6.7)
    nu = ones(Int64, length(z)) * (-kmax:kmax)'
    JK = z * ones(2kmax + 1)'
    Bjk = besselj.(α .+ nu, JK)
    nck = nck_mat(floor(Int64, 1.25 * kmax)) # nchoosek
    AA = [Bjk[:, kmax + 1] zeros(npts, kmax)]
    fact = 1
    for k in 1:kmax
        sgn = 1
        for l in 0:k
            AA[:, k + 1] = AA[:, k + 1] + sgn * nck[k, l + 1] * Bjk[:, kmax + 2 * l - k + 1]
            sgn = -sgn
        end
        fact = k * fact
        AA[:, k + 1] ./= 2^k * fact
    end
    # Evaluate Taylor series:
    Ja = zeros(npts)
    for k in 1:npts
        Ja[k] = dot(AA[k, :], H[:, k])
    end

    return Ja
end

function nck_mat(n::Integer)
    # almost triangular matrix storing n choose k
    M = zeros(Int64, n - 1, n)
    M[:, 1] .= 1
    M[1, 2] = 1
    for i in 2:(n - 1)
        for j in 2:(i + 1)
            M[i, j] = M[i - 1, j - 1] + M[i - 1, j]
        end
    end
    return M
end

# Fourth-order Jacobi root iteration and local Taylor propagation, independently
# implemented from Gil, Segura and Temme, arXiv:2008.08641, sections 2–4, and
# arXiv:2509.16716, section 2.3. No code from the accompanying software is used.
# U=(1-x)^((α+1)/2)(1+x)^((β+1)/2)P_n satisfies
# (1-x²)² U'' + (Ω+1)U=0. The iteration uses z=atanh(x), for which
# Ω=Ωmax-c(x-xe)². This vertex form avoids cancellation for large α=β.

function jacobi_gst_regime(n, α, β)
    return (50 <= n <= 100) || (2 <= n < 50 && min(α, β) < -0.9 && max(α, β) < 5) ||
        (max(α, β) >= 5 && (n >= 20 || α == β)) ||
        (n <= 100 && (min(α, β) + 1) / (n * (n + α + β + 1)) < eps(typeof(α)))
end

function jacobi_gst_start(n, α::T, β::T, x::T) where {T}
    # Form the positive shifted parameters before adding them. When both
    # parameters approach -1, α+β+2 and k+α-1 can lose their small values.
    xx = x
    a1, b1 = α + 1, β + 1
    shifted_sum = a1 + b1
    difference, parameter_sum = α - β, α + β
    p, dp = one(T), zero(T)
    q, dq = (difference + shifted_sum * xx) / 2, shifted_sum / 2
    changes = q < 0 ? 1 : 0
    lastsign = iszero(q) ? one(T) : sign(q)
    for k in 2:n
        offset = k - 2
        s, sm2 = 2(k - 1) + shifted_sum, 2offset + shifted_sum
        a = (s - 1) * (s * sm2 * xx + difference * parameter_sum)
        da = (s - 1) * s * sm2
        b = 2 * (offset + a1) * (offset + b1) * s
        c = 2k * (offset + shifted_sum) * sm2
        p, q, dp, dq = q, (a * q - b * p) / c,
            dq, (da * q + a * dq - b * dp) / c
        if !iszero(q)
            sign(q) != lastsign && (changes += 1)
            lastsign = sign(q)
        end
        scale = max(abs(p), abs(q), abs(dp), abs(dq))
        if scale > (T === Float64 ? 1.0e100 : sqrt(floatmax(T)))
            p /= scale
            q /= scale
            dp /= scale
            dq /= scale
        end
    end
    d = dq + (b1 / (2 * (1 + xx)) - a1 / (2 * (1 - xx))) * q
    scale = max(abs(q), abs(d))
    return q / scale, d / scale, changes
end

# A single local coordinate q represents x=origin+q. In the center origin=0;
# near an endpoint origin=±1, so the small endpoint distance is never recovered
# from a rounded output node. These are coordinate charts, not residual parts.
@inline function jacobi_gst_geometry(q, origin)
    x = origin + q
    left = origin == 1 ? -q : (1 - origin) - q
    right = origin == -1 ? q : (1 + origin) + q
    return x, left, right
end

function jacobi_gst_omega(q, origin, c, xe, ωmax, α, β, shifted = false)
    x, left, right = jacobi_gst_geometry(q, origin)
    main = c * left * right
    offset = (origin - xe) + q
    if main < ωmax
        if shifted
            # Form Ω+1 directly when Ω is close to -1 near an endpoint.
            return main + ((1 - α) * (1 + α) * right + (1 - β) * (1 + β) * left) / 2
        end
        return main - (α * α * right + β * β * left) / 2
    end
    return muladd(-c, offset * offset, ωmax + shifted)
end

# Taylor coefficients include powers of the actual representable step.
function jacobi_gst_taylor(q, origin, h, u, d, c, xe, ωmax, α, β)
    h == 0 && return u, d
    T = typeof(q)
    x, left, right = jacobi_gst_geometry(q, origin)
    s = left * right
    h2 = h * h
    q0 = s * s
    q1 = -4x * s * h
    q2 = (6x * x - 2) * h2
    q3 = 4x * h2 * h
    q4 = h2 * h2
    r0 = jacobi_gst_omega(q, origin, c, xe, ωmax, α, β, true) * h2
    r1 = -2c * ((origin - xe) + q) * h2 * h
    r2 = -c * h2 * h2
    am2, am1, a0, a1 = zero(u), zero(u), u, h * d
    value, deriv = u + a1, a1
    quiet = 0
    tolerance = eps(T) / 8
    for j in 0:max(100, 2precision(T))
        a2 = -muladd(
            q1 * j * (j + 1), a1,
            muladd(
                q2 * j * (j - 1) + r0, a0,
                muladd(
                    q3 * (j - 1) * (j - 2) + r1, am1,
                    (q4 * (j - 2) * (j - 3) + r2) * am2
                )
            )
        ) / (q0 * (j + 1) * (j + 2))
        value += a2
        deriv += (j + 2) * a2
        value_small = abs(a2) < tolerance * max(abs(value), abs(u), abs(h * d))
        deriv_small = abs((j + 2) * a2) < tolerance * max(abs(deriv), abs(h * d))
        if (value_small && deriv_small) || iszero(a2)
            quiet += 1
            quiet == 3 && return value, deriv / h
        else
            quiet = 0
        end
        am2, am1, a0, a1 = am1, a0, a1, a2
    end
    error("Jacobi Taylor series did not converge in $T; the type may have insufficient precision or exponent range.")
end

function jacobi_gst_propagate(q, origin, h, u, d, c, xe, ωmax, α, β)
    target = q + h
    _, left, right = jacobi_gst_geometry(target, origin)
    isfinite(target) && min(left, right) > 0 || error("Jacobi Taylor target is outside (-1, 1) in $(typeof(q)); the type may have insufficient precision or exponent range.")
    _, left, right = jacobi_gst_geometry(q, origin)
    distance = min(left, right)
    subdivisions = 0
    while abs(target - q) > distance / 2
        subdivisions += 1
        distance > 0 && subdivisions <= 100 || error("Jacobi Taylor propagation did not converge")
        next = q + copysign(distance / 2, h)
        next != q || error("Jacobi Taylor propagation made no progress")
        u, d = jacobi_gst_taylor(q, origin, next - q, u, d, c, xe, ωmax, α, β)
        q = next
        _, left, right = jacobi_gst_geometry(q, origin)
        distance = min(left, right)
    end
    u, d = jacobi_gst_taylor(q, origin, target - q, u, d, c, xe, ωmax, α, β)
    return target, u, d
end

function jacobi_gst_sweep!(nodes, logs, divisors, n, α, β, x, u, d, c, xe, ωmax, transform, reflected)
    T = typeof(x)
    origin = abs(x) <= T(0.5) ? zero(T) : copysign(one(T), x)
    q = x - origin
    count = length(nodes)
    for k in 1:count
        if k == count && α < T(-0.9)
            _, left, _ = jacobi_gst_geometry(q, origin)
            left, logs[k] = jacobi_gst_endpoint(n, α, β, left, u, d)
            nodes[k] = 1 - left
            divisors[k] = jacobi_gst_divisor(left, 2 - left, transform, reflected)
            break
        end
        for iter in 1:20
            x, left, right = jacobi_gst_geometry(q, origin)
            ω = jacobi_gst_omega(q, origin, c, xe, ωmax, α, β)
            rootω = sqrt(abs(ω))
            s = left * right
            ratio = rootω * u / muladd(s, d, x * u)
            if ω > 0
                angle = atan(ratio)
                if iter == 1 && (k > 1 || angle >= 0)
                    angle -= T(π)
                end
                step = -angle / rootω
            elseif iszero(ω)
                step = -u / muladd(s, d, x * u)
            else
                abs(ratio) < 1 || error("Jacobi iteration left the convergence interval in $T; the type may have insufficient precision or exponent range.")
                step = -atanh(ratio) / rootω
            end
            t = tanh(step)
            denominator = x < 0 && t > 0 ? (1 - t) + right * t : muladd(x, t, one(T))
            h = (s * t) / denominator
            previous = q
            q, u, d = jacobi_gst_propagate(q, origin, h, u, d, c, xe, ωmax, α, β)
            x, left, right = jacobi_gst_geometry(q, origin)
            if q == previous || abs(h) <= 2eps(T) * min(max(abs(x), 1 / sqrt(c)), left)
                break
            end
            neworigin = abs(x) <= T(0.5) ? zero(T) : copysign(one(T), x)
            if neworigin != origin
                q = (origin - neworigin) + q
                origin = neworigin
            end
            iter == 20 && error("Jacobi fourth-order iteration did not converge in $T.")
        end
        x, left, right = jacobi_gst_geometry(q, origin)
        nodes[k] = x
        if iszero(origin)
            logleft, logright = log1p(-q), log1p(q)
            logproduct = log1p(-q * q)
        else
            logleft, logright = log(left), log(right)
            logproduct = logleft + logright
        end
        a = min(α, β)
        logs[k] = a * logproduct + (α - a) * logleft + (β - a) * logright - 2log(abs(d))
        divisors[k] = jacobi_gst_divisor(left, right, transform, reflected)
    end
    return nothing
end

@inline jacobi_gst_divisor(left, right, ::Val{:none}, reflected) = one(left)
@inline jacobi_gst_divisor(left, right, ::Val{:radau}, reflected) = reflected ? left : right
@inline jacobi_gst_divisor(left, right, ::Val{:lobatto}, reflected) = left * right

function jacobi_gst(n::Integer, α::T, β::T, transform = Val(:none)) where {T <: AbstractFloat}
    s = α + β
    c = (T(n) + s / 2) * (T(n) + s / 2 + 1)
    xe = (β - α) * s / (4c)
    center = xe
    v = 4T(n) * (T(n) + s + 1) + 2s
    ωmax = (v / (4c)) * (v / 4 + α * β)
    if !(isfinite(c) && isfinite(center) && isfinite(ωmax))
        error("The Jacobi iteration coefficients are not representable in $T.")
    end
    x, logs, divisors = Vector{T}(undef, n), Vector{T}(undef, n), Vector{T}(undef, n)
    if α == β
        m = n ÷ 2
        u, d = isodd(n) ? (zero(T), one(T)) : (one(T), zero(T))
        # For Radau, retain the small reflected endpoint distance. The other
        # divisor is 2 minus this distance, so symmetry needs only one sweep.
        jacobi_gst_sweep!(
            @view(x[(n - m + 1):n]), @view(logs[(n - m + 1):n]), @view(divisors[(n - m + 1):n]), n,
            α, β, zero(T), u, d, c, zero(T), ωmax, transform, transform isa Val{:radau}
        )
        @inbounds for k in 1:m
            x[k], logs[k] = -x[n + 1 - k], logs[n + 1 - k]
            divisors[k] = divisors[n + 1 - k]
            if transform isa Val{:radau}
                divisors[n + 1 - k] = 2 - divisors[k]
            end
        end
        if isodd(n)
            x[m + 1], logs[m + 1] = zero(T), zero(T)
            divisors[m + 1] = one(T)
        end
    else
        u, d, right = jacobi_gst_start(n, α, β, center)
        left = n - right
        first_right = left + 1
        if iszero(u)
            x[left] = center
            logs[left] = α * log1p(-center) + β * log1p(center) - 2log(abs(d))
            divisors[left] = jacobi_gst_divisor(1 - center, 1 + center, transform, false)
            left -= 1
        end
        jacobi_gst_sweep!(
            @view(x[first_right:n]), @view(logs[first_right:n]), @view(divisors[first_right:n]), n,
            α, β, center, u, d, c, xe, ωmax, transform, false
        )
        jacobi_gst_sweep!(
            @view(x[left:-1:1]), @view(logs[left:-1:1]), @view(divisors[left:-1:1]), n,
            β, α, -center, u, -d, c, -xe, ωmax, transform, true
        )
        @inbounds for k in 1:left
            x[k] = -x[k]
        end
    end
    w = jacobi_gst_weights(logs, divisors, α, β)
    validate_rule(x, w)
    return x, w
end

# Normalize the Jacobi weights before applying Radau or Lobatto divisors. The logs
# retain weights that would underflow before multiplication by the total mass.
function jacobi_gst_weights(logs, divisors, α::T, β::T) where {T <: AbstractFloat}
    peak = maximum(logs)
    w = exp.(logs .- peak)
    mass = jacobi_gst_mass(α, β)
    if T !== Float64 && !(isfinite(mass) && mass > zero(T))
        error("The quadrature normalization is not representable in $T.")
    end
    scale = mass / sum(w)
    logscale = log(scale)
    for k in eachindex(w)
        # A normalized weight can underflow even when multiplication by a large
        # total mass (or the Radau divisor) makes it representable.
        w[k] = w[k] < floatmin(T) ? exp(logs[k] - peak + logscale - log(divisors[k])) :
            (w[k] * scale) / divisors[k]
    end
    return w
end

# The root tending to 1 as α tends to -1 lies beyond the turning point in
# atanh(x). Refine it in δ=(1-x)/2 with the terminating hypergeometric series.
# Match its derivative to the preceding Taylor value; only the final root
# requires this evaluation, and δ stays accurate when the output rounds to 1.
function jacobi_gst_endpoint(n, a, b, left, u, d)
    δprevious = left / 2
    f, df = jacobi_gst_hypergeometric(n, a, b, δprevious)
    logfactor = (a + 1) / 2 * log(2δprevious) + (b + 1) / 2 * log(2 * (1 - δprevious))
    if abs(d) > abs(u)
        derivative = ((b + 1) / (4 * (1 - δprevious)) - (a + 1) / (4δprevious)) * f - df / 2
        logscale = log(abs(d)) - logfactor - log(abs(derivative))
    else
        logscale = log(abs(u)) - logfactor - log(abs(f))
    end
    δ = (a + 1) / (n * (n + a + b + 1))
    converged = false
    for _ in 1:20
        f, df = jacobi_gst_hypergeometric(n, a, b, δ)
        step = f / df
        δ -= step
        if abs(step) <= 8eps(typeof(δ)) * δ
            converged = true
            break
        end
    end
    converged || error("Jacobi endpoint series iteration did not converge")
    _, df = jacobi_gst_hypergeometric(n, a, b, δ)
    return 2δ, -log(4δ * (1 - δ)) - 2log(abs(df / 2)) - 2logscale
end

function jacobi_gst_hypergeometric(n, a, b, δ)
    term, value, derivative = one(δ), one(δ), zero(δ)
    for k in 1:n
        term *= ((k - 1 - n) * (n + a + b + k) / ((a + k) * k)) * δ
        value += term
        derivative += k * term / δ
        abs(term) <= eps(typeof(δ)) * abs(value) && k > 5 && break
    end
    return value, derivative
end

# A precision-aware Stirling correction avoids subtracting large log-gamma
# values in the mass of fixed-precision rules with large Jacobi parameters.
function jacobi_gst_stirling(x::T) where {T <: AbstractFloat}
    coefficients = (
        1 // 12, -1 // 360, 1 // 1260, -1 // 1680,
        1 // 1188, -691 // 360360, 1 // 156, -3617 // 122400,
    )
    if precision(T) > 64
        coefficients = (
            coefficients..., 43867 // 244188, -174611 // 125400,
            77683 // 5796, -236364091 // 1506960, 657931 // 300,
            -3392780147 // 93960, 1723168255201 // 2492028,
            -7709321041217 // 505920,
        )
    end
    return inv(x) * evalpoly(inv(x)^2, map(T, coefficients))
end

function jacobi_gst_mass_generic(α::T, β::T) where {T <: AbstractFloat}
    bits = precision(T)
    # For the retained terms, the first omitted Stirling term is below the
    # working epsilon after shifting both gamma arguments to this threshold.
    shift = bits <= 64 ? cld(bits + 8, 17) : cld(bits + 40, 33)
    shift > 10 && return jacobimoment(α, β)
    a, b = α + 1, β + 1
    correction = zero(T)
    while a < max(8, 1 << shift)
        correction += log((a + b) / (2a))
        a += 1
    end
    while b < max(8, 1 << shift)
        correction += log((a + b) / (2b))
        b += 1
    end
    s = a + b
    r = (a - b) / s
    if abs(r) <= one(T) / 2
        logmass = (s - 1) / 2 * log1p(-r^2) + (a - b) * atanh(r) + log(T(2) * T(π) / s) / 2
    else
        hi, lo = max(a, b), min(a, b)
        q = lo / hi
        logq = log1p(q)
        log2 = log(T(2))
        logmass = hi * (log2 - logq) + lo * (log2 + log(q) - logq) +
            (log(T(π)) - log2 + logq - log(lo)) / 2
    end
    return exp(
        logmass + jacobi_gst_stirling(a) + jacobi_gst_stirling(b) -
            jacobi_gst_stirling(s) + correction
    )
end

# Shared mass policy for generic Jacobi arithmetic, including Radau's endpoint
# weight. Keep the cancellation-resistant Stirling formula for large parameters.
function jacobi_mass(α::T, β::T) where {T <: AbstractFloat}
    return T === BigFloat || max(α, β) < 5 ? jacobimoment(α, β) : jacobi_gst_mass_generic(α, β)
end

function jacobi_gst_mass(α::T, β::T) where {T <: AbstractFloat}
    if (iszero(α) && (iszero(β) || β == one(T))) || (α == one(T) && iszero(β))
        return T(2)
    elseif α == β == one(T)
        return T(4) / 3
    else
        return jacobi_mass(α, β)
    end
end

function jacobi_gst_mass(α::Float64, β::Float64)
    α == β && return beta(0.5, α + 1)
    a, b = α + 1, β + 1
    correction = 0.0
    while a < 8
        correction += log((a + b) / (2a))
        a += 1
    end
    while b < 8
        correction += log((a + b) / (2b))
        b += 1
    end
    s = a + b
    r = (a - b) / s
    logmass = (s - 1) / 2 * log1p(-r^2) + (a - b) * atanh(r) +
        log(2π / s) / 2 + jacobi_gst_stirling(a) + jacobi_gst_stirling(b) - jacobi_gst_stirling(s)
    return exp(logmass + correction)
end
