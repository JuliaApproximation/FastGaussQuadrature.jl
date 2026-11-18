gausslaguerre(::Type{T}, n::Integer) where {T} = gausslaguerre(n, T(0))
@doc raw"""
    gausslaguerre([T=Float64,] n::Integer) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Laguerre quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature) with type `T`.
The order `n` must be nonnegative. For a nonzero Laguerre parameter or the
`reduced` keyword, use `gausslaguerre(n, α; reduced=false)`.

```math
\int_{0}^{+\infty} f(x) e^{-x} dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gausslaguerre(3);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ 24
true
```
"""
gausslaguerre(n::Integer) = gausslaguerre(Float64, n)

@doc raw"""
    gausslaguerre(n::Integer, α::Real; reduced=false) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of generalized [Gauss-Laguerre quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature).
The parameters must satisfy `n ≥ 0` and `α > -1`. A floating-point `α`
determines the element type of both output vectors; for example,
`gausslaguerre(32, 0.5f0)` returns `Float32` nodes and weights.

```math
\int_{0}^{+\infty} f(x) x^\alpha e^{-x} dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gausslaguerre(3, 1.0);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ 120
true
```

When the working type is `Float64`, `n ≥ 128`, and `α < 5`, `reduced=true`
returns a prefix of the rule, stopping when weights fall below
`10floatmin(Float64)`. This can omit representable normal and subnormal weights.
Other cases return the full rule. The cutoff is not an error bound for an
arbitrary integrand.

The default algorithm supports compatible floating-point types. The following
alternative routines can also be invoked manually:

- `FastGaussQuadrature.gausslaguerre_rec`: Newton iteration using the recurrence
  for orders at most 127.
- `FastGaussQuadrature.gausslaguerre_asy`: asymptotic expansions tuned for `Float64`.

These alternatives have different precision and performance limitations from
the default dispatch; see their docstrings.
"""
gausslaguerre(n::Integer, α::Real; reduced = false) = gausslaguerre_rule(n, α; reduced)
function gausslaguerre(n::Integer, α::T; reduced = false)::Tuple{Vector{T}, Vector{T}} where {T <: AbstractFloat}
    return gausslaguerre_rule(n, α; reduced)
end

function gausslaguerre_rule(n::Integer, α::Real; reduced = false)
    if α ≤ -1
        throw(DomainError(α, "The Laguerre parameter α ≤ -1 corresponds to a nonintegrable weight function"))
    end
    if n < 0
        throw(DomainError(n, "gausslaguerre($n,$α) not defined: n must be positive."))
    end

    # Guess the numerical type from the supplied type of α
    α = float(α)
    T = typeof(α)
    if n == 0
        x, w = T[], T[]
        validate_rule(x, w)
        return x, w
    elseif precision(T) < precision(Float32)
        x, w = gausslaguerre(n, Float32(α); reduced)
        x, w = T.(x), T.(w)
        validate_rule(x, w)
        return x, w
    end
    if n == 1
        x, w = [1 + α], [gamma(1 + α)]
        validate_rule(x, w)
        return x, w
    elseif n == 2
        if iszero(α)
            s = sqrt(T(2))
            x, w = [s / (s + 1), 2 + s], [(1 + inv(s)) / 2, inv(2s * (s + 1))]
            validate_rule(x, w)
            return x, w
        end
        s = sqrt(α + 2)
        μ = gamma(α + 1)
        # Rationalize the small root and weight as α approaches -1.
        # The original formulas subtract nearly equal numbers in both places.
        x = [(α + 1) * (s / (s + 1)), α + 2 + s]
        w = [μ * ((1 + inv(s)) / 2), μ * ((α + 1) / (2s * (s + 1)))]
        validate_rule(x, w)
        x[1] > 0 || error("Laguerre nodes are not positive in $T; try a type with higher precision or wider exponent range.")
        return x, w
    elseif T != Float64 || n < 128 || α ≥ 5
        return gausslaguerre_gst(n, float(α))
    end
    # Preserve the established Float64 reduced-rule region and accurate
    # half-parameter endpoint expansions; all other nodes use native GST.
    x, w = gausslaguerre_gst(n, float(α); reduced)
    if abs(α) == 0.5
        # At the first few nodes the half-integer Bessel expansions retain
        # close to full relative precision in both nodes and weights.
        d = inv(4n + 2α + 2)
        for k in 1:(n ≥ 1000 ? 12 : 2)
            x[k], w[k], _ = gausslaguerre_asy_bessel(n, α, π * (α > 0 ? k : k - 0.5), d, -1, false)
        end
    end
    if reduced
        cutoff = findfirst(wi -> abs(wi) < 10floatmin(T), w)
        if cutoff !== nothing
            resize!(x, cutoff - 1)
            resize!(w, cutoff - 1)
        end
    end
    validate_rule(x, w)
    return x, w
end


"""
    FastGaussQuadrature.gausslaguerre_asy(n, α; reduced=false, T=max(1, ceil(Int, 50/log(n))), recompute=false)

Compute Gauss–Laguerre nodes and weights using explicit asymptotic expansions.
This routine is tuned for large-order `Float64` rules with fixed `α > -1`;
`float(α)` must be `Float64`. Its coefficients and initial approximations do
not support other working formats; use `gausslaguerre` for those formats.
Large `α^2/n` can reduce accuracy. Prefer `gausslaguerre` for automatic selection.

Optional keywords:

- `reduced`: stop at the first weight below `10floatmin(typeof(float(α)))`,
  omitting that node and the remaining tail, including representable weights.
- `T`: number of expansion terms, not a floating-point type. Set `T=-1` to
  choose the truncation adaptively from the term sizes.
- `recompute`: refine a node and weight by recurrence and Newton iteration when
  a heuristic error estimate exceeds a tolerance. This option is limited to
  `n ≤ 1000`, since refining a full rule this way can take quadratic work.
"""
function gausslaguerre_asy(
        n::Integer, α;
        reduced = false,
        T = max(1, ceil(Int, 50 / log(n))),  # Heuristic for number of terms
        recompute = false
    )

    recompute && n > 1000 && throw(DomainError(n, "Recurrence recomputation is limited to n ≤ 1000; use recompute=false or gausslaguerre for larger rules."))

    if α^2 / n > 1
        @warn "A large value of α may lead to inaccurate results."
    end

    ELT = typeof(float(α))
    ELT === Float64 || throw(ArgumentError("The asymptotic helper supports Float64 parameters only; use gausslaguerre for other floating-point types."))

    x = zeros(ELT, reduced ? 0 : n)
    w = zeros(ELT, reduced ? 0 : n)

    # The expansions are given in powers of 1/(4n+2α+2)
    d = one(ELT) / (4n + 2α + 2)

    # Heuristical indices for Bessel and Airy regions
    k_bessel = max(ceil(Int, sqrt(n)), 7)

    # The Bessel region
    # First compute the roots of the Bessel function of order α
    # The half-integer Bessel zeros are exact multiples of π.
    jak_vector = abs(α) == 0.5 ? ELT[] : approx_besselroots(α, k_bessel)

    bessel_wins = true
    k = 0
    while bessel_wins && k < n
        k += 1
        if k == 1 && α + 1 <= 1 / 8
            xk, wk, _ = laguerre_endpoint(n, float(α))
            if reduced
                if abs(wk) < 10floatmin(ELT)
                    validate_rule(x, w)
                    return x, w
                end
                push!(x, xk); push!(w, wk)
            else
                x[k] = xk; w[k] = wk
            end
            continue
        end
        # We iterate until the estimated error of the bulk expansion is smaller
        # than the one of the Bessel expansion
        jak = if α == 0.5
            ELT(π) * k
        elseif α == -0.5
            ELT(π) * (k - ELT(1) / 2)
        else
            (k < k_bessel) ? jak_vector[k] : McMahon(α, k)
        end

        xk, wk, δ_bessel = gausslaguerre_asy_bessel(n, α, jak, d, T, k ≤ 12)
        xkb, wkb, δ_bulk = gausslaguerre_asy_bulk(n, α, k, d, T)
        if δ_bulk < δ_bessel && !(abs(α) == 0.5 && n ≥ 10000 && k ≤ 20)
            bessel_wins = false
            xk = xkb
            wk = wkb
        end
        if recompute
            if min(δ_bessel, δ_bulk) > 1.0e-13
                xk_rec, wk_rec = gl_rec_newton(xk, n, α)
                if abs(xk_rec - xk) < 100 * min(δ_bessel, δ_bulk)
                    xk = xk_rec
                    wk = wk_rec
                end
            end
        end
        if reduced
            if abs(wk) < 10floatmin(ELT)
                validate_rule(x, w)
                return x, w
            else
                push!(x, xk); push!(w, wk)
            end
        else
            x[k] = xk; w[k] = wk
        end
    end

    # The bulk region
    # - First we go from where we left of to our heuristic
    while k < floor(Int, 0.9 * n) - 1
        k += 1
        xk, wk, δ_bulk = gausslaguerre_asy_bulk(n, α, k, d, T)
        if recompute
            if δ_bulk > 1.0e-13
                xk_rec, wk_rec = gl_rec_newton(xk, n, α)
                if abs(xk_rec - xk) < 100δ_bulk
                    xk = xk_rec
                    wk = wk_rec
                end
            end
        end
        if reduced
            if abs(wk) < 10floatmin(ELT)
                validate_rule(x, w)
                return x, w
            else
                push!(x, xk); push!(w, wk)
            end
        else
            x[k] = xk; w[k] = wk
        end
    end

    # - Then we compare to Airy until it wins, and then we switch to just Airy
    bulk_wins = true
    while bulk_wins && k < n
        k += 1
        xk, wk, δ_bulk = gausslaguerre_asy_bulk(n, α, k, d, T)
        xka, wka, δ_airy = gausslaguerre_asy_airy(n, α, k, d, T)
        if δ_airy < δ_bulk
            bulk_wins = false
            xk = xka
            wk = wka
        end
        if recompute
            if min(δ_airy, δ_bulk) > 1.0e-13
                xk_rec, wk_rec = gl_rec_newton(xk, n, α)
                if abs(xk_rec - xk) < 100 * min(δ_airy, δ_bulk)
                    xk = xk_rec
                    wk = wk_rec
                end
            end
        end
        if reduced
            if abs(wk) < 10floatmin(ELT)
                validate_rule(x, w)
                return x, w
            else
                push!(x, xk); push!(w, wk)
            end
        else
            x[k] = xk; w[k] = wk
        end
    end

    # The Airy region
    while k < n
        k += 1
        xk, wk, δ_airy = gausslaguerre_asy_airy(n, α, k, d, T)
        if recompute
            if δ_airy > 1.0e-13
                xk_rec, wk_rec = gl_rec_newton(xk, n, α)
                if abs(xk_rec - xk) < 100δ_airy
                    xk = xk_rec
                    wk = wk_rec
                end
            end
        end
        if reduced
            if abs(wk) < 10floatmin(ELT)
                validate_rule(x, w)
                return x, w
            else
                push!(x, xk); push!(w, wk)
            end
        else
            x[k] = xk; w[k] = wk
        end
    end

    validate_rule(x, w)
    all(t -> 0 < t ≤ 4 * n + 2 * α + 2, x) ||
        error("Laguerre nodes lie outside the expected positive interval.")
    return x, w
end

## Expansion coefficients
# These are explicit formulas of the coefficients, up to a simple postprocessing
# that is common to all factors and not included here (see below).
#
# General expressions are given in terms of α, more specific expressions
# follow for the special case α = 0.

## The bulk

# Bulk expansions use (1-t)⁻¹, (1-t)/t, and d=(4n+2α+2)⁻¹.
# Note: there is always one division by an integer, placed such that it preserves the type of `d`
function gl_bulk(t, d, α)
    α² = α * α
    tinv = inv(1 - t)
    _t = (1 - t) / t
    _t² = _t * _t
    d² = d * d
    d⁴ = d² * d²

    c1 = evalpoly(α², (-4, 12))
    x3 = -evalpoly(tinv, (c1, -4, 5)) * d / 12
    c1, c2 = evalpoly(α², (224, -960, 480)), evalpoly(α², (7, -30, 15))
    x5 = d² * d * _t / 720 * evalpoly(tinv, (c1, -48 * c2, -16, -576, 2814, -3815, 1600))
    c1, c2, c3 = evalpoly(α², (-285696, 1354752, -967680, 193536)), evalpoly(α², (-31, 147, -105, 21)), evalpoly(α², (-1346, 6405, -4620, 945))
    c4, c5 = evalpoly(α², (43, -126, 63)), evalpoly(α², (-221, -630, 315))
    x7 = -d⁴ * d / 181440 * _t² * evalpoly(tinv, (c1, -23040 * c2, 384 * c3, -320 * c4, 80 * c5, -1727136, 16131880, -48469876, 66424575, -43122800, 10797500))
    c1, c2, c3 = evalpoly(α², (210677760, -1028505600, 812851200, -232243200, 24883200)), evalpoly(α², (127, -620, 490, -140, 15)), evalpoly(α², (1193053, -5826660, 4613070, -1324260, 143325))
    c4, c5, c6 = evalpoly(α², (555239, -2716980, 2163630, -631260, 70875)), evalpoly(α², (-641, 2960, -2155, 450)), evalpoly(α², (-1598, 17685, -13905, 3375))
    c7, c8, c9 = evalpoly(α², (-7823, -9042, 4521)), evalpoly(α², (15948182, -206850, 103425)), evalpoly(α², (64957561, -24000, 12000))
    x9 = (d⁴ * d²) * d / 10886400 * (_t² * _t) * evalpoly(tinv, (c1, -5806080 * c2, 768 * c3, -768 * c4, 16128 * c5, -1792 * c6, 3360 * c7, -192 * c8, 672 * c9, -212307298152, 518401904799, -714465642135, 566519158800, -241928673000, 43222750000))

    w3 = d² / 6 * (2 * t + 3) / (t - 1)^3
    c1 = evalpoly(α², (7, -30, 15))
    w5 = d⁴ * _t² / 720 * evalpoly(tinv, (0, 0, 16 * c1, 32, 1712, -12408, 27517, -24860, 8000))
    c1, c2, c3 = evalpoly(α², (-31, 147, -105, 21)), evalpoly(α², (-416, 1995, -1470, 315)), evalpoly(α², (43, -126, 63))
    c4, c5 = evalpoly(α², (-89, -378, 189)), evalpoly(α², (53752, -630, 315))
    w7 = -(d⁴ * d²) * (_t² * _t) / 90720 * evalpoly(tinv, (0, 0, 2304 * c1, -384 * c2, 480 * c3, -320 * c4, 80 * c5, -50986344, 201908326, -386872990, 393326325, -204917300, 43190000))
    return (x3, x5, x7, x9), (w3, w5, w7)
end

function gl_bulk(t, d)
    tinv = inv(1 - t)
    _t = (1 - t) / t
    _t² = _t * _t
    d² = d * d
    d⁴ = d² * d²

    x3 = -d / 12 * evalpoly(tinv, (-4, -4, 5))
    x5 = d² * d * _t / 720 * evalpoly(tinv, (224, -336, -16, -576, 2814, -3815, 1600))
    x7 = -d⁴ * d / 181440 * _t² * evalpoly(tinv, (-285696, 714240, -516864, -13760, -17680, -1727136, 16131880, -48469876, 66424575, -43122800, 10797500))
    x9 = (d⁴ * d²) * d / 10886400 * (_t² * _t) * evalpoly(tinv, (210677760, -737372160, 916264704, -426423552, -10338048, 2863616, -26285280, -3062050944, 43651480992, -212307298152, 518401904799, -714465642135, 566519158800, -241928673000, 43222750000))

    w3 = d² / 6 * (2 * t + 3) / (t - 1)^3
    w5 = d⁴ / 720 * _t² * evalpoly(tinv, (0, 0, 112, 32, 1712, -12408, 27517, -24860, 8000))
    w7 = -(d⁴ * d²) * (_t² * _t) / 90720 * evalpoly(tinv, (0, 0, -71424, 159744, 20640, 28480, 4300160, -50986344, 201908326, -386872990, 393326325, -204917300, 43190000))
    return (x3, x5, x7, x9), (w3, w5, w7)
end

## The hard edge (Bessel region)
function gl_bessel(jak, d, α)
    jak² = jak * jak
    d² = d * d
    d⁴ = d² * d²
    α² = α * α

    c1 = evalpoly(α², (-2, 2))
    x3 = d² * evalpoly(jak², (c1, 1)) / 3
    w3 = d² * evalpoly(jak², (c1, 2)) / 3
    c1, c2 = evalpoly(α², (94, -140, 46)), evalpoly(α², (-19, 11))
    x5 = d⁴ * evalpoly(jak², (c1, 3 * c2, 11)) / 45
    w5 = d⁴ * evalpoly(jak², (c1, 6 * c2, 33)) / 45
    c1, c2, c3 = evalpoly(α², (-12077, 19887, -9303, 1493)), evalpoly(α², (14051, -10750, 2459)), evalpoly(α², (-181, 73))
    x7 = d⁴ * d² * evalpoly(jak², (4 * c1, 2 * c2, 36 * c3, 657)) / 2835
    w7 = 4 * d⁴ * d² * evalpoly(jak², (c1, c2, 27 * c3, 657)) / 2835
    c1, c2, = evalpoly(α², (6029959, -10087180, 5095482, -1146220, 107959)), evalpoly(α², (-2201939, 1678761, -507801, 63299))
    c3, c4 = evalpoly(α², (1456807, -729422, 125671)), evalpoly(α², (-2879, 887))
    x9 = (d⁴ * d⁴) * evalpoly(jak², (2 * c1, 3 * c2, c3, 60 * c4, 10644)) / 42525
    w9 = (d⁴ * d⁴) * evalpoly(jak², (2 * c1, 6 * c2, 3 * c3, 240 * c4, 53220)) / 42525
    return (x3, x5, x7, x9), (w3, w5, w7, w9)
end

# And for α = 0:
function gl_bessel(jak, d)
    jak² = jak * jak
    d² = d * d
    d⁴ = d² * d²
    d⁸ = d⁴ * d⁴

    x3 = d² * evalpoly(jak², (-2, 1)) / 3
    x5 = d⁴ * evalpoly(jak², (94, -57, 11)) / 45
    x7 = d⁴ * d² * evalpoly(jak², (-48308, 28102, -6516, 657)) / 2835
    x9 = d⁸ * evalpoly(jak², (12059918, -6605817, 1456807, -172740, 10644)) / 42525
    x11 = d⁸ * d² * evalpoly(jak², (-11427291076, 6028914206, -1248722004, 138902061, -9908262, 410649)) / 1403325

    w3 = 2 * d² * evalpoly(jak², (-1, 1)) / 3
    w5 = d⁴ * evalpoly(jak², (94, -114, 33)) / 45
    w7 = 4 * d⁴ * d² * evalpoly(jak², (-12077, 14051, -4887, 657)) / 2835
    w9 = d⁸ * evalpoly(jak², (12059918, -13211634, 4370421, -690960, 53220)) / 42525
    w11 = 2 * d⁸ * d² * evalpoly(jak², (-5713645538, 6028914206, -1873083006, 277804122, -24770655, 1231947)) / 1403325
    return (x3, x5, x7, x9, x11), (w3, w5, w7, w9, w11)
end

## The soft edge (Airy region)
function gl_airy(ak, d, α)
    T = eltype(α)
    cbrt_d = cbrt(d)
    ak² = ak * ak

    # Airy corrections in powers of d^(1/3), starting from x≈d⁻¹+a_k(4/d)^(1/3).
    x1 = 1 / d + ak * cbrt(T(4)) / cbrt_d
    x3 = ak² * cbrt(T(16)) / 5 * cbrt_d + (T(11) / 35 - α^2 - T(12) / 175 * (ak² * ak)) * d + (T(16) / 1575 * ak + T(92) / 7875 * (ak² * ak²)) * cbrt(T(4)) * d * (cbrt_d * cbrt_d)
    x5 = -(T(15152) / 3031875 * (ak² * ak²) * ak + T(1088) / 121275 * ak²) * cbrt(T(2)) * d * d * cbrt_d
    return x1, x3, x5
end

# Sum the first Q elements of vals, and return the absolute value of the next
# value in the list (or of the last value in the list)
function sum_explicit(vals, Q)
    z = zero(eltype(vals[1]))
    for q in min(Q, length(vals)):-1:1
        z += vals[q]
    end
    if Q < length(vals)
        return z, abs(vals[Q + 1])
    else
        return z, abs(vals[end])
    end
end

function sum_adaptive(vals)
    z = vals[1]
    i = 1
    while (i < length(vals)) && (abs(vals[i + 1]) < abs(vals[i]))
        i += 1
        z += vals[i]
    end
    return z, abs(vals[min(i + 1, length(vals))])
end

function gl_bulk_solve_t(n, k, d, α = zero(d))
    T = typeof(d)
    # With t=sin(θ/2)^2 the phase equation is θ+sin(θ)=phase.
    # Near t=1 solve for the complementary angle, η-sin(η)=π-phase,
    # forming either small phase directly rather than subtracting from π.
    complement = k > n ÷ 2
    phase = T(π) * ((complement ? 4(n - k) + 3 : 4k + 2α - 1) * d)
    angle = complement ? cbrt(6phase) : phase / 2
    for _ in 1:20
        s, c = sincos(angle)
        step = (angle + (complement ? -one(T) : one(T)) * s - phase) / (1 + (complement ? -one(T) : one(T)) * c)
        next = angle - step
        # Stop at the rounding floor in t; a tighter test can oscillate
        # between adjacent angles until the iteration limit is reached.
        if next == angle || abs(step) <= (complement ? eps(T) / abs(s) : 2eps(T) * angle)
            angle = next
            break
        end
        angle = next
    end
    return complement ? cos(angle / 2)^2 : sin(angle / 2)^2
end

function gausslaguerre_asy_bulk(n, α, k, d, T)
    if α == 0
        return gausslaguerre_asy0_bulk(n, k, d, T)
    end

    t = gl_bulk_solve_t(n, k, d, α)
    xs, ws = gl_bulk(t, d, α)

    xk, xdelta = (T > 0) ? sum_explicit(xs, (T - 1) >> 1) : sum_adaptive(xs)
    wk, wdelta = (T > 0) ? sum_explicit(ws, (T - 1) >> 1) : sum_adaptive(ws)

    xk += t / d

    wfactor = xk^α * exp(-xk) * 2π * sqrt(t / (1 - t))
    wk = wfactor * (1 + wk)
    wdelta *= wfactor

    return xk, wk, max(xdelta, wdelta)
end


function gausslaguerre_asy0_bulk(n, k, d, T)
    t = gl_bulk_solve_t(n, k, d)
    xs, ws = gl_bulk(t, d)


    xk, xdelta = (T > 0) ? sum_explicit(xs, (T - 1) >> 1) : sum_adaptive(xs)
    wk, wdelta = (T > 0) ? sum_explicit(ws, (T - 1) >> 1) : sum_adaptive(ws)

    xk += t / d

    wfactor = exp(-xk) * 2π * sqrt(t / (1 - t))
    wk = wfactor * (1 + wk)
    wdelta *= wfactor

    return xk, wk, max(xdelta, wdelta)
end


function gausslaguerre_asy_bessel(n, α, jak, d, T, refine = false)
    if α == 0
        return gausslaguerre_asy0_bessel(n, jak, d, T)
    end

    # At zeros of J_{±1/2}, J_{α-1}² = 2/(π*j).
    derivative = abs(α) == 0.5 ? sqrt(2 / (oftype(jak, π) * jak)) : besselj(α - 1, jak)
    if refine && abs(α) != 0.5
        # Piessens's first six zeros and low-index McMahon approximations
        # need refinement before the direct quadrature expansion. One Newton
        # correction suffices at these already accurate starting values.
        # Reuse the weight's Bessel value; its differential recurrence updates
        # it to the corrected argument with a negligible quadratic remainder.
        value = besselj(α, jak)
        step = value / (derivative - α / jak * value)
        corrected = jak - step
        derivative = muladd((α - 1) / jak * derivative - value, corrected - jak, derivative)
        jak = corrected
    end
    xs, ws = gl_bessel(jak, d, α)

    xk, xdelta = (T > 0) ? sum_explicit(xs, (T - 1) >> 1) : sum_adaptive(xs)
    wk, wdelta = (T > 0) ? sum_explicit(ws, (T - 1) >> 1) : sum_adaptive(ws)

    xk = (jak^2 * d) * (1 + xk)
    xdelta *= jak^2 * d

    wfactor = 4d * xk^α * exp(-xk) / derivative^2
    wk = wfactor * (1 + wk)
    wdelta *= wfactor

    return xk, wk, max(xdelta, wdelta)
end

function gausslaguerre_asy0_bessel(n, jak, d, T)
    xs, ws = gl_bessel(jak, d)

    xk, xdelta = (T > 0) ? sum_explicit(xs, (T - 1) >> 1) : sum_adaptive(xs)
    wk, wdelta = (T > 0) ? sum_explicit(ws, (T - 1) >> 1) : sum_adaptive(ws)

    xk = (jak^2 * d) * (1 + xk)
    xdelta *= jak^2 * d

    wfactor = 4d * exp(-xk) / besselj(-1, jak)^2
    wk = wfactor * (1 + wk)
    wdelta *= wfactor

    return xk, wk, max(xdelta, wdelta)
end

function compute_airy_root(n, k)
    index = n - k + 1
    if index ≤ 11
        return AIRY_ROOTS[index]
    else
        t = 3 * π / 2 * (index - 0.25)
        return -t^(2 / 3) * evalpoly(inv(t * t), (6967296, 725760, -967680, 6478500, -10856875)) / 6967296
    end
end

function gausslaguerre_asy_airy(n, α, k, d, T)
    ak = compute_airy_root(n, k)
    xs = gl_airy(ak, d, α)

    xk, xdelta = (T > 0) ? sum_explicit(xs, (T + 1) >> 1) : sum_adaptive(xs)

    wk = 4^(1 / 3) * xk^(α + 1 / 3) * exp(-xk) / (airyaiprime(ak))^2
    return xk, wk, max(xdelta, abs(wk))
end

########################## Routines for the forward recurrence ##########################

function gl_rec_newton(x0, n, α; maxiter = 20, computeweight = true)
    T = eltype(x0)
    step = x0
    iter = 0
    xk = x0

    xk_prev = xk
    pn_prev = floatmax(T)
    pn_deriv = zero(T)
    while (abs(step) > 40eps(T) * xk) && (iter < maxiter)
        iter += 1
        pn, pn_deriv = evalLaguerreRec(n, α, xk)
        if abs(pn) ≥ abs(pn_prev) * (1 - 50eps(T))
            # The function values do not decrease enough any more due to roundoff errors.
            xk = xk_prev # Set to the previous value and quit.
            break
        end
        step = pn / pn_deriv
        xk_prev = xk
        xk -= step
        pn_prev = pn
    end
    if (xk < 0) || (xk > 4n + 2α + 2) || (iter == maxiter)
        @warn "Newton method may not have converged in gausslaguerre_rec($n,$α)."
    end
    wk = oftype(xk, 0)
    if computeweight
        # Evaluate the derivative at the accepted node, including when Newton
        # reverted a step. The derivative from the loop belongs to an older node.
        _, pn_deriv = evalLaguerreRec(n, α, xk)
        # x*p_n' = n*p_n + sqrt(n*(n+α))*p_{n-1}; at a root
        # this removes the poorly conditioned evaluation of p_{n-1} near zero.
        wk = inv(pn_deriv) / xk / pn_deriv
    end

    return xk, wk
end

"""
    FastGaussQuadrature.gausslaguerre_rec(n, α; reduced=false)

Compute a Gauss–Laguerre rule of order `0 ≤ n ≤ 127` using polynomial recurrence
and Newton iteration. Larger orders throw `DomainError`; use `gausslaguerre`
for automatic selection of a linear-work large-order algorithm.
`float(α)` must be `Float64`; use `gausslaguerre` for other working formats.
This bounded routine takes quadratic work at bounded iteration count. Initial
guesses use `Float64` Bessel-root approximations and extrapolation; the routine
warns if Newton iteration may not have converged. Prefer `gausslaguerre` for
automatic algorithm selection and generic-precision rules.

With `reduced=true`, stop at the first weight below
`10floatmin(typeof(float(α)))` and omit the remaining tail. Otherwise return
all nodes, with weights beyond that cutoff set to zero after the first
below-cutoff weight.
"""
function gausslaguerre_rec(n, α; reduced = false)
    0 ≤ n ≤ 127 || throw(DomainError(n, "The recurrence helper supports 0 ≤ n ≤ 127; use gausslaguerre for larger rules."))
    T = typeof(float(α))
    T === Float64 || throw(ArgumentError("The recurrence helper supports Float64 parameters only; use gausslaguerre for other floating-point types."))

    w = zeros(T, reduced ? 0 : n)
    x = zeros(T, reduced ? 0 : n)

    # We compute up to 7 starting values for the Newton iterations
    # j_{α,k}²/(4n+2α+2) is a lower bound by DLMF 18.16.10.
    x_pre = T.(approx_besselroots(α, min(n, 7))) .^ 2 / (4n + 2α + 2)

    noUnderflow = true  # false after the first weight below the tail cutoff
    for k in 1:n
        # Use sextic extrapolation for a new initial guess
        xk = (k ≤ min(n, 7)) ? x_pre[k] : 7 * x[k - 1] - 21 * x[k - 2] + 35 * x[k - 3] - 35 * x[k - 4] + 21 * x[k - 5] - 7 * x[k - 6] + x[k - 7]

        xk, wk = gl_rec_newton(xk, n, α, maxiter = 20, computeweight = noUnderflow)
        if noUnderflow && abs(wk) < 10floatmin(T)
            noUnderflow = false
        end

        if reduced
            if !noUnderflow
                validate_rule(x, w)
                return x, w
            else
                push!(x, xk); push!(w, wk)
            end
        else
            x[k] = xk
            w[k] = wk
        end
    end

    validate_rule(x, w)
    return x, w
end


"""
Evaluate the orthonormal associated Laguerre polynomial with positive leading coefficient,
as well as its derivative, in the point x using the recurrence relation.
"""
function evalLaguerreRec(n, α, x)
    T = typeof(α)
    pnprev = zero(T)
    pn = 1 / sqrt(gamma(α + 1))
    pndprev = zero(T)
    pnd = zero(T)
    for k in 1:n
        pnold = pn
        pn = (x - 2 * k - α + 1) / sqrt(k * (α + k)) * pnold - sqrt((k - 1 + α) * (k - 1) / k / (k + α)) * pnprev
        pnprev = pnold
        pndold = pnd
        pnd = (pnold + (x - 2 * k - α + 1) * pndold) / sqrt(k * (α + k)) - sqrt((k - 1 + α) * (k - 1) / k / (α + k)) * pndprev
        pndprev = pndold
    end

    return pn, pnd
end

# Gil--Segura--Temme, Numer. Math. 143 (2019), §4, and
# SIAM J. Sci. Comput. 48 (2026), §3. With x=z², the scaled polynomial
# y=z^(α+1/2)exp(-z²/2)L_n^α(z²) satisfies
# z²y'' + (-z⁴+(4n+2α+2)z²+1/4-α²)y = 0.
# A continued fraction or centered recurrence supplies the starting values.
# Local Taylor series and fourth-order iteration then sweep in both directions
# from the maximum of the normal-form coefficient.

function laguerre_gst_start(n, α::T, z) where {T}
    x = z * z
    if α > n
        # A centered recurrence avoids a very long continued fraction when
        # the parameter is larger than the polynomial degree.
        offset = x - α
        p, q = one(T), one(T) - offset
        dp, dq = zero(T), -one(T)
        limit = sqrt(floatmax(T))
        for k in 1:(n - 1)
            nextq = ((2k + 1 - offset) * q - (k + α) * p) / (k + 1)
            nextdq = ((2k + 1 - offset) * dq - q - (k + α) * dp) / (k + 1)
            p, q, dp, dq = q, nextq, dq, nextdq
            scale = max(abs(p), abs(q), abs(dp), abs(dq))
            if scale > limit
                p, q, dp, dq = p / scale, q / scale, dp / scale, dq / scale
            end
        end
        derivative = ((one(T) / 2 - offset) / z) * q + 2z * dq
        scale = max(abs(q), abs(derivative))
        return q / scale, derivative / scale
    end
    tiny = T === Float64 ? T(1.0e-300) : sqrt(floatmin(T))
    # L_n^α/L_n^(α-1)=(n+α)/f, with
    # f=x+α - x(n+α+1)/(x+α+1 - x(n+α+2)/(...)).
    f = x + α
    abs(f) < tiny && (f = copysign(tiny, f))
    c, d = f, zero(T)
    for k in 1:100000
        d = muladd(-x * (n + α + k), d, x + α + k)
        c = (x + α + k) + (-x * (n + α + k)) / c
        abs(d) < tiny && (d = copysign(tiny, d))
        abs(c) < tiny && (c = copysign(tiny, c))
        d = inv(d)
        f *= c * d
        if abs(c * d - 1) ≤ eps(T)
            derivative = ((one(T) / 2) - α + 2f) / z - z
            scale = max(one(T), abs(derivative))
            return inv(scale), derivative / scale
        end
    end
    error("Laguerre continued fraction did not converge in $T; try a type with higher precision.")
end

function laguerre_gst_count(n, α::T, x) where {T}
    # Sturm count for the tridiagonal Laguerre matrix, without eigenvectors.
    pivot = α + 1 - x
    count = Int(pivot < 0)
    for k in 2:n
        iszero(pivot) && (pivot = -floatmin(T))
        pivot = 2k + α - 1 - x - ((k - 1) * (k - 1 + α)) / pivot
        count += pivot < 0
    end
    return count
end

# Stored reciprocals avoid a division per Taylor term in hardware arithmetic.
# The 107 entries cover k=0:106 at 53 bits; higher precisions divide natively.
const LAGUERRE_TAYLOR_RECIPROCALS = ntuple(k -> 1.0 / (k * (k + 1)), 107)

function laguerre_gst_taylor(z::T, y, derivative, h, ν, c) where {T}
    iszero(h) && return y, derivative
    # For α = ±1/2 the inverse-square term vanishes and the normal-form ODE
    # is exactly the Hermite equation. Use its shorter recurrence with ordinary
    # native sums for the scaled Taylor terms.
    iszero(c) && return hermite_gst_taylor(ν, z, h, y, derivative, false)
    # Propagate t_k=y^(k)(z)h^k/k! by coefficient matching in the ODE.
    # The scaled coefficients avoid forming very large derivatives.
    iz = inv(z)
    r, h2 = h * iz, h * h
    q0 = h2 * (ν - z * z + c * iz * iz)
    q1 = h2 * r * (2ν - 4z * z)
    q2 = h2 * r * r * (ν - 6z * z)
    q3, q4 = -4h2 * h2 * r, -h2 * h2 * r * r
    tm4 = tm3 = tm2 = tm1 = zero(T)
    t0, t1 = y, h * derivative
    value, scaled_derivative = y + t1, t1
    small_terms = 0
    for k in 0:max(100, 2precision(T))
        numerator = -muladd(
            2r * k * (k + 1), t1,
            muladd(
                r * r * k * (k - 1) + q0, t0,
                muladd(q1, tm1, muladd(q2, tm2, muladd(q3, tm3, q4 * tm4)))
            )
        )
        # Keep a single rounded division for the leading Taylor terms.
        t2 = k >= 4 && (T === Float64 || T === Float32) ? numerator * T(LAGUERRE_TAYLOR_RECIPROCALS[k + 1]) : numerator / ((k + 1) * (k + 2))
        value += t2
        scaled_derivative += (k + 2) * t2
        # Require four consecutive small terms for both y and h*y'.
        if abs(t2) <= (eps(T) / 32) * max(abs(value), abs(h * derivative)) &&
                abs((k + 2) * t2) <= (eps(T) / 32) * max(abs(scaled_derivative), abs(h * derivative))
            small_terms += 1
            small_terms ≥ 4 && return value, scaled_derivative / h
        else
            small_terms = 0
        end
        tm4, tm3, tm2, tm1, t0, t1 = tm3, tm2, tm1, t0, t1, t2
    end
    error("Laguerre Taylor series did not converge in $T; try a type with higher precision.")
end

function laguerre_gst_advance(z::T, y, derivative, target, ν, c) where {T}
    # The origin is a singularity; subdivide a step to stay safely inside
    # the Taylor series' convergence disk. At c=0 the ODE is entire, so no
    # subdivision is needed. Always use rounded displacements.
    isfinite(target) && target > 0 || error("Laguerre iteration left the positive axis in $T; try a type with higher precision.")
    while !iszero(c) && abs(target - z) > (T(7) / 20) * z
        next = z + copysign((T(7) / 20) * z, target - z)
        next != z || error("Laguerre Taylor step is not representable in $T; try a type with higher precision.")
        y, derivative = laguerre_gst_taylor(z, y, derivative, next - z, ν, c)
        z = next
    end
    y, derivative = laguerre_gst_taylor(z, y, derivative, target - z, ν, c)
    return target, y, derivative
end

function laguerre_gst_root(z::T, y, derivative, ν, c, direction, nextroot) where {T}
    a = iszero(c) ? ν - z * z : ν - z * z + c / (z * z)
    isfinite(a) && a > 0 || error("Laguerre iteration cannot resolve the next node in $T; try a type with higher precision.")
    frequency = sqrt(a)
    angle = atan(frequency * y / derivative)
    if nextroot || direction * angle ≥ 0
        angle -= direction * T(π)
    end
    target = z - angle / frequency
    for _ in 1:max(20, precision(T))
        z, y, derivative = laguerre_gst_advance(z, y, derivative, target, ν, c)
        a = iszero(c) ? ν - z * z : ν - z * z + c / (z * z)
        frequency = sqrt(abs(a))
        ratio = frequency * y / derivative
        if a < 0 && !(abs(ratio) < 1)
            error("Laguerre iteration cannot resolve the next node in $T; try a type with higher precision.")
        end
        step = iszero(frequency) ? y / derivative : (a > 0 ? atan(ratio) : atanh(ratio)) / frequency
        target = z - step
        # At half-integer α this is the Hermite iteration, whose accepted
        # fourth-order correction has relative error -step^4/6 to leading order.
        if abs(step) ≤ (iszero(c) ? sqrt(sqrt(eps(T))) / 2 : 2eps(T) * abs(z)) || target == z
            return laguerre_gst_advance(z, y, derivative, target, ν, c)
        end
    end
    error("Laguerre fourth-order iteration did not converge in $T; try a type with higher precision.")
end

function gausslaguerre_gst(n::Integer, α::T; reduced = false) where {T <: AbstractFloat}
    endpoint = α + one(T) <= one(T) / 8
    firstx, firstw, logp = endpoint ? laguerre_endpoint(n, α) : (zero(T), zero(T), zero(T))
    μ = gamma(α + 1)
    ν, c = 4n + 2α + 2, ((one(T) / 2) - α) * ((one(T) / 2) + α)
    isfinite(μ) && μ > 0 || error("Laguerre normalization is not representable in $T; try a type with a wider exponent range.")
    isfinite(ν) && isfinite(c) || error("Laguerre parameters are not representable in $T; try a type with a wider exponent range.")
    if abs(α) > one(T) / 2
        pivot = sqrt(-c)
        left = laguerre_gst_count(n, α, pivot)
    else
        pivot, left = zero(T), 0
    end
    if iszero(left)
        # GST (2019), equation (30), bounds the first zero for every α>-1.
        # Near |α|=1/2, the stationary point tends to zero; starting at least
        # this far from the origin avoids an ill-conditioned continued fraction.
        upper = (
            2n^2 + n * (α - 1) + 2 * (α + 1) +
                2 * (n - 1) * sqrt(n^2 + (n + 2) * (α + 1))
        ) / (n + 2)
        lower = ((α + 1) * (n * (α + 5) + 2 * (α - 1))) / ((n + 2) * upper)
        pivot = max(pivot, lower)
    end
    z0 = sqrt(pivot)
    y0, derivative0 = laguerre_gst_start(n, α, z0)
    # Positivity and polynomial exactness give wᵢ ≤ Γ(α+m+1)/xᵢ^m for
    # 0 ≤ m ≤ 2n-1. Stop beyond a point where this bound is below the
    # cutoff, rather than assuming that successive quadrature weights decrease.
    # The omitted mass is at most n times the cutoff, far below roundoff.
    stopx = if reduced
        logcutoff = log(10floatmin(T)) - log(T(4))
        m = min(2n - 1, max(1, ceil(Int, -logcutoff)))
        exp((loggamma(α + m + 1) - logcutoff) / m)
    else
        T(Inf)
    end
    x, weights = Vector{T}(undef, left), Vector{T}(undef, left)
    expected = reduced ? min(n, max(left, ceil(Int, 2sqrt(T(n) * stopx) / T(π)))) : n
    sizehint!(x, expected)
    sizehint!(weights, expected)
    for (direction, indices) in ((-1, left:-1:1), (1, (left + 1):n))
        z, y, derivative = z0, y0, derivative0
        for i in indices
            endpoint && i == 1 && break
            z, y, derivative = laguerre_gst_root(z, y, derivative, ν, c, direction, i != first(indices))
            # Evaluate the weight in the native z coordinate. Factoring the
            # difference of squares retains the scale near the starting point.
            square = z * z
            if direction > 0
                push!(x, square)
                push!(weights, zero(T))
            else
                x[i] = square
            end
            isfinite(x[i]) && x[i] > 0 || error("A Laguerre node is not representable in $T; try a type with higher precision.")
            distance = z - z0
            logratio = z < z0 / 2 ? log(z) - log(z0) : log1p(distance / z0)
            weights[i] = muladd(2α + 1, logratio, -distance * (z + z0)) - 2log(abs(derivative))
            reduced && direction > 0 && x[i] >= stopx && break
        end
    end
    if endpoint
        # The first weight can dominate the total mass near α = -1.  Anchor
        # the remaining weights at the second node instead of subtracting it.
        δ = α + one(T)
        term, p = one(T), one(T)
        for k in 1:(n - 1)
            term *= ((k - n) / T(k)) * (x[2] / (δ + k))
            p += term
            abs(term) <= eps(T) * abs(p) / 8 && break
        end
        anchor = gamma(one(T) + δ) * exp(-logp) / (n * x[2] * p^2)
        isfinite(anchor) && anchor > 0 || error("A Laguerre weight is not representable in $T; try a type with higher precision.")
        offset = log(anchor) - weights[2]
        @inbounds for i in 2:length(weights)
            weights[i] = exp(weights[i] + offset)
        end
        x[1], weights[1] = firstx, firstw
    else
        shift = maximum(weights)
        if μ <= one(T)
            # Since max(exp(logw-shift))=1, the normalization scale is at most
            # one. No underflowed exponential can be rescued into a normal
            # weight, so compute each exponential only once in this regime.
            weights .= exp.(weights .- shift)
            scale = μ / sum(weights)
            weights .*= scale
        else
            scale = μ / sum(t -> exp(t - shift), weights)
            @inbounds for i in eachindex(weights)
                exponent = weights[i] - shift
                # Multiplying in two pieces avoids underflow before a large gamma
                # normalization restores a representable ordinary weight.
                weights[i] = exponent > log(floatmin(T)) / 2 ? exp(exponent) * scale :
                    exp(exponent / 2) * scale * exp(exponent / 2)
            end
        end
    end
    validate_rule(x, weights)
    return x, weights
end

# Near α = -1, solve the normalized polynomial ₁F₁(-n; δ; x) = 0
# in t = n*x/δ.  Its first root is close to t = 1, avoiding the
# cancellation in a direct expansion of the node.  If d is the derivative
# with respect to t and p = ∏ₖ₌₁ⁿ⁻¹(1 + δ/k), the weight is
# gamma(δ)/(p*t*d^2).
function laguerre_endpoint(n::Integer, α::T) where {T <: AbstractFloat}
    δ = α + one(T)
    t = one(T)
    d = -one(T)
    for it in 1:max(10, 2ceil(Int, log2(precision(T))))
        term = -t
        f, d = one(T) - t, -one(T)
        for k in 2:n
            term *= ((k - 1 - n) / T(k)) * ((δ * t / n) / (δ + k - 1))
            f += term
            d += k * term / t
            abs(term) <= eps(T) * max(abs(f), abs(d)) / 8 && break
        end
        next = t - f / d
        (next == t || abs(next - t) <= eps(T) * abs(t)) && break
        t = next
        it == max(10, 2ceil(Int, log2(precision(T)))) && error("Laguerre endpoint iteration did not converge in $T; try a type with higher precision.")
    end
    logp = T === Float64 ? -logbeta(T(n), δ) - log(δ) :
        sum(k -> log1p(δ / k), 1:(n - 1); init = zero(T))
    firstx = δ / n * t
    firstw = gamma(δ) * exp(-logp) / (t * d^2)
    isfinite(firstx) && firstx > 0 && isfinite(firstw) && firstw > 0 ||
        error("The first Laguerre node or weight is not representable in $T; try a type with higher precision or wider exponent range.")
    return firstx, firstw, logp
end
