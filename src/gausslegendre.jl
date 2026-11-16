@doc raw"""
    gausslegendre([T=Float64,] n::Integer) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Legendre quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature) with type `T`.

For `n > 60` and `precision(T) <= precision(Float64)`, evaluate Bogaert's
expansion with a precision-dependent number of terms. Precisions below `Float32` compute in
`Float32` and convert the completed rule. Higher precisions use Jacobi iteration.
An error is thrown if the returned nodes are not distinct in `T`.

```math
\int_{-1}^{1} f(x) dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gausslegendre(3);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ 2/5
true
```
"""
@inline function gausslegendre(::Type{T}, n::Integer) where {T}
    # GAUSSLEGENDRE(n) COMPUTE THE GAUSS-LEGENDRE NODES AND WEIGHTS IN O(n) time.

    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    elseif n == 0
        x, w = T[], T[]
    elseif n == 1
        x, w = T[0], T[2]
    elseif precision(T) < precision(Float32)
        x, w = gausslegendre(Float32, n)
        x, w = T.(x), T.(w)
    elseif n == 2
        sqrt3 = sqrt(T(3))
        x, w = [T(-1) / sqrt3, 1 / sqrt3], T[1, 1]
    elseif n == 3
        sqrt3_5 = sqrt(T(3) / 5)
        x, w = [-sqrt3_5, T(0), sqrt3_5], [5, 8, 5] ./ T(9)
    elseif n == 4
        a = 2 * sqrt(T(6) / 5) / 7
        T3_7 = T(3) / 7
        sqrt30 = sqrt(T(30))
        x, w = (
            [-sqrt(T3_7 + a), -sqrt(T3_7 - a), sqrt(T3_7 - a), sqrt(T3_7 + a)],
            [
                (18 - sqrt30) / 36, (18 + sqrt30) / 36,
                (18 + sqrt30) / 36, (18 - sqrt30) / 36,
            ],
        )
    elseif n == 5
        b = 2 * sqrt(T(10) / 7)
        sqrt70 = sqrt(T(70))
        x, w = (
            [
                -sqrt(5 + b) / 3, -sqrt(5 - b) / 3, T(0),
                sqrt(5 - b) / 3, sqrt(5 + b) / 3,
            ],
            [
                (322 - 13 * sqrt70) / 900, (322 + 13 * sqrt70) / 900, T(128) / 225,
                (322 + 13 * sqrt70) / 900, (322 - 13 * sqrt70) / 900,
            ],
        )
    elseif n ≤ 60
        # NEWTON'S METHOD WITH THREE-TERM RECURRENCE:
        return rec(T, n)
    elseif precision(T) <= precision(Float64)
        # The expansion and coefficient tables are accurate to Float64 precision.
        return asy(n, T)
    else
        # FALL BACK TO GAUSS-JACOBI at higher precision.
        return jacobi_gst(n, zero(T), zero(T))
    end
    validate_rule(x, w)
    return x, w
end
gausslegendre(n::Integer) = gausslegendre(Float64, n)

function asy(n, ::Type{T} = Float64) where {T}
    # Bogaert's direct expansions, evaluated through the smooth G/H functions.
    # Keep the Bessel phase residual separate to preserve small central nodes.
    # The stored coefficients are accurate through Float64 precision. Convert
    # them once so every expansion operation uses the requested working type.
    x = Vector{T}(undef, n)
    w = Vector{T}(undef, n)
    g1coeffs, g2coeffs, g3coeffs, g4coeffs = map(T, LEGENDRE_G1), map(T, LEGENDRE_G2),
        map(T, LEGENDRE_G3), map(T, LEGENDRE_G4)
    h2coeffs, h3coeffs, h4coeffs = map(T, LEGENDRE_H2), map(T, LEGENDRE_H3), map(T, LEGENDRE_H4)
    pi = T(π)
    N = T(n) + T(0.5)
    v = inv(N)
    # For 0 <= a < 1.58, q <= (1.6/N)^2 and the coefficient 1-norms
    # of G2:G4 and H2:H4 are below 1/32. Choose the retained depth so
    # the additional omitted angle is at most eps(T)/(32*N*(1-q)).
    qmax = ((T(8) / 5) * v)^2
    q2 = qmax * qmax
    terms = if N * q2 <= eps(T)
        1
    elseif N * q2 * qmax <= eps(T)
        2
    elseif N * q2 * q2 <= eps(T)
        3
    else
        4
    end
    @inbounds for k in 1:((n + 1) >> 1)
        phase = pi * (T(k) - T(0.25))
        d, beta_correction = legendre_bessel(k, phase)
        a = v * (phase + d)
        central = k > n ÷ 4
        sb, cb = central ? sincospi(T(n + 1 - 2k) / T(2n + 1)) : sincos(a)
        vd = v * d
        sa = central ? muladd(cb, -vd^2 / 2, muladd(sb, vd * (1 - vd^2 / 6), cb)) : sb
        r = a / sa
        t = (a / T(1.58))^2
        g1 = evalpoly(t, g1coeffs)
        q = (v * v) * r^2
        g = g1
        h = zero(T)
        if terms == 4
            g = evalpoly(q, (g1, evalpoly(t, g2coeffs), evalpoly(t, g3coeffs), evalpoly(t, g4coeffs)))
            h = evalpoly(q, (evalpoly(t, h2coeffs), evalpoly(t, h3coeffs), evalpoly(t, h4coeffs)))
        elseif terms == 3
            g = evalpoly(q, (g1, evalpoly(t, g2coeffs), evalpoly(t, g3coeffs)))
            h = muladd(q, evalpoly(t, h3coeffs), evalpoly(t, h2coeffs))
        elseif terms == 2
            g = muladd(q, evalpoly(t, g2coeffs), g1)
            h = evalpoly(t, h2coeffs)
        end
        correction = a / r * q * g
        if central
            delta = muladd(v, d, correction)
            x[n + 1 - k] = muladd(sb, -delta^2 / 2, muladd(-cb, delta * (1 - delta^2 / 6), sb))
        else
            # cos(a + correction), with a tiny-angle rotation. At n >= 61
            # the omitted fourth-order term is below Float64 rounding error.
            x[n + 1 - k] = muladd(cb, -correction^2 / 2, muladd(-sa, correction * (1 - correction^2 / 6), cb))
        end
        # W1 = 1/8 + F1/a. Retain it at every order, including n > 850000.
        # w_k = π v sin(a) / (beta * (1 + q² H + v²(r G₁ + 1/8))).
        x[k] = -x[n + 1 - k]
        c = muladd(q^2, h, (v * v) * muladd(r, g1, T(0.125)))
        if k <= 20
            # Combine π*j/beta = 2/J₁(j)^2 before rounding. Retain sinc(a)-1
            # and c separately so rounding near 1 does not erase corrections.
            base = T(LEGENDRE_BESSEL_WEIGHT_FACTORS[k]) / (N * N)
            adjustment = (legendre_sincm1(a) - c) / (1 + c)
            w[k] = w[n + 1 - k] = muladd(adjustment, base, base)
        else
            # beta*(1+c) = 1+delta; evaluate division as a relative correction.
            delta = muladd(beta_correction, c, beta_correction + c)
            base = (pi * sa) / N
            w[k] = w[n + 1 - k] = muladd(-delta / (1 + delta), base, base)
        end
    end
    @inbounds isodd(n) && (x[(n + 1) >> 1] = zero(T))
    validate_rule(x, w)
    return x, w
end

@inline function legendre_bessel(k, phase::T) where {T}
    if k <= 20
        return T(LEGENDRE_BESSEL_RESIDUALS[k]), T(LEGENDRE_BESSEL_BETA_CORRECTIONS[k])
    end
    # McMahon residual and (π/2)j*J1(j)^2 - 1, retaining the small correction.
    # DLMF 10.17.3 and 10.21.19.
    return (
        inv(phase) * evalpoly(inv(phase)^2, map(T, (1 / 8, -31 / 384, 3779 / 15360, -6277237 / 3440640, 2092163573 / 82575360))),
        inv(phase)^2 * evalpoly(inv(phase)^2, map(T, (1 / 8, -29 / 128, 3599 / 3072, -6102707 / 491520, 2058116637 / 9175040))),
    )
end

function legendre_sincm1(a::T) where {T}
    # sin(a)/a - 1 for the first 20 roots and n>=61, where 0<a<1.01.
    # The alternating terms decrease; the first omitted term bounds the tail.
    z = a * a
    term = -z / 6
    total = term
    for k in 2:precision(T)
        term *= -z / ((2k) * (2k + 1))
        abs(term) <= eps(T) * abs(total) && return total
        total += term
    end
    error("Legendre sinc series did not converge at the requested precision")
end

function legpts_nodes(n, a::AbstractVector{T}) where {T}
    # Asymptotic seeds, evaluated at no more than 53 bits, for 6 ≤ n ≤ 60.
    nodes = cot.(a)
    vn² = (1 / (T(n) + T(0.5)))^2
    c2 = map(T, (2595 / 15360, 6350 / 15360, 3779 / 15360))
    c31, c32, c33, c4 = T(31 / 1024), T(11 / 1024), T(25 / 3072), T(1073 / 5120)
    @inbounds for i in 1:length(a)
        u = nodes[i]
        u² = u^2
        ai = a[i]
        ai² = ai * ai
        ai³ = ai² * ai
        node = ai + (u - 1 / ai) / 8 * vn²
        v1 = (6 * (1 + u²) / ai + 25 / ai³ - u * muladd(31, u², 33)) / 384
        v2 = u * evalpoly(u², c2)
        v3 = (1 + u²) * (
            -muladd(c31, u², c32) / ai +
                u / 512 / ai² + -c33 / ai³
        )
        v4 = (v2 - c4 / (ai² * ai³) + v3)
        node = muladd(v1, vn² * vn², node)
        node = muladd(v4, (vn² * vn²) * vn², node)
        nodes[i] = node
    end
    @inbounds for jj in 1:length(a)
        nodes[jj] = cos(nodes[jj])
    end

    return nodes
end

function rec(::Type{T}, n) where {T}
    # COMPUTE GAUSS-LEGENDRE NODES AND WEIGHTS USING NEWTON'S METHOD.
    # THREE-TERM RECURENCE IS USED FOR EVALUATION. COMPLEXITY O(n^2).

    # Seeds use no wider arithmetic than T; Newton recovers its precision.
    x = leg_initial_guess(n, T)

    # Perform Newton to find zeros of Legendre polynomial:
    PP1, PP2 = similar(x), similar(x)

    if T === Float64
        # The accurate seeds need two Newton steps.
        innerRec!(PP1, PP2, n, x)
        @inbounds for i in eachindex(x)
            x[i] -= PP1[i] / PP2[i]
        end
        innerRec!(PP1, PP2, n, x)
        @inbounds for i in eachindex(x)
            step = PP1[i] / PP2[i]
            oldx = x[i]
            x[i] -= step
            # Transport (1-x²)P'² through the final tiny step, using the
            # Legendre ODE. Preserve the unrounded root correction in the
            # weight as well as in P'; omitted terms are O(step²).
            factor = muladd(-oldx, oldx, 1.0)
            factor = muladd(-2oldx, step, factor)
            factor = muladd(2n * (n + 1), step^2, factor)
            PP2[i] = 2 / (factor * PP2[i]^2)
        end
    else
        # Refine to the requested precision, then recompute the derivative.
        num_iters = 20 + ndigits(max(precision(T) - 1, 1); base = 2)
        for iteration in 1:num_iters
            innerRec!(PP1, PP2, n, x)
            converged = true
            @inbounds for i in eachindex(x)
                step = PP1[i] / PP2[i]
                x[i] -= step
                converged &= abs(step) <= 8eps(T) * max(abs(x[i]), inv(T(n)))
            end
            if converged
                innerRec!(PP1, PP2, n, x)
                break
            elseif iteration == num_iters
                error("Legendre Newton iteration did not converge at the requested precision")
            end
        end
        @inbounds for i in eachindex(x)
            # Divide in two stages so squaring the derivative cannot overflow
            # in a narrow-exponent type when the final weight is representable.
            # w_i = 2 / ((1 - x_i²) P_n'(x_i)²).
            PP2[i] = (T(2) / PP2[i]) / (((1 - x[i]) * (1 + x[i])) * PP2[i])
        end
    end

    # Both paths have replaced the derivatives in PP2 with weights.
    # Use symmetry to get the other Legendre nodes and weights:
    m = length(x) # Save the computed half-length before resizing x.
    resize!(x, n)
    resize!(PP2, n)
    @inbounds for i in 1:(m - 1)
        x[n + 1 - i] = -x[i]
        PP2[n + 1 - i] = PP2[i]
    end
    @inbounds isodd(n) && (x[(n + 1) ÷ 2] = zero(T))
    validate_rule(x, PP2)
    return x, PP2
end

@inline function innerRec!(myPm1, myPPm1, n, x)
    # EVALUATE LEGENDRE AND ITS DERIVATIVE USING THREE-TERM RECURRENCE RELATION.
    T = eltype(x)
    @inbounds for j in 1:length(x)
        xj = x[j]
        Pm2 = T(1)
        Pm1 = xj
        PPm1 = T(1)
        PPm2 = T(0)
        for k in 1:(n - 1)
            Pm2, Pm1 = Pm1, muladd((2 * k + 1) * Pm1, xj, - k * Pm2) / (k + 1)
            PPm2, PPm1 = PPm1, (
                    (2 * k + 1) * muladd(xj, PPm1, Pm2) -
                    k * PPm2
                ) / (k + 1)
        end
        myPm1[j] = Pm1
        myPPm1[j] = PPm1
    end
    return myPm1, myPPm1
end

@inline function leg_initial_guess(n, ::Type{T} = Float64) where {T}
    # Seeds for the first n÷2+1 roots, used by Newton for 6 ≤ n ≤ 60.
    # Accurate Float64 seeds keep the two-step Float64 iteration inexpensive;
    # higher precisions refine these lower-precision guesses in T.
    W = precision(T) < precision(Float64) ? T : Float64
    a = besselZeroRoots((n ÷ 2) + 1, W)
    rmul!(a, 1 / (W(n) + W(0.5)))
    x = legpts_nodes(n, a)
    rmul!(x, -one(W))
    return convert(Vector{T}, x)
end

# Verify with `julia --startup-file=no test/generate_legendre.jl`.
# Generated Bogaert coefficient tables; do not edit by hand.
# Smooth Bogaert coefficients in t = (a / 1.58)^2, 0 <= a <= 1.58.
const LEGENDRE_G1 = (
    -0.041666666666666664,
    0.010401666666666561,
    -0.0009273828809501513,
    4.2872567093440924e-5,
    -1.2162167086534596e-6,
    2.335495444630439e-8,
    -3.2373173382003673e-10,
    3.2876434565730634e-12,
)
const LEGENDRE_G2 = (
    0.008159722217286646,
    -0.005218031372275117,
    0.0017581556101258153,
    -0.0003940721599245086,
    6.290073552040636e-5,
    -7.2963446374149055e-6,
    5.323157582940204e-7,
)
const LEGENDRE_G3 = (
    -0.004160117999150595,
    0.0032114635558404834,
    -0.0015636946199803392,
    0.0006383497815630678,
    -0.0001960876851496558,
    3.225685804318767e-5,
)
const LEGENDRE_G4 = (
    0.0043270099115829365,
    -0.0028146344099921217,
    0.0013847291495654196,
    -0.0008100693263105142,
    0.000244753385829156,
)
const LEGENDRE_H2 = (
    -0.011111111111298923,
    0.0067143033750598794,
    -0.002538281803920811,
    0.0007249403085730113,
    -0.00014830475932460258,
    2.0534934070615628e-5,
    -1.7613825390116413e-6,
    6.836291635211902e-8,
)
const LEGENDRE_H3 = (
    0.006569662960802897,
    -0.00023646090297784855,
    -0.0006614590245693976,
    -0.0006391238107272089,
    0.0007234263089345481,
    -0.00029955185370806817,
    5.3050457982907744e-5,
)
const LEGENDRE_H4 = (
    -0.007569603223222529,
    -0.007504919110026087,
    0.005692867027199045,
    0.0016201340527773997,
    -0.002525587721852953,
    0.0007324450147449329,
)
const LEGENDRE_BESSEL_RESIDUALS = (
    0.04863106750342784,
    0.022290966504172484,
    0.014348115539080811,
    0.010561988052556969,
    0.008352603936268065,
    0.006906209769611422,
    0.005886218148154599,
    0.005128465428405139,
    0.004543413129563959,
    0.004078095931491043,
    0.003699187483291371,
    0.003384673983973428,
    0.0031194313583755044,
    0.0028927263170733285,
    0.0026967312123637515,
    0.0025256033585736677,
    0.002374893485959285,
    0.002241153801149329,
    0.0021216712723189117,
    0.0020142818287534232,
)
const LEGENDRE_BESSEL_BETA_CORRECTIONS = (
    0.018087218694297996,
    0.00392013621674922,
    0.001636551843562515,
    0.0008892796349325885,
    0.000556868684819165,
    0.0003809718021230442,
    0.0002768652596050394,
    0.00021022693776474386,
    0.00016502823272164556,
    0.00013297372137155157,
    0.0001094222755141307,
    9.16133256630134e-5,
    7.782167548719328e-5,
    6.692432209045039e-5,
    5.816481661508719e-5,
    5.101855968357026e-5,
    4.511249113423512e-5,
    4.017545054632166e-5,
    3.6006519266516735e-5,
    3.2454268447740936e-5,
)
const LEGENDRE_BESSEL_WEIGHT_FACTORS = (
    7.420761371418964,
    17.27411993534628,
    27.14206863490362,
    37.01128458651283,
    46.88075495998108,
    56.750301539491645,
    66.61987707748293,
    76.48946551165082,
    86.35906037972798,
    96.22865874398447,
    106.09825913898929,
    115.96786077832688,
    125.83746321430178,
    135.7070661793183,
    145.5766695067482,
    155.44627308915338,
    165.31587685514634,
    175.18548075603482,
    185.0550847578332,
    194.9246888363317,
)
