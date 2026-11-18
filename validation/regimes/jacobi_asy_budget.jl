module JacobiAsymptoticBudget
using SpecialFunctions

# Reference-only finite engineering model. No production call or measured error
# enters this file. See jacobi-asy-rationale.md for scope and assumptions.

function bernoulli_factors(count)
    b = Rational{BigInt}[1]
    for m in 1:2count
        push!(b, -sum(binomial(big(m + 1), k) * b[k + 1] for k in 0:(m - 1)) / (m + 1))
    end
    return [(-1)^k * (2k - 1) * b[2k + 1] / factorial(big(2k)) for k in 1:count]
end
const F_FACTORS = bernoulli_factors(33)

function product_coefficient(a, b, k)
    return sum(a[j + 1] * b[k - j + 1] for j in max(0, k + 1 - length(b)):min(k, length(a) - 1))
end

function boundary_next(a, b, theta)
    # f=g'/2 is analytic for |theta|<pi. q=theta^2, A0=1.
    # Four recurrences produce the first omitted B3 and A4 independently of
    # the production degree-nine Chebyshev coefficient approximation.
    aa, bb = BigFloat(1) / 4 - a * a, BigFloat(1) / 4 - b * b
    f = [BigFloat(F_FACTORS[k]) * (aa + bb * (big(4)^k - 1)) for k in 1:length(F_FACTORS)]
    av = zeros(BigFloat, length(f)); av[1] = 1
    bv = BigFloat[]
    for _ in 0:3
        count = length(av) - 1
        bv = [
            -(k + 1) * (1 + (1 + 2a) / (2k + 1)) * av[k + 2] +
                product_coefficient(f, av, k) / (2 * (2k + 1)) for k in 0:(count - 1)
        ]
        next = zeros(BigFloat, count)
        for k in 1:(count - 1)
            next[k + 1] = (k - a) * bv[k + 1] - product_coefficient(f, bv, k - 1) / (4k)
        end
        av = next
    end
    q = theta * theta
    # Factor2 reserves the remaining analytic power tail. theta<=1/2 and
    # 29+ coefficients put this well beyond the degree-nine production data.
    bn = 2sum(abs(bv[k]) * q^(k - 1) for k in eachindex(bv))
    an = 2sum(abs(av[k]) * q^(k - 1) for k in eachindex(av))
    return bn, an
end

function hahn_tail(n, a, b, theta)
    nu = n + (a + b + 1) / 2
    pa, pb = ones(BigFloat, 23), ones(BigFloat, 23)
    s, c = 2sin(theta / 2), 2cos(theta / 2)
    for k in 1:22
        pa[k + 1] = pa[k] * abs((k - BigFloat(1) / 2)^2 - a * a) / (k * s)
        pb[k + 1] = pb[k] * abs((k - BigFloat(1) / 2)^2 - b * b) / (k * c)
    end
    denominator = one(nu); envelope = BigFloat[]
    for m in 0:22
        m > 0 && (denominator *= 2nu + m)
        push!(envelope, sum(pa[l + 1] * pb[m - l + 1] for l in 0:m) / denominator)
    end
    ratio = iszero(envelope[22]) ? zero(nu) : envelope[23] / envelope[22]
    # The first possible signed stopping term has degree11. Cover all terms
    # that might then be omitted, including incomplete degrees18 and19.
    tail = 2sum(envelope[12:21]) + (ratio < 1 ? 2envelope[22] / (1 - ratio) : BigFloat(Inf))
    return tail, ratio
end

function prefactor(n, a, b, theta, boundary)
    nu = n + (a + b + 1) / 2
    logden = (a + BigFloat(1) / 2) * log(sin(theta / 2)) + (b + BigFloat(1) / 2) * log(cos(theta / 2))
    if boundary
        return exp(loggamma(n + a + 1) - loggamma(BigFloat(n + 1)) - a * log(nu) - logden) * sqrt(theta / 2)
    end
    return exp(
        2nu * log(BigFloat(2)) + loggamma(n + a + 1) + loggamma(n + b + 1) -
            loggamma(2n + a + b + 2) - log(BigFloat(pi)) - logden
    )
end

function polynomial_tail(n, a, b, theta, boundary)
    if boundary
        bn, an = boundary_next(a, b, theta)
        nu = n + (a + b + 1) / 2
        ratio = (1 + abs(a) + abs(b))^2 / nu^2
        # The B3 term is theta/nu^7; only in the fixed boundary chart is it
        # also O(nu^-8). Keep both terms explicit. Eight is the frozen finite
        # next-term reserve, not a uniform asymptotic remainder theorem.
        tail = 8 * (1 + abs(a)) * (theta * bn / nu^7 + an / nu^8) / (1 - ratio)
    else
        tail, ratio = hahn_tail(n, a, b, theta)
    end
    return prefactor(n, a, b, theta, boundary) * tail, ratio
end

function constant_error(n, a, b, e, kweight)
    # Actual weight constant uses a degree20 hypergeometric sum, with a
    # next-term stopping test. Bound its remaining full-depth term separately.
    term = one(e)
    for m in 0:20
        term *= -(m + a) * (m + b) / ((m + 1) * (n - m))
    end
    ratio = abs((21 + a) * (21 + b) / (22 * (n - 21)))
    scale = kweight / exp((a + b + 1) * log(BigFloat(2)))
    tail = ratio < 1 ? abs(term) / (abs(scale) * (1 - ratio)) : BigFloat(Inf)
    return 8e * (1 + abs(a + b + 1) + abs(a * b)) + e / 10 + 2tail, ratio
end

"""Finite Float64 ASY budget at an independently indexed interior root.

`n,i` are the underlying Jacobi degree/index, excluding fixed endpoints.
`w` is the final consumer weight: Jacobi, Radau, or Lobatto according to
`transform`. All numeric arguments are exact represented parameters/reference
values, evaluated at >=181 BigFloat bits. Extra Lobatto8e/32e diagnostics must
be checked independently. Unknown consumers/transforms are rejected.
"""
function point_budget(n, a, b, bits, i, x, w; transform = :none, consumer = :jacobi)
    bits == 53 || error("This fixed-depth ASY model is Float64-only")
    precision(BigFloat) >= bits + 128 || error("Reference precision required")
    transform in (:none, :radau, :lobatto) || error("Unknown weight transform")
    consumer in (:jacobi, :radau, :lobatto) || error("Unknown ASY consumer")
    a, b, x, w = BigFloat(a), BigFloat(b), BigFloat(x), BigFloat(w)
    ordinary = consumer == :jacobi && transform == :none && n >= 768 &&
        min(a, b) >= BigFloat(-0.99) && max(a, b) < 5
    standard_radau = consumer == :radau && transform == :radau && n >= 101 && iszero(a) && b == 1
    generalized_radau = consumer == :radau && transform == :radau && n >= 500 && -1 < a < 5 && 0 < b < 5
    lobatto = consumer == :lobatto && transform == :lobatto && n >= 768 && a == b == 1
    ordinary || standard_radau || generalized_radau || lobatto || error("Outside actual finite ASY selection domain")
    1 <= i <= n && abs(x) < 1 && w > 0 || error("Invalid interior reference")
    e = BigFloat(2)^(1 - bits)
    s = sqrt((1 - x) * (1 + x)); nu = n + (a + b + 1) / 2
    reflected = x < 0
    aa, bb = reflected ? (b, a) : (a, b)
    theta = acos(abs(x)); boundary = min(i, n + 1 - i) <= 10
    (!boundary || theta <= BigFloat(1) / 2) || error("Boundary coefficient chart exceeded")
    divisor = transform == :none ? one(e) : transform == :radau ? 1 + x : s * s
    baseweight = w * divisor
    kw = exp(
        (a + b + 1) * log(BigFloat(2)) + loggamma(n + a + 1) + loggamma(n + b + 1) -
            loggamma(BigFloat(n + 1)) - loggamma(n + a + b + 1)
    )
    derivative = sqrt(kw / baseweight)
    pn, r1 = polynomial_tail(n, aa, bb, theta, boundary)
    pm, r2 = polynomial_tail(n - 1, aa, bb, theta, boundary)
    xlocal = abs(x); degree = 2n + a + b
    dtail = (abs(n * (aa - bb - degree * xlocal)) * pn + 2 * (n + aa) * (n + bb) * pm) / (degree * s)
    angle_tail = pn / derivative
    # A coherent phase perturbation shifts the root by O(e theta). Residual
    # amplitude/derivative cancellation has its own location-dependent chi.
    d = 1 + a * a + b * b
    chi = d / ((boundary ? nu * nu : nu) * s * s)
    logscale = 1 + abs(a) + abs(b) + abs((aa + BigFloat(1) / 2) * log(sin(theta / 2))) +
        abs((bb + BigFloat(1) / 2) * log(cos(theta / 2))) + log(nu)
    angle_round = e * (8theta + 8 * (1 + d + chi) / nu)
    angle = angle_round + angle_tail
    node = 2e * abs(x) + s * angle + angle * angle / 2
    a == b && isodd(n) && i == (n + 1) ÷ 2 && (node = zero(e))
    sensitivity = abs(a - b + (a + b + 1) * x) / s
    derivative_error = dtail / derivative + 8e * (1 + d + chi + logscale) + sensitivity * angle
    transform_sensitivity = transform == :none ? zero(e) : transform == :radau ? s / (1 + x) : 2abs(x) / s
    transform_error = transform == :none ? zero(e) : 4e + transform_sensitivity * angle
    mass_error, r3 = constant_error(n, a, b, e, kw)
    weight = (1 + mass_error) * (1 + 8e) / ((1 - derivative_error)^2 * (1 - transform_error)) - 1
    conclusive = all(isfinite, (node, weight)) && node >= 0 && 0 <= weight < BigFloat(1) / 8 &&
        max(r1, r2, r3) < BigFloat(1) / 2 && derivative_error < BigFloat(1) / 8 &&
        transform_error < BigFloat(1) / 8 && nu * angle < BigFloat(1) / 8
    return (;
        node_absolute = node, weight_relative = weight, conclusive,
        regime = boundary ? "jacobi_ASY_boundary" : "jacobi_ASY_interior",
        node_tail = s * angle_tail + angle_tail^2 / 2, derivative_tail = dtail / derivative,
        angle_round, chi, derivative_error, transform_error, mass_error,
        extra_node_absolute = consumer == :lobatto ? 8e : missing,
        extra_weight_relative = consumer == :lobatto ? 32e : missing,
    )
end
end
