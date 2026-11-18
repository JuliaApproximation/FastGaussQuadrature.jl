module LaguerreAsymptoticBudget
using SpecialFunctions
include("laguerre_coefficients.jl")
include("formulas_rec.jl")

# Record the algorithm selected by the existing Float64 helper. These values
# select the formula/depth; no measured quadrature error enters the allowance.
depth(v) = begin
    k = 1
    while k < length(v)&&abs(v[k + 1]) < abs(v[k])
        k += 1
    end
    k
end
function trace(mod, n, a = 0.0; hybrid = false)
    d = 1.0 / (4n + 2a + 2);out = NamedTuple[]
    count = hybrid ? (n >= 1000 ? 12 : 2) : n
    kb = max(ceil(Int, sqrt(n)), 7)
    roots = abs(a) == 0.5 ? Float64[] : mod.approx_besselroots(a, kb)
    state = :bessel
    for i in 1:count
        candidates = NamedTuple[]
        if state == :bessel
            j = a == 0.5 ? pi * i : a == -0.5 ? pi * (i - 0.5) : i < kb ? roots[i] : mod.McMahon(a, i)
            x, w, delta = mod.gausslaguerre_asy_bessel(n, a, j, d, -1, i <= 12)
            xs, ws = iszero(a) ? mod.gl_bessel(j, d) : mod.gl_bessel(j, d, a)
            push!(candidates, (; chart = :bessel, x, w, delta, j, t = NaN, ak = NaN, qx = depth(xs), qw = depth(ws)))
        end
        if !hybrid&&state != :airy
            x, w, delta = mod.gausslaguerre_asy_bulk(n, a, i, d, -1)
            t = mod.gl_bulk_solve_t(n, i, d, a);xs, ws = iszero(a) ? mod.gl_bulk(t, d) : mod.gl_bulk(t, d, a)
            push!(candidates, (; chart = :bulk, x, w, delta, j = NaN, t, ak = NaN, qx = depth(xs), qw = depth(ws)))
        end
        if !hybrid&&(state == :airy||(state == :bulk&&i >= floor(Int, 0.9n)))
            x, w, delta = mod.gausslaguerre_asy_airy(n, a, i, d, -1);ak = mod.compute_airy_root(n, i)
            xs = mod.gl_airy(ak, d, a)
            push!(candidates, (; chart = :airy, x, w, delta, j = NaN, t = NaN, ak, qx = depth(xs), qw = 0))
        end
        chosen = candidates[argmin(getfield.(candidates, :delta))]
        chosen.chart == :bulk&&(state = :bulk);chosen.chart == :airy&&(state = :airy)
        recomputed = false
        if !hybrid&&n <= 1000&&chosen.delta > 1.0e-13
            xx, _ = mod.gl_rec_newton(chosen.x, n, a)
            recomputed = abs(xx - chosen.x) < 100chosen.delta
        end
        push!(out, merge(chosen, (; recomputed)))
    end
    return out
end

function bulk_phase(n, i, d, a)
    # Monotone bisection of theta+sin(theta)=pi*(4i+2a-1)*d.
    lo, hi = zero(d), BigFloat(pi);target = BigFloat(pi) * (4i + 2a - 1) * d
    for _ in 1:(precision(BigFloat) + 8)
        mid = (lo + hi) / 2
        mid + sin(mid) < target ? (lo = mid) : (hi = mid)
    end
    return sin((lo + hi) / 4)^2
end
function bessel_root(i, a, seed)
    abs(a) == BigFloat(0.5)&&return BigFloat(pi) * (a > 0 ? i : i - BigFloat(0.5))
    j = BigFloat(seed)
    for _ in 1:12
        value = besselj(a, j);step = value / (besselj(a - 1, j) - a / j * value);j -= step
        abs(step) <= 8eps(BigFloat) * abs(j)&&return j
    end
    error("Bessel constant reference did not converge")
end
function airy_root(seed)
    a = BigFloat(seed)
    for _ in 1:12
        # MPFR supplies Ai but this SpecialFunctions version has no BigFloat
        # Ai' method. A reference-only centered difference retains about2p/3
        # derivative bits; Newton is stopped by the full-precision Ai residual.
        h = BigFloat(2)^(-div(precision(BigFloat), 3)) * max(one(a), abs(a))
        derivative = (airyai(a + h) - airyai(a - h)) / (2h)
        step = airyai(a) / derivative;a -= step
        abs(step) <= 8eps(BigFloat) * abs(a)&&return a
    end
    error("Airy constant reference did not converge")
end
tail(v, q) = abs(v[min(q + 1, length(v))])

function point_budget(n, a, bits, i, x, w, selected; minimum_normal = BigFloat(floatmin(Float64)))
    bits == 53&&n >= 128&&a in (0, BigFloat(0.5), BigFloat(-0.5))||error("Outside finite Float64 Laguerre expansion scope")
    precision(BigFloat) >= bits + 128||error("Reference precision required")
    a, x, w = BigFloat(a), BigFloat(x), BigFloat(w);e = BigFloat(2)^(1 - bits);d = inv(BigFloat(4n) + 2a + 2)
    if selected.recomputed
        return RemainingBudgets.point_budget(:laguerre, n, a, 0, bits, i, x, w; route = :recompute)
    end
    chart = selected.chart;constant_node = zero(e);constant_weight = zero(e)
    if chart == :bessel
        j = bessel_root(i, a, selected.j)
        xs, ws = iszero(a) ? gl_bessel(j, d) : gl_bessel(j, d, a)
        scale = j * j * d
        tx = 8scale * tail(xs, selected.qx)
        tw = 8tail(ws, selected.qw) / abs(1 + sum(ws[1:selected.qw]))
        dj = abs(BigFloat(selected.j) - j)
        # Root-dependent prefactor: d(j^2*d)/dj=2jd; at a Bessel zero,
        # J_(a-1)'/J_(a-1)=(a-1)/j. Expansion factors are a small correction.
        constant_node = 4abs(j * d) * dj
        constant_weight = 4abs(a - 1) * dj / abs(j)
        roundnode = 8e * x + 32e * scale * sum(abs, xs)
        roundweight = 8e * (1 + x + abs(a * log(x))) + 32e * sum(abs, ws)
    elseif chart == :bulk
        t = bulk_phase(n, i, d, a);xs, ws = iszero(a) ? gl_bulk(t, d) : gl_bulk(t, d, a)
        tx = 8tail(xs, selected.qx);tw = 8tail(ws, selected.qw) / abs(1 + sum(ws[1:selected.qw]))
        roundnode = 8e * x + 32e * sum(abs, xs)
        roundweight = 8e * (1 + x + abs(a * log(x)) + inv(1 - t)) + 32e * sum(abs, ws)
    else
        ak = airy_root(selected.ak);xs = gl_airy(ak, d, a)
        tx = 8tail(xs, selected.qx)
        # Only a leading Airy weight is supplied. Its first missing relative
        # correction has d^(2/3) scale; this is deliberately not called native
        # relative accuracy. Public normal weights use bounded recomputation.
        tw = 8cbrt(d)^2 * (1 + abs(ak)^2)
        constant_node = 4cbrt(4 / d) * abs(BigFloat(selected.ak) - ak)
        roundnode = 8e * x + 32e * sum(abs, xs[2:end])
        roundweight = 8e * (1 + x + abs((a + BigFloat(1) / 3) * log(x)))
    end
    node = roundnode + tx + constant_node
    weight = roundweight + tw + constant_weight + abs(a / x - 1) * node
    conclusive = all(isfinite, (node, weight))&&node >= 0&&weight >= 0&&(w < minimum_normal||weight < BigFloat(1) / 8)
    return (;
        node_absolute = node, weight_relative = weight, conclusive, regime = "laguerre_ASY_$(chart)",
        node_tail = tx, weight_tail = tw, constant_node, constant_weight,
    )
end
end
