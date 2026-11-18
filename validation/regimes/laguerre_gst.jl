module NativeLaguerreBudget
using SpecialFunctions

# Frozen equations in budget-v1.md/json. This module uses reference arithmetic,
# never production intermediate values. It does not run quadrature.

mass_reserve(s) = 8 * (2 + abs(loggamma(s)) + s * abs(log(s)))
coordinate(c, x) = 16x + 16c.n * sqrt(x / c.nu)
local_reserve(c, x, w) = 16 * (c.n + 1) + abs((c.a + 1 / BigFloat(2)) / x - 1) * coordinate(c, x) +
    8 * (1 + abs(c.a + 1 / BigFloat(2)) * abs(log(x / c.pivot)) + abs(x - c.pivot) + abs(log(w) - c.logmass))

function startup(n, a)
    abs(a) > 1 / BigFloat(2) && return sqrt(a^2 - 1 / BigFloat(4))
    upper = (2n^2 + n * (a - 1) + 2 * (a + 1) + 2 * (n - 1) * sqrt(n^2 + (n + 2) * (a + 1))) / (n + 2)
    return (a + 1) * (n * (a + 5) + 2 * (a - 1)) / ((n + 2) * upper)
end

"""
Finite 1F1(-degree; b; x), differentiated with respect to `variable` when x is
proportional to that variable. The positive magnitude ratios decrease with k.
An omitted tail is bounded geometrically for the value and both derivatives;
absolute sums include the remainder bound. Failure to resolve the tail raises
an error for the assessor to retain as an inconclusive model evaluation.
"""
function finite_series(degree, b, x, variable = x)
    degree >= 0 && b > 0 && x > 0 && variable > 0 || error("Invalid finite-series arguments")
    term = one(x)
    value, derivative, second = one(x), zero(x), zero(x)
    absolute, absolute_derivative = one(x), zero(x)
    maxterms = min(degree, max(128, 4precision(BigFloat)))
    for k in 1:maxterms
        term *= ((k - 1 - degree) / BigFloat(k)) * (x / (b + k - 1))
        first_term, second_term = k * term / variable, k * (k - 1) * term / variable^2
        value += term
        derivative += first_term
        second += second_term
        absolute += abs(term)
        absolute_derivative += abs(first_term)
        k == degree && return (; value, derivative, second, absolute, absolute_derivative, terms = k, tail = zero(x))
        ratio = (degree - k) * x / ((k + 1) * (b + k))
        if ratio < 1 / BigFloat(2)
            nextterm = abs(term) * ratio
            denominator = 1 - ratio
            tail = nextterm / denominator
            tail_first = nextterm / variable * ((k + 1) / denominator + ratio / denominator^2)
            tail_second = nextterm / variable^2 * (
                k * (k + 1) / denominator + (2k + 1) * ratio / denominator^2 +
                    ratio * (1 + ratio) / denominator^3
            )
            if tail <= eps(BigFloat) * absolute &&
                    tail_first <= eps(BigFloat) * max(one(x), abs(derivative)) &&
                    tail_second <= eps(BigFloat) * max(one(x), abs(second))
                return (;
                    value, derivative, second, absolute = absolute + tail,
                    absolute_derivative = absolute_derivative + tail_first, terms = k, tail,
                )
            end
        end
    end
    degree == 0 && return (; value, derivative, second, absolute, absolute_derivative, terms = 0, tail = zero(x))
    error("Reference finite series did not resolve its remainder")
end

function context(
        n::Integer, a, bits; reduced = false, minimum_normal, float64 = false,
        x1 = nothing, w1 = nothing, x2 = nothing, w2 = nothing
    )
    n >= 3 && bits >= 24 || error("Native GST budget requires n>=3 and p>=24")
    a = BigFloat(a)
    delta, degree, e = a + 1, BigFloat(n), BigFloat(2)^(1 - bits)
    delta > 0 || error("Nonintegrable Laguerre input")
    nu, pivot = 4degree + 2a + 2, startup(degree, a)
    logmass = loggamma(delta)
    mass = exp(logmass)
    all(isfinite, (nu, pivot, logmass, mass)) && min(nu, pivot, mass) > 0 || error("Unavailable finite reference scale")
    base = (;
        n = degree, order = n, a, bits, e, delta, nu, pivot, logmass, mass,
        endpoint = delta <= 1 / BigFloat(8),
    )
    if !base.endpoint
        lower = delta / degree
        upper = 2degree + a - 1 + 2sqrt((degree - 1) * (degree - 1 + a))
        magnitude = abs(a + 1 / BigFloat(2))
        Kglobal = 16 * (magnitude + upper) + 16degree / sqrt(nu) * (magnitude / sqrt(lower) + sqrt(upper))
        Gglobal = magnitude * max(abs(log(lower / pivot)), abs(log(upper / pivot))) + max(abs(lower - pivot), abs(upper - pivot))
        mean = 16 * (degree + 1) + Kglobal + 8 * (1 + Gglobal + log(degree))
        denominator = 1 - (degree - 1) * e / 2
        tau = reduced ? min(one(e), degree * 10 * BigFloat(minimum_normal) / mass) : zero(e)
        sum_error = denominator > 0 ? (degree - 1) * e / (2denominator) : BigFloat(Inf)
        tail_error = tau < 1 ? tau / (1 - tau) : BigFloat(Inf)
        return merge(
            base, (;
                lower, upper, mean, denominator, tau, sum_error, tail_error,
                mass_error = mass_reserve(delta), context_conclusive = denominator > 0 && tau < 1,
            )
        )
    end
    all(v -> v !== nothing && isfinite(v) && v > 0, (x1, w1, x2, w2)) || error("Endpoint model needs two independently valid reference pairs")
    x1, w1, x2, w2 = BigFloat.((x1, w1, x2, w2))
    t = degree * x1 / delta
    first_series = finite_series(n, delta, x1, t)
    Fprime, Fsecond = first_series.derivative, first_series.second
    !iszero(Fprime) || error("Unresolved first-root derivative")
    Ct = 8 * (degree + 1) * first_series.absolute / abs(Fprime) + 8abs(t)
    ED = 8 * (degree + 1) * first_series.absolute_derivative
    logp = sum(k -> log1p(delta / k), 1:(n - 1); init = zero(e))
    LP = float64 ? 8 * (1 + abs(logbeta(degree, delta)) + abs(log(delta)) + abs(logp)) :
        8 * (degree + abs(logp))
    first_node = e * (delta * Ct / degree + 4x1)
    first_relative = e * (mass_reserve(delta) + LP + 2ED / abs(Fprime) + Ct * abs(1 / t + 2Fsecond / Fprime) + 8)
    anchor_series = finite_series(n - 1, delta + 1, x2)
    P, Pprime = anchor_series.value, anchor_series.derivative
    !iszero(P) || error("Unresolved second-root anchor series")
    chi2 = anchor_series.absolute / abs(P)
    anchor = mass_reserve(delta + 1) + LP + 16 * (degree + 1) * chi2 +
        abs(1 / x2 + 2Pprime / P) * (coordinate(base, x2) + x2) + 8
    anchor_local = local_reserve(base, x2, w2)
    return merge(
        base, (;
            first_node, first_relative, anchor, anchor_local, chi2, Ct, ED, LP, logp,
            first_series_terms = first_series.terms, anchor_series_terms = anchor_series.terms,
            context_conclusive = all(isfinite, (first_node, first_relative, anchor, anchor_local)),
        )
    )
end

function point_budget(c, index, x, weight)
    1 <= index <= c.order && x > 0 && weight > 0 || error("Invalid independent reference pair")
    if c.endpoint && index == 1
        node, relative = c.first_node, c.first_relative
        model = "GST_endpoint_series"
    else
        node = c.e * (coordinate(c, x) + x)
        local_error = local_reserve(c, x, weight)
        if c.endpoint
            relative = c.e * (c.anchor + local_error + c.anchor_local + 8)
            model = "GST_second_root_anchor"
        else
            relative = c.e * (c.mass_error + local_error + c.mean + 16) + c.sum_error + c.tail_error
            model = "GST_mass_normalization"
        end
    end
    conclusive = c.context_conclusive && all(isfinite, (node, relative)) && relative < 1 / BigFloat(8)
    return (; node_absolute = node, weight_relative = relative, conclusive, model)
end
end
