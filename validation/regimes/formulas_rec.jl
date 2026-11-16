module RemainingBudgets
using SpecialFunctions

"""Finite reference-only budgets; evaluate with at least p+128 BigFloat bits.

Arguments a,b,x,w are independent reference quantities for exactly represented
inputs. This module does not inspect computed errors or production output.
Previously frozen GST/Legendre/Lobatto models remain separate and unchanged.
"""
function result(node, weight, regime; extra_node = missing, extra_weight = missing)
    conclusive = isfinite(node) && node >= 0 && isfinite(weight) && 0 <= weight < BigFloat(1) / 8
    return (;
        node_absolute = node, weight_relative = weight, conclusive, regime,
        extra_node_absolute = extra_node, extra_weight_relative = extra_weight,
    )
end

function mass_error(a, b; recurrence = false)
    if recurrence
        # Direct gamma ratio used by the bounded Jacobi recurrence.
        A, B, S = a + 2, b + 2, a + b + 2
        return 8 * (
            1 + abs((a + b + 1) * log(BigFloat(2))) + abs(loggamma(A)) + abs(loggamma(B)) +
                abs(loggamma(S)) + abs(log(a + 1)) + abs(log(b + 1)) +
                A * abs(digamma(A)) + B * abs(digamma(B)) + S * abs(digamma(S))
        )
    end
    # Retain the previously frozen log-mass coefficient8.
    return 8 * (
        1 + abs((a + b + 1) * log(BigFloat(2))) + abs(loggamma(a + 1)) +
            abs(loggamma(b + 1)) + abs(loggamma(a + b + 2))
    )
end

gamma_error(a) = 8 * (1 + abs(loggamma(a + 1)) + (a + 1) * abs(digamma(a + 1)))

function jacobi_rec(n, a, b, bits, x)
    e = BigFloat(2)^(1 - bits); u = e / 2; nu = n + (a + b + 1) / 2
    X = 8max(abs(x), inv(nu)); s = (1 - x) * (1 + x)
    K = 2abs(a - b + (a + b + 1) * x) / s
    local_error = 16(n + 1) + K * X + 2 / s + 8
    # P'(±1)/P(±1) and positivity of reciprocal endpoint distances.
    D = n * (n + a + b + 1)
    lo = -1 + 2(b + 1) / D; hi = 1 - 2(a + 1) / D
    xmax = max(abs(lo), abs(hi)); smin = min(1 - lo^2, 1 - hi^2)
    Kmax = 2 * (abs(a - b) + abs(a + b + 1) * xmax) / smin
    global_error = 16(n + 1) + Kmax * 8max(xmax, inv(nu)) + 2 / smin + 8
    summation = (n - 1) * u / (1 - (n - 1) * u)
    weight = e * (local_error + global_error + mass_error(a, b; recurrence = true) + 8) + summation
    node = iszero(a - b)&&isodd(n)&&iszero(x) ? zero(e) : e * X
    return result(node, weight, "jacobi_REC")
end

function point_budget(family, n, a, b, bits, i, x, w; format = "Float64", normalize = false, route = :auto)
    precision(BigFloat) >= bits + 128 || error("Budget requires independent reference precision")
    bits >= 24 || error("Assess completed Float32 calculation before floor conversion")
    family = Symbol(family);a, b, x, w = BigFloat(a), BigFloat(b), BigFloat(x), BigFloat(w)
    e = BigFloat(2)^(1 - bits);u = e / 2
    if family in (:chebyshev1, :chebyshev2, :chebyshev3, :chebyshev4)
        kind = findfirst(==(family), (:chebyshev1, :chebyshev2, :chebyshev3, :chebyshev4))
        q = max(abs(x), inv(BigFloat(n)))
        node = kind <= 2&&isodd(n)&&i == (n + 1) ÷ 2 ? zero(e) : 6e * q
        if kind == 1
            weight = 4e
        elseif kind == 2
            weight = 12e
        else
            j = kind == 4 ? n + 1 - i : i
            z = kind == 4 ? -x : x
            weight = j <= div(2n + 1, 8) ? 12e : e * (8 + 6q / (1 + z))
        end
        return result(node, weight, "chebyshev_$(kind)_formula")
    elseif family == :legendre && n <= 5
        node = iszero(x) ? zero(e) : 8e * max(abs(x), inv(BigFloat(n)))
        weight = n <= 2 ? zero(e) : 16e
        return result(node, weight, "legendre_explicit")
    elseif family == :hermite && n <= 1
        return result(zero(e), normalize ? zero(e) : 4e, "hermite_explicit")
    elseif family == :jacobi && n == 1
        return result(4e * abs(x), e * mass_error(a, b), "jacobi_explicit")
    elseif family == :laguerre && n <= 2
        return result((n == 1 ? 2e : 8e) * abs(x), e * (gamma_error(a) + (n == 1 ? 0 : 16)), "laguerre_explicit")
    elseif family == :radau && n <= 2
        if iszero(a)&&iszero(b)
            return result(i == 1 ? zero(e) : 2e * abs(x), zero(e), "radau_standard_explicit")
        end
        return result(i == 1 ? zero(e) : 8e, e * (mass_error(a, b) + (n == 1 ? 0 : 16)), "radau_positive_explicit")
    elseif family == :lobatto && n <= 3
        return result(zero(e), n == 2 ? zero(e) : 4e, "lobatto_explicit")
    elseif family == :jacobi && (route == :rec || (route == :auto && n <= (format == "Float64" ? 24 : 20) && min(a, b) >= -BigFloat(1) / 2 && max(a, b) <= 2))
        return jacobi_rec(n, a, b, bits, x)
    elseif family == :laguerre && route in (:rec, :recompute)
        n <= (route == :rec ? 127 : 1000) || error("Outside bounded recurrence scope")
        #40eps*x is the actual Newton trigger;8eps*x reserves finite
        #recurrence residual/subtraction error on this assessed bounded regime.
        X = 48x
        K = abs(2 - (2a + 1) / x)
        amplitude = 8(n + 1) * (1 + abs(a))
        weight = e * (2amplitude + K * X + gamma_error(a) + 16)
        return result(e * X, weight, route == :rec ? "laguerre_REC" : "laguerre_bounded_recompute")
    end
    error("Use the previously frozen model or an explicit asymptotic context for this route")
end

end
