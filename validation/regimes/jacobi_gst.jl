module NativeJacobiBudget

using SpecialFunctions

export context, point_budget, radau_endpoint_budget

"""
    context(m, a, b, bits; lower, upper, ordinary_lower, ordinary_upper,
            lower_series, upper_series)

Evaluate the frozen native-GST engineering envelope in the caller's reference
precision. Parameters are the exact represented inputs, with Radau's mathematical
beta+1 shift performed in reference precision. Bounds must enclose the indicated
reference roots; this routine does not verify the brackets or run quadrature.
The returned object contains reference quantities only, never production state.
"""
function context(
        m, a, b, bits; lower, upper, ordinary_lower = nothing,
        ordinary_upper = nothing, lower_series = false, upper_series = false
    )
    m >= 2 || throw(ArgumentError("GST budget requires degree at least two"))
    bits >= 24 || throw(ArgumentError("Apply the Float32 floor before this budget"))
    -1 < lower <= upper < 1 || throw(ArgumentError("Root envelope must lie inside (-1,1)"))
    a, b = BigFloat(a), BigFloat(b)
    e = BigFloat(2)^(1 - bits)
    nu = m + (a + b + 1) / 2
    nu > 0 || throw(ArgumentError("Nonpositive oscillation scale"))
    z = min(a, b)
    H0 = max(abs(log((1 - lower) * (1 + lower))), abs(log((1 - upper) * (1 + upper))))
    Hm = max(abs(log1p(-lower)), abs(log1p(-upper)))
    Hp = max(abs(log1p(lower)), abs(log1p(upper)))
    G_global = abs(z) * H0 + abs(a - z) * Hm + abs(b - z) * Hp
    KX_global = BigFloat(0)
    if ordinary_lower !== nothing
        ordinary_upper !== nothing || throw(ArgumentError("Incomplete ordinary-root envelope"))
        -1 < ordinary_lower <= ordinary_upper < 1 || throw(ArgumentError("Invalid ordinary-root envelope"))
        KX_global = 8 * (abs(a) + abs(b)) + 8 * (m / nu) * (
            abs(a) * sqrt((1 + ordinary_upper) / (1 - ordinary_upper)) +
                abs(b) * sqrt((1 - ordinary_lower) / (1 + ordinary_lower))
        )
    end
    if lower_series || upper_series
        KX_global = max(KX_global, 32 * (abs(a) + abs(b)))
    end
    mean_over_e = 16 * (m + 1) + KX_global + 8 * (1 + G_global + log(BigFloat(m)))
    logmass = (a + b + 1) * log(BigFloat(2)) + loggamma(a + 1) +
        loggamma(b + 1) - loggamma(a + b + 2)
    mass_over_e = 8 * (
        1 + abs((a + b + 1) * log(BigFloat(2))) +
            abs(loggamma(a + 1)) + abs(loggamma(b + 1)) + abs(loggamma(a + b + 2))
    )
    summation_denominator = 1 - (m - 1) * e / 2
    sum_error = summation_denominator > 0 ?
        (m - 1) * e / (2 * summation_denominator) + 8e : BigFloat(Inf)
    return (;
        m, a, b, bits, e, nu, z, lower_series, upper_series,
        logmass, mean_over_e, mass_over_e, sum_error, G_global, KX_global,
        summation_denominator,
    )
end

"""Return pointwise frozen allowances; `v` is the reference Jacobi weight."""
function point_budget(c, index, x, v; radau = false)
    1 <= index <= c.m || throw(ArgumentError("Invalid Jacobi root index"))
    -1 < x < 1 && v > 0 || throw(ArgumentError("Invalid reference node or weight"))
    left, right = 1 - x, 1 + x
    series = (index == 1 && c.lower_series) || (index == c.m && c.upper_series)
    C = if index == 1 && c.lower_series
        32right
    elseif index == c.m && c.upper_series
        32left
    else
        8 * min(left, right, max(abs(x), 1 / c.nu)) +
            8 * (c.m / c.nu) * sqrt(left * right)
    end
    K_J = abs(-c.a / left + c.b / right)
    K = radau ? abs(-c.a / left + (c.b - 1) / right) : K_J
    G = abs(c.z * log(left * right)) + abs(c.a - c.z) * abs(log(left)) +
        abs(c.b - c.z) * abs(log(right))
    logp = log(v) - c.logmass
    # The references can produce a tiny positive log(p) by final rounding only.
    logp <= 16eps(BigFloat) || throw(ArgumentError("A normalized reference weight exceeds one"))
    local_over_e = 16 * (c.m + 1) + K * C + 8 * (1 + G + abs(logp))
    node_absolute = c.e * (C + abs(x))
    weight_relative = c.e * (c.mass_over_e + local_over_e + c.mean_over_e) +
        c.sum_error + (radau ? 8c.e : zero(c.e))
    conclusive = c.summation_denominator > 0 && isfinite(weight_relative) && weight_relative < BigFloat(1) / 8
    return (;
        node_absolute, weight_relative, conclusive, series, C, K, K_J,
        G, logp, local_over_e, mean_over_e = c.mean_over_e,
        mass_over_e = c.mass_over_e, sum_error = c.sum_error,
        KX_global = c.KX_global, G_global = c.G_global,
    )
end

"""Separate fixed-left-endpoint Radau relative-weight allowance."""
function radau_endpoint_budget(n, a, b, bits)
    e = BigFloat(2)^(1 - bits)
    allowance = 8e * (
        n + 1 + abs((a + b + 1) * log(BigFloat(2))) +
            abs(loggamma(b + 1)) + abs(loggamma(b + 2)) + abs(loggamma(n + a)) +
            abs(loggamma(BigFloat(n))) + abs(loggamma(n + b + 1)) + abs(loggamma(n + a + b + 1))
    )
    return (;
        node_absolute = BigFloat(0), weight_relative = allowance,
        conclusive = isfinite(allowance) && allowance < BigFloat(1) / 8,
    )
end

end
