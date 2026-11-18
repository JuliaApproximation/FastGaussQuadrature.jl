module HermiteGSTBudget
function context(n, bits, extreme)
    n > 20 && bits >= 24 || error("GST model scope")
    e = BigFloat(2)^(1 - bits); u = e / 2; λ = BigFloat(2n + 1)
    R = abs(BigFloat(extreme)); amin = λ - R^2; amin > 0 || error("Root beyond turning point")
    m = n ÷ 2
    # Sturm comparison: outward root count k <= 1+sqrt(λ)*abs(x)/pi.
    kmean = 1 + sqrt(λ / 2) / BigFloat(π)
    meanL = 16 * (kmean + 1) + 8 + 8 + 16 * (sqrt(λ) / (2BigFloat(π)) + inv(sqrt(BigFloat(2)))) / sqrt(amin) + BigFloat(0.5)
    dmax = e * (8R + 8m / sqrt(amin))
    meanQ = 2 * (λ + 2R * dmax + dmax^2) * dmax^2
    ρ = e * meanL + meanQ
    sumden = 1 - (n - 1) * u; γ = (n - 1) * u / sumden
    η = ρ + γ * (1 + ρ)
    return (; n, bits, e, u, λ, R, amin, meanL, meanQ, η, sumden)
end
function allowance(c, i, x; view = :physical)
    x = BigFloat(x);a = c.λ - x^2; a > 0 || error("Invalid normal-form coefficient")
    k = cld(abs(2i - c.n - 1), 2)
    X = 8abs(x) + 8k / sqrt(a); d = c.e * X
    Q = 2 * (a + 2abs(x) * d + d^2) * d^2
    L = 16 * (k + 1) + 8 + (view == :scaled ? zero(x) : 2abs(x) * X + x^2)
    raw = c.e * L + Q
    weight = (1 + raw) * (1 + 4c.e) / (1 - c.η) - 1
    node = iszero(x) ? zero(x) : d
    if view == :normalized
        δ = 2c.u + c.u^2
        node = sqrt(BigFloat(2)) * (node + δ * (abs(x) + node))
        weight += 2c.e
    end
    phase = d * sqrt(a + 2abs(x) * d + d^2)
    conclusive = c.sumden > 0&&c.η < BigFloat(1) / 8&&phase < BigFloat(1) / 8&&isfinite(weight)&&0 < weight < BigFloat(1) / 8
    return (;
        node_absolute = node, weight_relative = weight, conclusive, steps = k, coordinate = X,
        stationary = Q, local_coefficient = L, normalization = c.η, phase,
    )
end
end
