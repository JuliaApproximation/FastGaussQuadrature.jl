module SmallHermiteBudget

# Reference-only model: run in the independently chosen high-precision scope.
# Inputs are physical Hermite reference nodes and weights, before normalization.
function context(n::Integer, nodes, weights, bits::Integer)
    2 <= n <= 20 || error("Small Hermite recurrence budget is restricted to n2:20")
    bits >= 24 || error("Assess the Float32 working rule separately below 24 bits")
    precision(BigFloat) >= bits + 128 || error("Insufficient reference-model precision")
    length(nodes) == length(weights) == n || error("Incomplete reference rule")
    all(isfinite, nodes) && all(w -> isfinite(w) && w > 0, weights) ||
        error("Invalid independent Hermite references")
    x, w = BigFloat.(nodes), BigFloat.(weights)
    e = BigFloat(2)^(1 - bits)
    u = e / 2
    X = [4max(one(t), abs(t)) for t in x]
    physical = [4n + 8 + 4abs(x[i]) * X[i] for i in eachindex(x)]
    scaled = [4n + 8 + 2abs(x[i]) * X[i] + x[i]^2 for i in eachindex(x)]
    mean = sum(w .* physical) / sqrt(BigFloat(π))
    denominator = 1 - (n - 1) * u
    gamma = (n - 1) * u / denominator
    return (; n, bits, x, w, e, u, X, physical, scaled, mean, gamma, denominator)
end

function allowance(c, i::Integer; normalize::Bool = false, scaled::Bool = false)
    normalize && scaled && error("Normalized and scaled outputs are separate views")
    x, e = c.x[i], c.e
    node = iszero(x) ? zero(x) : e * c.X[i]
    local_weight = scaled ? c.scaled[i] : c.physical[i]
    weight = e * (local_weight + c.mean + 8) + c.gamma
    if normalize
        delta = 2c.u + c.u^2
        node = sqrt(BigFloat(2)) * (node + delta * (abs(x) + node))
        weight += 2e
    end
    conclusive = c.denominator > 0 && isfinite(node) && node >= 0 &&
        isfinite(weight) && 0 < weight < BigFloat(1) / 8
    regime = normalize ? "small_recurrence_normalized" :
        scaled ? "small_recurrence_scaled" : "small_recurrence_physical"
    return (;
        node_absolute = node, weight_relative = weight, conclusive, regime,
        exact_zero_required = iszero(x),
    )
end

end
