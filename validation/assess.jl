module QuadratureAssessment

using Dates, SHA, Quadmath, SpecialFunctions, DelimitedFiles, FastGaussQuadrature
include("regimes/jacobi_gst.jl")
include("regimes/laguerre_gst.jl")
include("regimes/hermite_rec.jl")
include("regimes/hermite_gst.jl")
include("regimes/formulas_rec.jl")
include("regimes/jacobi_asy_budget.jl")
include("regimes/laguerre_asy_budget.jl")
include(joinpath(@__DIR__, "..", "test", "quadrature_reference.jl"))

# This adapter composes previously frozen models without using measured errors.
# Half-parameter Laguerre supplies selected formula/depth and native Bessel-seed
# metadata; its error allowances use independently evaluated reference formulas.
function records(path)
    data, header = readdlm(path, ',', String; header = true, quotes = true)
    names = vec(header)
    length(unique(names)) == length(names) || error("Duplicate CSV header")
    size(data, 2) == length(names) || error("CSV field count mismatch")
    return [Dict(names .=> vec(data[i, :])) for i in axes(data, 1)]
end
csvfield(v) = "\"" * replace(string(v), "\"" => "\"\"") * "\""
csvrow(io, row) = println(io, join(csvfield.(row), ','))
key(r) = (r["implementation"], r["format"], r["case_id"])
num(r, k) = parse(BigFloat, r[k])
int(r, k) = parse(Int, r[k])
flag(r, k) = r[k] == "true"
function format_type(s)
    s == "Float32" && return Float32
    s == "Float64" && return Float64
    s == "Float128" && return Float128
    occursin(r"^BigFloat[0-9]+$", s) && return BigFloat
    error("Unsupported assessment format: $s")
end
function represented(r, k)
    T = format_type(r["format"])
    value = setprecision(BigFloat, int(r, "output_bits")) do
        # cases.jl stores Float64 parameter literals, which the measurement
        # converts to T. Preserve that binary value before reference lifting.
        # Output strings instead describe a value already rounded in T.
        T(k in ("a", "b") ? parse(Float64, r[k]) : parse(BigFloat, r[k]))
    end
    return BigFloat(value)
end
ratio(e, b) = b > 0 ? e / b : iszero(e) ? zero(e) : BigFloat(Inf)
result(node, weight, regime; conclusive = true) = (;
    node_absolute = node, weight_relative = weight, regime,
    conclusive = conclusive&&isfinite(node)&&node >= 0&&isfinite(weight)&&0 <= weight < BigFloat(1) / 8,
)

function jacobi_route(n, a, b, f)
    f == "Float64"&&n >= 768&&min(a, b) >= BigFloat(-0.99)&&max(a, b) < 5&&return :ASY
    n <= (f == "Float64" ? 24 : 20)&&min(a, b) >= -BigFloat(1) / 2&&max(a, b) <= 2&&return :REC
    return :GST
end
function route(r)
    n, p = int(r, "n"), int(r, "output_bits");a, b = represented(r, "a"), represented(r, "b")
    f, family = r["format"], Symbol(r["family"])
    n > 0 && p >= 24 || error("Assess empty rules and Float32-floor conversion separately")
    T = format_type(f)
    (T === BigFloat ? parse(Int, f[9:end]) == p : precision(T) == p) || error("Format/precision mismatch")
    family in (:jacobi, :laguerre, :radau_jacobi) || (iszero(a) && iszero(b)) ||
        error("Unused parameters must be zero for an unparameterized rule")
    startswith(string(family), "chebyshev")&&return :formula
    family == :legendre&&return n <= 5 ? :formula : n <= 60 ? :legendre_REC : p <= 53 ? :legendre_ASY : :GST
    family == :hermite&&return n <= 1 ? :formula : n <= 20 ? :hermite_REC : :hermite_GST
    if family == :jacobi
        ((iszero(a) && iszero(b)) || (abs(a) == abs(b) == BigFloat(1) / 2)) &&
            error("Assess Jacobi's Legendre/Chebyshev special pairs with their dedicated family")
        return n == 1 ? :formula : jacobi_route(n, a, b, f)
    end
    if family == :laguerre
        n <= 2&&return :formula
        return :laguerre_GST
    elseif family == :lobatto
        n <= 3&&return :formula
        f == "Float64"&&n >= 770&&return :ASY
        return n <= (f == "Float64" ? 26 : 22) ? :REC : :GST
    elseif family == :radau
        n <= 2&&return :formula
        f == "Float64"&&n > 101&&return :ASY
        (!(f in ("Float32", "Float64"))||(f != "Float64"&&n > 101))&&return :GST
        return jacobi_route(n - 1, zero(a), one(b), f) == :REC ? :REC : :GST_postdivide
    elseif family == :radau_jacobi
        n <= 2&&return :formula
        # Non-Float64 generalized rules use transformed GST. For Float64,
        # restrict this adapter to the declared matrix rather than assigning
        # its small postdivision and asymptotic routes a GST budget.
        f != "Float64" && return :GST
        large = a == 10 && b == 20 && n in (20, 21, 100, 257)
        near = n == 100 && a == BigFloat(-0.99) && b == BigFloat(0.2)
        (large || near) && return :GST
        error("Generalized Radau route outside the assessed matrix")
    end
    error("Unrecognized family")
end

function postdivide_allowance(rule, row, base, x, weight, regime)
    n, bits = int(rule, "n"), int(rule, "output_bits")
    rule["family"] == "radau" && rule["format"] in ("Float32", "Float64") && n <= 101 ||
        error("Postdivision adapter is restricted to standard Radau through n=101")
    iszero(represented(rule, "a")) && iszero(represented(rule, "b")) || error("Nonstandard postdivision parameters")
    T = format_type(rule["format"])
    u = BigFloat(2)^(-bits)
    eta = base.node_absolute / (1 + x)
    normal_min, normal_max = BigFloat(floatmin(T)), BigFloat(floatmax(T))
    xa, wa = represented(row, "node_actual"), represented(row, "weight_actual")
    actual_denominator = BigFloat(
        setprecision(BigFloat, bits) do
            one(T) + T(xa)
        end
    )
    reference_numerator = (1 + x) * weight
    implied_lower = wa * actual_denominator / (1 + u)
    implied_upper = wa * actual_denominator / (1 - u)
    quotient_lower = weight * (1 - base.weight_relative) / ((1 + eta) * (1 + u))
    quotient_upper = eta < 1 ? weight * (1 + base.weight_relative) / ((1 - eta) * (1 - u)) : BigFloat(Inf)
    normal_regime = normal_min <= quotient_lower <= quotient_upper <= normal_max &&
        normal_min <= weight <= normal_max && normal_min <= wa <= normal_max &&
        normal_min <= reference_numerator <= normal_max &&
        normal_min <= 1 + x <= normal_max && normal_min <= actual_denominator <= normal_max &&
        normal_min <= implied_lower <= implied_upper <= normal_max &&
        eta < 1 && (1 + x) * (1 - eta) * (1 - u) >= normal_min
    relative = eta < 1 ? (1 + base.weight_relative) * (1 + u) / ((1 - eta) * (1 - u)) - 1 : BigFloat(Inf)
    return result(base.node_absolute, relative, regime; conclusive = base.conclusive && normal_regime)
end

function root_bracket(r, m, a, b, index)
    x = num(r, "node_reference")
    radius = min((1 - abs(x)) / 4, BigFloat(2)^(-int(r, "reference_bits") + 32))
    lo, hi = x - radius, x + radius
    QuadratureReference.roots_below(:jacobi, m, a, b, lo) == index - 1||error("Lower Sturm envelope")
    QuadratureReference.roots_below(:jacobi, m, a, b, hi) == index||error("Upper Sturm envelope")
    return (; lo, hi)
end
function context(rule, rows, method)
    n, p = int(rule, "n"), int(rule, "output_bits");a, b = represented(rule, "a"), represented(rule, "b")
    f = Symbol(rule["family"]);e = BigFloat(2)^(1 - p)
    byindex = Dict(int(r, "index") => r for r in rows)
    if method == :hermite_REC
        scale = flag(rule, "normalize") ? sqrt(BigFloat(2)) : one(e)
        mass = flag(rule, "normalize") ? sqrt(BigFloat(π)) : one(e)
        return SmallHermiteBudget.context(
            n, [num(byindex[i], "node_reference") / scale for i in 1:n],
            [num(byindex[i], "weight_reference") * mass for i in 1:n], p
        )
    elseif method == :hermite_GST
        scale = flag(rule, "normalize") ? sqrt(BigFloat(2)) : one(e)
        return HermiteGSTBudget.context(n, p, abs(num(byindex[n], "node_reference")) / scale)
    elseif method == :laguerre_GST
        T = format_type(rule["format"])
        gst = NativeLaguerreBudget.context(
            n, a, p; reduced = flag(rule, "reduced"),
            minimum_normal = BigFloat(floatmin(T)), float64 = T === Float64,
            x1 = num(byindex[1], "node_reference"), w1 = num(byindex[1], "weight_reference"),
            x2 = num(byindex[2], "node_reference"), w2 = num(byindex[2], "weight_reference")
        )
        # The first two/twelve half-parameter outputs replace GST values with
        # the existing Bessel expansion. Assess that actual algorithm instead.
        hybrid = T === Float64 && n >= 128 && abs(a) == BigFloat(1) / 2
        trace = hybrid ? LaguerreAsymptoticBudget.trace(FastGaussQuadrature, n, Float64(a); hybrid = true) : nothing
        return (; gst, trace)
    elseif method in (:laguerre_ASY, :laguerre_recompute)
        return LaguerreAsymptoticBudget.trace(FastGaussQuadrature, n, Float64(a))
    elseif method in (:GST, :GST_postdivide, :REC)
        offset = f in (:radau, :radau_jacobi, :lobatto) ? 1 : 0
        m = f == :lobatto ? n - 2 : n - offset
        a, b = f == :lobatto ? (one(e), one(e)) : f in (:radau, :radau_jacobi) ? (a, b + 1) : (a, b)
        lower = root_bracket(byindex[1 + offset], m, a, b, 1).lo
        upper = root_bracket(byindex[m + offset], m, a, b, m).hi
        if method == :REC
            rmax = max(abs(lower), abs(upper));smin = 1 - rmax^2
            denominator = 1 - (m - 1) * e / 2
            return (;
                m, a, b, e, rmax, smin, qmax = max(rmax, 1 / (m + BigFloat(1.5))),
                kappamax = 2rmax / smin, denominator, gamma_sum = (m - 1) * e / (2denominator),
            )
        end
        threshold = setprecision(BigFloat, p) do;
            BigFloat(-0.9)
        end
        ls, us = b < threshold, a < threshold
        ol = root_bracket(byindex[1 + Int(ls) + offset], m, a, b, 1 + Int(ls)).lo
        ou = root_bracket(byindex[m - Int(us) + offset], m, a, b, m - Int(us)).hi
        return NativeJacobiBudget.context(
            m, a, b, p; lower, upper, ordinary_lower = ol, ordinary_upper = ou,
            lower_series = ls, upper_series = us
        )
    end
    return nothing
end

function allowance(rule, r, c, method)
    n, p, i = int(rule, "n"), int(rule, "output_bits"), int(r, "index")
    a, b = represented(rule, "a"), represented(rule, "b");x, w = num(r, "node_reference"), num(r, "weight_reference")
    f = Symbol(rule["family"]);e = BigFloat(2)^(1 - p)
    if f == :lobatto && i in (1, n)
        return result(zero(e), 4e, "Lobatto_endpoint")
    elseif f in (:radau, :radau_jacobi)&&i == 1&&n > 2
        z = NativeJacobiBudget.radau_endpoint_budget(n, a, b, p)
        return merge(z, (; regime = "Radau_endpoint"))
    elseif method == :formula
        aliases = Dict(:chebyshevt => :chebyshev1, :chebyshevu => :chebyshev2, :chebyshevv => :chebyshev3, :chebyshevw => :chebyshev4, :radau_jacobi => :radau)
        return RemainingBudgets.point_budget(get(aliases, f, f), n, a, b, p, i, x, w; format = rule["format"], normalize = flag(rule, "normalize"))
    elseif method == :legendre_REC
        q = max(abs(x), inv(BigFloat(n)));kappa = 2abs(x) / (1 - x^2)
        return result(4e * q, e * (8n + 8 + 4kappa * q), "Legendre_REC")
    elseif method == :legendre_ASY
        return result(4e * max(abs(x), inv(BigFloat(n))), 8e, "Legendre_ASY")
    elseif method == :hermite_REC
        return SmallHermiteBudget.allowance(c, i; normalize = flag(rule, "normalize"))
    elseif method == :hermite_GST
        normalized = flag(rule, "normalize")
        z = HermiteGSTBudget.allowance(c, i, normalized ? x / sqrt(BigFloat(2)) : x; view = normalized ? :normalized : :physical)
        return merge(z, (; regime = "Hermite_GST"))
    elseif method == :laguerre_REC
        return RemainingBudgets.point_budget(:laguerre, n, a, b, p, i, x, w; route = :rec)
    elseif method in (:laguerre_ASY, :laguerre_recompute)
        return LaguerreAsymptoticBudget.point_budget(n, a, p, i, x, w, c[i])
    elseif method == :laguerre_GST
        if c.trace !== nothing && i <= length(c.trace)
            return LaguerreAsymptoticBudget.point_budget(n, a, p, i, x, w, c.trace[i])
        end
        return merge(NativeLaguerreBudget.point_budget(c.gst, i, x, w), (; regime = "Laguerre_GST"))
    elseif method == :ASY
        offset = f in (:radau, :radau_jacobi, :lobatto) ? 1 : 0
        m = f == :lobatto ? n - 2 : n - offset
        aa, bb = f == :lobatto ? (one(e), one(e)) : f in (:radau, :radau_jacobi) ? (a, b + 1) : (a, b)
        transform = f == :lobatto ? :lobatto : f in (:radau, :radau_jacobi) ? :radau : :none
        consumer = f in (:radau, :radau_jacobi) ? :radau : f == :lobatto ? :lobatto : :jacobi
        return JacobiAsymptoticBudget.point_budget(m, aa, bb, p, i - offset, x, w; transform, consumer)
    elseif method == :REC
        if f == :lobatto
            s = 1 - x^2;kappa = 2abs(x) / s;q = max(abs(x), 1 / (c.m + BigFloat(1.5)));X = 8q;A = 8(c.m + 1)
            raw = 2A + 3c.kappamax * (8c.qmax) + 2 / c.smin + 8
            relative = e * (2A + 2kappa * X + 2 / s + 8 + raw + 16) + c.gamma_sum
            return result(e * X, relative, "Lobatto_REC"; conclusive = c.denominator > 0)
        end
        z = RemainingBudgets.jacobi_rec(c.m, c.a, c.b, p, x)
        f == :jacobi&&return z
        return postdivide_allowance(rule, r, z, x, w, "Radau_REC_postdivide")
    elseif method in (:GST, :GST_postdivide)
        offset = f in (:radau, :radau_jacobi, :lobatto) ? 1 : 0
        v = f == :lobatto ? (1 - x^2) * w : f in (:radau, :radau_jacobi) ? (1 + x) * w : w
        transformed = f in (:radau, :radau_jacobi)&&method != :GST_postdivide
        z = NativeJacobiBudget.point_budget(c, i - offset, x, v; radau = transformed)
        relative = z.weight_relative
        if f == :lobatto
            relative += e * (2abs(x) / (1 - x^2) * z.C + 8)
        elseif method == :GST_postdivide
            return postdivide_allowance(rule, r, z, x, w, "Radau_GST_postdivide")
        end
        return result(z.node_absolute, relative, string(f, "_", method); conclusive = z.conclusive)
    end
    error("Unsupported assessment route: $method")
end

function overflow_certificate(rule)
    # If the average full-rule weight exceeds floatmax, at least one required
    # weight is unrepresentable. A normalization overflow alone is insufficient.
    n = int(rule, "n")
    rule["status"] == "exception" && n > 0 && !flag(rule, "reduced") ||
        return (; proven = false, logmean = "", loglimit = "", message = "No full-rule overflow certificate")
    a, b = represented(rule, "a"), represented(rule, "b")
    f = Symbol(rule["family"])
    f in (:jacobi, :radau_jacobi) && min(a, b) <= -1 && error("Nonintegrable input")
    f == :laguerre && a <= -1 && error("Nonintegrable input")
    logmass = if f in (:legendre, :radau, :lobatto)
        log(BigFloat(2))
    elseif f == :hermite
        flag(rule, "normalize") ? zero(a) : log(BigFloat(π)) / 2
    elseif f in (:jacobi, :radau_jacobi)
        (a + b + 1) * log(BigFloat(2)) + loggamma(a + 1) + loggamma(b + 1) - loggamma(a + b + 2)
    elseif f == :laguerre
        loggamma(a + 1)
    elseif f in (:chebyshevt, :chebyshevv, :chebyshevw)
        log(BigFloat(π))
    elseif f == :chebyshevu
        log(BigFloat(π) / 2)
    else
        error("Unknown mass for representation certificate")
    end
    logmean = logmass - log(BigFloat(n))
    loglimit = log(BigFloat(floatmax(format_type(rule["format"]))))
    proven = isfinite(logmean) && isfinite(loglimit) && logmean > loglimit
    return (; proven, logmean, loglimit, message = proven ? "Independent mass/n exceeds floatmax" : "Required weights are not proved unrepresentable")
end

function assess(rulepaths, samplepaths, output)
    ruleoutput = splitext(output)[1] * "_rules.csv"
    (isfile(output) || isfile(ruleoutput)) && error("Refusing existing assessment")
    rules = vcat(records.(rulepaths)...);points = vcat(records.(samplepaths)...)
    rulekeys = Set(key.(rules))
    length(rulekeys) == length(rules) || error("Duplicate rule identities")
    all(r -> key(r) in rulekeys, points) || error("Sample without a corresponding rule")
    groups = Dict(k => [r for r in points if key(r) == k] for k in rulekeys)
    for rule in rules
        rows = groups[key(rule)]
        length(rows) == int(rule, "sample_count") || error("Sample count mismatch for $(key(rule))")
        length(Set(int(r, "index") for r in rows)) == length(rows) || error("Duplicate sample indices")
        all(r -> 1 <= int(r, "index") <= int(rule, "returned_nodes"), rows) || error("Invalid sample index")
        rule["status"] == "ok" && isempty(rows) && error("Valid rule has no reference samples")
        for row in rows, field in ("family", "n", "a", "b", "reduced", "normalize", "output_bits")
            row[field] == rule[field] || error("Sample/rule identity mismatch: $field")
        end
    end
    candidates = filter(r -> r["implementation"] == "candidate", rules)
    isempty(candidates) && error("No candidate rules to assess")
    open(output, "w") do io
        open(ruleoutput, "w") do rio
            println(io, "implementation,format,case_id,index,regime,status,node_error,node_allowance,node_fraction,weight_error,weight_absolute_allowance,weight_fraction,weight_class")
            println(rio, "implementation,format,case_id,raw_status,status,sample_count,assessed_points,point_passes,point_failures,point_inconclusive,log_mass_per_node,log_floatmax,message")
            for rule in candidates
                setprecision(BigFloat, int(rule, "output_bits") + 256) do
                    rows = groups[key(rule)]
                    if rule["status"] != "ok"
                        certificate = overflow_certificate(rule)
                        status = certificate.proven ? "proved_representation_limit" : "fail"
                        csvrow(
                            rio, (
                                key(rule)..., rule["status"], status, length(rows), 0, 0, 0, 0,
                                certificate.logmean, certificate.loglimit, certificate.message,
                            )
                        )
                        return
                    end
                    all(
                        k -> flag(rule, k), (
                            "type_ok", "working_precision_preserved", "caller_precision_preserved",
                            "finite", "strict_order", "nonnegative",
                        )
                    ) || error("Contradictory valid-rule metadata")
                    nodes, weights, n = int(rule, "returned_nodes"), int(rule, "returned_weights"), int(rule, "n")
                    eligible_reduced = flag(rule, "reduced") && rule["family"] == "laguerre" &&
                        rule["format"] == "Float64" && n >= 128 && represented(rule, "a") < 5
                    nodes == weights && (eligible_reduced ? 0 < nodes <= n : nodes == n) || error("Invalid rule length")
                    int(rule, "reference_successes") == length(rows) &&
                        flag(rule, "reference_stable") && flag(rule, "root_indices_ok") || error("Incomplete rule references")
                    method = route(rule);c = context(rule, rows, method)
                    statuses = String[]
                    for r in rows
                        r["sample_status"] == "ok" && flag(r, "reference_stable") && flag(r, "root_index_ok") || error("Invalid reference or sample")
                        z = allowance(rule, r, c, method)
                        dx, dw = num(r, "node_abs"), num(r, "weight_abs")
                        all(isfinite, (dx, dw)) && min(dx, dw) >= 0 || error("Invalid error metrics")
                        reference_weight = num(r, "weight_reference")
                        reference_weight > 0 || error("Invalid positive-weight reference")
                        r["reference_weight_class"] == "normal" && iszero(represented(r, "weight_actual")) &&
                            error("A required normal weight was lost")
                        wa = reference_weight * z.weight_relative
                        if r["reference_weight_class"] in ("subnormal", "true_underflow", "zero")
                            T = format_type(r["format"])
                            wa += 2BigFloat(nextfloat(zero(T)))
                        elseif r["reference_weight_class"] != "normal"
                            error("Unsupported weight representation class")
                        end
                        node_allowance = z.node_absolute
                        # The new tail model supplements the earlier stricter
                        # Lobatto diagnostic; it does not replace that check.
                        extra_node = get(z, :extra_node_absolute, missing)
                        extra_weight = get(z, :extra_weight_relative, missing)
                        !ismissing(extra_node) && (node_allowance = min(node_allowance, extra_node))
                        !ismissing(extra_weight) && (wa = min(wa, reference_weight * extra_weight))
                        xr, wr = ratio(dx, node_allowance), ratio(dw, wa)
                        status = !z.conclusive ? "inconclusive" : max(xr, wr) <= 1 ? "pass" : "fail"
                        push!(statuses, status)
                        println(
                            io, join(
                                (
                                    key(rule)..., r["index"], z.regime, status,
                                    dx, node_allowance, xr, dw, wa, wr, r["reference_weight_class"],
                                ), ','
                            )
                        )
                    end
                    failed, pending = count(==("fail"), statuses), count(==("inconclusive"), statuses)
                    status = failed > 0 ? "fail" : pending > 0 ? "inconclusive" : "pass"
                    csvrow(
                        rio, (
                            key(rule)..., rule["status"], status, length(rows), length(statuses),
                            count(==("pass"), statuses), failed, pending, "", "", "",
                        )
                    )
                end
            end
        end
    end
    return (; points = output, rules = ruleoutput)
end
end
