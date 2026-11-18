module ValidationAccuracy

using SpecialFunctions, Dates, SHA, Quadmath
include("cases.jl")

# Additional reference probes supplement the fixed performance matrix.
const DIAGNOSTIC_CASES = [
    (; id = "jacobi_issue48", family = :jacobi, n = 4097, a = 100.0, b = 0.2, reduced = false, normalize = false),
    (; id = "jacobi_endpoint1024", family = :jacobi, n = 1024, a = 0.25, b = 0.0, reduced = false, normalize = false),
    ((; id = string(f, "_257"), family = f, n = 257, a = 0.0, b = 0.0, reduced = false, normalize = false) for f in (:chebyshevt, :chebyshevu, :chebyshevv, :chebyshevw))...,
]
const GENERIC_CASE_IDS = Set(
    [
        "legendre_6", "legendre_60", "legendre_61", "legendre_1000",
        "hermite_20", "hermite_21", "hermite_201", "hermite_201_normalized",
        "jacobi_24_quarter_zero", "jacobi_25_quarter_zero", "jacobi_100_near_endpoint",
        "jacobi_100_large_parameters", "jacobi_endpoint1024",
        "laguerre_14_zero_full", "laguerre_15_zero_full", "laguerre_39_44p8_full", "laguerre_128_point2_full",
        "radau_3", "radau_102", "radau_jacobi_100_large_parameters", "lobatto_4", "lobatto_103",
        "chebyshevt_3", "chebyshevu_3", "chebyshevv_3", "chebyshevw_3",
        "chebyshevt_257", "chebyshevu_257", "chebyshevv_257", "chebyshevw_257",
    ]
)

include(joinpath(@__DIR__, "..", "test", "quadrature_reference.jl"))


# This is an evidence driver, not production code or a routine test dependency.
# All references use the exact binary parameters supplied to each implementation.

function csvfield(value)
    text = string(value)
    return occursin(r"[,\"\n\r]", text) ? "\"" * replace(text, "\"" => "\"\"") * "\"" : text
end
csvrow(io, values) = println(io, join(csvfield.(values), ','))

function load_implementations(directory)
    return map(("baseline", "candidate")) do label
        wrapper = Module(gensym(:AccuracySnapshot))
        Base.include(wrapper, joinpath(directory, label, "src", "FastGaussQuadrature.jl"))
        (label, Base.invokelatest(getfield, wrapper, :FastGaussQuadrature))
    end
end

function call_rule(mod, c, T)
    a, b = T(c.a), T(c.b)
    name = c.family == :radau_jacobi ? :gaussradau : Symbol("gauss", c.family)
    f = getfield(mod, name)
    if c.family in (:jacobi, :radau_jacobi)
        return Base.invokelatest(f, c.n, a, b)
    elseif c.family == :laguerre
        return Base.invokelatest(f, c.n, a; reduced = c.reduced)
    elseif c.family == :hermite
        return Base.invokelatest(f, T, c.n; normalize = c.normalize)
    end
    return Base.invokelatest(f, T, c.n)
end

function default_formats()
    formats = [
        (label = "Float32", T = Float32, bits = 24),
        (label = "Float64", T = Float64, bits = 53),
        (label = "BigFloat128", T = BigFloat, bits = 128),
        (label = "BigFloat256", T = BigFloat, bits = 256),
    ]
    quadmath_error = ""
    try
        precision(Quadmath.Float128) == 113 || error("Float128 does not have 113 significand bits")
        insert!(formats, 3, (label = "Float128", T = Quadmath.Float128, bits = 113))
    catch err
        err isa InterruptException && rethrow()
        quadmath_error = sprint(showerror, err)
    end
    return formats, quadmath_error
end

function sample_indices(x, w, T, bits; exhaustive = bits <= 53 ? 128 : 40)
    n = length(x)
    n <= exhaustive && return collect(1:n)
    indices = [
        1, 2, 3, n ÷ 8, n ÷ 4, n ÷ 3, n ÷ 2 - 1, n ÷ 2,
        n ÷ 2 + 1, n ÷ 2 + 2, 2n ÷ 3, 3n ÷ 4, 7n ÷ 8, n - 2, n - 1, n,
    ]
    # Include both sides of each normal/subnormal/zero transition.
    previous_class = ""
    for i in eachindex(w)
        current_class = value_class(w[i], T)
        if current_class != previous_class
            append!(indices, (i - 1, i, i + 1))
        end
        previous_class = current_class
    end
    return sort!(unique!(clamp.(indices, 1, n)))
end

function value_class(value, T)
    !isfinite(value) && return isnan(value) ? "nan" : "infinite"
    iszero(value) && return "zero"
    value < 0 && return "negative"
    return abs(value) < floatmin(T) ? "subnormal" : "normal"
end

function reference_parameters(c, T, bits)
    # Round inputs at the output precision before entering reference arithmetic.
    a, b = setprecision(BigFloat, bits) do
        T(c.a), T(c.b)
    end
    return BigFloat(a), BigFloat(b)
end

function chebyshev_reference(c, index)
    n, i, pi = BigFloat(c.n), BigFloat(index), BigFloat(π)
    if c.family == :chebyshevt
        angle = (2i - 1) * pi / (2n)
        x = -cos(angle)
        isodd(c.n) && index == (c.n + 1) ÷ 2 && (x = zero(x))
        return x, pi / n
    elseif c.family == :chebyshevu
        angle = i * pi / (n + 1)
        x = -cos(angle)
        isodd(c.n) && index == (c.n + 1) ÷ 2 && (x = zero(x))
        return x, 2pi / (n + 1) * sin(angle)^2 / 2
    elseif c.family == :chebyshevv
        angle = 2i * pi / (2n + 1)
        return -cos(angle), 4pi / (2n + 1) * sin(angle / 2)^2
    else
        angle = (2i - 1) * pi / (2n + 1)
        return -cos(angle), 4pi / (2n + 1) * cos(angle / 2)^2
    end
end

function radau_endpoint(n, a, b)
    # The product mass * Π k/(k+b+1) * (k+a)/(k+a+b+1)
    # reduces to this gamma ratio. Positive gamma arguments for a,b > -1.
    nn = BigFloat(n)
    return exp(
        (a + b + 1) * log(BigFloat(2)) +
            loggamma(b + 1) + loggamma(b + 2) + loggamma(nn + a) +
            loggamma(nn) - loggamma(nn + b + 1) - loggamma(nn + a + b + 1)
    )
end

function reference(c, index, seed, T, bits; output_bits = bits)
    reference_bits = bits + 128
    return setprecision(BigFloat, reference_bits) do
        a, b = reference_parameters(c, T, output_bits)
        if startswith(string(c.family), "chebyshev")
            x, w = chebyshev_reference(c, index)
            return (; x, w, indexed = true, root_index = index, reference_bits)
        end
        family, degree, root_index = :jacobi, c.n, index
        if c.family == :legendre
            a = b = BigFloat(0)
        elseif c.family in (:radau, :radau_jacobi)
            c.family == :radau && (a = b = BigFloat(0))
            if index == 1
                return (
                    x = BigFloat(-1), w = radau_endpoint(c.n, a, b),
                    indexed = true, root_index = 0, reference_bits,
                )
            end
            degree, root_index = c.n - 1, index - 1
            b += 1
        elseif c.family == :lobatto
            if index == 1 || index == c.n
                return (
                    x = BigFloat(index == 1 ? -1 : 1),
                    w = BigFloat(2) / (BigFloat(c.n) * (c.n - 1)),
                    indexed = true, root_index = 0, reference_bits,
                )
            end
            degree, root_index, a, b = c.n - 2, index - 1, BigFloat(1), BigFloat(1)
        elseif c.family in (:laguerre, :hermite)
            family = c.family
        end
        reference_seed = BigFloat(seed)
        c.family == :hermite && c.normalize && (reference_seed /= sqrt(BigFloat(2)))
        # Symmetry gives an exact central zero, independently of the supplied seed.
        if isodd(degree) && root_index == (degree + 1) ÷ 2 &&
                (family == :hermite || (family == :jacobi && a == b))
            reference_seed = BigFloat(0)
        end
        x, w = QuadratureReference.nodeweight(family, degree, a, b, reference_seed; bits)
        radius = BigFloat(2)^(-reference_bits ÷ 2) * max(one(x), abs(x))
        below = QuadratureReference.roots_below(family, degree, a, b, x - radius)
        above = QuadratureReference.roots_below(family, degree, a, b, x + radius)
        indexed = below == root_index - 1 && above == root_index
        if c.family == :radau
            p = QuadratureReference.polynomials(:jacobi, c.n - 1, BigFloat(0), BigFloat(0), x)[end]
            w = (1 - x) / (BigFloat(c.n)^2 * p^2)
        elseif c.family == :radau_jacobi
            w /= 1 + x
        elseif c.family == :lobatto
            p = QuadratureReference.polynomials(:jacobi, c.n - 1, BigFloat(0), BigFloat(0), x)[end]
            w = 2 / (BigFloat(c.n) * (c.n - 1) * p^2)
        elseif c.family == :hermite && c.normalize
            x *= sqrt(BigFloat(2))
            w /= sqrt(BigFloat(π))
        end
        return (; x, w, indexed, root_index, reference_bits)
    end
end

function target_round(value, T, bits)
    return T === BigFloat ? setprecision(() -> BigFloat(value), BigFloat, bits) : T(value)
end

function error_metrics(value, expected, T, bits)
    actual = BigFloat(value)
    absolute = abs(actual - expected)
    relative = iszero(expected) ? BigFloat(NaN) : absolute / abs(expected)
    rounded = target_round(expected, T, bits)
    # Use spacing adjacent to the rounded reference toward the observed value.
    # Absolute error is retained for zero and very small nodes.
    spacing = if isfinite(rounded)
        direction = actual - BigFloat(rounded)
        # If already correctly rounded, use the interval containing the exact
        # reference. The two intervals differ at a power-of-two boundary.
        iszero(direction) && (direction = expected - BigFloat(rounded))
        neighbor = setprecision(BigFloat, bits) do
            direction >= 0 ? nextfloat(rounded) : prevfloat(rounded)
        end
        abs(BigFloat(neighbor) - BigFloat(rounded))
    else
        BigFloat(NaN)
    end
    ulp = iszero(absolute) ? BigFloat(0) :
        (isfinite(spacing) && spacing > 0 ? absolute / spacing : BigFloat(NaN))
    return (; absolute, relative, ulp, rounded)
end

function rule_mass(c, T, bits)
    a, b = reference_parameters(c, T, bits)
    c.family in (:legendre, :radau, :lobatto) && return BigFloat(2)
    c.family == :hermite && return c.normalize ? BigFloat(1) : sqrt(BigFloat(π))
    c.family == :laguerre && return gamma(a + 1)
    c.family == :chebyshevu && return BigFloat(π) / 2
    startswith(string(c.family), "chebyshev") && return BigFloat(π)
    return exp(
        (a + b + 1) * log(BigFloat(2)) +
            loggamma(a + 1) + loggamma(b + 1) - loggamma(a + b + 2)
    )
end

function expected_weight_class(expected, T, bits)
    rounded = target_round(expected, T, bits)
    iszero(rounded) && expected > 0 && return "true_underflow"
    !isfinite(rounded) && return "overflow"
    return value_class(rounded, T)
end

const SAMPLE_HEADER = (
    "implementation", "case_id", "family", "n", "a", "b", "reduced", "normalize",
    "format", "output_bits", "index", "reference_bits", "sample_status",
    "node_actual", "node_reference", "node_abs", "node_rel", "node_ulp",
    "weight_actual", "weight_reference", "weight_abs", "weight_rel", "weight_ulp",
    "weight_class", "reference_weight_class", "reference_node_abs_change",
    "reference_weight_rel_change", "reference_stable", "root_index_ok", "message",
)

const RULE_HEADER = (
    "implementation", "case_id", "family", "n", "a", "b", "reduced", "normalize",
    "format", "output_bits", "status", "returned_nodes", "returned_weights",
    "node_type", "weight_type", "type_ok", "working_precision_preserved",
    "caller_precision_preserved", "finite", "strict_order", "nonnegative",
    "sample_count", "reference_successes", "reference_stable", "root_indices_ok",
    "max_node_abs", "max_node_ulp", "max_weight_rel", "max_weight_ulp",
    "max_normal_weight_ulp", "max_subnormal_weight_ulp", "mass_abs_error",
    "symmetry_abs_error", "message",
)

function case_prefix(label, c, format)
    return (
        label, c.id, c.family, c.n, c.a, c.b, c.reduced, c.normalize,
        format.label, format.bits,
    )
end

function audit_case(ruleio, sampleio, label, mod, c, format; max_reference_order)
    T, bits = format.T, format.bits
    prefix = case_prefix(label, c, format)
    x = w = nothing
    caller_precision = precision(BigFloat)
    working_precision_preserved = true
    try
        x, w = setprecision(BigFloat, bits) do
            try
                call_rule(mod, c, T)
            finally
                working_precision_preserved = precision(BigFloat) == bits
            end
        end
    catch err
        err isa InterruptException && rethrow()
        message = sprint(showerror, err)
        csvrow(
            ruleio, (
                prefix..., "exception", 0, 0, "", "", false,
                working_precision_preserved, precision(BigFloat) == caller_precision,
                false, false, false, 0, 0, false, false, "", "", "", "", "", "", "", "", message,
            )
        )
        return
    end
    caller_precision_preserved = precision(BigFloat) == caller_precision
    return setprecision(BigFloat, bits + 192) do
        count = min(length(x), length(w))
        matching_lengths = length(x) == length(w)
        expected_length = c.reduced && c.family == :laguerre ?
            (c.n == 0 ? count == 0 : 0 < count <= c.n) : count == c.n
        type_ok = eltype(x) === T && eltype(w) === T &&
            (T !== BigFloat || (all(z -> precision(z) == bits, x) && all(z -> precision(z) == bits, w)))
        finite = all(isfinite, x) && all(isfinite, w)
        ordered = all(i -> x[i] < x[i + 1], 1:max(0, length(x) - 1))
        nonnegative = all(z -> z >= 0, w)
        positive_weight = c.n == 0 || any(z -> z > 0, w)
        status = matching_lengths && expected_length && finite && ordered &&
            nonnegative && positive_weight && type_ok && working_precision_preserved &&
            caller_precision_preserved ? "ok" : "structural_failure"
        indices = matching_lengths ? sample_indices(x, w, T, bits) : Int[]
        c.n > max_reference_order && (status = status == "ok" ? "reference_skipped" : status)
        (c.n > max_reference_order || !finite) && empty!(indices)
        maximum_node_abs = maximum_node_ulp = maximum_weight_rel = maximum_weight_ulp = BigFloat(0)
        maximum_normal_ulp = maximum_subnormal_ulp = BigFloat(0)
        references_ok, stable_all, index_all = 0, true, true
        normal_weight_lost = false
        messages = String[]
        for index in indices
            try
                first_reference = reference(c, index, x[index], T, bits)
                second_reference = reference(c, index, x[index], T, bits + 64; output_bits = bits)
                xr, wr = second_reference.x, second_reference.w
                dx = abs(first_reference.x - xr)
                dw = abs(first_reference.w - wr) / abs(wr)
                stable = dx <= BigFloat(2)^(-bits - 32) * max(one(xr), abs(xr)) &&
                    dw <= BigFloat(2)^(-bits - 32)
                indexed = first_reference.indexed && second_reference.indexed
                xm = error_metrics(x[index], xr, T, bits)
                wm = error_metrics(w[index], wr, T, bits)
                weight_class = value_class(w[index], T)
                reference_class = expected_weight_class(wr, T, bits)
                normal_weight_lost |= reference_class == "normal" && iszero(w[index])
                sample_status = reference_class == "normal" && iszero(w[index]) ? "normal_weight_loss" : !stable ? "unstable_reference" :
                    !indexed ? "wrong_root_index" : "ok"
                csvrow(
                    sampleio, (
                        prefix..., index, second_reference.reference_bits,
                        sample_status, x[index], xr, xm.absolute, xm.relative, xm.ulp,
                        w[index], wr, wm.absolute, wm.relative, wm.ulp, weight_class,
                        reference_class, dx, dw, stable, indexed, "",
                    )
                )
                maximum_node_abs = max(maximum_node_abs, xm.absolute)
                isfinite(xm.ulp) && (maximum_node_ulp = max(maximum_node_ulp, xm.ulp))
                isfinite(wm.relative) && (maximum_weight_rel = max(maximum_weight_rel, wm.relative))
                isfinite(wm.ulp) && (maximum_weight_ulp = max(maximum_weight_ulp, wm.ulp))
                if reference_class == "normal" && isfinite(wm.ulp)
                    maximum_normal_ulp = max(maximum_normal_ulp, wm.ulp)
                elseif reference_class == "subnormal" && isfinite(wm.ulp)
                    maximum_subnormal_ulp = max(maximum_subnormal_ulp, wm.ulp)
                end
                references_ok += 1
                stable_all &= stable
                index_all &= indexed
            catch err
                err isa InterruptException && rethrow()
                message = sprint(showerror, err)
                push!(messages, "index $(index): $(message)")
                csvrow(
                    sampleio, (
                        prefix..., index, bits + 192, "reference_exception",
                        x[index], "", "", "", "", w[index], "", "", "", "",
                        value_class(w[index], T), "", "", "", false, false, message,
                    )
                )
                stable_all = index_all = false
            end
        end
        if !isempty(indices) && (references_ok != length(indices) || !stable_all || !index_all)
            status = status == "ok" ? "reference_failure" : status
        end
        normal_weight_lost && (status = "normal_weight_loss")
        mass_error = finite && c.n > 0 ? abs(sum(BigFloat, w) - rule_mass(c, T, bits)) : BigFloat(NaN)
        symmetric = c.family in (:legendre, :hermite, :lobatto, :chebyshevt, :chebyshevu) ||
            (c.family == :jacobi && T(c.a) == T(c.b))
        symmetry_error = symmetric && count > 0 && finite ? maximum(
                abs(BigFloat(x[i]) + BigFloat(x[end + 1 - i])) for i in eachindex(x)
            ) : BigFloat(NaN)
        no_references = references_ok == 0 && c.n != 0
        values = no_references ? ("", "", "", "", "", "") :
            (
                maximum_node_abs, maximum_node_ulp, maximum_weight_rel, maximum_weight_ulp,
                maximum_normal_ulp, maximum_subnormal_ulp,
            )
        csvrow(
            ruleio, (
                prefix..., status, length(x), length(w), eltype(x), eltype(w),
                type_ok, working_precision_preserved, caller_precision_preserved,
                finite, ordered, nonnegative, length(indices), references_ok,
                !no_references && stable_all, !no_references && index_all, values...,
                mass_error, symmetry_error, join(messages, " | "),
            )
        )
    end
end

"""
    run(evidence_dir; implementations=nothing, formats=nothing,
        include_diagnostics=true, include_generic=true, max_reference_order=100000,
        selected_case_ids=nothing, output_prefix="accuracy")

Save rule-level and sample-level CSVs. Float64 uses all timing rows and optional
diagnostics; other formats use GENERIC_CASE_IDS. Exceptions and explicit reference
skips are preserved. No low-precision-floor policy assertion is made by this
driver; the separate Legendre/contract checks cover that policy.
"""
function run(
        evidence_dir; implementations = nothing, formats = nothing,
        include_diagnostics = true, include_generic = true, max_reference_order = 100000,
        selected_case_ids = nothing, output_prefix = "accuracy"
    )
    implementations === nothing && (implementations = load_implementations(evidence_dir))
    quadmath_error = ""
    if formats === nothing
        formats, quadmath_error = default_formats()
    end
    cases = vcat(TIMING_CASES, DIAGNOSTIC_CASES)
    lookup = Dict(c.id => c for c in cases)
    all_cases = include_diagnostics ? cases : TIMING_CASES
    selected = selected_case_ids === nothing ? nothing : Set(selected_case_ids)
    any(isfile(joinpath(evidence_dir, output_prefix * suffix)) for suffix in ("_rules.csv", "_samples.csv", "_configuration.txt")) && error("Refusing existing accuracy outputs")
    rulepath = joinpath(evidence_dir, output_prefix * "_rules.csv")
    samplepath = joinpath(evidence_dir, output_prefix * "_samples.csv")
    open(joinpath(evidence_dir, output_prefix * "_configuration.txt"), "w") do io
        println(io, "Started UTC: ", now(UTC))
        println(io, "Julia: ", VERSION)
        println(io, "Project: ", Base.active_project())
        println(io, "Driver SHA256: ", bytes2hex(sha256(read(@__FILE__))))
        for (label, mod) in implementations
            source = string(which(mod.gausslegendre, (Type{Float64}, Int)).file)
            for name in sort(readdir(dirname(source); join = true))
                endswith(name, ".jl") || continue
                println(io, label, " source ", name, " SHA256 ", bytes2hex(sha256(read(name))))
            end
        end
        println(io, "formats: ", [(f.label, f.bits) for f in formats])
        println(io, "include_diagnostics: ", include_diagnostics)
        println(io, "include_generic: ", include_generic)
        println(io, "max_reference_order: ", max_reference_order)
        println(io, "selected_case_ids: ", selected_case_ids)
        println(io, "Quadmath availability error: ", quadmath_error)
        println(io, "Reference precision: output bits +128, repeated at +192")
        println(io, "All-node threshold: n<=128 at <=53 bits; n<=40 otherwise")
        println(io, "Remaining indices: endpoints, center, fixed fractions, weight-class transitions")
        println(io, "ULP convention: adjacent target-format spacing toward actual value from rounded reference")
        println(io, "Rule status ok means structural/reference checks; accuracy acceptance requires justified native-precision regime budgets")
    end
    open(rulepath, "w") do ruleio
        open(samplepath, "w") do sampleio
            csvrow(ruleio, RULE_HEADER)
            csvrow(sampleio, SAMPLE_HEADER)
            for format in formats
                format.T !== Float64 && !include_generic && continue
                format_cases = format.T === Float64 ? all_cases :
                    [lookup[id] for id in GENERIC_CASE_IDS]
                for (label, mod) in implementations, c in format_cases
                    selected !== nothing && !(c.id in selected) && continue
                    println("accuracy ", label, " ", format.label, " ", c.id)
                    flush(stdout)
                    Base.invokelatest(audit_case, ruleio, sampleio, label, mod, c, format; max_reference_order)
                    flush(ruleio)
                    flush(sampleio)
                end
            end
        end
    end
    open(joinpath(evidence_dir, output_prefix * "_configuration.txt"), "a") do io
        println(io, "Completed UTC: ", now(UTC))
    end
    return (; rulepath, samplepath, quadmath_error)
end

end
