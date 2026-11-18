module ValidationTimings
using BenchmarkTools, Statistics, LinearAlgebra, Logging
include("cases.jl")

csvrow(io, row) = println(io, join(("\"" * replace(string(v), "\"" => "\"\"") * "\"" for v in row), ','))

function rule_closure(mod, c, T)
    n, a, b, reduced, normalize = c.n, T(c.a), T(c.b), c.reduced, c.normalize
    name = c.family == :radau_jacobi ? :gaussradau : Symbol("gauss", c.family)
    f = getfield(mod, name)
    if c.family in (:jacobi, :radau_jacobi)
        return () -> f(n, a, b)
    elseif c.family == :laguerre
        return () -> f(n, a; reduced)
    elseif c.family == :hermite
        return () -> f(T, n; normalize)
    else
        return () -> f(T, n)
    end
end

function run(
        implementations; evidence_dir = @__DIR__, cases = TIMING_CASES,
        T = Float64, bits = precision(T), batches = 3, samples = 30, tag = "timing_initial"
    )
    BLAS.set_num_threads(1)
    any(isfile(joinpath(evidence_dir, tag * suffix)) for suffix in ("_raw.csv", "_summary.csv", "_status.csv")) && error("Refusing to overwrite timing outputs")
    open(joinpath(evidence_dir, tag * "_raw.csv"), "w") do raw
        open(joinpath(evidence_dir, tag * "_summary.csv"), "w") do summary
            open(joinpath(evidence_dir, tag * "_status.csv"), "w") do status
                csvrow(raw, ("id", "implementation", "type", "bits", "batch", "sample", "time_ns", "gc_ns", "bytes", "allocations"))
                csvrow(summary, ("id", "implementation", "type", "bits", "samples", "median_ns", "minimum_ns", "bytes", "allocations"))
                csvrow(status, ("id", "implementation", "type", "bits", "status", "returned_length", "normal", "subnormal", "zero", "message"))
                setprecision(BigFloat, bits) do
                    for c in cases
                        ready = Pair{String, Any}[]
                        for (label, mod) in implementations
                            try
                                f = rule_closure(mod, c, T)
                                x, w = f()
                                valid = length(x) == length(w) && (c.reduced ? length(x) <= c.n : length(x) == c.n) &&
                                    all(isfinite, x) && all(isfinite, w) && all(>=(0), w) && all(>(0), diff(x))
                                normal = count(z -> abs(z) >= floatmin(T), w)
                                subnormal = count(z -> 0 < abs(z) < floatmin(T), w)
                                zeros = count(iszero, w)
                                csvrow(status, (c.id, label, T, bits, valid ? "structurally_valid" : "invalid_output", length(x), normal, subnormal, zeros, "Accuracy is assessed separately"))
                                f(); f() # exclude compilation and warm both implementations before alternating
                                push!(ready, label => f)
                            catch err
                                err isa InterruptException && rethrow()
                                csvrow(status, (c.id, label, T, bits, "exception", 0, 0, 0, 0, sprint(showerror, err)))
                            end
                        end
                        timings = Dict(label => Float64[] for (label, _) in ready)
                        allocations = Dict{String, Tuple{Int, Int}}()
                        for batch in 1:batches
                            order = isodd(batch) ? ready : reverse(ready)
                            for (label, f) in order
                                trial = with_logger(NullLogger()) do
                                    @benchmark $f() samples = samples evals = 1 seconds = 30
                                end
                                append!(timings[label], trial.times)
                                allocations[label] = (trial.memory, trial.allocs)
                                for i in eachindex(trial.times)
                                    csvrow(raw, (c.id, label, T, bits, batch, i, trial.times[i], trial.gctimes[i], trial.memory, trial.allocs))
                                end
                            end
                        end
                        for (label, _) in ready
                            times = timings[label]
                            csvrow(summary, (c.id, label, T, bits, length(times), median(times), minimum(times), allocations[label]...))
                        end
                        flush(raw); flush(summary); flush(status)
                    end
                end
            end
        end
    end
    return nothing
end
end
