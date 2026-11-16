module QuadratureMeasurements
using FastGaussQuadrature, Quadmath, LinearAlgebra, Dates, SHA
include("benchmarks.jl")
const B = ValidationTimings
function load_snapshot(path)
    wrapper = Module(gensym(:TimingSource))
    Base.include(wrapper, joinpath(path, "src", "FastGaussQuadrature.jl"))
    return Base.invokelatest(getfield, wrapper, :FastGaussQuadrature)
end
case(family, n; a = 0.0, b = 0.0, reduced = false, normalize = false) =
    (; id = string(family, "_", n, "_", a, "_", b, "_", reduced, "_", normalize), family, n, a, b, reduced, normalize)
const NATIVE_CASES = [
    case(:legendre, 61), case(:legendre, 257), case(:hermite, 20), case(:hermite, 257), case(:hermite, 257; normalize = true),
    case(:jacobi, 25; a = 0.25, b = 0.0), case(:jacobi, 257; a = 0.25, b = 0.125),
    case(:laguerre, 128; a = 0.2), case(:laguerre, 257; a = 0.2),
    case(:radau, 257), case(:radau_jacobi, 257; a = 10.0, b = 20.0), case(:lobatto, 257),
    (case(f, 257) for f in (:chebyshevt, :chebyshevu, :chebyshevv, :chebyshevw))...,
]
const SCALING_CASES = [
    case(f, n; a = a, b = b, reduced = reduced) for (f, a, b, reduced) in
        (
            (:legendre, 0.0, 0.0, false), (:hermite, 0.0, 0.0, false), (:jacobi, 10.0, 20.0, false),
            (:jacobi, 0.25, 0.125, false), (:laguerre, 0.0, 0.0, false), (:laguerre, 0.2, 0.0, false),
            (:laguerre, 0.2, 0.0, true), (:radau, 0.0, 0.0, false), (:radau_jacobi, 10.0, 20.0, false),
            (:lobatto, 0.0, 0.0, false), (:chebyshevt, 0.0, 0.0, false), (:chebyshevu, 0.0, 0.0, false),
            (:chebyshevv, 0.0, 0.0, false), (:chebyshevw, 0.0, 0.0, false),
        ) for n in (2048, 4096, 8192, 16384)
]
function configuration(path)
    return open(path, "w") do io
        println(io, "Started UTC ", now(UTC), " Julia ", VERSION, " project ", Base.active_project())
        println(io, "Julia threads ", Threads.nthreads(), " BLAS threads ", BLAS.get_num_threads())
        println(io, "Float64: exact agreed60 cases; eight alternating batches ×30 samples; evals1; warmed pooled medians")
        println(io, "Native diagnostics: sixteen cases, Float32/Float128/BigFloat128/256; three alternating batches ×10samples")
        println(io, "Scaling: fourteen family/regime/semantics combinations ×four doubling orders; three batches ×10samples")
        println(io, "Compilation excluded by warmup; invalid baseline rules remain visible and are not cap passes")
        for file in sort(readdir(dirname(pathof(FastGaussQuadrature)); join = true))
            endswith(file, ".jl")&&println(io, file, " SHA256 ", bytes2hex(sha256(read(file))))
        end
        println(io, "Driver SHA256 ", bytes2hex(sha256(read(@__FILE__))))
    end
end
function run(baseline_path, output_dir; tag = "final")
    BLAS.set_num_threads(1)
    mkpath(output_dir)
    master = load_snapshot(baseline_path)
    config = joinpath(output_dir, tag * "-timing-configuration.txt");isfile(config)&&error("Existing configuration")
    configuration(config)
    open(config, "a") do io
        println(io, "Baseline source: ", abspath(baseline_path))
        for file in sort(readdir(joinpath(baseline_path, "src"); join = true))
            endswith(file, ".jl")&&println(io, "Baseline ", basename(file), " SHA256 ", bytes2hex(sha256(read(file))))
        end
        for mod in (FastGaussQuadrature, Quadmath, B.BenchmarkTools)
            println(io, nameof(mod), " version ", Base.pkgversion(mod))
        end
    end
    Base.invokelatest(B.run, ["master" => master, "candidate" => FastGaussQuadrature]; evidence_dir = output_dir, batches = 8, samples = 30, tag = tag * "_float64")
    for (T, bits, label) in ((Float32, 24, "float32"), (Float128, 113, "float128"), (BigFloat, 128, "bigfloat128"), (BigFloat, 256, "bigfloat256"))
        Base.invokelatest(B.run, ["master" => master, "candidate" => FastGaussQuadrature]; evidence_dir = output_dir, cases = NATIVE_CASES, T, bits, batches = 3, samples = 10, tag = tag * "_" * label)
    end
    Base.invokelatest(B.run, ["candidate" => FastGaussQuadrature]; evidence_dir = output_dir, cases = SCALING_CASES, batches = 3, samples = 10, tag = tag * "_scaling")
    return open(config, "a") do io
        println(io, "Completed UTC ", now(UTC))
    end
end
end
