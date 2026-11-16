# Keep the package project active and put this directory on JULIA_LOAD_PATH.
# Include this file, then call QuadratureQuality.run(). No source is modified.
# The pinned tooling environment requires Julia >= 1.12; package tests retain
# their separate Julia >= 1.10 compatibility and do not depend on these tools.
module QuadratureQuality

using FastGaussQuadrature, Quadmath
using Aqua, CodeComplexity, ExplicitImports, JET, Runic
using Dates, InteractiveUtils, LinearAlgebra, SHA

chebyshev32() = (gausschebyshevt(Float32, 31), gausschebyshevu(Float32, 31), gausschebyshevv(Float32, 31), gausschebyshevw(Float32, 31))
legendre_rec32() = gausslegendre(Float32, 60)
legendre_asy32() = gausslegendre(Float32, 61)
legendre_gst128() = gausslegendre(Quadmath.Float128, 61)
hermite_rec32() = gausshermite(Float32, 20; normalize = true)
hermite_gst32() = gausshermite(Float32, 21)
hermite_gst64() = gausshermite(201)
hermite_gst128() = gausshermite(Quadmath.Float128, 21)
laguerre_small32() = gausslaguerre(14, -0.99f0)
laguerre_rec64() = FastGaussQuadrature.gausslaguerre_rec(42, 0.0)
laguerre_gst64() = gausslaguerre(39, 44.8)
laguerre_gst128() = gausslaguerre(32, Quadmath.Float128(0.2))
laguerre_asy64() = FastGaussQuadrature.gausslaguerre_asy(10000, 0.0; reduced = true)
jacobi_gst32() = gaussjacobi(101, 0.25f0, 0.125f0)
jacobi_gst64() = gaussjacobi(101, 0.25, 0.125)
jacobi_asy64() = gaussjacobi(4096, 0.25, 0.0)
radau_gst32() = gaussradau(100, 10.0f0, 20.0f0)
radau_gst64() = gaussradau(21, 10.0, 20.0)
radau_edge64() = gaussradau(501, nextfloat(-1.0), 0.2)
radau_asy64() = gaussradau(10000)
lobatto_gst64() = gausslobatto(103)
lobatto_asy64() = gausslobatto(10000)
bessel_roots64() = approx_besselroots(0.25, 10)

const JET_CASES = (
    chebyshev32, legendre_rec32, legendre_asy32, legendre_gst128,
    hermite_rec32, hermite_gst32, hermite_gst64, hermite_gst128,
    laguerre_small32, laguerre_rec64, laguerre_gst64, laguerre_gst128,
    laguerre_asy64, jacobi_gst32, jacobi_gst64, jacobi_asy64,
    radau_gst32, radau_gst64, radau_edge64, radau_asy64, lobatto_gst64,
    lobatto_asy64, bessel_roots64,
)

function environment(io)
    println(io, "Recorded UTC: ", Dates.now(Dates.UTC))
    versioninfo(io)
    println(io, "Active project: ", Base.active_project())
    println(io, "LOAD_PATH: ", LOAD_PATH)
    println(io, "Julia threads: ", Threads.nthreads(), "; BLAS threads: ", BLAS.get_num_threads())
    for mod in (FastGaussQuadrature, Quadmath, Aqua, CodeComplexity, ExplicitImports, JET, Runic)
        println(io, nameof(mod), " ", Base.pkgversion(mod), " ", pathof(mod))
    end
    return nothing
end

function format_report(io, files)
    differences = String[]
    for file in files
        source = read(file, String)
        formatted = Runic.format_string(source; filemode = true)
        source == formatted || push!(differences, file)
        println(
            io, source == formatted ? "FORMAT MATCH " : "FORMAT DIFF ", file,
            " SHA256=", bytes2hex(sha256(source))
        )
    end
    isempty(differences) || error("Runic differences: $(join(differences, ", "))")
    return nothing
end

function imports_report(io)
    ExplicitImports.print_explicit_imports(io, FastGaussQuadrature; show_locations = true)
    ExplicitImports.check_no_implicit_imports(FastGaussQuadrature)
    ExplicitImports.check_no_stale_explicit_imports(FastGaussQuadrature)
    ExplicitImports.check_all_explicit_imports_are_public(FastGaussQuadrature)
    ExplicitImports.check_all_qualified_accesses_are_public(FastGaussQuadrature)
    return nothing
end

function jet_report(io)
    JET.JET_AVAILABLE || error("JET is unavailable for this Julia version; this check has not passed.")
    count = 0
    for f in JET_CASES
        println(io, "\nJET: ", nameof(f))
        result = JET.report_call(f, Tuple{})
        show(io, MIME"text/plain"(), result)
        println(io)
        count += length(JET.get_reports(result))
    end
    iszero(count) || error("JET reported $count possible errors; inspect all retained reports.")
    return nothing
end

function complexity_report(io, root)
    # Descriptive evidence, without a numeric pass threshold. Review complete
    # algorithms and their helper chains, not only individual function scores.
    for metric in (CodeComplexity.CognitiveComplexity(), CodeComplexity.CyclomaticComplexity())
        total = 0
        for file in sort(readdir(joinpath(root, "src"); join = true))
            endswith(file, ".jl") || continue
            result = CodeComplexity.measure_file(metric, file)
            total += result.total_value
            println(io, nameof(typeof(metric)), " ", basename(file), ": ", result.total_value)
            for f in result.functions
                println(io, "  ", f.name, " line ", f.line, " score ", f.value)
            end
        end
        println(io, "All production source ", nameof(typeof(metric)), ": ", total)
    end
    return nothing
end

function run(;
        io = stdout, root = dirname(dirname(pathof(FastGaussQuadrature))),
        files = [
            joinpath(root, folder, name) for folder in ("src", "test")
                for name in readdir(joinpath(root, folder)) if endswith(name, ".jl")
        ],
        checks = (:format, :imports, :aqua, :jet, :complexity)
    )
    realpath(joinpath(root, "Project.toml")) == realpath(Base.active_project()) ||
        error("Keep the package under validation as the active project.")
    environment(io)
    failures = Symbol[]
    for check in checks
        println(io, "\n=== ", check, " ===")
        try
            if check === :format
                format_report(io, files)
            elseif check === :imports
                imports_report(io)
            elseif check === :aqua
                Aqua.test_all(FastGaussQuadrature)
            elseif check === :jet
                jet_report(io)
            elseif check === :complexity
                complexity_report(io, root)
            else
                error("Unknown check: $check")
            end
            println(io, "PASS ", check)
        catch err
            err isa InterruptException && rethrow()
            push!(failures, check)
            showerror(io, err, catch_backtrace())
            println(io)
        end
        flush(io)
    end
    println(io, "Completed UTC: ", Dates.now(Dates.UTC))
    isempty(failures) || error("Quality checks failed: $(join(failures, ", "))")
    return nothing
end

end
