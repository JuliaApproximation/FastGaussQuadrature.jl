# Check the committed Legendre tables with:
#     julia --startup-file=no test/generate_legendre.jl
# Pass an output path to write the generated coefficient section instead.
using LinearAlgebra

# Reproduce the Legendre tables with BigFloat arithmetic. Bogaert (2014),
# DOI 10.1137/140954969, equations (2.26–27), (3.10–11), (4.6),
# and Appendices A.2–3 and B.2–3. Byte-identical output was verified with
# Julia 1.13.3; other versions and platforms may differ.
# Use check_legendre_coefficients() to compare against the source, or write
# the generated section to an IO stream with generate_legendre_coefficients(io).
const LEGENDRE_COEFFICIENT_MARKER = "# Generated Bogaert coefficient tables; do not edit by hand."

function _legendre_coefficients(a)
    u = cot(a)
    q = u^2
    f1 = (u - 1 / a) / 8
    f2 = (6 * (1 + q) / a + 25 / a^3 - u * (31q + 33)) / 384
    f3 = u * (2595 + 6350q + 3779q^2) / 15360 - 1073 / (5120a^5) +
        (1 + q) * (-(31q + 11) / (1024a) + u / (512a^2) - 25 / (3072a^3))
    f4 = -u * (2407755 + 10808595q + 14682157q^2 + 6277237q^3) / 3440640 +
        375733 / (229376a^7) + (1 + q) * (
        (519 + 3810q + 3779q^2) / (24576a) -
            u * (21 + 31q) / (4096a^2) + (279 + 787q) / (49152a^3) -
            25u / (12288a^4) + 1073 / (40960a^5)
    )
    w2 = (81 / a^4 - 31u / a^3 - (3 - 6q) / a^2 + 6u / a - (27 + 84q + 56q^2)) / 384
    w3 = 153 // 1024 + 295q / 256 + 187q^2 / 96 + 151q^3 / 160 -
        u * (65 // 1024 + 119q / 768 + 35q^2 / 384) / a +
        (5 // 512 + 15q / 512 + 7q^2 / 384) / a^2 +
        u * (q / 512 - 13 // 1536) / a^3 + (53 // 3072 - 7q / 384) / a^4 +
        3749u / (15360a^5) - 1125 / (1024a^6)
    w4 = -21429 // 32768 - 36941q / 4096 - 27351q^2 / 1024 - 669667q^3 / 23040 -
        3626248438009q^4 / 338228674560 +
        u * (2513 // 8192 + 8639q / 6144 + 7393q^2 / 3840 + 997510355q^3 / 1207959552) / a -
        (371 // 16384 + 1483q / 8192 + 1909q^2 / 6144 + 1837891769q^3 / 12079595520) / a^2 +
        u * (675 // 16384 + 1849q / 18432 + 355532953q^2 / 6039797760) / a^3 -
        (1565 // 98304 + 1183q / 24576 + 147456121q^2 / 4831838208) / a^4 +
        u * (6823 // 245760 - 19906471q / 6039797760) / a^5 +
        (149694043q / 2415919104 - 156817 // 1474560) / a^6 -
        76749336551u / (42278584320a^7) + 568995840001 / (48318382080a^8)
    s = sin(a) / a
    return [s^(2m - 1) * f / a for (m, f) in enumerate((f1, f2, f3, f4))],
        [s^(2m) * w for (m, w) in zip(2:4, (w2, w3, w4))]
end

# Independent convergent Taylor series (not Float64 Bessel tables).
function _legendre_bessel01(x)
    t0, t1 = one(x), x / 2
    j0, j1 = t0, t1
    for k in 1:1000
        t0 *= -x^2 / (4k^2)
        t1 *= -x^2 / (4k * (k + 1))
        j0 += t0
        j1 += t1
        max(abs(t0), abs(t1)) < eps(x) / 100 && return j0, j1
    end
    error("Bessel series did not converge")
end

function generate_legendre_coefficients(io::IO; bits = 384, degrees = (7, 6, 5, 4, 7, 6, 5))
    # (j/phase) times (π*phase*J1(j)^2/2), using McMahon and the
    # J1 expansion already present in the baseline (DLMF 10.17.3, 10.21.19).
    root_ratio = (1 // 1, 1 // 8, -31 // 384, 3779 // 15360, -6277237 // 3440640, 2092163573 // 82575360)
    bessel_ratio = (1 // 1, 0 // 1, -7 // 48, 151 // 160, -172913 // 16128, 461797 // 2304)
    combined = [sum(root_ratio[k + 1] * bessel_ratio[j - k + 1] for k in 0:j) for j in 0:5]
    @assert combined == [1 // 1, 1 // 8, -29 // 128, 3599 // 3072, -6102707 // 491520, 2058116637 // 9175040]
    return setprecision(bits) do
        hi = BigFloat("1.58")
        coefs = map(enumerate(degrees)) do (j, degree)
            ts = [(1 + cospi(BigFloat(2k - 1) / (2degree + 2))) / 2 for k in 1:(degree + 1)]
            values = [vcat(_legendre_coefficients(hi * sqrt(t))...)[j] for t in ts]
            vandermonde = [t^k for t in ts, k in 0:degree]
            vandermonde \ values
        end
        # Recheck the interpolant on a separate dense grid, including a tiny angle.
        errors = zeros(BigFloat, length(coefs))
        rounded_errors = zeros(BigFloat, length(coefs))
        for k in 0:4096
            a = max(hi * k / 4096, BigFloat("1e-12"))
            ys = vcat(_legendre_coefficients(a)...)
            for j in eachindex(coefs)
                errors[j] = max(errors[j], abs(evalpoly((a / hi)^2, coefs[j]) - ys[j]))
                rounded_errors[j] = max(rounded_errors[j], abs(BigFloat(evalpoly((Float64(a) / 1.58)^2, Float64.(coefs[j]))) - ys[j]))
            end
        end
        println("grid absolute interpolation errors = ", Float64.(errors))
        println("grid absolute rounded evaluation errors = ", Float64.(rounded_errors))
        residuals, beta_corrections, weight_factors = BigFloat[], BigFloat[], BigFloat[]
        for k in 1:20
            phase = big(π) * (k - 1 // 4)
            z = phase + inv(8phase)
            for _ in 1:20
                j0, j1 = _legendre_bessel01(z)
                dz = j0 / j1
                z += dz
                abs(dz) < eps(z) * 8 && break
            end
            j0, j1 = _legendre_bessel01(z)
            @assert abs(j0) < big"1e-80"
            push!(residuals, z - phase)
            push!(beta_corrections, big(π) / 2 * z * j1^2 - 1)
            push!(weight_factors, 2 / j1^2)
        end
        println(io, LEGENDRE_COEFFICIENT_MARKER)
        println(io, "# Smooth Bogaert coefficients in t = (a / 1.58)^2, 0 <= a <= 1.58.")
        for (name, coefficients) in zip(("G1", "G2", "G3", "G4", "H2", "H3", "H4"), coefs)
            println(io, "const LEGENDRE_$name = (")
            for c in coefficients
                println(io, "    ", repr(Float64(c)), ",")
            end
            println(io, ")")
        end
        for (name, vals) in (
                ("LEGENDRE_BESSEL_RESIDUALS", residuals),
                ("LEGENDRE_BESSEL_BETA_CORRECTIONS", beta_corrections),
                ("LEGENDRE_BESSEL_WEIGHT_FACTORS", weight_factors),
            )
            println(io, "const $name = (")
            for c in vals
                println(io, "    ", repr(Float64(c)), ",")
            end
            println(io, ")")
        end
    end
end

function check_legendre_coefficients(path = joinpath(@__DIR__, "..", "src", "gausslegendre.jl"))
    sections = split(read(path, String), LEGENDRE_COEFFICIENT_MARKER; limit = 2)
    length(sections) == 2 || error("Generated coefficient marker not found in $path")
    sprint(generate_legendre_coefficients) == LEGENDRE_COEFFICIENT_MARKER * sections[2] ||
        error("Generated coefficients differ from $path")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    if isempty(ARGS)
        check_legendre_coefficients()
    else
        open(generate_legendre_coefficients, only(ARGS), "w")
    end
end
