@testset "Gauss–Hermite" begin
    @testset "Check error" begin
        @test_throws DomainError gausshermite(-1)
    end

    @testset "Small n case" begin
        @test gausshermite(0) == (Float64[], Float64[])
        @testset "Type: $T" for T in (Float16, Float32, Float64, BigFloat)
            @test gausshermite(T, 0) == (T[], T[])
            @test gausshermite(T, 1; normalize = true) == (T[0], T[1])
            @test gausshermite(T, 1) == (T[0], T[sqrt(T(π))])
            @test FastGaussQuadrature.unweightedgausshermite(T, 0) == (T[], T[])
            @test_throws DomainError FastGaussQuadrature.unweightedgausshermite(T, -1)
        end
        x, w = gausshermite(1)
        @test x ≈ [0.0]
        @test w ≈ [√π]
    end

    tol = 1.0e-14

    @testset "Golub-Welsch" begin
        n = 18
        x, w = gausshermite(n)
        @test isa(x, Vector{Float64})
        @test isa(w, Vector{Float64})
        @test (length(x) == n && length(w) == n)
        @test (dot(w, x) < tol && abs(dot(w, x .^ 2) - sqrt(π) / 2) < tol)
    end

    @testset "Recurrence" begin
        n = 42
        x, w = gausshermite(n)
        @test isa(x, Vector{Float64})
        @test isa(w, Vector{Float64})
        @test (length(x) == n && length(w) == n)
        @test (dot(w, x) < tol && abs(dot(w, x .^ 2) - sqrt(π) / 2) < tol)
        @test abs(x[37] - 5.660357581283058) < tol
        @test abs(w[17] - 0.032202101288908) < tol

        f = x -> 1 + x + x^2 + x^3
        @test w'f.(x) ≈ 2.6586807763789717
    end

    @testset "Asymptotics"  begin
        n = 251
        x, w = gausshermite(n)
        @test isa(x, Vector{Float64})
        @test isa(w, Vector{Float64})
        @test (length(x) == n && length(w) == n)
        @test (dot(w, x) < tol && abs(dot(w, x .^ 2) - sqrt(π) / 2) < 300 * tol)
        @test abs(x[37] - -13.292221459334638) < tol
        @test abs(w[123] - 0.117419270715955) < 2 * tol

        n = 500
        x, w = gausshermite(n)
        @test isa(x, Vector{Float64})
        @test isa(w, Vector{Float64})
        @test (length(x) == n && length(w) == n)
        @test (dot(w, x) < tol && abs(dot(w, x .^ 2) - sqrt(π) / 2) < 300 * tol)
    end

    @testset "Recurrence" begin
        x = 0.1
        He0 = 1 # Probabilists He_n(x) = H_n(x/sqrt(2))/2^(n/2)
        He1 = x; He2 = x^2 - 1; He3 = x * (x^2 - 3); He4 = 3 - 6x^2 + x^4
        @test FastGaussQuadrature.hermpoly_rec(1, x)[1] ≈ He1 * exp(-x^2 / 4)
        @test FastGaussQuadrature.hermpoly_rec(2, x)[1] ≈ He2 * exp(-x^2 / 4) / sqrt(2)
        @test FastGaussQuadrature.hermpoly_rec(3, x)[1] ≈ He3 * exp(-x^2 / 4) / sqrt(2 * 3)
        @test FastGaussQuadrature.hermpoly_rec(4, x)[1] ≈ He4 * exp(-x^2 / 4) / sqrt(2 * 3 * 4)

        @test FastGaussQuadrature.hermpoly_rec(1, 60)[1] ≈ 0.0
        @test FastGaussQuadrature.hermpoly_rec(20^2, 20)[1] ≈ 0.20019391063012504
        @test FastGaussQuadrature.hermpoly_rec(60^2, 60)[1] ≈ 0.07918022667865038
        @test FastGaussQuadrature.hermpoly_rec(60^2, 100)[1] ≈ -0.14191347555044895
        @test FastGaussQuadrature.hermpoly_rec(60^2, 130)[1] ≈ 1.7334093012299562e-73
        @test FastGaussQuadrature.hermpoly_rec(60^2, 200)[1] ≈ 0.0

        @test FastGaussQuadrature.hermpoly_rec(1:9, x) ≈ [FastGaussQuadrature.hermpoly_rec(n, x)[1] for n in 1:9]
        @test FastGaussQuadrature.hermpoly_rec(0:(20^2), 20)[end] ≈ FastGaussQuadrature.hermpoly_rec(20^2, 20)[1]
    end

    @testset "Type inference" begin
        @inferred gausshermite(Float32, 10)
        @inferred gausshermite(BigFloat, 10)
        @inferred gausshermite(BigFloat, 250)
    end

    @testset "Float16 output from Float32 arithmetic" begin
        # The scaled helper has its own conversion path; its tail weights
        # remain representable after the physical Float16 weights underflow.
        x32, w32 = FastGaussQuadrature.unweightedgausshermite(Float32, 21)
        x, w = FastGaussQuadrature.unweightedgausshermite(Float16, 21)
        @test x isa Vector{Float16} && w isa Vector{Float16}
        @test (x, w) == (Float16.(x32), Float16.(w32))
        xr, wr = QuadratureReference.nodeweight(:hermite, 21, 0, 0, x[1]; bits = 128)
        setprecision(BigFloat, 256) do
            delta = BigFloat(2)^(-128) * abs(xr)
            @test QuadratureReference.roots_below(:hermite, 21, BigFloat(0), BigFloat(0), xr - delta) == 0
            @test QuadratureReference.roots_below(:hermite, 21, BigFloat(0), BigFloat(0), xr + delta) == 1
        end
        @test x[1] == Float16(xr)
        @test w[1] == Float16(wr * exp(xr^2))
        for normalize in (false, true)
            x, w = gausshermite(Float16, 10; normalize)
            seed = normalize ? Float64(x[2]) / sqrt(2) : Float64(x[2])
            xr, wr = QuadratureReference.nodeweight(:hermite, 10, 0, 0, seed; bits = 128)
            @test x[2] == Float16(normalize ? sqrt(big(2)) * xr : xr)
            @test w[2] == Float16(normalize ? wr / sqrt(big(π)) : wr)
        end
        @testset "n: $n" for n in (2, 3, 20, 128)
            x, w = gausshermite(Float16, n)
            xr, wr = gausshermite(BigFloat, n)
            @test eltype(x) === eltype(w) === Float16
            @test x == -reverse(x) && w == reverse(w)
            @test isapprox(x, xr; atol = 24eps(Float16), rtol = 24eps(Float16))
            @test all(wr[i] < floatmin(Float16) || isapprox(w[i], wr[i]; rtol = 24eps(Float16)) for i in eachindex(w))
        end
        @test_throws ErrorException gausshermite(Float16, 5000)
        # Float32 normalization retains distinct nodes where native Float16
        # intermediate rounding previously merged them.
        xn, wn = gausshermite(Float16, 2350; normalize = true)
        @test all(>(0), diff(xn)) && all(isfinite, xn) && all(isfinite, wn)
        @test sum(Float64.(wn)) ≈ 1 rtol = 2eps(Float16)
        @test dot(Float64.(wn), Float64.(xn) .^ 2) ≈ 1 rtol = 2eps(Float16)
        for i in (1, 1175, 2350)
            xr, wr = QuadratureReference.nodeweight(:hermite, 2350, 0, 0, Float64(xn[i]) / sqrt(2); bits = 128)
            @test Float64(xn[i]) ≈ sqrt(big(2)) * xr rtol = eps(Float16)
            @test wn[i] == Float16(wr / sqrt(big(π)))
        end
        setprecision(BigFloat, 4) do
            x, w = gausshermite(BigFloat, 30)
            @test all(>(0), diff(x))
            @test all(t -> precision(t) == 4, x) && all(t -> precision(t) == 4, w)
            for i in (1, 15, 30)
                xr, wr = QuadratureReference.nodeweight(:hermite, 30, 0, 0, x[i]; bits = 128)
                @test x[i] == BigFloat(xr; precision = 4)
                @test w[i] == BigFloat(wr; precision = 4)
            end
        end
    end

    @testset "Normalize" begin
        for t in 1:6
            x, w = gausshermite(t; normalize = true)
            N = 2t - 1
            v = 1.0
            for n in 2:2:N
                v *= (n - 1)
                @test v ≈ dot(x .^ n, w)
            end
            for n in 1:2:N
                @test abs(0.0 - dot(x .^ n, w)) < 10 * tol
            end
        end
    end

    @testset "All" begin
        for n in 2:220
            x, w = gausshermite(n)
            @test isa(x, Vector{Float64})
            @test isa(w, Vector{Float64})
            @test (length(x) == n && length(w) == n)
            # The tol should be large for large n.
            @test (dot(w, x) < tol && abs(dot(w, x .^ 2) - sqrt(π) / 2) < 1000 * tol)
        end
    end

    @testset "Converged scaled weights and symmetry" begin
        # Independent 256-bit classical Hermite recurrence/Newton reference.
        x, w = gausshermite(42)
        @test isapprox(w[1], 6.1678589258107145e-31; rtol = 8.0e-14)
        @test x == -reverse(x)
        @test w == reverse(w)
    end

    @testset "Gil–Segura–Temme iteration" begin
        # Independent 384-bit classical-polynomial recurrence/Newton values.
        @testset "n: $n" for (n, reference) in (
                (
                    201, (
                        (1, -19.38970039958089, 3.1562517583901727e-164),
                        (102, 0.15649522104787772, 0.15271218047550833),
                    ),
                ),
                (251, ((1, -21.767664215263075, 9.505980773681771e-207),)),
                (10000, ((5001, 0.01110692968010143, 0.022211119186859443),)),
            )
            x, w = gausshermite(n)
            @testset "i: $i" for (i, node, weight) in reference
                @test isapprox(x[i], node; rtol = 4eps(Float64))
                @test isapprox(w[i], weight; rtol = 2.0e-13)
            end
            @test issorted(x)
            @test all(isfinite, x)
            @test all(≥(0), w)
            @test x == -reverse(x)
            @test w == reverse(w)
        end
        # Check high-precision tail and central weights against the independent
        # classical recurrence at two working precisions.
        @testset "bits: $bits, n: $n" for (bits, n) in ((128, 100), (256, 201))
            setprecision(bits) do
                x, w = gausshermite(BigFloat, n)
                for i in (1, n ÷ 2 + 2)
                    xr, wr = QuadratureReference.nodeweight(:hermite, n, 0, 0, x[i]; bits)
                    @test isapprox(x[i], xr; rtol = 20eps(BigFloat))
                    @test isapprox(w[i], wr; rtol = 1000eps(BigFloat))
                end
            end
        end
        x, scaled = FastGaussQuadrature.unweightedgausshermite(1000)
        xr, wr = QuadratureReference.nodeweight(:hermite, 1000, 0, 0, x[end]; bits = 256)
        @test iszero(last(gausshermite(1000)[2]))
        @test isapprox(scaled[end], wr * exp(xr^2); rtol = 2.0e-12)
    end

    @testset "Recurrence and GST across precisions" begin
        # Independent polynomial references cover both sides of the shared
        # cutoff, including Float32's squared-derivative overflow at n=21.
        @testset "Type: $T, n: $n" for T in (Float32, Float64, BigFloat), n in (2, 3, 19, 20, 21, 22)
            x, w = @inferred gausshermite(T, n)
            @test x isa Vector{T} && w isa Vector{T} && length(x) == length(w) == n
            @test x == -reverse(x) && w == reverse(w)
            @test all(isfinite, x) && all(isfinite, w) && all(>(0), w)
            @test isapprox(sum(w), sqrt(T(π)); rtol = 8n * eps(T))
            @test isapprox(dot(w, x .^ 2), sqrt(T(π)) / 2; rtol = 8n * eps(T))
            for i in (1, cld(n, 2))
                xr, wr = QuadratureReference.nodeweight(:hermite, n, 0, 0, x[i]; bits = max(128, precision(T)))
                # Root rounding affects exp(-x²), with relative error O(n*eps).
                @test isapprox(x[i], xr; atol = 4eps(T), rtol = 4eps(T))
                @test isapprox(w[i], wr; rtol = 8n * eps(T))
            end
        end
    end

    @testset "Transform" begin
        n = 500
        x, w = FastGaussQuadrature.unweightedgausshermite(n)
        @test w[1] ≠ 0

        V = Array{Float64}(undef, n, n)
        for k in 1:n
            V[k, :] = FastGaussQuadrature.hermpoly_rec(0:(n - 1), sqrt(2) * x[k])
        end

        f = x -> first(FastGaussQuadrature.hermpoly_rec(1, sqrt(2) * x))
        @test V' * (w .* f.(x)) ≈ [0; sqrt(π); zeros(n - 2)]
        f = x -> first(FastGaussQuadrature.hermpoly_rec(2, sqrt(2) * x))
        @test V' * (w .* f.(x)) ≈ [0; 0; sqrt(π); zeros(n - 3)]
        f = x -> first(FastGaussQuadrature.hermpoly_rec(3, sqrt(2) * x))
        @test V' * (w .* f.(x)) ≈ [0; 0; 0; sqrt(π); zeros(n - 4)]
    end
end

@testset "Small Hermite native precision and normalization" begin
    previousbits = precision(BigFloat)
    setprecision(BigFloat, 32) do
        n = 20
        x, w = gausshermite(BigFloat, n; normalize = true)
        @test all(v -> precision(v) == 32, x) && all(v -> precision(v) == 32, w)
        @test precision(BigFloat) == 32
        e = eps(BigFloat)
        for i in (1, n ÷ 2)
            xr, wr = setprecision(BigFloat, 256) do
                seed = BigFloat(x[i]) / sqrt(big(2))
                r, v = QuadratureReference.nodeweight(:hermite, n, 0, 0, seed; bits = 128)
                sqrt(big(2)) * r, v / sqrt(big(π))
            end
            @test abs(x[i] - xr) <= 8e * max(one(xr), abs(xr))
            @test abs(w[i] - wr) <= 8n * e * wr
        end
    end
    @test precision(BigFloat) == previousbits
end

@testset "Hermite physical-weight transport" begin
    x, w = gausshermite(Float32, 128)
    xs, ws = FastGaussQuadrature.unweightedgausshermite(Float32, 128)
    @test x == xs
    xr, wr = QuadratureReference.nodeweight(:hermite, 128, 0, 0, x[22]; bits = 128)
    # A normal tail weight: preserve the native scaled amplitude while removing
    # the avoidable sensitivity to rounding its output node before exp(-x^2).
    @test isapprox(w[22], wr; rtol = 64eps(Float32))
    @test isapprox(ws[22], wr * exp(xr^2); rtol = 64eps(Float32))
end
