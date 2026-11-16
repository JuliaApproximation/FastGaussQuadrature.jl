@testset "Gauss–Legendre" begin

    # Check error
    @test_throws DomainError gausslegendre(-1)

    # Check all special cases
    @test gausslegendre(0) == (Float64[], Float64[])
    for ν in 0:10
        (x, w) = gausslegendre(ν)
        for i in 0:(ν - 1)
            r = 2i
            @test abs(dot(w, (x .^ r)) - 2 / (r + 1)) < 10 * eps()
        end
    end

    tol = 1.0e-14
    n = 42
    x, w = @inferred gausslegendre(n)
    @test length(x) == n && length(w) == n
    @test ≈(x[37], 0.910959724904127; atol = tol)
    @test ≈(w[37], 0.030479240699603; atol = tol)

    @test dot(w, (x .^ 2)) ≈ 2 / 3
    @test dot(w, exp.(x)) ≈ exp(1) - exp(-1)

    # Test a larger n (using ASY)
    n = 251
    x, w = gausslegendre(n)
    @test all(length(x) == n) && all(length(w) == n)
    @test abs(x[37] + 0.896467746955729) < tol
    @test abs(w[37] - 0.005535005742012) < tol

    @test dot(w, (x .^ 2)) ≈ 2 / 3
    @test dot(w, exp.(x)) ≈ exp(1) - exp(-1)

    x, w = gausslegendre(1013)
    @test norm(x[2] - -0.99998516758611, Inf) < tol
    @test norm(x[13] - -0.999218995240887, Inf) < tol
    @test norm(w[2] - 1.681691163200592e-5, Inf) < tol
    @test norm(w[13] - 1.224755309137936e-4, Inf) < tol

    @test dot(w, (x .^ 2)) ≈ 2 / 3
    @test dot(w, exp.(x)) ≈ exp(1) - exp(-1)

    x, w = gausslegendre(10013)
    @test norm(x[2] - -0.999999848054223, Inf) < tol
    @test norm(x[13] - -0.999991998242661, Inf) < tol
    @test norm(w[2] - 1.722757320118474e-7, Inf) < tol
    @test norm(w[13] - 1.25498054003247e-6, Inf) < tol

    @test dot(w, (x .^ 2)) ≈ 2 / 3
    @test dot(w, exp.(x)) ≈ exp(1) - exp(-1)
end

@testset "Small Legendre rules at generic precision" begin
    @testset "$T, $bits bits, n: $n" for (T, bits, n) in (
            (Float32, 24, 6), (Float32, 24, 33), (Float32, 24, 60),
            (BigFloat, 128, 6), (BigFloat, 256, 33), (BigFloat, 1024, 60),
        )
        setprecision(BigFloat, bits) do
            x, w = @inferred gausslegendre(T, n)
            @test eltype(x) == eltype(w) == T
            @test all(z -> precision(z) == bits, x)
            @test all(z -> precision(z) == bits, w)
            @test x == -reverse(x) && w == reverse(w)
            @test all(>(0), diff(x)) && all(>(0), w)
            @testset "i: $i" for i in (1, 2, n ÷ 2)
                xr, wr = QuadratureReference.nodeweight(:jacobi, n, 0, 0, x[i]; bits = max(bits, 128))
                @test isapprox(x[i], xr; atol = 0, rtol = 8eps(T))
                # Endpoint weight sensitivity grows with n² in this recurrence.
                @test isapprox(w[i], wr; atol = 0, rtol = 8n^2 * eps(T))
            end
        end
    end
end

@testset "Legendre central phases and large-order weights" begin
    # Independent 192-bit polynomial/derivative recurrence and converged Newton
    # references.
    x, w = gausslegendre(61)
    @test x == -reverse(x)
    @test w == reverse(w)
    @test iszero(x[31])
    @test x[30] ≈ -0.05105890670797435 rtol = 5eps(Float64)
    @test w[30] ≈ 0.051014487038697265 rtol = 5eps(Float64)

    x, w = gausslegendre(1000000)
    @test x[500000] ≈ -1.5707955413962836e-6 rtol = 5eps(Float64)
    @test w[500000] ≈ 3.1415910827899833e-6 rtol = 5eps(Float64)
    @test w[1] ≈ 7.420753950655386e-12 rtol = 5eps(Float64)
end

@testset "Large Legendre rules in Float32" begin
    @testset "n: $n" for n in (61, 62, 1000)
        x, w = @inferred gausslegendre(Float32, n)
        @test eltype(x) == eltype(w) == Float32
        @test x == -reverse(x)
        @test w == reverse(w)
        @test all(isfinite, x) && all(>(0), diff(x))
        @test all(isfinite, w) && all(>(0), w)
        @test sum(w) ≈ 2 rtol = 10eps(Float32)
        @testset "i: $i" for i in (1, 2, cld(n, 2))
            xr, wr = QuadratureReference.nodeweight(:jacobi, n, 0, 0, x[i]; bits = 128)
            @test isapprox(x[i], xr; atol = 0, rtol = eps(Float32))
            @test isapprox(w[i], wr; atol = 0, rtol = eps(Float32))
        end
    end
    # Rounding cannot represent distinct endpoint nodes at sufficiently high order.
    @test_throws ErrorException gausslegendre(Float32, 100000)
end

@testset "Legendre asymptotics at low BigFloat precision" begin
    previousbits = precision(BigFloat)
    @testset "$bits bits, n: $n" for bits in (24, 53), n in (61, 1000)
        setprecision(BigFloat, bits) do
            x, w = @inferred gausslegendre(BigFloat, n)
            @test eltype(x) == eltype(w) == BigFloat
            @test all(z -> precision(z) == bits, x)
            @test all(z -> precision(z) == bits, w)
            @test x == -reverse(x)
            @test w == reverse(w)
            @test all(>(0), diff(x)) && all(>(0), w)
            @testset "i: $i" for i in (1, 2, cld(n, 2))
                xr, wr = QuadratureReference.nodeweight(:jacobi, n, 0, 0, x[i]; bits = 128)
                @test isapprox(x[i], xr; atol = 0, rtol = 4eps(BigFloat))
                @test isapprox(w[i], wr; atol = 0, rtol = 4eps(BigFloat))
            end
            @test precision(BigFloat) == bits
        end
    end
    @test precision(BigFloat) == previousbits

    # Just above the cutoff, the central node needs more than Float64 accuracy.
    setprecision(BigFloat, 54) do
        x, _ = gausslegendre(BigFloat, 1000)
        xr, _ = QuadratureReference.nodeweight(:jacobi, 1000, 0, 0, x[500]; bits = 128)
        @test isapprox(x[500], xr; atol = 0, rtol = eps(BigFloat))
    end
end

@testset "Legendre rules in Float16" begin
    # Wider arithmetic preserves endpoint weights before rounding the rule.
    @testset "n: $n" for n in (6, 33, 60, 200)
        x, w = @inferred gausslegendre(Float16, n)
        xr, wr = gausslegendre(n)
        @test eltype(x) == eltype(w) == Float16
        @test all(isfinite, x) && all(>(0), diff(x))
        @test all(isfinite, w) && all(>(0), w)
        @test maximum(abs.(Float64.(x) .- xr)) < 2eps(Float16)
        @test maximum(abs.(Float64.(w) ./ wr .- 1)) < 4eps(Float16)
    end
    @test_throws ErrorException gausslegendre(Float16, 226)
end
