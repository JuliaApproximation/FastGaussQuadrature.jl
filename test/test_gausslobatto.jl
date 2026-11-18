@testset "Gauss–Lobatto" begin
    # Check error
    @test_throws DomainError gausslobatto(-1)
    @test_throws DomainError gausslobatto(0)
    @test_throws DomainError gausslobatto(1)

    # Type stability
    @inferred gausslobatto(Float32, 10)
    @inferred gausslobatto(BigFloat, 10)

    n = 2
    x, w = gausslobatto(n)
    @test x[1] ≈ -1.0
    @test x[2] ≈ 1.0
    @test w[1] ≈ 1.0
    @test w[2] ≈ 1.0

    n = 3
    x, w = gausslobatto(n)
    @test x[1] ≈ -1.0
    @test abs(x[2]) < 1.0e-15
    @test x[3] ≈ 1.0
    @test w[1] ≈ 1 / 3
    @test w[2] ≈ 4 / 3
    @test w[3] ≈ 1 / 3

    tol = 1.0e-14
    n = 42
    x, w = gausslobatto(n)
    @test (length(x) == n) && (length(w) == n)
    @test abs(x[37] - 0.922259214258616) < tol
    @test abs(w[37] - 0.029306411216166) < tol
    @test (x[1] == -1 && x[n] == 1)

    @test dot(w, (x .^ 2)) ≈ 2 / 3
    @test dot(w, exp.(x)) ≈ exp(1) - exp(-1)

    @testset "endpoint weights retain their angles" begin
        # Independent 256-bit Jacobi-root/Legendre-weight references.
        @testset "n: $n" for (n, expected) in (
                (103, 0.0011732757161776752),
                (10000, 1.2330538839742197e-7),
            )
            x, w = gausslobatto(n)
            @test x[1] == -1 && x[end] == 1
            @test x == -reverse(x)
            @test w == reverse(w)
            @test w[2] ≈ expected rtol = 1.0e-14
        end
    end
end

@testset "Low-precision Lobatto rules" begin
    n = 128
    x, w = @inferred gausslobatto(Float16, n)
    @test eltype(x) === eltype(w) === Float16
    @test x[1] == -1 && x[end] == 1 && all(>(0), diff(x))
    @test x == -reverse(x) && w == reverse(w)
    @test w[1] == Float16(BigFloat(2) / (n * (n - 1)))
    for i in (2, cld(n, 2), n - 1)
        xr, wr = QuadratureReference.nodeweight(:jacobi, n - 2, 1, 1, x[i]; bits = 128)
        @test isapprox(x[i], xr; atol = eps(Float16), rtol = 0)
        @test isapprox(w[i], wr / ((1 - xr) * (1 + xr)); atol = eps(zero(Float16)), rtol = 2eps(Float16))
    end
    @test_throws ErrorException gausslobatto(Float16, 512)
end
