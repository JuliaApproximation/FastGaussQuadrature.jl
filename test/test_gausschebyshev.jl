@testset "Gauss–Chebyshev" begin

    @testset "Check error" begin
        @test_throws DomainError gausschebyshev1(-1)
        @test_throws DomainError gausschebyshev2(-1)
        @test_throws DomainError gausschebyshev3(-1)
        @test_throws DomainError gausschebyshev4(-1)
    end

    n = 10

    @testset "x.^2" begin
        x, w = gausschebyshev1(n)
        @test dot(x.^2,w) ≈ π/2
        x, w = gausschebyshev2(n)
        @test dot(x.^2,w) ≈ π/8
        x, w = gausschebyshev3(n)
        @test dot(x.^2,w) ≈ π/2
        x, w = gausschebyshev4(n)
        @test dot(x.^2,w) ≈ π/2
    end

    @testset "x^3" begin
        x, w = gausschebyshev1(n)
        @test abs(dot(x.^3,w)) < 1e-15
        x, w = gausschebyshev2(n)
        @test abs(dot(x.^3,w)) < 1e-15
        x, w = gausschebyshev3(n)
        @test dot(x.^3,w) ≈ 3*π/8
        x, w = gausschebyshev4(n)
        @test dot(x.^3,w) ≈ -3*π/8
    end

    @testset "deprecated" begin
        n = 42
        @test gausschebyshev1(n) == gausschebyshev(n, 1)
        @test gausschebyshev2(n) == gausschebyshev(n, 2)
        @test gausschebyshev3(n) == gausschebyshev(n, 3)
        @test gausschebyshev4(n) == gausschebyshev(n, 4)
        @test_throws ArgumentError gausschebyshev(0,5)
    end
end
