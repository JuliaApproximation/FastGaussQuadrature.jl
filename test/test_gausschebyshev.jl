@testset "Gauss–Chebyshev" begin

    @testset "Check error" begin
        @test_throws DomainError gausschebyshevT(-1)
        @test_throws DomainError gausschebyshevU(-1)
        @test_throws DomainError gausschebyshevV(-1)
        @test_throws DomainError gausschebyshevW(-1)
    end

    n = 10

    @testset "x.^2" begin
        x, w = gausschebyshevT(n)
        @test dot(x.^2,w) ≈ π/2
        x, w = gausschebyshevU(n)
        @test dot(x.^2,w) ≈ π/8
        x, w = gausschebyshevV(n)
        @test dot(x.^2,w) ≈ π/2
        x, w = gausschebyshevW(n)
        @test dot(x.^2,w) ≈ π/2
    end

    @testset "x^3" begin
        x, w = gausschebyshevT(n)
        @test abs(dot(x.^3,w)) < 1e-15
        x, w = gausschebyshevU(n)
        @test abs(dot(x.^3,w)) < 1e-15
        x, w = gausschebyshevV(n)
        @test dot(x.^3,w) ≈ 3*π/8
        x, w = gausschebyshevW(n)
        @test dot(x.^3,w) ≈ -3*π/8
    end

    @testset "deprecated" begin
        n = 42
        @test gausschebyshevT(n) == gausschebyshev(n, 1)
        @test gausschebyshevU(n) == gausschebyshev(n, 2)
        @test gausschebyshevV(n) == gausschebyshev(n, 3)
        @test gausschebyshevW(n) == gausschebyshev(n, 4)
        @test_throws ArgumentError gausschebyshev(0,5)
    end
end
