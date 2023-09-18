@testset "Gauss–Chebyshev" begin

    @testset "Check error" begin
        @test_throws DomainError gausschebyshevt(-1)
        @test_throws DomainError gausschebyshevu(-1)
        @test_throws DomainError gausschebyshevv(-1)
        @test_throws DomainError gausschebyshevw(-1)
    end

    n = 10

    @testset "x.^2" begin
        x, w = gausschebyshevt(n)
        @test dot(x.^2,w) ≈ π/2
        x, w = gausschebyshevu(n)
        @test dot(x.^2,w) ≈ π/8
        x, w = gausschebyshevv(n)
        @test dot(x.^2,w) ≈ π/2
        x, w = gausschebyshevw(n)
        @test dot(x.^2,w) ≈ π/2
    end

    @testset "x^3" begin
        x, w = gausschebyshevt(n)
        @test abs(dot(x.^3,w)) < 1e-15
        x, w = gausschebyshevu(n)
        @test abs(dot(x.^3,w)) < 1e-15
        x, w = gausschebyshevv(n)
        @test dot(x.^3,w) ≈ 3*π/8
        x, w = gausschebyshevw(n)
        @test dot(x.^3,w) ≈ -3*π/8
    end

    @testset "deprecated" begin
        n = 42
        @test gausschebyshevt(n) == gausschebyshev(n, 1)
        @test gausschebyshevu(n) == gausschebyshev(n, 2)
        @test gausschebyshevv(n) == gausschebyshev(n, 3)
        @test gausschebyshevw(n) == gausschebyshev(n, 4)
        @test_throws ArgumentError gausschebyshev(0,5)
    end
end
