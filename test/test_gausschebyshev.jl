@testset "Gauss–Chebyshev" begin

    @testset "Check error" begin
        @test_throws DomainError gausschebyshev(-1,1)
        @test_throws DomainError gausschebyshev(-1,2)
        @test_throws DomainError gausschebyshev(-1,3)
        @test_throws DomainError gausschebyshev(-1,4)
        @test_throws ArgumentError gausschebyshev(0,5)
    end

    n = 10

    @testset "x.^2" begin
        x, w = gausschebyshev(n,1)
        @test dot(x.^2,w) ≈ π/2
        x, w = gausschebyshev(n,2)
        @test dot(x.^2,w) ≈ π/8
        x, w = gausschebyshev(n,3)
        @test dot(x.^2,w) ≈ π/2
        x, w = gausschebyshev(n,4)
        @test dot(x.^2,w) ≈ π/2
    end

    @testset "x^3" begin
        x, w = gausschebyshev(n,1)
        @test abs(dot(x.^3,w))<1e-15
        x, w = gausschebyshev(n,2)
        @test abs(dot(x.^3,w))<1e-15
        x, w = gausschebyshev(n,3)
        @test dot(x.^3,w) ≈ 3*π/8
        x, w = gausschebyshev(n,4)
        @test dot(x.^3,w) ≈ -3*π/8
    end
end