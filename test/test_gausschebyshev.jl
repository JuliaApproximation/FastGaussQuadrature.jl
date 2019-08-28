@testset "Gauss–Chebyshev" begin

    n = 10

    @testset "x.^2" begin
        x, w = gausschebyshev(n,1)
        @test dot(x.^2,w) ≈ pi/2
        x, w = gausschebyshev(n,2)
        @test dot(x.^2,w) ≈ pi/8
        x, w = gausschebyshev(n,3)
        @test dot(x.^2,w) ≈ pi/2
        x, w = gausschebyshev(n,4)
        @test dot(x.^2,w) ≈ pi/2
    end

    @testset "x^3" begin
        x, w = gausschebyshev(n,1)
        @test abs(dot(x.^3,w))<1e-15
        x, w = gausschebyshev(n,2)
        @test abs(dot(x.^3,w))<1e-15
        x, w = gausschebyshev(n,3)
        @test dot(x.^3,w) ≈ 3*pi/8
        x, w = gausschebyshev(n,4)
        @test dot(x.^3,w) ≈ -3*pi/8
    end
end