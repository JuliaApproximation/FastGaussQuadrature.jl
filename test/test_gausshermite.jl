using FastGaussQuadrature

@testset "Gauss–Hermite" begin
    tol = 1e-14

    @testset "Golub-Welsch" begin
        n = 18;
        x,w = gausshermite( n )
        @test isa(x,Vector{Float64})
        @test isa(w,Vector{Float64})
        @test (length(x) == n && length(w) == n)
        @test (dot(w,x) < tol && abs(dot(w,x.^2) - sqrt(pi)/2) < tol)
    end

    @testset "Recurrence" begin
        n = 42
        x,w = gausshermite( n )
        @test isa(x,Vector{Float64})
        @test isa(w,Vector{Float64})
        @test (length(x) == n && length(w) == n)
        @test (dot(w,x) < tol && abs(dot(w,x.^2) - sqrt(pi)/2) < tol)
        @test abs(x[37] - 5.660357581283058) < tol
        @test abs(w[17] - 0.032202101288908) < tol

        f = x -> 1 + x + x^2 + x^3
        @test w'f.(x) ≈ 2.6586807763789717
    end
    
    @testset "Asymptotics"  begin
        n = 251
        x,w = gausshermite( n )
        @test isa(x,Vector{Float64})
        @test isa(w,Vector{Float64})
        @test (length(x) == n && length(w) == n)
        @test (dot(w,x) < tol && abs(dot(w,x.^2) - sqrt(pi)/2) < 300*tol)
        @test abs(x[37] - -13.292221459334638) < tol
        @test abs(w[123] - 0.117419270715955) < 2*tol

        n = 500
        x,w = gausshermite( n )
        @test isa(x,Vector{Float64})
        @test isa(w,Vector{Float64})
        @test (length(x) == n && length(w) == n)
        @test (dot(w,x) < tol && abs(dot(w,x.^2) - sqrt(pi)/2) < 300*tol)
    end

    @testset "Recurrence" begin
        x = 0.1
        @test FastGaussQuadrature.hermpoly_rec(1,x)[1] ≈ (2x*exp(-x^2/4))/2
        FastGaussQuadrature.hermpoly_rec(2,x)
    end

    @testset "Unweighted" begin
        n = 500
        x,w = FastGaussQuadrature.unweightedgausshermite( n )
        @test w[1] ≠ 0
    end
end