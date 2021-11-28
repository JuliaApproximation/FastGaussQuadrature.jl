using FastGaussQuadrature

@testset "Gauss–Hermite" begin
    @testset "Check error" begin
        @test_throws DomainError gausshermite(-1)
    end
    
    @testset "Small n case" begin
        @test gausshermite(0) == (Float64[], Float64[])
        x, w = gausshermite(1)
        @test x ≈ [0.0]
        @test w ≈ [√π]
    end

    tol = 1e-14

    @testset "Golub-Welsch" begin
        n = 18;
        x,w = gausshermite(n)
        @test isa(x,Vector{Float64})
        @test isa(w,Vector{Float64})
        @test (length(x) == n && length(w) == n)
        @test (dot(w,x) < tol && abs(dot(w,x.^2) - sqrt(π)/2) < tol)
    end

    @testset "Recurrence" begin
        n = 42
        x,w = gausshermite(n)
        @test isa(x,Vector{Float64})
        @test isa(w,Vector{Float64})
        @test (length(x) == n && length(w) == n)
        @test (dot(w,x) < tol && abs(dot(w,x.^2) - sqrt(π)/2) < tol)
        @test abs(x[37] - 5.660357581283058) < tol
        @test abs(w[17] - 0.032202101288908) < tol

        f = x -> 1 + x + x^2 + x^3
        @test w'f.(x) ≈ 2.6586807763789717
    end

    @testset "Asymptotics"  begin
        n = 251
        x,w = gausshermite(n)
        @test isa(x,Vector{Float64})
        @test isa(w,Vector{Float64})
        @test (length(x) == n && length(w) == n)
        @test (dot(w,x) < tol && abs(dot(w,x.^2) - sqrt(π)/2) < 300*tol)
        @test abs(x[37] - -13.292221459334638) < tol
        @test abs(w[123] - 0.117419270715955) < 2*tol

        n = 500
        x,w = gausshermite(n)
        @test isa(x,Vector{Float64})
        @test isa(w,Vector{Float64})
        @test (length(x) == n && length(w) == n)
        @test (dot(w,x) < tol && abs(dot(w,x.^2) - sqrt(π)/2) < 300*tol)
    end

    @testset "Recurrence" begin
        x = 0.1 
        He0 = 1 # Probabilists He_n(x) = H_n(x/sqrt(2))/2^(n/2)
        He1 = x; He2 = x^2-1; He3 = x*(x^2-3); He4 = 3 - 6x^2 + x^4
        @test FastGaussQuadrature.hermpoly_rec(1,x)[1] ≈ He1*exp(-x^2/4)
        @test FastGaussQuadrature.hermpoly_rec(2,x)[1] ≈ He2*exp(-x^2/4)/sqrt(2)
        @test FastGaussQuadrature.hermpoly_rec(3,x)[1] ≈ He3*exp(-x^2/4)/sqrt(2*3)
        @test FastGaussQuadrature.hermpoly_rec(4,x)[1] ≈ He4*exp(-x^2/4)/sqrt(2*3*4)
        
        @test FastGaussQuadrature.hermpoly_rec(1,60)[1] ≈ 0.0
        @test FastGaussQuadrature.hermpoly_rec(20^2,20)[1] ≈ 0.20019391063012504
        @test FastGaussQuadrature.hermpoly_rec(60^2,60)[1] ≈ 0.07918022667865038
        @test FastGaussQuadrature.hermpoly_rec(60^2,100)[1] ≈ -0.14191347555044895
        @test FastGaussQuadrature.hermpoly_rec(60^2,130)[1] ≈ 1.7334093012299562E-73
        @test FastGaussQuadrature.hermpoly_rec(60^2,200)[1] ≈ 0.0

        @test FastGaussQuadrature.hermpoly_rec(1:9,x) ≈ [FastGaussQuadrature.hermpoly_rec(n,x)[1] for n=1:9]
        @test FastGaussQuadrature.hermpoly_rec(0:20^2,20)[end] ≈ FastGaussQuadrature.hermpoly_rec(20^2,20)[1]
    end

    @testset "All" begin
        for n in 2:220
            x,w = gausshermite(n)
            @test isa(x,Vector{Float64})
            @test isa(w,Vector{Float64})
            @test (length(x) == n && length(w) == n)
            # The tol should be large for large n.
            @test (dot(w,x) < tol && abs(dot(w,x.^2) - sqrt(π)/2) < 1000*tol)
        end
    end

    @testset "Transform" begin
        n = 500
        x,w = FastGaussQuadrature.unweightedgausshermite(n)
        @test w[1] ≠ 0

        V = Array{Float64}(undef,n,n)
        for k=1:n
            V[k,:] = FastGaussQuadrature.hermpoly_rec(0:n-1, sqrt(2)*x[k])
        end

        f = x -> first(FastGaussQuadrature.hermpoly_rec(1, sqrt(2)*x))
        @test V' * (w.* f.(x)) ≈ [0; sqrt(π); zeros(n-2)]
        f = x -> first(FastGaussQuadrature.hermpoly_rec(2, sqrt(2)*x))
        @test V' * (w.* f.(x)) ≈ [0; 0; sqrt(π); zeros(n-3)]
        f = x -> first(FastGaussQuadrature.hermpoly_rec(3, sqrt(2)*x))
        @test V' * (w.* f.(x)) ≈ [0; 0; 0; sqrt(π); zeros(n-4)]
    end
end
