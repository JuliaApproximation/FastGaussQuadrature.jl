using FastGaussQuadrature, Test

@testset "Gauss–Radau" begin
    @testset "Legendre" begin
        n = 1
        x, w = gaussradau(n)
        @test x[1] ≈ -1.
        @test w[1] ≈ 2.

        n = 2
        x, w = gaussradau(n)
        @test x[1] ≈ -1.
        @test x[2] ≈ 1/3
        @test w[1] ≈ .5
        @test w[2] ≈ 1.5

        tol = 1e-14
        n = 42
        ntests = 4
        x, w = gaussradau(n)
        @test length(x) == n && length(w) == n
        @test abs(x[37] - 0.908847278001044) < tol
        @test abs(w[37] - 0.031190846817016) < tol
        @test x[1] == -1

        @test dot( w,exp.(x)) ≈ exp(1)-exp(-1)
    end
    @testset "Jacobi" begin
        a,b = 0.1, 0.2
        # compare with Gauss–Jacobi
        for n = 1:5
            x, w = gaussradau(n, a, b)
            x̃, w̃ = gaussjacobi(n, a, b)
            @test sum(w) ≈ sum(w̃)
            @test dot(w,x) ≈ dot(w̃,x̃)
            @test dot(w,x.^2) ≈ dot(w̃,x̃.^2)
        end
    end
end