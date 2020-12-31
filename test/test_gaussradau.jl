using FastGaussQuadrature, Test

@testset "Gauss–Radau" begin
    @testset "Check error" begin
        @test_throws DomainError gaussradau(0)
        @test_throws DomainError gaussradau(0,0,0)
    end

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
        for (a,b) in ((0.1, 0.2), (0,-0.5))
            # compare with Gauss–Jacobi
            for n = 1:5
                x, w = gaussradau(n, a, b)
                x̃, w̃ = gaussjacobi(n, a, b)

                @test length(x) == n

                @test sum(w) ≈ sum(w̃)
                for j = 1:2n-2
                    @test dot(w,x.^j) ≈ dot(w̃,x̃.^j)
                end
            end
        end

        # compare with Gauss-Radau-Legendre
        @test all(gaussradau(2, 0, 0) .≈ gaussradau(2))
        @test all(gaussradau(3, 0, 0) .≈ gaussradau(3))
    end
end