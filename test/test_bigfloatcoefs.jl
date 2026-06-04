using LinearAlgebra, GenericLinearAlgebra

@testset "Gauss–Legendre" begin
    @testset "quad: $quad" for quad in (gausslegendre, gaussradau)
        @testset "n: $n" for n in (2, 3, 4, 5, 6, 7, 8, 9, 40, 100)
            # check BigFloat
            x, w = @inferred quad(n)
            xbig, wbig = @inferred quad(BigFloat, n)
            @test length(xbig) == n && length(wbig) == n
            @test ≈(x, xbig; atol = 1.0e-14)
            @test ≈(w, wbig; atol = 1.0e-14)
            # tests bigfloat coefs by integrating against random polynomial
            poly_coefs = rand(BigFloat, 2n - 1)
            poly_coefs = FastGaussQuadrature.jacobi_gw(2n - 1, big(0.0), big(0.0))[1]
            test_fun(x) = evalpoly(x, poly_coefs)
            # analytic integral of f
            integral_testfun(x) = x * evalpoly(x, [p / i for (i, p) in enumerate(poly_coefs)])
            @test ≈(dot(wbig, test_fun.(xbig)), integral_testfun(1) - integral_testfun(-1); atol = 1.0e-70)
        end
    end
end

@testset "Gauss–Hermite" begin
    # n values cover: special cases (1), GW (≤20), rec (21–200), asy→rec for BigFloat (>200)
    @testset "n: $n" for n in (1, 2, 10, 21, 100, 250)
        x, w = @inferred gausshermite(n)
        xbig, wbig = @inferred gausshermite(BigFloat, n)
        @test length(xbig) == n && length(wbig) == n
        @test ≈(x, xbig; atol = n ≤ 200 ? 1.0e-14 : 1.0e-13)
        @test ≈(w, wbig; atol = n ≤ 200 ? 1.0e-14 : 1.0e-12)
        # ∫ x^2 exp(-x^2) dx = sqrt(π)/2, exact for n ≥ 2
        n ≥ 2 && @test ≈(dot(wbig, xbig .^ 2), sqrt(big(π)) / 2; atol = n == 21 ? 1.0e-50 : 1.0e-70)
    end
end

@testset "Gauss–Laguerre" begin
    # n values cover: special cases (1,2), GW (<15), rec (15–127), asy→rec for BigFloat (≥128)
    @testset "n: $n" for n in (1, 2, 3, 10, 20, 130)
        x, w = @inferred gausslaguerre(n)
        xbig, wbig = @inferred gausslaguerre(BigFloat, n)
        @test length(xbig) == n && length(wbig) == n
        @test ≈(x, xbig; atol = n < 128 ? 1.0e-13 : 1.0e-12)
        @test ≈(w, wbig; atol = 1.0e-13)
        # ∫ x exp(-x) dx = 1, exact for n ≥ 2
        n ≥ 2 && @test ≈(dot(wbig, xbig), big(1); atol = 1.0e-70)
    end

    # non-zero α: type inferred from α, or set explicitly via T
    @testset "α ≠ 0, n: $n" for n in (3, 10, 20, 130)
        α = 0.5
        x, w = @inferred gausslaguerre(n, α)
        xbig, wbig = @inferred gausslaguerre(n, big(α))
        @test ≈(x, xbig; atol = 1.0e-12)
        @test ≈(w, wbig; atol = 1.0e-14)
        # ∫ x * x^α exp(-x) dx = Γ(α+2), exact for n ≥ 1
        @test ≈(dot(wbig, xbig), gamma(big(α) + 2); atol = 1.0e-70)
    end
end

@testset "Gauss–Lobatto" begin
    @testset "n: $n" for n in (3, 4, 5, 6, 7, 8, 9, 40, 100)
        # check BigFloat
        x, w = @inferred gausslobatto(n)
        xbig, wbig = @inferred gausslobatto(BigFloat, n)
        @test length(xbig) == n && length(wbig) == n
        @test ≈(x, xbig; atol = 1.0e-14)
        @test ≈(w, wbig; atol = 1.0e-14)
        # tests bigfloat coefs by integrating against a polynomial of degree 2n-3
        # (exact for Lobatto with n points, which integrates degree 2n-3 exactly)
        poly_coefs = FastGaussQuadrature.jacobi_gw(2n - 3, big(0.0), big(0.0))[1]
        test_fun(x) = evalpoly(x, poly_coefs)
        integral_testfun(x) = x * evalpoly(x, [p / i for (i, p) in enumerate(poly_coefs)])
        @test ≈(dot(wbig, test_fun.(xbig)), integral_testfun(1) - integral_testfun(-1); atol = 1.0e-70)
    end
end
