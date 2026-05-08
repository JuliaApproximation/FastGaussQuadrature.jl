using LinearAlgebra, GenericLinearAlgebra

@testset "Gauss–Legendre" begin
    @testset "quad: $quad" for quad in (gausslegendre, gaussradau)
        @testset "n: $n" for n in (2, 3, 4, 5, 6, 7, 8, 9, 40, 100)
            # check BigFloat
            x, w = quad(n)
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
