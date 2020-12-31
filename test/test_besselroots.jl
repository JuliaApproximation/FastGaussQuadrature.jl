import SpecialFunctions
import SpecialFunctions: besselj

@testset "Bessel Roots" begin
    @test_throws DomainError besselroots(0.0, -1)
 
    tol = 1e-11

    # Check if besselj(nu, besselroots(nu, n) ) is small:

    for nu = 0.:0.1:5.
        n = 10
        @test norm( besselj.(nu, besselroots(nu, n) ), Inf ) < tol
    end
end
