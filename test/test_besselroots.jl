import SpecialFunctions
import SpecialFunctions: besselj

@testset "Bessel Roots" begin
    @test_throws DomainError approx_besselroots(0.0, -1)
 
    tol = 1e-11

    # Check if besselj(ν, approx_besselroots(ν, n) ) is small:

    for ν = 0.:0.1:5.
        n = 10
        @test norm( besselj.(ν, approx_besselroots(ν, n) ), Inf ) < tol
    end
end
