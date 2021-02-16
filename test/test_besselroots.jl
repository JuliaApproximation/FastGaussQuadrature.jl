import SpecialFunctions
import SpecialFunctions: besselj

@testset "Bessel Roots" begin
    @test_throws DomainError approx_besselroots(0.0, -1)

    msg = "`besselroots` has been renamed to `approx_besselroots`, and `besselroots` will be removed in the next major release."
    @test_logs (:warn, msg) besselroots(0.2,3)

    # Check if besselj(ν, approx_besselroots(ν, n) ) is small:
    tol = 1e-11

    for ν = 0.:0.1:5.
        n = 10
        @test norm( besselj.(ν, approx_besselroots(ν, n) ), Inf ) < tol
    end
end
