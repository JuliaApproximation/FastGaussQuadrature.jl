import SpecialFunctions
import SpecialFunctions: besselj

@testset "Bessel Roots" begin
    tol = 1e-11

    # Check if besselj(nu, besselroots(nu, n) ) is small:

    for nu = 0.:0.1:5.
        let n = 10
            @test norm( besselj.(nu, besselroots(nu, n) ), Inf ) < tol
        end
    end
end