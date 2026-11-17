@testset "Bessel Roots" begin
    @test_throws DomainError approx_besselroots(0.0, -1)

    @test_deprecated besselroots(0.2, 3)

    # Check if besselj(ν, approx_besselroots(ν, n) ) is small:
    tol = 1.0e-11

    for ν in 0.0:0.1:5.0
        n = 10
        @test norm(besselj.(ν, approx_besselroots(ν, n)), Inf) < tol
    end

    @testset "Real parameters and asymptotic zeros" begin
        # Cross the J₀ table boundary and the ν > 5 approximation branch.
        for ν in (0, 6)
            roots = @inferred approx_besselroots(ν, 21)
            @test roots isa Vector{Float64} && length(roots) == 21
            @test all(>(0), roots) && all(>(0), diff(roots))
            @test abs(besselj(ν, roots[end])) < tol
        end
        # J₁/₂(x) = sqrt(2/(π*x))*sin(x), so its kth positive zero is kπ.
        @test FastGaussQuadrature.McMahon(1 // 2, 7) == 7π
    end
end
