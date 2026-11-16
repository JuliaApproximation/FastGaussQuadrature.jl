@testset "Legendre expansion working precision" begin
    # Both the tabulated and asymptotic Bessel factors must retain Float32.
    for k in (2, 21)
        phase = Float32(π) * (Float32(k) - 0.25f0)
        d, beta_correction = @inferred FastGaussQuadrature.legendre_bessel(k, phase)
        @test d isa Float32 && beta_correction isa Float32
    end

    previousbits = precision(BigFloat)
    for bits in (24, 53)
        setprecision(BigFloat, bits) do
            x, w = gausslegendre(BigFloat, 61)
            @test precision(BigFloat) == bits
            @test all(t -> precision(t) == bits, x)
            @test all(t -> precision(t) == bits, w)
        end
    end
    @test precision(BigFloat) == previousbits
end

@testset "Legendre Float32 floor across the recurrence cutoff" begin
    for T in (Float16, BigFloat), n in (6, 60, 61, 200)
        setprecision(BigFloat, 16) do
            x32, w32 = gausslegendre(Float32, n)
            x, w = gausslegendre(T, n)
            @test x == T.(x32)
            @test w == T.(w32)
            @test precision(BigFloat) == 16
        end
    end
end
@testset "Legendre seed working precision" begin
    for n in (6, 7, 60)
        seeds = @inferred FastGaussQuadrature.leg_initial_guess(n, Float32)
        @test eltype(seeds) == Float32
        @test length(seeds) == n ÷ 2 + 1
        @test all(isfinite, seeds) && all(>(0), diff(seeds))
    end
    previousbits = precision(BigFloat)
    for bits in (24, 40, 52)
        setprecision(BigFloat, bits) do
            seeds = @inferred FastGaussQuadrature.leg_initial_guess(60, BigFloat)
            @test eltype(seeds) == BigFloat
            @test all(t -> precision(t) == bits, seeds)
            @test precision(BigFloat) == bits
        end
    end
    @test precision(BigFloat) == previousbits
    # Include both tabulated Bessel zeros and their asymptotic continuation.
    @test eltype(@inferred FastGaussQuadrature.besselZeroRoots(31, Float32)) == Float32
    phases = FastGaussQuadrature.besselZeroRoots(31, Float32) ./ 60.5f0
    @test eltype(@inferred FastGaussQuadrature.legpts_nodes(60, phases)) == Float32
end

@testset "Legendre precision-dependent expansion depth" begin
    # Exercise different truncation transitions in Float32 and Float64. The
    # central node needs an absolute allowance proportional to its 1/n scale.
    for (T, n) in ((Float32, 379), (Float32, 380), (Float64, 294), (Float64, 295))
        x, w = gausslegendre(T, n)
        for i in (20, n ÷ 2)
            xr, wr = setprecision(BigFloat, 256) do
                seed = -cospi((BigFloat(i) - 1 // 4) / (BigFloat(n) + 1 // 2))
                QuadratureReference.nodeweight(:jacobi, n, 0, 0, seed; bits = 128)
            end
            @test abs(x[i] - xr) <= 4eps(T) * max(abs(xr), inv(BigFloat(n)))
            @test abs(w[i] - wr) <= 8eps(T) * wr
        end
    end
end
