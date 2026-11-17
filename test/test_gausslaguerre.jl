# Test gausslaguerre.jl

# Test integration
tol = 1.0e-13

# Evaluate the exact integral of x^α * p(x) *exp(-x) on the positive halfline,
# where p(x) = sum( coef[i+1]*x^i, i=0..degree(p)) is a polynomial given by its
# monomial coefficients
function exact_integral_poly(α, coef)
    # do calculations in BigFloat precision
    z = big(0.0)
    for i in 1:length(coef)
        z += gamma(big(α) + big(i)) * big(coef[i])
    end
    return Float64(z)
end

function polyval(x, coef)
    z = zeros(size(x))
    d = length(coef) - 1
    for i in 1:length(x)
        z[i] = sum(coef .* x[i] .^ (0:d))
    end
    return z
end

# Evaluate the exact integral of cos(a*x) * x^α * exp(-x) on the positive halfline
# - For α = 0
exact_integral_cos1(a) = 1 / (a^2 + 1)
# - For α = 1/2
exact_integral_cos2(a) = sqrt(π) * cos(3 / 2 * atan(a)) / 2 / (a^2 + 1)^(3 / 4)

Random.seed!(0)

@testset "Gauss–Laguerre" begin
    # Check error
    @test_throws DomainError gausslaguerre(1, -1.4)
    @test_throws DomainError gausslaguerre(-1)

    # Type inference
    @inferred gausslaguerre(Float32, 10)
    @inferred gausslaguerre(BigFloat, 10)
    @inferred gausslaguerre(BigFloat, 130)

    α = 0.5
    @inferred gausslaguerre(10, Float32(α))
    @inferred gausslaguerre(10, big(α))
    @inferred gausslaguerre(130, big(α))

    # Check optional argument
    for n in 0:20
        @test gausslaguerre(n) == gausslaguerre(n, 0) == gausslaguerre(n, 0.0)
    end
    for n in 30:10:200
        @test gausslaguerre(n) == gausslaguerre(n, 0) == gausslaguerre(n, 0.0)
    end

    ##########
    # Test the special cases
    ##########

    x, w = gausslaguerre(0)
    @test x == Float64[]
    @test w == Float64[]

    x, w = gausslaguerre(1, 0.4)
    @test x == [1.4]
    @test w == [gamma(1.4)]

    for α in (1, 1 // 2)
        x, w = gausslaguerre(1, α)
        @test x isa Vector{Float64} && w isa Vector{Float64}
        @test (x, w) == gausslaguerre(1, float(α))
    end

    x, w = gausslaguerre(2, -0.37)
    @test x ≈ [0.35328546651962944; 2.90671453348037]
    @test w ≈ [1.269857265034167; 0.15433992954023468]

    @testset "Float16 rules" begin
        for n in 0:2, α in (Float16(0), Float16(0.5))
            x, w = gausslaguerre(n, α)
            @test eltype(x) == eltype(w) == Float16
            xr, wr = gausslaguerre(n, Float64(α))
            @test all(isapprox.(x, xr; rtol = 4eps(Float16)))
            @test all(isapprox.(w, wr; rtol = 4eps(Float16)))
        end
        for (n, a) in ((3, 0), (3, 0.5), (8, -0.99), (32, 0.5), (128, 0), (128, -0.99))
            α = Float16(a)
            x, w = gausslaguerre(n, α)
            @test eltype(x) == eltype(w) == Float16
            @test all(isfinite, x) && all(>(0), diff(x))
            @test all(isfinite, w) && all(>=(0), w)
            for i in unique((1, cld(n, 2), n))
                xr, wr = QuadratureReference.nodeweight(:laguerre, n, α, 0, x[i]; bits = 128)
                @test isapprox(x[i], xr; rtol = 2eps(Float16), atol = eps(zero(Float16)))
                @test isapprox(w[i], wr; rtol = 2eps(Float16), atol = eps(zero(Float16)))
            end
        end
        # Validate the converted rule as well as its Float32 working values.
        @test_throws ErrorException gausslaguerre(1, Float16(10))
    end

    @testset "Two-point rule near the parameter boundary" begin
        @testset "T: $T" for T in (Float32, Float64, BigFloat)
            α = nextfloat(-one(T))
            x, w = gausslaguerre(2, α)
            @test all(>(0), x)
            @test all(isfinite, w)
            @test all(>(0), w)
            @test isapprox(sum(w), gamma(α + 1); rtol = 8eps(T))
            @test isapprox(dot(w, x), gamma(α + 2); rtol = 8eps(T))
            @test isapprox(x[1] / (α + 1), T(1) / 2; rtol = 8eps(T))
        end
    end

    @testset "Small rule near the Laguerre parameter boundary" begin
        α = nextfloat(-1.0)
        x, _ = gausslaguerre(14, α)
        @test all(>(0), x)
        # For δ=α+1, n*x₁/δ = 1 + (n-1)δ/(2n) + O(δ²).
        @test x[1] / (α + 1) ≈ 1 / 14 rtol = 8eps(Float64)
    end

    @testset "Small Laguerre roots" begin
        # Independent roots catch the small-node loss in the bounded recurrence.
        for (n, a, i) in ((31, -0.99, 1), (127, -0.75, 1), (31, 0.2, 1))
            x, _ = gausslaguerre(n, a)
            xr, _ = QuadratureReference.nodeweight(:laguerre, n, a, 0, x[i]; bits = 128)
            @test abs(BigFloat(x[i]) - xr) <= 48eps(Float64) * abs(xr)
        end
    end

    @testset "Full-rule upper Laguerre node" begin
        n, i = 1001, 998
        x, _ = gausslaguerre(n, 0.0)
        xr, _ = QuadratureReference.nodeweight(:laguerre, n, 0, 0, x[i]; bits = 128)
        # Frozen native GST node envelope, including the final square rounding.
        allowance = eps(Float64) * (17xr + 16n * sqrt(xr / (4n + 2)))
        @test abs(BigFloat(x[i]) - xr) <= allowance
    end

    @testset "Float64 rules near the Laguerre parameter boundary" begin
        # Independent 512-bit classical Laguerre polynomial references. The
        # parameters are the exact values of the corresponding Float64 inputs.
        @testset "n: $n" for (n, reference) in (
                (
                    128, (
                        (1, 7.851191089816142e-5, 95.59891032780894),
                        (2, 0.02895527853616579, 1.565850545585225),
                        (6, 0.5310611974292319, 0.22304104755100823),
                        (7, 0.7531111239342755, 0.15049886314476896),
                        (12, 2.4442295469824766, 0.015595426392628143),
                    ),
                ),
                (
                    1000, (
                        (1, 1.0049866900161312e-5, 93.6601611162678),
                        (2, 0.003706333359401857, 1.5731853741246702),
                        (6, 0.06795503336156661, 0.3470821120707875),
                        (7, 0.09635501708952898, 0.28419686703769),
                        (12, 0.31238217820973163, 0.1285968688848128),
                    ),
                ),
            )
            x, w = gausslaguerre(n, -0.99)
            @test issorted(x) && all(>(0), x)
            @test all(isfinite, x) && all(isfinite, w) && all(≥(0), w)
            @testset "i: $i" for (i, node, weight) in reference
                @test x[i] ≈ node rtol = 5.0e-13
                @test w[i] ≈ weight rtol = 5.0e-13
                @test w[i] > 0
            end
        end
        @testset "α: $α" for (α, node, weight) in (
                (-0.99999999, 7.812500078013256e-11, 9.999999498325518e7),
                (nextfloat(-1.0), 8.673617379884035e-19, 9.007199254740987e15),
            )
            x, w = gausslaguerre(128, α)
            @test issorted(x) && all(>(0), x)
            @test all(isfinite, w) && all(≥(0), w)
            @test x[1] ≈ node rtol = 5.0e-13
            @test w[1] ≈ weight rtol = 5.0e-13
            @test w[1] > 0
        end
    end

    @testset "Derivative at the accepted node" begin
        # Independent 256-bit classical Laguerre recurrence/Newton reference.
        x, w = FastGaussQuadrature.gausslaguerre_rec(15, 0.0)
        @test isapprox(w[end], 1.6005949062111332e-20; rtol = 5.0e-14)
        @test FastGaussQuadrature.gausslaguerre_rec(15, 0.0; reduced = true) == (x, w)
        # The inverse derivative squared can overflow before division by x,
        # even though the weight itself is representable.
        xlarge, wlarge = gausslaguerre(2, 170.0)
        _, wrec = FastGaussQuadrature.gl_rec_newton(xlarge[1], 2, 170.0)
        @test isfinite(wrec)
        @test isapprox(wrec, wlarge[1]; rtol = 2.0e-13)
    end

    @testset "Native startup near the half parameter" begin
        setprecision(BigFloat, 128) do
            a = prevfloat(-BigFloat(0.5))
            x, w = gausslaguerre(129, a)
            # Independent recurrence references and Sturm indices detect the
            # loss from a starting point that collapses towards the origin.
            QuadratureReference.check(:laguerre, 129, a, 0, x, w; indices = (1, 65, 129))
        end
    end

    @testset "Large-parameter rules" begin
        @testset "α: $α" for α in (nextfloat(5.0), 50.0, 169.0)
            x, w = gausslaguerre(15, α)
            @test issorted(x) && all(>(0), x)
            @test all(isfinite, w) && all(>(0), w)
            @test isapprox(sum(w) / gamma(α + 1), 1; rtol = 2.0e-13)
            @test isapprox(dot(w / gamma(α + 1), x), α + 1; rtol = 2.0e-13)
        end
        # A mass check alone misses the old skipped-root failure at α=50.
        x, w = gausslaguerre(128, 50.0)
        @test all(>(0), diff(x))
        @test all(isfinite, x) && all(isfinite, w)
        xref, wref = setprecision(BigFloat, 256) do
            a = BigFloat(50)
            references = [QuadratureReference.nodeweight(:laguerre, 128, a, 0, seed; bits = 128) for seed in x]
            refs = first.(references)
            finer = [QuadratureReference.nodeweight(:laguerre, 128, a, 0, seed; bits = 192)[1] for seed in x]
            @test all(isapprox.(refs, finer; rtol = BigFloat(2)^(-192), atol = 0))
            # Independent recurrence/Sturm counts detect skipped or duplicated
            # roots even though Newton starts from the computed nodes.
            indexed = all(eachindex(refs)) do i
                radius = BigFloat(2)^(-128) * max(one(a), abs(refs[i]))
                QuadratureReference.roots_below(:laguerre, 128, a, 0, refs[i] - radius) == i - 1 &&
                    QuadratureReference.roots_below(:laguerre, 128, a, 0, refs[i] + radius) == i
            end
            @test indexed
            refs, last.(references)
        end
        @test isapprox(x, xref; rtol = 2.0e-13)
        for k in (1, 2, 64, 128)
            xr, wr = xref[k], wref[k]
            @test isapprox(x[k], xr; rtol = 8.0e-14)
            if wr ≥ floatmin(Float64)
                @test isapprox(w[k], wr; rtol = 1.0e-12)
                @test w[k] > 0
            end
        end
        # Full rules retain finite nodes after weights have reached the cutoff.
        x, w = gausslaguerre(1000, 10.0)
        @test all(isfinite, x) && all(>(0), diff(x))
        @test all(isfinite, w) && all(≥(0), w) && iszero(w[end])
        for k in (1, 2, 500, 1000)
            xr, wr = QuadratureReference.nodeweight(:laguerre, 1000, 10.0, 0, x[k]; bits = 128)
            @test isapprox(x[k], xr; rtol = 8.0e-14)
            if wr ≥ floatmin(Float64)
                @test isapprox(w[k], wr; rtol = 1.0e-12)
                @test w[k] > 0
            end
        end
    end

    @testset "Asymptotic complementary bulk phase" begin
        for (n, a, k) in ((10000, 0.0, 100), (10000, 0.5, 100), (1000, 0.2, 100))
            x, w, _ = FastGaussQuadrature.gausslaguerre_asy_bulk(n, a, k, inv(4n + 2a + 2), -1)
            xr, wr = QuadratureReference.nodeweight(:laguerre, n, a, 0, x; bits = 128)
            @test isapprox(x, xr; rtol = 4.0e-15)
            @test isapprox(w, wr; rtol = 2.0e-14)
        end
    end

    @testset "Manual asymptotic rules and reduced tails" begin
        # Explicit/adaptive/default depths, both half-parameters and alpha near
        # -1. At n=190/200 the normal Airy tail needs recomputation;
        # the unrecomputed large-order case has an underflowed Airy tail.
        for (n, a, options) in (
                (128, 0.0, (; T = 3, recompute = true)),
                (128, -0.99, (; T = -1, recompute = true)),
                (190, -0.5, (; recompute = true)),
                (200, 0.5, (; T = -1, recompute = true)),
                (1000, 0.2, (;)),
            )
            x, w = FastGaussQuadrature.gausslaguerre_asy(n, a; options...)
            xr, wr = FastGaussQuadrature.gausslaguerre_asy(n, a; options..., reduced = true)
            @test length(x) == length(w) == n
            @test all(>(0), diff(x)) && all(t -> 0 < t <= 4n + 2a + 2, x)
            @test all(isfinite, w) && all(>=(0), w)
            cutoff = findfirst(<(10floatmin(Float64)), w)
            retained = isnothing(cutoff) ? n : cutoff - 1
            @test xr == x[1:retained] && wr == w[1:retained]

            # Independent references check Bessel, bulk and Airy nodes, plus
            # both sides of an omitted tail. These selected-point regressions
            # use the frozen bounded-REC 48eps node allowance and the existing
            # 1e-12 normal-weight threshold. Recompute is conditional: these
            # checks do not certify every node or the whole manual ASY domain.
            for i in unique((1, n ÷ 2, n, retained, min(retained + 1, n)))
                node, weight = QuadratureReference.nodeweight(:laguerre, n, a, 0, x[i]; bits = 128)
                setprecision(BigFloat, 256) do
                    delta = BigFloat(2)^(-128) * max(one(node), abs(node))
                    @test QuadratureReference.roots_below(:laguerre, n, BigFloat(a), BigFloat(0), node - delta) == i - 1
                    @test QuadratureReference.roots_below(:laguerre, n, BigFloat(a), BigFloat(0), node + delta) == i
                end
                @test isapprox(x[i], node; rtol = 48eps(Float64))
                if weight >= floatmin(Float64)
                    @test isapprox(w[i], weight; rtol = 1.0e-12)
                else
                    @test abs(BigFloat(w[i]) - weight) <= 2BigFloat(nextfloat(0.0))
                end
                i == retained && (@test weight >= 10floatmin(Float64))
                i == retained + 1 && (@test weight < 10floatmin(Float64))
            end
        end
    end

    @testset "GST bulk nodes and weights" begin
        for n in (128, 251)
            x, w = gausslaguerre(n, 0.0)
            for k in (n ÷ 2 + 1, 9n ÷ 10)
                xr, wr = QuadratureReference.nodeweight(:laguerre, n, 0, 0, x[k]; bits = 128)
                @test isapprox(x[k], xr; rtol = 4.0e-15)
                @test isapprox(w[k], wr; rtol = 4.0e-13)
            end
        end
    end

    @testset "Gil–Segura–Temme iteration" begin
        # Classical Laguerre polynomials provide an independent reference for
        # the transformed-ODE iteration and its scaled derivative weights.
        @testset "n: $n, α: $a" for (n, a) in (
                (32, -0.5), (32, -0.03741494404284086), (127, 0.5),
                (15, 5.0), (251, 100.0),
            )
            x, w = gausslaguerre(n, a)
            @test all(>(0), diff(x))
            @test all(isfinite, x) && all(isfinite, w)
            for k in unique([1, 2, n ÷ 2, n])
                xr, wr = QuadratureReference.nodeweight(:laguerre, n, a, 0, x[k]; bits = 128)
                @test isapprox(x[k], xr; rtol = 8.0e-14)
                if wr ≥ floatmin(Float64)
                    @test isapprox(w[k], wr; rtol = 1.0e-12)
                    @test w[k] > 0
                end
            end
        end
        # Direct asymptotics need more accurate Bessel zeros than the public
        # twelve-digit approximations, for both full and reduced rules.
        full = gausslaguerre(10000, 0.5)
        indices = (1, 6, 7, 12)
        references = [QuadratureReference.nodeweight(:laguerre, 10000, 0.5, 0, full[1][k]; bits = 128) for k in indices]
        @testset "reduced: $reduced" for reduced in (false, true)
            x, w = reduced ? gausslaguerre(10000, 0.5; reduced) : full
            for (k, (xr, wr)) in zip(indices, references)
                @test isapprox(x[k], xr; rtol = 4.0e-15)
                @test isapprox(w[k], wr; rtol = 8.0e-15)
            end
        end
        @testset "Reduced GST tail and normalization" begin
            for (n, a) in ((1000, -0.99), (1000, 4.9), (10000, -0.5), (10000, 0.2))
                xf, wf = gausslaguerre(n, a)
                x, w = gausslaguerre(n, a; reduced = true)
                cutoff = findfirst(<(10floatmin(Float64)), wf)
                @test length(x) == cutoff - 1
                @test x == xf[1:length(x)]
                @test isapprox(w, wf[1:length(w)]; rtol = 8eps(Float64), atol = 0)
                for i in (1, 100, length(x))
                    xr, wr = QuadratureReference.nodeweight(:laguerre, n, a, 0, x[i]; bits = 128)
                    @test isapprox(x[i], xr; rtol = 4.0e-15)
                    @test isapprox(w[i], wr; rtol = 1.0e-12)
                end
            end
        end
    end


    ##########
    n = 42
    ##########
    α = 0.0
    coef = rand(17)
    Z = exact_integral_poly(α, coef)

    x, w = gausslaguerre(n, α)
    @test isa(x, Vector{Float64})
    @test isa(w, Vector{Float64})
    @test abs(sum(w) - 1) < tol
    @test abs(dot(w, x) - 1) < tol
    @test abs(x[37] - 98.388267163326702) < tol
    @test abs(w[7] - 0.055372813167092) < tol
    @test abs(dot(w, polyval(x, coef)) - Z) / abs(Z) < tol

    x_rec, w_rec = FastGaussQuadrature.gausslaguerre_rec(n, α)
    @test abs(x_rec[37] - 98.388267163326702) < tol
    @test abs(w_rec[7] - 0.055372813167092) < tol
    @test abs(dot(w_rec, polyval(x_rec, coef)) - Z) / abs(Z) < tol
    @test_throws DomainError FastGaussQuadrature.gausslaguerre_rec(128, α)
    @test_throws DomainError FastGaussQuadrature.gausslaguerre_asy(1001, α; recompute = true)
    for helper in (FastGaussQuadrature.gausslaguerre_rec, FastGaussQuadrature.gausslaguerre_asy)
        @test_throws ArgumentError helper(42, Float32(α))
        @test_throws ArgumentError helper(42, BigFloat(α))
    end


    α = 0.5
    coef = rand(17)
    Z = exact_integral_poly(α, coef)
    x, w = gausslaguerre(n, α)
    @test abs(dot(w, polyval(x, coef)) - Z) / abs(Z) < tol


    ##########
    n = 251
    ##########
    α = 0.0
    x, w = gausslaguerre(n, α)
    a = 4
    Z = exact_integral_cos1(a)
    @test isa(x, Vector{Float64})
    @test isa(w, Vector{Float64})
    @test abs(x[37] - 13.309000189442097) < tol
    @test abs(w[3] - 0.050091759039996) < tol
    @test abs(dot(w, cos.(a * x)) - Z) < tol

    α = 0.5
    x, w = gausslaguerre(n, α)
    @test abs(x[44] - 19.095577327730616) < tol
    @test abs(w[18] - 0.026245779174690266) < tol
    a = 4
    Z = exact_integral_cos2(a)
    @test abs(dot(w, cos.(a * x)) - Z) < tol


    ############
    n = 350000
    ############
    α = 0.0
    a = 50
    Z = exact_integral_cos1(a)
    x, w = gausslaguerre(n, α)
    @test isa(x, Vector{Float64})
    @test isa(w, Vector{Float64})
    @test abs(dot(w, cos.(a * x)) - Z) < tol

    α = 0.44
    x, w = gausslaguerre(n, α; reduced = true)
    @test w[end] < 100 * floatmin(Float64)
end
