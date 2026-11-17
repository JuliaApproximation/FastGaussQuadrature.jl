@testset "Gauss–Jacobi" begin
    @testset "check error" begin
        @test_throws DomainError gaussjacobi(-1, 0.1, 0.1)
        @test_throws DomainError gaussjacobi(-1, -1.1, 0.1)
        @test_throws DomainError gaussjacobi(-1, 0.1, -1.1)
        @test_throws DomainError gaussjacobi(2, -1.0, 0.0)
        @test_throws DomainError gaussjacobi(2, 0.0, -1.0)
        # Finite parameters can exceed the working coefficient or weight range.
        @test_throws ErrorException gaussjacobi(2, floatmax(Float32), 0.0f0)
        @test_throws ErrorException gaussjacobi(2, 200.0f0, 0.0f0)
    end

    @testset "native endpoint coordinates" begin
        for T in (Float32, Float64)
            q = -eps(T) / 4
            x, left, right = @inferred FastGaussQuadrature.jacobi_gst_geometry(q, one(T))
            @test x == one(T)
            @test left == -q > 0
            @test FastGaussQuadrature.jacobi_gst_divisor(left, right, Val(:radau), true) == left
            # A step below coordinate resolution must not advance the ODE state.
            q = T(0.125)
            next, u, d = @inferred FastGaussQuadrature.jacobi_gst_propagate(
                q, zero(T), eps(q) / 4, one(T), one(T), T(42), zero(T), T(42), zero(T), zero(T)
            )
            @test (next, u, d) == (q, one(T), one(T))
        end
    end

    tol = 1.0e-14
    tol_x = 4.0e-15
    tol_w = 1.0e-13

    n_range = vcat(0:10, 20:10:100, 120:20, 300)

    @testset "compatible with Gauss-Legendre" begin
        α₀ = 0.0
        β₀ = 0.0
        α₊ = nextfloat(α₀)
        β₊ = nextfloat(β₀)
        α₋ = prevfloat(α₀)
        β₋ = prevfloat(β₀)
        @testset "n: $n" for n in n_range
            x, w = gausslegendre(n)
            x₀₀, w₀₀ = gaussjacobi(n, α₀, β₀)
            x₀₊, w₀₊ = gaussjacobi(n, α₀, β₊)
            x₀₋, w₀₋ = gaussjacobi(n, α₀, β₋)
            x₊₀, w₊₀ = gaussjacobi(n, α₊, β₀)
            x₊₊, w₊₊ = gaussjacobi(n, α₊, β₊)
            x₊₋, w₊₋ = gaussjacobi(n, α₊, β₋)
            x₋₀, w₋₀ = gaussjacobi(n, α₋, β₀)
            x₋₊, w₋₊ = gaussjacobi(n, α₋, β₊)
            x₋₋, w₋₋ = gaussjacobi(n, α₋, β₋)
            for _w in (w₀₀, w₀₊, w₀₋, w₊₀, w₊₊, w₊₋, w₋₀, w₋₊, w₋₋)
                @test norm(w - _w) < tol_w
            end
            for _x in (x₀₀, x₀₊, x₀₋, x₊₀, x₊₊, x₊₋, x₋₀, x₋₊, x₋₋)
                @test norm(x - _x) < tol_x
            end
        end
    end

    @testset "compatible with Gauss-Chebyshev of the 1st kind" begin
        α₀ = -1 / 2
        β₀ = -1 / 2
        α₊ = nextfloat(α₀)
        β₊ = nextfloat(β₀)
        α₋ = prevfloat(α₀)
        β₋ = prevfloat(β₀)
        @testset "n: $n" for n in n_range
            x, w = @test_deprecated gausschebyshev(n, 1)
            x₀₀, w₀₀ = gaussjacobi(n, α₀, β₀)
            x₀₊, w₀₊ = gaussjacobi(n, α₀, β₊)
            x₀₋, w₀₋ = gaussjacobi(n, α₀, β₋)
            x₊₀, w₊₀ = gaussjacobi(n, α₊, β₀)
            x₊₊, w₊₊ = gaussjacobi(n, α₊, β₊)
            x₊₋, w₊₋ = gaussjacobi(n, α₊, β₋)
            x₋₀, w₋₀ = gaussjacobi(n, α₋, β₀)
            x₋₊, w₋₊ = gaussjacobi(n, α₋, β₊)
            x₋₋, w₋₋ = gaussjacobi(n, α₋, β₋)
            for _w in (w₀₀, w₀₊, w₀₋, w₊₀, w₊₊, w₊₋, w₋₀, w₋₊, w₋₋)
                @test norm(w - _w) < tol_w
            end
            for _x in (x₀₀, x₀₊, x₀₋, x₊₀, x₊₊, x₊₋, x₋₀, x₋₊, x₋₋)
                @test norm(x - _x) < tol_x
            end
        end
    end

    @testset "compatible with Gauss-Chebyshev of the 2nd kind" begin
        α₀ = 1 / 2
        β₀ = 1 / 2
        α₊ = nextfloat(α₀)
        β₊ = nextfloat(β₀)
        α₋ = prevfloat(α₀)
        β₋ = prevfloat(β₀)
        @testset "n: $n" for n in n_range
            x, w = @test_deprecated gausschebyshev(n, 2)
            x₀₀, w₀₀ = gaussjacobi(n, α₀, β₀)
            x₀₊, w₀₊ = gaussjacobi(n, α₀, β₊)
            x₀₋, w₀₋ = gaussjacobi(n, α₀, β₋)
            x₊₀, w₊₀ = gaussjacobi(n, α₊, β₀)
            x₊₊, w₊₊ = gaussjacobi(n, α₊, β₊)
            x₊₋, w₊₋ = gaussjacobi(n, α₊, β₋)
            x₋₀, w₋₀ = gaussjacobi(n, α₋, β₀)
            x₋₊, w₋₊ = gaussjacobi(n, α₋, β₊)
            x₋₋, w₋₋ = gaussjacobi(n, α₋, β₋)
            for _w in (w₀₀, w₀₊, w₀₋, w₊₀, w₊₊, w₊₋, w₋₀, w₋₊, w₋₋)
                @test norm(w - _w) < tol_w
            end
            for _x in (x₀₀, x₀₊, x₀₋, x₊₀, x₊₊, x₊₋, x₋₀, x₋₊, x₋₋)
                @test norm(x - _x) < tol_x
            end
        end
    end

    @testset "compatible with Gauss-Chebyshev of the 3rd kind" begin
        α₀ = -1 / 2
        β₀ = 1 / 2
        α₊ = nextfloat(α₀)
        β₊ = nextfloat(β₀)
        α₋ = prevfloat(α₀)
        β₋ = prevfloat(β₀)
        @testset "n: $n" for n in n_range
            x, w = @test_deprecated gausschebyshev(n, 3)
            x₀₀, w₀₀ = gaussjacobi(n, α₀, β₀)
            x₀₊, w₀₊ = gaussjacobi(n, α₀, β₊)
            x₀₋, w₀₋ = gaussjacobi(n, α₀, β₋)
            x₊₀, w₊₀ = gaussjacobi(n, α₊, β₀)
            x₊₊, w₊₊ = gaussjacobi(n, α₊, β₊)
            x₊₋, w₊₋ = gaussjacobi(n, α₊, β₋)
            x₋₀, w₋₀ = gaussjacobi(n, α₋, β₀)
            x₋₊, w₋₊ = gaussjacobi(n, α₋, β₊)
            x₋₋, w₋₋ = gaussjacobi(n, α₋, β₋)
            for _w in (w₀₀, w₀₊, w₀₋, w₊₀, w₊₊, w₊₋, w₋₀, w₋₊, w₋₋)
                @test norm(w - _w) < tol_w
            end
            for _x in (x₀₀, x₀₊, x₀₋, x₊₀, x₊₊, x₊₋, x₋₀, x₋₊, x₋₋)
                @test norm(x - _x) < tol_x
            end
        end
    end

    @testset "compatible with Gauss-Chebyshev of the 4th kind" begin
        α₀ = 1 / 2
        β₀ = -1 / 2
        α₊ = nextfloat(α₀)
        β₊ = nextfloat(β₀)
        α₋ = prevfloat(α₀)
        β₋ = prevfloat(β₀)
        @testset "n: $n" for n in n_range
            x, w = @test_deprecated gausschebyshev(n, 4)
            x₀₀, w₀₀ = gaussjacobi(n, α₀, β₀)
            x₀₊, w₀₊ = gaussjacobi(n, α₀, β₊)
            x₀₋, w₀₋ = gaussjacobi(n, α₀, β₋)
            x₊₀, w₊₀ = gaussjacobi(n, α₊, β₀)
            x₊₊, w₊₊ = gaussjacobi(n, α₊, β₊)
            x₊₋, w₊₋ = gaussjacobi(n, α₊, β₋)
            x₋₀, w₋₀ = gaussjacobi(n, α₋, β₀)
            x₋₊, w₋₊ = gaussjacobi(n, α₋, β₊)
            x₋₋, w₋₋ = gaussjacobi(n, α₋, β₋)
            for _w in (w₀₀, w₀₊, w₀₋, w₊₀, w₊₊, w₊₋, w₋₀, w₋₊, w₋₋)
                @test norm(w - _w) < tol_w
            end
            for _x in (x₀₀, x₀₊, x₀₋, x₊₀, x₊₊, x₊₋, x₋₀, x₋₊, x₋₋)
                @test norm(x - _x) < tol_x
            end
        end
    end

    @testset "a small n" begin
        n = 42
        a = -0.1
        b = 0.3
        x, w = gaussjacobi(n, a, b)
        @test all(length(x) == n) && all(length(w) == n)
        @test abs(x[37] - 0.912883347814032) < tol
        @test abs(w[37] - 0.046661910947553) < tol

        @test dot(w, exp.(x)) ≈ 2.7568520361985516
    end
    @testset "a larger n (using GST)" begin
        a = -0.7
        b = 1.3
        n = 251
        x, w = gaussjacobi(n, a, b)
        @test all(length(x) == n) && all(length(w) == n)
        @test abs(x[37] + 0.893103435898983) < tol
        @test abs(w[37] - 1.962534523788093e-4) < tol

        @test dot(w, exp.(x)) ≈ 16.722957039404044
    end
    @testset "n = 1:" begin
        a = 1.0
        b = 2.0
        x, w = gaussjacobi(1, a, b)
        @test abs(x[1] - (b - a) / (a + b + 2)) < tol
        @test abs(w[1] - 2^(a + b + 1) * beta(a + 1, b + 1)) < tol

        @test dot(w, fill(1, length(x))) ≈ 1.3333333333333333
    end
    @testset "0.9, -0.1" begin
        x, w = gaussjacobi(1013, 0.9, -0.1)
        @test abs(x[2] + 0.999986012231899) < tol
        @test abs(x[13] + 0.999225722939832) < tol
        @test abs(w[2] - 9.314674169892358e-5) < tol
        @test abs(w[13] - 4.654651764553262e-4) < tol

        @test dot(w, exp.(x)) ≈ 1.6915068974063106

        x, w = gaussjacobi(10013, 0.9, -0.1)
        @test abs(x[2] + 0.999999856605293) < tol
        @test abs(x[13] + 0.999992061552711) < tol
        @test abs(w[2] - 1.509654630405615e-6) < tol
        @test abs(w[13] - 7.548275262993863e-6) < tol

        @test dot(w, exp.(x)) ≈ 1.6915068974063106
    end
    @testset "bug where row vectors were returned" begin
        x, w = gaussjacobi(10013, 0.9, -0.1)
        @test isa(x, Vector{Float64})
        @test isa(w, Vector{Float64})
    end
    @testset "last alpha and beta parameters:" begin
        x, w = gaussjacobi(100, 19.0, 21.0)
        @test abs(x[87] - 0.83221144617604) < tol
        @test abs(w[50] - 0.026363584978877305) < tol # new test without condition sum(w)==1
    end
    @testset "for small alpha and beta:" begin
        x, w = gaussjacobi(10000, 0.1, 0.2)
        @test abs(x[1] - -0.999999963363548) < tol
        @test abs(w[500] - 2.183393039546711e-5) < tol
    end
    @testset " bug where row vectors were returned" begin
        x, w = gaussjacobi(20, 11.0, 0.0)
        @test isa(x, Vector{Float64})
        @test isa(w, Vector{Float64})
    end
    @testset "missing promotion with integers (#117)" begin
        x1, w1 = gaussjacobi(500, 0, -1 / 2)
        x2, w2 = gaussjacobi(500, 0.0, -1 / 2)
        @test x1 == x2
        @test w1 == w2
    end

    @testset "jacobi_rec against a precomputed rule in BigFloat" begin
        x, w = FastGaussQuadrature.jacobi_rec(10, big(0), big(0))
        @test abs(x[6] - big(1488743389816312108848260011297199846175648594206916957079892535159036173556674) / big(10)^79) < 1.0e-70
        @test abs(w[4] - big(2692667193099963550912269215694693528597599384608837958005632762421534323191819) / big(10)^79) < 1.0e-70

        x, w = FastGaussQuadrature.jacobi_rec(10, big(2 // 10), big(-1 // 30))
        @test abs(x[7] - big(4146701176053244724446765533149623814305606864483746748725151214196203624388371) / big(10)^79) < 1.0e-70
        @test abs(w[3] - big(2482452398859023537636458227096017466104571146686194593251519240011574719667464) / big(10)^79) < 1.0e-70
    end

    @testset "Float32 parameters with n > 100" begin
        x, w = gaussjacobi(160, Float32(1), Float32(1))
        @test eltype(x) == Float32
        @test eltype(w) == Float32
    end

    @testset "boundary asymptotics were previously not implemented (#58,#130)" begin
        # check example from issue #58
        _, w = gaussjacobi(2^10, 0.25, 0)
        @test abs(w[end] - 3.607554904604311e-7) < tol

        # check example from issue #130
        _, w = gaussjacobi(2^16, -0.9, 0)
        @test abs(w[end] - 1.223027243345722) < tol
        @test abs(sum(w) - 2^(-0.9 + 1) / (-0.9 + 1)) < tol

        # double check NaN fixing in barycentric interpolation implementation
        x = range(0, stop = 1, length = 5)
        N = 10
        fvals = ones(N)
        xk = 0.5 * (sin.(pi * (-(N - 1):2:(N - 1)) / (2 * (N - 1))) .+ 1)
        vk = [0.5; ones(N - 1)]
        vk[2:2:end] .= -1
        vk[end] *= 0.5
        @test FastGaussQuadrature.bary(x, fvals, xk, vk) ≈ ones(5) rtol = 1.0e-16
    end

    @testset "symmetric recurrence, iterative and asymptotic rules" begin
        # One odd/even pair for each recurrence, GST and asymptotic route.
        @testset "n: $n, α: $a" for (n, a) in (
                (20, 0.1), (21, 0.1),
                (99, 4.9), (100, 4.9),
                (768, 0.1), (769, 0.1),
            )
            x, w = gaussjacobi(n, a, a)
            @test x == -reverse(x)
            @test w == reverse(w)
            @test !isodd(n) || iszero(x[cld(n, 2)])
        end
        @testset "T: $T" for T in (Float32, BigFloat)
            x, w = @inferred gaussjacobi(21, T(1), T(1))
            @test eltype(x) == eltype(w) == T
            @test x == -reverse(x)
            @test w == reverse(w)
            @test iszero(x[11])
        end
    end

    @testset "accuracy across routing boundaries" begin
        # Check classical polynomial roots and derivative weights rather than
        # comparing kernels. The small cases check every root; the large cases
        # include both sides of the asymptotic boundary/interior splice.
        @testset "T: $T, n: $n, α: $a, β: $b" for (T, n, a, b) in (
                (Float64, 24, 2.0, -0.5), (Float64, 25, 2.0, -0.5),
                (Float64, 24, 0.0, 2.0),
                (Float32, 20, 2.0f0, -0.5f0), (Float32, 21, 2.0f0, -0.5f0),
                (Float64, 20, prevfloat(-0.5), 0.2),
                (Float64, 20, 0.2, nextfloat(2.0)),
                (Float64, 767, 4.9, -0.99), (Float64, 768, 4.9, -0.99),
                (Float64, 769, 4.9, -0.99),
                (Float64, 768, 4.9, prevfloat(-0.99)),
                (Float64, 768, 5.0, -0.99),
            )
            x, w = @inferred gaussjacobi(n, a, b)
            @test length(x) == length(w) == n
            @test eltype(x) === eltype(w) === T
            @test all(isfinite, x) && all(>(0), diff(x))
            @test all(isfinite, w) && all(>(0), w)
            indices = n <= 25 ? eachindex(x) : (1, 10, 11, cld(n, 2), n - 10, n - 9, n)
            for i in indices
                xr, wr = QuadratureReference.nodeweight(:jacobi, n, a, b, x[i]; bits = 256)
                @test isapprox(x[i], xr; atol = 4eps(T), rtol = 0)
                # These long GST cases previously used twofold amplitudes.
                # Apply the native order/parameter budget; retain tighter checks
                # for the unchanged recurrence and asymptotic cases.
                native_gst = T === Float64 && n > 100 &&
                    !(n >= 768 && min(a, b) >= -0.99 && max(a, b) < 5)
                weight_budget = native_gst ? QuadratureReference.native_jacobi_weight_budget(n, a, b, wr) : 128eps(T)
                @test isapprox(w[i], wr; atol = 0, rtol = weight_budget)
                δ = minimum((xr + 1, 1 - xr, i == 1 ? 2 : xr - x[i - 1], i == n ? 2 : x[i + 1] - xr)) / 4
                @test QuadratureReference.roots_below(:jacobi, n, BigFloat(a), BigFloat(b), xr - δ) == i - 1
                @test QuadratureReference.roots_below(:jacobi, n, BigFloat(a), BigFloat(b), xr + δ) == i
            end
        end

        @testset "BigFloat: $bits bits, n: $n" for bits in (128, 256), n in (20, 21)
            setprecision(bits) do
                a, b = BigFloat(2), BigFloat(-0.5)
                x, w = @inferred gaussjacobi(n, a, b)
                QuadratureReference.check(:jacobi, n, a, b, x, w)
            end
        end
        # The explicit asymptotic adapter promotes rational/integer parameters.
        @test FastGaussQuadrature.jacobi_asy(768, 1 // 4, 0) == gaussjacobi(768, 0.25, 0.0)
    end

    @testset "symmetric mass without premature overflow" begin
        a = 1.0e6
        x, w = gaussjacobi(1, a, a)
        # Duplication reduces the exact integral to B(1/2, a + 1).
        reference = setprecision(256) do
            beta(BigFloat(1) / 2, BigFloat(a) + 1)
        end
        @test iszero(x[1])
        @test w[1] ≈ reference rtol = 2.0e-15
        x, w = gaussjacobi(20, 1000.0, 1000.0)
        @test sum(w) ≈ beta(big"0.5", big"1001") rtol = 3.0e-15
    end

    @testset "Generic high-precision Jacobi mass" begin
        setprecision(BigFloat, 128) do
            a, b = BigFloat(63), BigFloat(127)
            mass = FastGaussQuadrature.jacobi_gst_mass_generic(a, b)
            reference = setprecision(BigFloat, 256) do
                aa, bb = BigFloat(a), BigFloat(b)
                2^(aa + bb + 1) * beta(aa + 1, bb + 1)
            end
            # Exponentiating a rounded logarithm amplifies error by its size.
            @test isapprox(mass, reference; rtol = 8eps(BigFloat) * (1 + abs(log(reference))))
        end
    end

    @testset "small asymmetric rules with large parameters" begin
        # Nearly equal large parameters expose cancellation in the old weight
        # normalization. Check every root against classical polynomials.
        @testset "T: $T, n: $n, α: $a, β: $b" for T in (Float32, Float64), (n, α, β) in (
                    (2, 1.0e6, 1.0e6 + 1), (3, 5.0, -0.99), (19, 5.0, 0.2),
                    (3, 50.0, 0.0),
                ), (a, b) in ((T(α), T(β)), (T(β), T(α)))
            x, w = @inferred gaussjacobi(n, a, b)
            @test length(x) == length(w) == n
            @test eltype(x) === eltype(w) === T
            @test all(isfinite, x) && all(>(0), diff(x))
            @test all(isfinite, w) && all(>(0), w)
            for i in eachindex(x)
                xr, wr = QuadratureReference.nodeweight(:jacobi, n, a, b, x[i]; bits = 256)
                @test isapprox(x[i], xr; atol = 4eps(T), rtol = 0)
                @test isapprox(w[i], wr; atol = 0, rtol = 64eps(T))
                δ = minimum((xr + 1, 1 - xr, i == 1 ? 2 : xr - x[i - 1], i == n ? 2 : x[i + 1] - xr)) / 4
                @test QuadratureReference.roots_below(:jacobi, n, BigFloat(a), BigFloat(b), xr - δ) == i - 1
                @test QuadratureReference.roots_below(:jacobi, n, BigFloat(a), BigFloat(b), xr + δ) == i
            end
        end
    end

    @testset "large-order endpoint conditioning" begin
        # ASY is accurate at -0.99; closer to -1, its boundary weights can
        # lose relative accuracy even when the rounded nodes remain accurate.
        n = 1024
        @testset "α: $a, β: $b" for parameter in (-0.99, -1 + 1.0e-8, nextfloat(-1.0)),
                (a, b) in ((parameter, 0.2), (0.2, parameter))
            x, w = gaussjacobi(n, a, b)
            @test all(isfinite, x) && all(>(0), diff(x))
            @test all(isfinite, w) && all(>(0), w)
            for i in (1, n ÷ 2, n)
                xr, wr = QuadratureReference.nodeweight(:jacobi, n, a, b, x[i]; bits = 128)
                @test isapprox(x[i], xr; atol = 2eps(Float64), rtol = 0)
                # Near -1 the native sweep's scale matching also accumulates
                # roundoff; the independently frozen GST budget includes it.
                weight_budget = parameter >= -0.99 ? 128eps(Float64) :
                    QuadratureReference.native_jacobi_weight_budget(n, a, b, wr)
                @test isapprox(w[i], wr; atol = 0, rtol = weight_budget)
            end
        end
    end

    @testset "long fourth-order sweeps" begin
        # Independent classical-polynomial references computed at 256 bits
        # and checked at 384 bits, using the exact Float64 parameters.
        @testset "n: $n, α: $a, β: $b" for (n, a, b, indices, reference) in (
                (257, 0.1, 0.2, [1, 129, 257], [1.9086855427474173e-5, 0.012193607107039362, 5.077203525306096e-5]),
                (4096, -0.7, 1.3, [1, 2048, 4096], [3.548063785105645e-15, 0.0007668406121600569, 0.08488120866553461]),
                (65536, -0.9, 0.0, [1, 32768, 65536], [9.25893310303011e-10, 4.7936294313505936e-5, 1.223027243345723]),
            )
            x, w = FastGaussQuadrature.jacobi_gst(n, a, b)
            @test length(x) == length(w) == n
            @test all(isfinite, x) && all(>(0), diff(x))
            @test all(isfinite, w) && all(>(0), w)
            for (i, expected) in zip(indices, reference)
                # The reference literals were independently verified within
                # one epsilon; native Taylor error has an order-dependent bound.
                allowance = QuadratureReference.native_jacobi_weight_budget(n, a, b, expected) + eps(Float64)
                @test w[i] ≈ expected rtol = allowance atol = 0
            end
            @test sum(w) ≈ FastGaussQuadrature.jacobimoment(a, b) rtol = 16eps(Float64)
            xr, wr = FastGaussQuadrature.jacobi_gst(n, b, a)
            @test xr ≈ -reverse(x) atol = eps(Float64) rtol = 0
            @test all(isapprox.(wr, reverse(w); rtol = 32eps(Float64), atol = 0))
            if n >= 768
                _, selected_weights = gaussjacobi(n, a, b)
                for (i, expected) in zip(indices, reference)
                    @test selected_weights[i] ≈ expected rtol = 32eps(Float64) atol = 0
                end
            end
        end
    end

    @testset "fourth-order iteration for large parameters" begin
        # Independent 256-bit Jacobi-polynomial references, recomputed at 384 bits.
        # Eigenvector weights previously lost these representable tail weights.
        x, w = gaussjacobi(100, 1000.0, 1000.0)
        @test all(>(0), w)
        @test x == -reverse(x)
        @test w == reverse(w)
        @test x[1] ≈ -0.39672435238005427 atol = 2.0e-16
        @test w[1] ≈ 7.747127990992642e-77 rtol = 3.0e-13
        x, w = gaussjacobi(1000, 1000.0, 1000.0)
        @test w[85] > 0
        @test w[85] ≈ 1.4094433524956e-310 rtol = 2.0e-12

        # The iteration does not need the old full-eigensystem order limit.
        x, w = gaussjacobi(4001, 5, 5)
        @test length(x) == length(w) == 4001
        @test issorted(x)
        @test all(isfinite, x)
        @test all(>(0), w)
        @test iszero(x[2001])
        @test x[1] ≈ -0.9999976034533193 atol = 2.0e-16
        @test w[1] ≈ 5.232463578418606e-33 rtol = 5.0e-13

        x, w = gaussjacobi(100, 1000.0, 20.0)
        @test x[50] ≈ -0.8535886827952174 atol = 2.0e-16
        @test w[50] ≈ 1.0848302167401123e+249 rtol = 5.0e-13
        xr, wr = gaussjacobi(100, 20.0, 1000.0)
        @test xr == -reverse(x)
        @test wr ≈ reverse(w) rtol = 3.0e-14

        # The starting recurrence exceeds 1e100 and must rescale its state.
        x, w = gaussjacobi(512, 1000.0, 20.0)
        for i in (1, 256, 512)
            xr, wr = QuadratureReference.nodeweight(:jacobi, 512, 1000.0, 20.0, x[i]; bits = 256)
            setprecision(BigFloat, 256) do
                delta = BigFloat(2)^(-128) * max(one(xr), abs(xr))
                @test QuadratureReference.roots_below(:jacobi, 512, BigFloat(1000), BigFloat(20), xr - delta) == i - 1
                @test QuadratureReference.roots_below(:jacobi, 512, BigFloat(1000), BigFloat(20), xr + delta) == i
            end
            @test isapprox(x[i], xr; atol = 4eps(Float64), rtol = 0)
            allowance = QuadratureReference.native_jacobi_weight_budget(512, 1000.0, 20.0, wr)
            @test isapprox(w[i], wr; atol = 0, rtol = allowance)
        end
    end

    @testset "endpoint distances below node resolution" begin
        x, w = gaussjacobi(20, 0.2, nextfloat(-1.0))
        @test x[1] == -1.0
        @test all(isfinite, w)
        @test all(>(0), w)
        @test w[1] ≈ 1.0346554967051492e+16 rtol = 3.0e-14
        @test w[2] ≈ 1.91991395416194 rtol = 3.0e-14
        @test issorted(x)
        xr, wr = gaussjacobi(20, nextfloat(-1.0), 0.2)
        @test xr ≈ -reverse(x) rtol = 2.0e-15
        @test wr ≈ reverse(w) rtol = 3.0e-14
    end

    @testset "weights near the Jacobi parameter boundary" begin
        # Independent 384-bit classical Jacobi recurrence references, with the
        # Float64 parameters converted exactly to BigFloat.
        x, w = gaussjacobi(20, -0.999, 0.2)
        @test issorted(x) && all(isfinite, x)
        @test all(isfinite, w) && all(>(0), w)
        @test x[1] ≈ -0.9909409680801299 rtol = 5.0e-13
        @test w[1] ≈ 0.004068480305281893 rtol = 5.0e-13
        @test x[10] ≈ -0.030710555300320148 rtol = 5.0e-13
        @test w[10] ≈ 0.1506752468841809 rtol = 5.0e-13
        @test x[20] ≈ 0.9999950472826589 rtol = 5.0e-13
        @test w[20] ≈ 1143.0118108108413 rtol = 5.0e-13
    end
end

@testset "Low-precision Jacobi rules" begin
    # Endpoint weights are compared before rounding the reference nodes. A
    # Float16 node can be accurate while 1-x has lost several significant bits.
    b = Float16(0.2)
    for (n, a) in ((8, Float16(0.1)), (32, nextfloat(-one(Float16))), (128, Float16(5)))
        x, w = @inferred gaussjacobi(n, a, b)
        @test eltype(x) === eltype(w) === Float16
        @test all(isfinite, x) && all(>(0), diff(x))
        @test all(isfinite, w) && all(>=(0), w)
        for i in (1, cld(n, 2), n)
            xr, wr = QuadratureReference.nodeweight(:jacobi, n, a, b, x[i]; bits = 128)
            @test isapprox(x[i], xr; atol = eps(Float16), rtol = 0)
            @test isapprox(w[i], wr; atol = eps(zero(Float16)), rtol = 2eps(Float16))
        end
    end
    @test_throws ErrorException gaussjacobi(512, Float16(0.1), Float16(0.2))
end
