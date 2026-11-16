@testset "Generic precision quadrature" begin
    @testset "BigFloat rules" begin
        bits = 256
        previousbits = precision(BigFloat)
        setprecision(bits) do
            a, b = BigFloat(-3) / 4, BigFloat(5) / 4
            x, w = gaussjacobi(17, a, b)
            QuadratureReference.check(:jacobi, 17, a, b, x, w; indices = (1, 9, 17))

            a = BigFloat(1) / 4
            x, w = gausslaguerre(17, a)
            QuadratureReference.check(:laguerre, 17, a, zero(a), x, w; indices = (1, 9, 17))
            @test gausslaguerre(17, a; reduced = true) == (x, w)

            x, w = gausshermite(BigFloat, 3)
            QuadratureReference.check(:hermite, 3, big(0), big(0), x, w)
            @test x == -reverse(x)
            @test w == reverse(w)
            xn, wn = gausshermite(BigFloat, 3; normalize = true)
            @test xn ≈ sqrt(big(2)) .* x
            @test wn ≈ w ./ sqrt(big(π))

            # The endpoint rules inherit the higher-precision Jacobi iteration.
            tolerance = BigFloat(2)^(-bits + 16)
            @testset "quad: $quad, n: $n" for (quad, n, aa, bb) in ((gaussradau, 17, 0, 1), (gausslobatto, 17, 1, 1))
                x, w = quad(BigFloat, n)
                @test eltype(x) == eltype(w) == BigFloat
                @test x[1] == -1
                @test all(z -> precision(z) == bits, x)
                @test all(z -> precision(z) == bits, w)
                lastinterior = quad == gaussradau ? n : n - 1
                @testset "i: $i" for i in (2, 9, lastinterior)
                    xr, wr = QuadratureReference.nodeweight(:jacobi, lastinterior - 1, aa, bb, x[i]; bits)
                    wr /= quad == gaussradau ? 1 + xr : (1 - xr) * (1 + xr)
                    @test isapprox(x[i], xr; atol = tolerance, rtol = tolerance)
                    @test isapprox(w[i], wr; atol = 0, rtol = tolerance)
                end
                if quad == gaussradau
                    @test w[1] == BigFloat(2) / n^2
                else
                    @test x[end] == 1
                    @test x == -reverse(x)
                    @test w == reverse(w)
                    @test w[1] == BigFloat(2) / (n * (n - 1))
                end
            end
            @test precision(BigFloat) == bits
        end
        @test precision(BigFloat) == previousbits
    end

    @testset "Higher precision" begin
        # One case beyond the usual precisions detects fixed-precision constants
        # or iteration limits without repeating the full rule matrix.
        setprecision(768) do
            a, b = BigFloat(-3) / 4, BigFloat(5) / 4
            x, w = gaussjacobi(17, a, b)
            QuadratureReference.check(:jacobi, 17, a, b, x, w; indices = (1, 9, 17))
        end
    end

    @testset "High-precision small Jacobi rule" begin
        # This difficult parameter pair uses GST even at small order; its
        # stopping rules must scale with the requested precision.
        setprecision(8192) do
            a, b = BigFloat(499) / 100, -BigFloat(9) / 10
            x, w = gaussjacobi(2, a, b)
            QuadratureReference.check(:jacobi, 2, a, b, x, w)
        end
    end

    @testset "Large Legendre rules at generic precision" begin
        # Cover the first GST order, even/odd orders and larger rules at chosen
        # precisions. Independent Jacobi references also check root indices.
        @testset "$bits bits, n: $n" for (bits, n) in ((64, 61), (128, 128), (256, 512))
            setprecision(bits) do
                x, w = gausslegendre(BigFloat, n)
                QuadratureReference.check(:jacobi, n, big(0), big(0), x, w; indices = (1, 2, n ÷ 2, n))
                @test x == -reverse(x)
                @test w == reverse(w)
            end
        end
    end

    @testset "Failure cleanup" begin
        for n in (0, -1)
            @test_throws DomainError gaussradau(BigFloat, n)
        end
        setprecision(BigFloat, 4) do
            @test_throws ErrorException gausslaguerre(64, BigFloat(0))
            @test_throws ErrorException gausslaguerre(2, BigFloat(2)^10)
        end
    end

    @testset "Fixed precision Jacobi iteration" begin
        n, a, b = 101, Float32(1) / 4, Float32(1) / 8
        previousbits = precision(BigFloat)
        x, w = gaussjacobi(n, a, b)
        @test precision(BigFloat) == previousbits
        @test eltype(x) == eltype(w) == Float32
        @test issorted(x) && all(>(0), diff(x))
        @test all(isfinite, w) && all(>(0), w)
        @testset "i: $i" for i in (1, cld(n, 2), n)
            xr, wr = QuadratureReference.nodeweight(:jacobi, n, a, b, x[i]; bits = 128)
            @test isapprox(x[i], xr; atol = 256eps(Float32), rtol = 256eps(Float32))
            @test isapprox(w[i], wr; atol = 0, rtol = 256eps(Float32))
        end
    end

    @testset "Parameters near the integrability boundary" begin
        setprecision(256) do
            bits = precision(BigFloat)
            n, a = 8, nextfloat(-one(BigFloat))
            tolerance = BigFloat(2)^(-bits + 16)
            @testset "$family" for (family, b, rule) in (
                    (:laguerre, big(0), () -> gausslaguerre(n, a)),
                    (:jacobi, BigFloat(1) / 4, () -> gaussjacobi(n, a, BigFloat(1) / 4)),
                    (:jacobi, nextfloat(a), () -> gaussjacobi(n, a, nextfloat(a))),
                )
                x, w = rule()
                @test precision(BigFloat) == bits
                @test all(z -> precision(z) == bits, x)
                @test all(z -> precision(z) == bits, w)
                @test all(isfinite, x) && all(isfinite, w) && all(>(0), w)
                @testset "i: $i" for i in (1, 4, n)
                    # Resolving α+1 as well as the requested digits needs
                    # extra precision in the independent classical recurrence.
                    xr, wr = QuadratureReference.nodeweight(family, n, a, b, x[i]; bits = 2bits)
                    @test isapprox(x[i], xr; atol = 0, rtol = tolerance)
                    @test isapprox(w[i], wr; atol = 0, rtol = tolerance)
                end
            end
            # Forming α+β+2 loses a significant bit when both parameters
            # approach -1, even in the one-node closed form.
            x, w = gaussjacobi(1, a, nextfloat(a))
            xr, wr = QuadratureReference.nodeweight(:jacobi, 1, a, nextfloat(a), x[1]; bits = 2bits)
            @test isapprox(x[1], xr; atol = 0, rtol = tolerance)
            @test isapprox(w[1], wr; atol = 0, rtol = tolerance)
        end
    end

    @testset "Large parameters and outer weights" begin
        setprecision(256) do
            # The old Laguerre seed/extrapolation path could return a first
            # node near 5.53688 instead of 5.35263, despite an accurate mass.
            x, w = gausslaguerre(128, BigFloat(50))
            QuadratureReference.check(:laguerre, 128, big(50), big(0), x, w; indices = (1, 2, 64, 127, 128))
            # A centered recurrence initializes the iteration when α > n.
            x, w = gausslaguerre(32, BigFloat(1000))
            QuadratureReference.check(:laguerre, 32, big(1000), big(0), x, w; indices = (1, 32))
            @testset "α: $a, β: $b" for (a, b) in ((big(35), big(80)), (big(1000), big(1000)))
                x, w = gaussjacobi(21, a, b)
                QuadratureReference.check(:jacobi, 21, a, b, x, w; indices = (1, 2, 11, 20, 21))
                if a == b
                    @test x == -reverse(x)
                    @test w == reverse(w)
                    @test iszero(x[11])
                end
            end
        end
    end

    @testset "Float16 reductions and tiny weights" begin
        @test gaussjacobi(8, Float16(0), Float16(0)) == gausslegendre(Float16, 8)
        @test gaussjacobi(8, Float16(-0.5), Float16(-0.5)) == gausschebyshevt(Float16, 8)
        @test_throws ErrorException gaussradau(512, Float16(0.1), Float16(0.2))
        x, w = gaussjacobi(64, Float16(5), Float16(0.2))
        # Keep an exact rounding check where relative error is limited by the
        # subnormal spacing; ordinary low-precision cases live with each rule.
        _, wr = QuadratureReference.nodeweight(:jacobi, 64, Float16(5), Float16(0.2), x[end]; bits = 256)
        @test wr < floatmin(Float16)
        @test w[end] == Float16(wr)
    end

    @testset "Mixed parameter types" begin
        setprecision(128) do
            # n > 100 previously selected the Float64 asymptotic wrapper when
            # only the second parameter was BigFloat, causing recursion.
            expected = gaussjacobi(101, big"0.25", big"1.25")
            @test gaussjacobi(101, 0.25, big"1.25") == expected
            @test gaussjacobi(101, big"0.25", 1.25) == expected
            @test eltype(expected[1]) == eltype(expected[2]) == BigFloat
        end
    end

    @testset "Generalized Radau" begin
        setprecision(256) do
            n, a, b = 19, BigFloat(21) / 4, -BigFloat(1) / 2
            x, w = gaussradau(n, a, b)
            tolerance = BigFloat(2)^(-precision(BigFloat) + 16)
            @test eltype(x) == eltype(w) == BigFloat
            @test x[1] == -1
            @test all(z -> precision(z) == precision(BigFloat), x)
            @test all(z -> precision(z) == precision(BigFloat), w)
            @test issorted(x) && all(>(0), diff(x))
            @test all(isfinite, w) && all(>(0), w)
            edge = setprecision(384) do
                aa, bb = BigFloat(a), BigFloat(b)
                2^(aa + bb + 1) * (bb + 1) * beta(BigFloat(n), bb + 1) * beta(n + aa, bb + 1)
            end
            @test isapprox(w[1], edge; atol = 0, rtol = tolerance)
            @testset "i: $i" for i in (2, 3, 10, 18, 19)
                xr, wr = QuadratureReference.nodeweight(:jacobi, n - 1, a, b + 1, x[i]; bits = 256)
                @test isapprox(x[i], xr; atol = tolerance, rtol = tolerance)
                @test isapprox(w[i], wr / (1 + xr); atol = 0, rtol = tolerance)
            end
            @test gaussradau(n, Float64(a), b) == (x, w)
        end
    end
end
