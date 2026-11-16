@testset "Gauss–Chebyshev" begin

    @testset "Check error" begin
        @test_throws DomainError gausschebyshevt(-1)
        @test_throws DomainError gausschebyshevu(-1)
        @test_throws DomainError gausschebyshevv(-1)
        @test_throws DomainError gausschebyshevw(-1)
    end

    n = 10

    @testset "x.^2" begin
        x, w = gausschebyshevt(n)
        @test dot(x .^ 2, w) ≈ π / 2
        x, w = gausschebyshevu(n)
        @test dot(x .^ 2, w) ≈ π / 8
        x, w = gausschebyshevv(n)
        @test dot(x .^ 2, w) ≈ π / 2
        x, w = gausschebyshevw(n)
        @test dot(x .^ 2, w) ≈ π / 2
    end

    @testset "x^3" begin
        x, w = gausschebyshevt(n)
        @test abs(dot(x .^ 3, w)) < 1.0e-15
        x, w = gausschebyshevu(n)
        @test abs(dot(x .^ 3, w)) < 1.0e-15
        x, w = gausschebyshevv(n)
        @test dot(x .^ 3, w) ≈ 3 * π / 8
        x, w = gausschebyshevw(n)
        @test dot(x .^ 3, w) ≈ -3 * π / 8
    end

    @testset "deprecated" begin
        n = 42
        @test gausschebyshevt(n) == @test_deprecated gausschebyshev(n, 1)
        @test gausschebyshevu(n) == @test_deprecated gausschebyshev(n, 2)
        @test gausschebyshevv(n) == @test_deprecated gausschebyshev(n, 3)
        @test gausschebyshevw(n) == @test_deprecated gausschebyshev(n, 4)
        @test_throws ArgumentError gausschebyshev(0, 5)
    end

    rules = (gausschebyshevt, gausschebyshevu, gausschebyshevv, gausschebyshevw)

    @testset "Inputs and compatibility" begin
        for f in rules
            x, w = @inferred f(Float64, 0)
            @test isempty(x) && isempty(w) && eltype(x) == eltype(w) == Float64
            @test f(Float16, 0) == (Float16[], Float16[])
            # A signed phase must not wrap when the order is unsigned.
            @test f(UInt128(8)) == f(8)
        end
    end

    @testset "Weighted moments" begin
        # Exact integrals of x^k times each weight, divided by π.
        moments = (
            (1, 0, 1 // 2, 0, 3 // 8, 0, 5 // 16, 0),
            (1 // 2, 0, 1 // 8, 0, 1 // 16, 0, 5 // 128, 0),
            (1, 1 // 2, 1 // 2, 3 // 8, 3 // 8, 5 // 16, 5 // 16, 35 // 128),
            (1, -1 // 2, 1 // 2, -3 // 8, 3 // 8, -5 // 16, 5 // 16, -35 // 128),
        )
        # An n-point Gaussian rule integrates degrees through 2n - 1.
        @testset "kind: $kind, n: $n" for (kind, f) in enumerate(rules), n in (1, 4)
            x, w = f(n)
            degrees = n == 1 ? (0, 1) : (2, 3, 6, 7)
            for k in degrees
                expected = moments[kind][k + 1] * π
                @test isapprox(dot(x .^ k, w), expected; atol = 16eps(Float64) * π, rtol = 16eps(Float64))
            end
        end
    end

    @testset "Symmetry and reflection" begin
        @testset "n: $n" for n in (8, 9)
            for f in (gausschebyshevt, gausschebyshevu)
                x, w = @inferred f(Float64, n)
                @test x == -reverse(x)
                @test w == reverse(w)
                @test issorted(x) && all(>(0), w)
                isodd(n) && @test iszero(x[div(n + 1, 2)])
            end
            xv, wv = @inferred gausschebyshevv(Float64, n)
            xw, ww = @inferred gausschebyshevw(Float64, n)
            @test issorted(xv) && all(>(0), wv)
            @test xv == -reverse(xw)
            @test wv == reverse(ww)
        end
    end

    # Independent high-precision formulas in terms of the original angles.
    function reference(kind, n, i)
        k = BigFloat(n + 1 - i)
        theta = if kind == 1
            big(π) * (2k - 1) / (2n)
        elseif kind == 2
            big(π) * k / (n + 1)
        elseif kind == 3
            big(π) * (2k - 1) / (2n + 1)
        else
            big(π) * 2k / (2n + 1)
        end
        x = cos(theta)
        w = if kind == 1
            big(π) / n
        elseif kind == 2
            big(π) / (n + 1) * sin(theta)^2
        elseif kind == 3
            4big(π) / (2n + 1) * cos(theta / 2)^2
        else
            4big(π) / (2n + 1) * sin(theta / 2)^2
        end
        return x, w
    end

    @testset "Phase and weight accuracy" begin
        setprecision(BigFloat, 256) do
            # Endpoints, both sides of the formula switch, and near-zero nodes.
            # The fourth kind is also covered by the reflection checks above.
            @testset "T: $T, n: $n" for (T, n) in ((Float32, 1000), (Float64, 1001))
                @testset "kind: $kind" for (kind, f) in enumerate(rules[1:3])
                    x, w = f(T, n)
                    for i in (1, div(n, 4), div(n, 4) + 1, div(n, 2))
                        xr, wr = reference(kind, n, i)
                        @test abs(x[i] - xr) ≤ 3eps(T) * abs(xr)
                        # First-kind weights are constant across the rule.
                        if kind != 1 || i == 1
                            @test abs(w[i] - wr) ≤ 6eps(T) * wr
                        end
                    end
                end
            end
        end
        # Compare 128-bit results with independently evaluated 256-bit formulas.
        setprecision(BigFloat, 128) do
            unit = eps(BigFloat)
            for (kind, f) in enumerate(rules)
                x, w = f(BigFloat, 100)
                @test all(z -> precision(z) == 128, x) && all(z -> precision(z) == 128, w)
                setprecision(BigFloat, 256) do
                    xr, _ = reference(kind, 100, 50)
                    i = kind == 4 ? 100 : 1
                    _, wr = reference(kind, 100, i)
                    @test abs(x[50] - xr) ≤ 4unit * abs(xr)
                    @test abs(w[i] - wr) ≤ 6unit * wr
                end
            end
        end
    end

    @testset "Limited output precision" begin
        # Validation rejects nearby nodes rounded to the same output value.
        for f in (gausschebyshevu, gausschebyshevv)
            @test_throws ErrorException f(Float16, 256)
        end
        @test_throws ErrorException gausschebyshevv(Float16, 1000)
        setprecision(BigFloat, 4) do
            @test_throws ErrorException gausschebyshevu(BigFloat, 256)
        end
    end
end
