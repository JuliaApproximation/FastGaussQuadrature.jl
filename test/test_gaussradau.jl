@testset "Gauss–Radau" begin
    @testset "Check error" begin
        @test_throws DomainError gaussradau(0)
        @test_throws DomainError gaussradau(0, 0, 0)
    end

    @testset "Legendre" begin
        n = 1
        x, w = gaussradau(n)
        @test x[1] ≈ -1.0
        @test w[1] ≈ 2.0

        n = 2
        x, w = gaussradau(n)
        @test x[1] ≈ -1.0
        @test x[2] ≈ 1 / 3
        @test w[1] ≈ 0.5
        @test w[2] ≈ 1.5

        tol = 1.0e-14
        n = 42
        ntests = 4
        x, w = gaussradau(n)
        @test length(x) == n && length(w) == n
        @test abs(x[37] - 0.908847278001044) < tol
        @test abs(w[37] - 0.031190846817016) < tol
        @test x[1] == -1

        @test dot(w, exp.(x)) ≈ exp(1) - exp(-1)
    end
    @testset "Jacobi" begin
        for (a, b) in ((0.1, 0.2), (0, -0.5))
            # compare with Gauss–Jacobi
            for n in 1:5
                x, w = gaussradau(n, a, b)
                x̃, w̃ = gaussjacobi(n, a, b)

                @test length(x) == n

                @test sum(w) ≈ sum(w̃)
                for j in 1:(2n - 2)
                    @test dot(w, x .^ j) ≈ dot(w̃, x̃ .^ j)
                end
            end
        end

        # compare with Gauss-Radau-Legendre
        @test all(gaussradau(2, 0, 0) .≈ gaussradau(2))
        @test all(gaussradau(3, 0, 0) .≈ gaussradau(3))
    end

    @testset "endpoint weights retain their angles" begin
        # Independent 256-bit Jacobi-root/Legendre-weight references.
        @testset "n: $n" for (n, expected) in (
                (102, 0.0011847757070271125),
                (10000, 1.2329305785888396e-7),
            )
            x, w = gaussradau(n)
            @test x[1] == -1
            @test w[1] == 2 / n^2
            @test w[2] ≈ expected rtol = 1.0e-14
        end
    end

    @testset "Generalized Jacobi reduction" begin
        # Independent 384-bit Jacobi recurrence/Newton references with shifted
        # parameters, followed by w_Radau = w_Jacobi/(1+x).
        @testset "n: $n, α: $a, β: $b" for (n, a, b, node, weight) in (
                (9, 0.1, 0.2, -0.9012718224204942, 0.10328170379512788),
                (21, 10.0, 20.0, -0.7281206083106859, 1.154626687715557e-10),
                (100, 1000.0, 1000.0, -0.3943472906912523, 7.223159788983759e-76),
                (1000, 10.0, 20.0, -0.9996593387282864, 5.940833881112758e-71),
                (500, 0.1, 0.2, -0.9999664138596807, 7.271743612332555e-6),
                (501, 0.1, 0.2, -0.9999665477615041, 7.2369683663365145e-6),
            )
            x, w = gaussradau(n, a, b)
            @test x[1] == -1
            @test length(x) == length(w) == n
            @test issorted(x)
            @test all(isfinite, x)
            @test all(>(0), w)
            @test isapprox(x[2], node; rtol = 4eps(Float64))
            @test isapprox(w[2], weight; rtol = 2.0e-13)
            mass, endpoint = setprecision(256) do
                aa, bb = BigFloat(a), BigFloat(b)
                2^(aa + bb + 1) * beta(aa + 1, bb + 1),
                    2^(aa + bb + 1) * (bb + 1) * beta(BigFloat(n), bb + 1) * beta(n + aa, bb + 1)
            end
            @test isapprox(w[1], endpoint; rtol = 2.0e-13)
            @test isapprox(sum(w), mass; rtol = 2.0e-13)
            s, difference = a + b + 2, b - a
            @test abs(dot(w, x) - mass * difference / s) ≤ 2.0e-13mass
            @test isapprox(dot(w, x .^ 2), mass * (difference^2 + s) / (s * (s + 1)); rtol = 2.0e-13)
        end
        @testset "n: $n" for n in (102, 500)
            @test gaussradau(n, 0.0, 0.0)[2][1] == 2 / n^2
        end
        @testset "T: $T" for T in (Float32, BigFloat)
            x, w = gaussradau(5, T(0.1), T(0.2))
            @test eltype(x) == eltype(w) == T
            @test x[1] == -1
            @test issorted(x)
        end
        @test_throws DomainError gaussradau(5, -1.0, 0.0)
        @test_throws DomainError gaussradau(5, 0.0, -1.0)
    end

    @testset "Float32 normal tail weights" begin
        # Independent high-precision references place weights 2:100 in the
        # normal range and the fixed-endpoint weight in the subnormal range.
        # Eigenvector weights previously lost the seven normal entries 2:8.
        x, w = @inferred gaussradau(100, 10.0f0, 20.0f0)
        @test eltype(x) === eltype(w) === Float32
        @test length(x) == length(w) == 100
        @test 0 < w[1] < floatmin(Float32)
        @test all(>=(floatmin(Float32)), @view(w[2:end]))
    end

    @testset "Symmetric Jacobi interiors" begin
        # α=β+1 makes the interior Jacobi rule symmetric, while the Radau
        # weights remain asymmetric. Preserve the small left-end divisor.
        @testset "n: $n, α: $a" for (n, a) in (
                (21, 1.0), (22, 5.0), (102, 2.0^-53),
            )
            x, w = gaussradau(n, a, a - 1)
            @test x[1] == -1
            @test x[2:end] == -reverse(x[2:end])
            @test all(isfinite, w) && all(>(0), w)
            for i in (2, n ÷ 2 + 1, n)
                xr, wr = QuadratureReference.nodeweight(:jacobi, n - 1, a, a, x[i]; bits = 256)
                @test isapprox(x[i], xr; atol = 3eps(Float64), rtol = 0)
                @test isapprox(w[i], wr / (1 + xr); atol = 0, rtol = 2.0e-13)
            end
        end
    end

end

@testset "Radau near the upper-endpoint parameter limit" begin
    n, a, b = 501, nextfloat(-1.0), 0.2
    x, w = gaussradau(n, a, b)
    setprecision(BigFloat, 384) do
        aa, bb = BigFloat(a), BigFloat(b) + 1
        xr, vr = QuadratureReference.nodeweight(:jacobi, n - 1, aa, bb, x[end]; bits = 256)
        # The last node rounds to one; the normal weight must still retain
        # native accuracy. The former ASY weight was about four times too big.
        radius = (1 - xr) / 4
        @test QuadratureReference.roots_below(:jacobi, n - 1, aa, bb, xr - radius) == n - 2 &&
            QuadratureReference.roots_below(:jacobi, n - 1, aa, bb, xr + radius) == n - 1
        tolerance = QuadratureReference.native_jacobi_weight_budget(n - 1, a, bb, vr) + 8eps(Float64)
        @test isapprox(w[end], vr / (1 + xr); atol = 0, rtol = tolerance)
    end
end

@testset "Radau across Jacobi routing boundaries" begin
    # The unparameterized rule transforms Jacobi weights after computing its
    # interior nodes. Check both sides of each working-type REC/GST cutoff.
    for (T, orders) in ((Float32, (21, 22)), (Float64, (25, 26))), n in orders
        x, w = gaussradau(T, n)
        for i in (2, n)
            xr, wr = QuadratureReference.nodeweight(:jacobi, n - 1, 0, 1, x[i]; bits = 256)
            @test isapprox(x[i], xr; atol = 4eps(T), rtol = 0)
            @test isapprox(w[i], wr / (1 + xr); atol = 0, rtol = 128eps(T))
        end
    end
end

@testset "Low-precision Radau rules" begin
    for (n, a, b) in ((8, Float16(0), Float16(0)), (128, Float16(0.1), Float16(0.2)))
        x, w = @inferred gaussradau(n, a, b)
        @test eltype(x) === eltype(w) === Float16
        @test x[1] == -1 && all(>(0), diff(x))
        @test all(isfinite, w) && all(>=(0), w)
        if iszero(a) && iszero(b)
            @test gaussradau(Float16, n) == (x, w)
        end
        for i in (2, cld(n, 2), n)
            xr, wr = QuadratureReference.nodeweight(:jacobi, n - 1, a, BigFloat(b) + 1, x[i]; bits = 128)
            @test isapprox(x[i], xr; atol = eps(Float16), rtol = 0)
            @test isapprox(w[i], wr / (1 + xr); atol = eps(zero(Float16)), rtol = 2eps(Float16))
        end
        setprecision(256) do
            aa, bb = BigFloat(a), BigFloat(b)
            wr = 2^(aa + bb + 1) * (bb + 1) * beta(BigFloat(n), bb + 1) * beta(n + aa, bb + 1)
            @test isapprox(w[1], wr; atol = eps(zero(Float16)), rtol = 2eps(Float16))
        end
    end
    @test_throws ErrorException gaussradau(Float16, 512)
end
