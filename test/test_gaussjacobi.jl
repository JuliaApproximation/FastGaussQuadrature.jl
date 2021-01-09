@testset "Gauss–Jacobi" begin
    @testset "check error" begin
        @test_throws DomainError gaussjacobi(-1, 0.1, 0.1)
        @test_throws DomainError gaussjacobi(-1, -1.1, 0.1)
        @test_throws DomainError gaussjacobi(-1, 0.1, -1.1)
    end

    tol = 1e-14
    tol_x = 4e-15
    tol_w = 1e-13

    n_range = vcat(0:10, 20:10:100, 120:20,300)

    @testset "compatible with Gauss-Legendre" begin
        α₀ = 0.0
        β₀ = 0.0
        α₊ = nextfloat(α₀)
        β₊ = nextfloat(β₀)
        α₋ = prevfloat(α₀)
        β₋ = prevfloat(β₀)
        for n in n_range
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
        α₀ = -1/2
        β₀ = -1/2
        α₊ = nextfloat(α₀)
        β₊ = nextfloat(β₀)
        α₋ = prevfloat(α₀)
        β₋ = prevfloat(β₀)
        for n in n_range
            x, w = gausschebyshev(n, 1)
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
        α₀ = 1/2
        β₀ = 1/2
        α₊ = nextfloat(α₀)
        β₊ = nextfloat(β₀)
        α₋ = prevfloat(α₀)
        β₋ = prevfloat(β₀)
        for n in n_range
            x, w = gausschebyshev(n, 2)
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
        α₀ = -1/2
        β₀ = 1/2
        α₊ = nextfloat(α₀)
        β₊ = nextfloat(β₀)
        α₋ = prevfloat(α₀)
        β₋ = prevfloat(β₀)
        for n in n_range
            x, w = gausschebyshev(n, 3)
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
        α₀ = 1/2
        β₀ = -1/2
        α₊ = nextfloat(α₀)
        β₊ = nextfloat(β₀)
        α₋ = prevfloat(α₀)
        β₋ = prevfloat(β₀)
        for n in n_range
            x, w = gausschebyshev(n, 4)
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
        a = -.1
        b = .3
        x, w = gaussjacobi(n, a, b)
        @test all(length(x) == n) && all(length(w) == n)
        @test abs(x[37] - 0.912883347814032) < tol
        @test abs(w[37] - 0.046661910947553) < tol

        @test dot(w,exp.(x)) ≈ 2.7568520361985516
    end
    @testset "a larger n (using ASY)" begin
        a = -.7
        b = 1.3
        n = 251
        x, w = gaussjacobi(n, a, b)
        @test all(length(x) == n) && all(length(w) == n)
        @test abs(x[37] + 0.893103435898983) < tol
        @test abs(w[37] - 1.962534523788093e-04) < tol

        @test dot(w,exp.(x)) ≈ 16.722957039404044
    end
    @testset "n = 1:" begin
        a = 1.0
        b = 2.0
        x, w = gaussjacobi(1, a, b)
        @test abs(x[1] - (b - a) / (a + b + 2)) < tol
        @test abs(w[1] - 2^(a + b + 1) * beta(a + 1, b + 1)) < tol

        @test dot(w,fill(1,length(x))) ≈ 1.3333333333333333
    end
    @testset "0.9, -0.1" begin
        x, w = gaussjacobi(1013, .9, -.1)
        @test abs(x[2] + 0.999986012231899) < tol
        @test abs(x[13] + 0.999225722939832) < tol
        @test abs(w[2] - 9.314674169892358e-05) < tol
        @test abs(w[13] - 4.654651764553262e-04) < tol

        @test dot(w,exp.(x)) ≈ 1.6915068974063106

        x, w = gaussjacobi(10013, .9, -.1)
        @test abs(x[2]  +0.999999856605293) < tol
        @test abs(x[13]  + 0.999992061552711) < tol
        @test abs(w[2] - 1.509654630405615e-06) < tol
        @test abs(w[13] - 7.548275262993863e-06) < tol

        @test dot(w,exp.(x)) ≈ 1.6915068974063106
    end
    @testset "bug where row vectors were returned" begin
        x, w = gaussjacobi(10013, .9, -.1)
        @test isa(x,Vector{Float64})
        @test isa(w,Vector{Float64})
    end
    @testset "last alpha and beta parameters:" begin
        x, w = gaussjacobi(100, 19., 21.)
        @test abs(x[87] - 0.832211446176040) < tol
        @test abs(w[50] - 0.026363584978877305) < tol # new test without condition sum(w)==1
    end
    @testset "for small alpha and beta:" begin
        x, w = gaussjacobi(10000, .1, .2)
        @test abs(x[1] - -0.999999963363548) < tol
        @test abs(w[500] - 2.183393039546711e-05) < tol
    end
    @testset " bug where row vectors were returned" begin
        x,w=gaussjacobi(20,11.,0.)
        @test isa(x,Vector{Float64})
        @test isa(w,Vector{Float64})
    end

    @testset "jacobi_rec against a precomputed rule in BigFloat" begin
        x,w = FastGaussQuadrature.jacobi_rec(10, big(0), big(0))
        @test abs(x[6] - big(1488743389816312108848260011297199846175648594206916957079892535159036173556674)/big(10)^79) < 1e-70
        @test abs(w[4] - big(2692667193099963550912269215694693528597599384608837958005632762421534323191819)/big(10)^79) < 1e-70

        x,w = FastGaussQuadrature.jacobi_rec(10, big(2//10), big(-1//30))
        @test abs(x[7] - big(4146701176053244724446765533149623814305606864483746748725151214196203624388371)/big(10)^79) < 1e-70
        @test abs(w[3] - big(2482452398859023537636458227096017466104571146686194593251519240011574719667464)/big(10)^79) < 1e-70
    end
end