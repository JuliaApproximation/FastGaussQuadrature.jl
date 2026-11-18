@testset "Quadrature rule validation" begin
    validate_rule = FastGaussQuadrature.validate_rule

    @testset "Valid rules" begin
        @test validate_rule(Float64[], Float64[]) === nothing
        @test validate_rule([0.0], [2.0]) === nothing
        # The shared validator applies to unbounded and mixed-type rules too.
        @test validate_rule([-10.0, 0.0, 20.0], Float32[0, 1, 0]) === nothing
        @test validate_rule([typemin(Int), typemax(Int)], [1, 1]) === nothing
        @test validate_rule([0.0, 1.0], [nextfloat(0.0), 0.0]) === nothing

        x, w = [-2.0, 0.0, 2.0], [0.25, 0.5, 0.25]
        x_before, w_before = copy(x), copy(w)
        @test validate_rule(x, w) === nothing
        @test x == x_before && w == w_before
        @test validate_rule(@view(x[1:2]), @view(w[1:2])) === nothing
    end

    @testset "Invalid rules" begin
        @test_throws DimensionMismatch validate_rule([0.0], Float64[])
        @test_throws DimensionMismatch validate_rule(Float64[], [1.0])
        for x in ([NaN], [Inf], [-Inf], [0.0, NaN], [-Inf, 0.0])
            @test_throws ErrorException validate_rule(x, ones(length(x)))
        end
        for x in ([0.0, 0.0], [1.0, -1.0], [-0.0, 0.0], [0.0, 2.0, 1.0])
            @test_throws ErrorException validate_rule(x, ones(length(x)))
        end
        for w in ([NaN], [Inf], [-Inf], [-1.0], [0.0], [-0.0])
            @test_throws ErrorException validate_rule([0.0], w)
        end
        @test_throws ErrorException validate_rule([0.0, 1.0], [1.0, -eps()])
        @test_throws ErrorException validate_rule([0.0, 1.0], [0.0, 0.0])
    end
end
