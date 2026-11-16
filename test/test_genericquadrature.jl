module MinimumPrecisionTest
    # The output type needs no elementary or special functions, and widening
    # must not influence the package's fixed Float32 working-precision policy.
    struct SmallFloat <: AbstractFloat
        value::Float16
    end
    SmallFloat(x::Real) = SmallFloat(Float16(x))
    Base.float(x::SmallFloat) = x
    Base.precision(::Type{SmallFloat}) = 11
    Base.widen(::Type{SmallFloat}) = error("widen must not be called")
    Base.Float32(x::SmallFloat) = Float32(x.value)
    Base.isfinite(x::SmallFloat) = isfinite(x.value)
    Base.zero(::Type{SmallFloat}) = SmallFloat(0)
    Base.:<(x::SmallFloat, y::SmallFloat) = x.value < y.value
    Base.:<=(x::SmallFloat, y::SmallFloat) = x.value <= y.value
    Base.:<(x::Integer, y::SmallFloat) = x < y.value
    Base.:<(x::SmallFloat, y::Integer) = x.value < y
    Base.:<=(x::Integer, y::SmallFloat) = x <= y.value
    Base.:<=(x::SmallFloat, y::Integer) = x.value <= y
end

@testset "Minimum Float32 working precision" begin
    rules = (
        W -> gausslegendre(W, 16),
        W -> gausshermite(W, 16),
        W -> gausshermite(W, 16; normalize = true),
        W -> gausslaguerre(16, W(0.25)),
        W -> gausslaguerre(16, W(0.25); reduced = true),
        W -> gaussjacobi(16, W(0.25), W(0.125)),
        W -> gaussradau(W, 16),
        W -> gaussradau(16, W(0.25), W(0.125)),
        W -> gausslobatto(W, 16),
        W -> gausschebyshevt(W, 16),
        W -> gausschebyshevu(W, 16),
        W -> gausschebyshevv(W, 16),
        W -> gausschebyshevw(W, 16),
    )
    @testset "Output type: $T" for T in (Float16, MinimumPrecisionTest.SmallFloat)
        for rule in rules
            x, w = rule(T)
            xr, wr = rule(Float32)
            @test eltype(x) == eltype(w) == T
            @test Float32.(x) == Float32.(T.(xr))
            @test Float32.(w) == Float32.(T.(wr))
        end
    end

    # Exercise the low end and the boundary immediately below native precision.
    # Intermediate floor conversions use the same path; Legendre separately
    # covers 16-bit outputs across its recurrence/expansion cutoff.
    @testset "BigFloat: $bits bits" for bits in (8, 23)
        setprecision(BigFloat, bits) do
            for rule in rules
                x, w = rule(BigFloat)
                xr, wr = rule(Float32)
                @test x == BigFloat.(xr)
                @test w == BigFloat.(wr)
                @test all(t -> precision(t) == bits, x)
                @test all(t -> precision(t) == bits, w)
            end
            @test gaussjacobi(16, Float16(0.25), BigFloat(0.125)) ==
                gaussjacobi(16, BigFloat(0.25), BigFloat(0.125))
            @test gaussradau(16, Float16(0.25), BigFloat(0.125)) ==
                gaussradau(16, BigFloat(0.25), BigFloat(0.125))
            # Float32 working range applies even when BigFloat could represent
            # the final weight. Validation rejects the working-type overflow.
            @test_throws ErrorException gausslaguerre(1, BigFloat(64))
            @test precision(BigFloat) == bits
        end
    end
    @test_throws ErrorException gausslaguerre(1, Float32(64))

    @testset "24-bit BigFloat keeps its own arithmetic and range" begin
        setprecision(BigFloat, 24) do
            x, w = gausslaguerre(1, BigFloat(64))
            @test x == [BigFloat(65)]
            @test isfinite(w[1]) && w[1] > floatmax(Float32)
            @test precision(x[1]) == precision(w[1]) == 24
        end
    end

    @testset "Fourth-kind reflection precedes directed rounding" begin
        setprecision(BigFloat, 8) do
            setrounding(BigFloat, RoundDown) do
                x, w = gausschebyshevw(BigFloat, 16)
                xr, wr = gausschebyshevw(Float32, 16)
                @test x == BigFloat.(xr)
                @test w == BigFloat.(wr)
            end
        end
    end

    @testset "Numerical accuracy and output conversion" begin
        # Independent moments exercise numerical accuracy after conversion.
        x, w = gausslegendre(Float16, 16)
        @test sum(Float64, w) ≈ 2 atol = 2eps(Float16)
        @test dot(Float64.(w), Float64.(x) .^ 2) ≈ 2 / 3 atol = 2eps(Float16)
        # Validate after conversion: merged nodes and overflowing weights fail.
        @test_throws ErrorException gausslegendre(Float16, 256)
        x, w = gausslaguerre(Float16, 32)
        @test any(iszero, w) && any(>(0), w)
        @test sum(Float64, w) ≈ 1 atol = 2eps(Float16)
    end
end
