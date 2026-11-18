using FastGaussQuadrature
using Test, Aqua, LinearAlgebra, Random, SpecialFunctions

include("quadrature_reference.jl")

Aqua.test_all(FastGaussQuadrature)

@testset "FastGaussQuadrature.jl" begin
    include("test_validate_rule.jl")
    include("test_genericquadrature.jl")
    include("test_genericprecision.jl")
    include("test_gausschebyshev.jl")
    include("test_gausslegendre.jl")
    include("test_legendre_nativeprecision.jl")
    include("test_gaussjacobi.jl")
    include("test_gaussradau.jl")
    include("test_gausslobatto.jl")
    include("test_besselroots.jl")
    include("test_gausslaguerre.jl")
    include("test_gausshermite.jl")

    include("test_bigfloatcoefs.jl")
    # Default rules and higher-precision references need no generic eigensolver.
    @test !Base.root_module_exists(Base.PkgId(Base.UUID("14197337-ba66-59df-a3e3-ca00e7dcff7a"), "GenericLinearAlgebra"))
end
