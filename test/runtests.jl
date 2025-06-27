using FastGaussQuadrature
using Test, Aqua, LinearAlgebra, Random, SpecialFunctions

Aqua.test_all(FastGaussQuadrature;deps_compat = (; check_extras=false))

@testset "FastGaussQuadrature.jl" begin
    include("test_gausschebyshev.jl")
    include("test_gausslegendre.jl")
    include("test_gaussjacobi.jl")
    include("test_gaussradau.jl")
    include("test_gausslobatto.jl")
    include("test_besselroots.jl")
    include("test_gausslaguerre.jl")
    include("test_gausshermite.jl")
end
