using FastGaussQuadrature
using Test, LinearAlgebra, Random, SpecialFunctions

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
