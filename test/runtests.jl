using FastGaussQuadrature, Compat
using Compat.Test, Compat.LinearAlgebra, Compat.Random, SpecialFunctions

if VERSION < v"0.7-"
    const floatmin = realmin
    const floatmax = realmax
end


println("Chebyshev tests")
include("test_gausschebyshev.jl")
println("Legendre tests")
include("test_gausslegendre.jl")
println("Jacobi tests")
include("test_gaussjacobi.jl")
println("Radau tests")
include("test_gaussradau.jl")
println("Lobatto tests")
include("test_gausslobatto.jl")
println("Bessel roots tests")
include("test_besselroots.jl")
println("Laguerre tests")
include("test_gausslaguerre.jl")
println("Hermite tests")
include("test_gausshermite.jl")
