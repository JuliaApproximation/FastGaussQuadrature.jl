module FastGauss
   using Base

export GaussLegendre 
export GaussChebyshev 

include("GaussLegendre.jl")
include("GaussChebyshev.jl")
end
