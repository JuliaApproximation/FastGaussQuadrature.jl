module FastGauss
   using Base

export GaussLegendre 
export GaussChebyshev 
export GaussJacobi
export GaussLobatto
export GaussRadau

include("GaussLegendre.jl")
include("GaussChebyshev.jl")
include("GaussJacobi.jl")
include("GaussLobatto.jl")
include("GaussRadau.jl")

end
