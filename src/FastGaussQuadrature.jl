module FastGaussQuadrature
   using Base

export gausslegendre 
export gausschebyshev 
export gaussjacobi
export gausslobatto
export gaussradau

include("gausslegendre.jl")
include("gausschebyshev.jl")
include("gaussjacobi.jl")
include("gausslobatto.jl")
include("gaussradau.jl")

end
