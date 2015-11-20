module FastGaussQuadrature

using Compat

export gausslegendre
export gausschebyshev
export gausshermite
export gaussjacobi
export gausslobatto
export gaussradau

include("gausslegendre.jl")
include("gausschebyshev.jl")
include("gausshermite.jl")
include("gaussjacobi.jl")
include("gausslobatto.jl")
include("gaussradau.jl")

end
