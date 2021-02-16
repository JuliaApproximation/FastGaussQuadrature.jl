module FastGaussQuadrature

using LinearAlgebra
using SpecialFunctions
using StaticArrays

export gausslegendre
export gausschebyshev
export gausslaguerre
export gausshermite
export gaussjacobi
export gausslobatto
export gaussradau
export approx_besselroots
export besselroots

import SpecialFunctions: besselj, airyai, airyaiprime

include("constants.jl")
include("gausslegendre.jl")
include("gausschebyshev.jl")
include("gausslaguerre.jl")
include("gausshermite.jl")
include("gaussjacobi.jl")
include("gausslobatto.jl")
include("gaussradau.jl")
include("besselroots.jl")

end
