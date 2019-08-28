module FastGaussQuadrature

using SpecialFunctions, LinearAlgebra
cumprod(A::AbstractArray) = Base.cumprod(A, dims=1)
cumprod(A::AbstractArray, d::Int) = Base.cumprod(A, dims=d)
sum(A::AbstractArray, n::Int) = Base.sum(A, dims=n)
sum(A) = Base.sum(A)
flipdim(A, d) = reverse(A, dims=d)



export gausslegendre
export gausschebyshev
export gausslaguerre
export gausshermite
export gaussjacobi
export gausslobatto
export gaussradau
export besselroots

import SpecialFunctions: besselj, airyai, airyaiprime

include("gausslegendre.jl")
include("gausschebyshev.jl")
include("gausslaguerre.jl")
include("gausshermite.jl")
include("gaussjacobi.jl")
include("gausslobatto.jl")
include("gaussradau.jl")
include("besselroots.jl")

end
