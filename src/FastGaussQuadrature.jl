module FastGaussQuadrature

using LinearAlgebra: I, dot, rmul!
using SpecialFunctions: airyaiprime, besselj, beta, gamma, logbeta, loggamma
using StaticArrays: @SMatrix, @SVector, SVector, setindex

export gausslegendre
export gausschebyshev
export gausschebyshevt
export gausschebyshevu
export gausschebyshevv
export gausschebyshevw
export gausslaguerre
export gausshermite
export gaussjacobi
export gausslobatto
export gaussradau
export approx_besselroots
export besselroots

"""
    validate_rule(x, w) -> nothing

Check that the node and weight vectors have equal lengths, the nodes are finite
and strictly increasing, and the weights are finite and nonnegative. A nonempty
rule must have at least one positive weight; individual zero weights are allowed
because tail weights can underflow. Empty rules are valid.

Throw `DimensionMismatch` for unequal lengths and an error for invalid nodes or
weights. The vectors are not mutated. This checks representation and ordering,
not the integration domain or the accuracy of the rule's moments.
"""
function validate_rule(x::AbstractVector, w::AbstractVector)
    length(x) == length(w) || throw(DimensionMismatch("Quadrature nodes and weights must have equal lengths."))
    all(isfinite, x) || error("Quadrature nodes must be finite; use a type with a wider exponent range.")
    for i in firstindex(x):(lastindex(x) - 1)
        x[i] < x[i + 1] || error("Quadrature nodes must be strictly increasing; use a type with greater precision.")
    end
    all(isfinite, w) && all(>=(0), w) ||
        error("Quadrature weights must be finite and nonnegative; use a type with greater precision or a wider exponent range.")
    isempty(w) || any(>(0), w) || error("A nonempty quadrature rule must have at least one positive weight.")
    return nothing
end

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
