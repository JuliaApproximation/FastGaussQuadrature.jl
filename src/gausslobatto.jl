@doc raw"""
    gausslobatto([T=Float64,] n::Integer) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Lobatto quadrature](https://mathworld.wolfram.com/LobattoQuadrature.html) with type `T`.

```math
\int_{-1}^{1} f(x) dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gausslobatto(4);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ 2/5
true
```

Note that the both ends of nodes are fixed at -1 and 1.

```jldoctest
julia> x, w = gausslobatto(4);

julia> x[1], x[end]
(-1.0, 1.0)
```
"""
function gausslobatto(::Type{T}, n::Integer) where {T}
    # Gauss-Legendre-Lobatto Quadrature Nodes and Weights
    if n ≤ 1
        throw(DomainError(n, "Lobatto undefined for n ≤ 1."))
    elseif n == 2
        return T[-1, 1], T[1, 1]
    elseif n == 3
        return T[-1, 0, 1], [T(1) / 3, T(4) / 3, T(1) / 3]
    else
        # Compute via GaussJacobi:
        x, w = gaussjacobi(n - 2, T(1), T(1))
        @inbounds for i in 1:length(x)
            w[i] = w[i] / (1 - x[i]^2)
        end
        pushfirst!(x, T(-1))
        push!(x, T(1))
        pushfirst!(w, T(2) / (n * (n - 1)))
        push!(w, T(2) / (n * (n - 1)))
        return x, w
    end
end
gausslobatto(n::Integer) = gausslobatto(Float64, n)
