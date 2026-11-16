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

The order must satisfy `n ≥ 2`. The first and last nodes are fixed at -1 and 1.

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
        x, w = T[-1, 1], T[1, 1]
    elseif precision(T) < precision(Float32)
        x, w = gausslobatto(Float32, n)
        x, w = T.(x), T.(w)
    elseif n == 3
        x, w = T[-1, 0, 1], [T(1) / 3, T(4) / 3, T(1) / 3]
    else
        # Divide the Jacobi(1,1) weights by 1-x² before discarding the
        # accurate endpoint distances or angles. Keep Jacobi's route cutoffs.
        if T === Float64 && n >= 770
            x, w = jacobi_asy(n - 2, T(1), T(1), Val(:lobatto))
        elseif n <= (T === Float64 ? 26 : 22)
            x, w = gaussjacobi(n - 2, T(1), T(1))
            @inbounds for i in eachindex(w)
                w[i] /= (1 - x[i]) * (1 + x[i])
            end
        else
            x, w = jacobi_gst(n - 2, T(1), T(1), Val(:lobatto))
        end
        pushfirst!(x, T(-1))
        push!(x, T(1))
        pushfirst!(w, T(2) / (n * (n - 1)))
        push!(w, T(2) / (n * (n - 1)))
    end
    validate_rule(x, w)
    return x, w
end
gausslobatto(n::Integer) = gausslobatto(Float64, n)
