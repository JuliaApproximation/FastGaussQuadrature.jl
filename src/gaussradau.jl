@doc raw"""
    gaussradau(n::Integer) -> Tuple{Vector{Float64},Vector{Float64}}

Return nodes and weights of [Gauss-Radau quadrature](https://mathworld.wolfram.com/RadauQuadrature.html).

```math
\int_{-1}^{1} f(x) dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gaussradau(3);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ 2/5
true
```

Note that the first node is fixed at -1.

```jldoctest
julia> x, w = gaussradau(3);

julia> x[1]
-1.0
```
"""
function gaussradau(n::Integer, T::Type=Float64)
    a = b = zero(T)
    # RADAUPTS   Gauss-Legendre-Radau Quadrature Nodes and Weights
    if n == 1
        return T[-1], T[2]
    elseif n == 2
        return [-1, one(T)/3], [one(T)/2, convert(T,3)/2]
    else
        # Compute via GaussJacobi:
        x, w = gaussjacobi(n - 1, a, b+1)
        @inbounds for i in 1:length(w)
            w[i] = w[i] / (1 + x[i])
        end
        pushfirst!(x, -1)
        pushfirst!(w, convert(T, 2) / n^2)
        return x, w
    end
end

function gaussradau(n::Integer, α, β)
    if n ≤ 0
        throw(DomainError(n, "Input N must be a positive integer"))
    end
    m = n - 1
    s = α + β
    T = float(eltype(s))
    μ = jacobimoment(α, β)
    n == 0 && return T[], T[]
    n == 1 && return [-one(T)], [μ]
    J = jacobi_jacobimatrix(n, α, β)
    aᴿ = -1 + 2m*convert(T,m+α)/((2m+s)*(2m+s+1))
    J.dv[end] = aᴿ
    x, V = eigen(J)
    w = V[1,:].^2 .* μ
    x[1] = -1 # fix rounding from eigen computation
    return x, w
end
