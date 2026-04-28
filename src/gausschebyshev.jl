@doc raw"""
    gausschebyshevt([T=Float64,] n::Integer) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Chebyshev quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) of the 1st kind and type `T`.

```math
\int_{-1}^{1} \frac{f(x)}{\sqrt{1-x^2}} dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gausschebyshevt(3);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ 3π/8
true
```
"""
function gausschebyshevt(::Type{T}, n::Integer) where T
    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    end
    return [cospi(T(2 * k - 1) / (2 * n)) for k in n:-1:1], fill(T(π) / n, n)
end
gausschebyshevt(n::Integer) = gausschebyshevt(Float64, n)

@doc raw"""
    gausschebyshevu([T=Float64,] n::Integer) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Chebyshev quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) of the 2nd kind and type `T`.

```math
\int_{-1}^{1} f(x)\sqrt{1-x^2} dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gausschebyshevu(3);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ π/16
true
```
"""
function gausschebyshevu(::Type{T}, n::Integer) where T
    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    end
    return [cospi(T(k) / (n + 1)) for k in n:-1:1], [T(π) / (n + 1) * sinpi(T(k) / (n + 1) )^2 for k in n:-1:1]
end
gausschebyshevu(n::Integer) = gausschebyshevu(Float64, n)

@doc raw"""
    gausschebyshevv([T=Float64,] n::Integer) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Chebyshev quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) of the 3rd kind and type T.

```math
\int_{-1}^{1} f(x)\sqrt{\frac{1+x}{1-x}} dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gausschebyshevv(3);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ 3π/8
true
```
"""
function gausschebyshevv(::Type{T}, n::Integer) where T
    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    end
    return [cospi(T(2k - 1) / (2n + 1)) for k in n:-1:1], [4T(π) / (2n + 1) * cospi(T(2k - 1) / (4n + 2))^2 for k in n:-1:1]
end
gausschebyshevv(n::Integer) = gausschebyshevv(Float64, n)

@doc raw"""
    gausschebyshevw(T=Float64, n::Integer) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Chebyshev quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) of the 4th kind and type T.

```math
\int_{-1}^{1} f(x)\sqrt{\frac{1-x}{1+x}} dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gausschebyshevw(3);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ 3π/8
true
```
"""
function gausschebyshevw(::Type{T}, n::Integer) where T
    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    end
    return [cospi(T(2k) / (2n + 1)) for k in n:-1:1], [4T(π) / (2n + 1) * sinpi(k / T(2n + 1))^2 for k in n:-1:1]
end
gausschebyshevw(n::Integer) = gausschebyshevw(Float64, n)

function gausschebyshev(n::Integer, kind::Integer = 1)
    # GAUSS-CHEBYSHEV NODES AND WEIGHTS.

    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    end

    # Use known explicit formulas. Complexity O(n).
    if kind == 1
        # Gauss-Chebyshevt quadrature, i.e., w(x) = 1/sqrt(1-x^2)
        Base.depwarn("`gausschebyshev(n, 1)` is deprecated and will be removed in the next breaking release. Please use `gausschebyshevt(n)` instead.", :gausschebyshev)
        return gausschebyshevt(n)
    elseif kind == 2
        # Gauss-Chebyshevu quadrature, i.e., w(x) = sqrt(1-x^2)
        Base.depwarn("`gausschebyshev(n, 2)` is deprecated and will be removed in the next breaking release. Please use `gausschebyshevu(n)` instead.", :gausschebyshev)
        return gausschebyshevu(n)
    elseif kind == 3
        # Gauss-Chebyshevv quadrature, i.e., w(x) = sqrt((1+x)/(1-x))
        Base.depwarn("`gausschebyshev(n, 3)` is deprecated and will be removed in the next breaking release. Please use `gausschebyshevv(n)` instead.", :gausschebyshev)
        return gausschebyshevv(n)
    elseif kind == 4
        # Gauss-Chebyshevw quadrature, i.e., w(x) = sqrt((1-x)/(1+x))
        Base.depwarn("`gausschebyshev(n, 4)` is deprecated and will be removed in the next breaking release. Please use `gausschebyshevw(n)` instead.", :gausschebyshev)
        return gausschebyshevw(n)
    else
        throw(ArgumentError("Chebyshev kind should be 1, 2, 3, or 4"))
    end
end
