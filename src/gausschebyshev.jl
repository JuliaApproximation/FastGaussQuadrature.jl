@doc raw"""
    gausschebyshevt(n::Integer) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Chebyshev quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) of the 1st kind.

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
function gausschebyshevt(n::Integer)
    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    end
    return [cos((2 * k - 1) * π / (2 * n)) for k = n:-1:1], fill(π / n, n)
end

@doc raw"""
    gausschebyshevu(n::Integer) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Chebyshev quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) of the 2nd kind.

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
function gausschebyshevu(n::Integer)
    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    end
    return [cos(k * π / (n + 1)) for k = n:-1:1], [π/(n + 1) * sin(k / (n + 1) * π)^2 for k = n:-1:1]
end

@doc raw"""
    gausschebyshevv(n::Integer) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Chebyshev quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) of the 3rd kind.

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
function gausschebyshevv(n::Integer)
    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    end
    return [cos((k - .5) * π / (n + .5)) for k = n:-1:1], [2π / (n + .5) * cos((k - .5) * π / (2 * (n + .5)))^2 for k = n:-1:1]
end

@doc raw"""
    gausschebyshevw(n::Integer) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Chebyshev quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) of the 4th kind.

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
function gausschebyshevw(n::Integer)
    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    end
    return [cos(k * π / (n + .5)) for k = n:-1:1], [2π / (n + .5) * sin(k * π / (2 * (n + .5)))^2 for k = n:-1:1]
end

function gausschebyshev(n::Integer, kind::Integer=1)
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
