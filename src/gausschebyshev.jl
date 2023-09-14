@doc raw"""
    gausschebyshev1(n::Integer) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Chebyshev quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) of the 1st kind.

```math
\int_{-1}^{1} \frac{f(x)}{\sqrt{1-x^2}} dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gausschebyshev1(3);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ 3π/8
true
```
"""
function gausschebyshev1(n::Integer)
    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    end
    return [cos((2 * k - 1) * π / (2 * n)) for k = n:-1:1], fill(π / n, n)
end

@doc raw"""
    gausschebyshev2(n::Integer) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Chebyshev quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) of the 2nd kind.

```math
\int_{-1}^{1} f(x)\sqrt{1-x^2} dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gausschebyshev2(3);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ π/16
true
```
"""
function gausschebyshev2(n::Integer)
    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    end
    return [cos(k * π / (n + 1)) for k = n:-1:1], [π/(n + 1) * sin(k / (n + 1) * π)^2 for k = n:-1:1]
end

@doc raw"""
    gausschebyshev3(n::Integer) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Chebyshev quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) of the 3rd kind.

```math
\int_{-1}^{1} f(x)\sqrt{\frac{1+x}{1-x}} dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gausschebyshev3(3);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ 3π/8
true
```
"""
function gausschebyshev3(n::Integer)
    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    end
    return [cos((k - .5) * π / (n + .5)) for k = n:-1:1], [2π / (n + .5) * cos((k - .5) * π / (2 * (n + .5)))^2 for k = n:-1:1]
end

@doc raw"""
    gausschebyshev4(n::Integer) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Chebyshev quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) of the 4th kind.

```math
\int_{-1}^{1} f(x)\sqrt{\frac{1-x}{1+x}} dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gausschebyshev4(3);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ 3π/8
true
```
"""
function gausschebyshev4(n::Integer)
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
        # Gauss-ChebyshevT quadrature, i.e., w(x) = 1/sqrt(1-x^2)
        Base.depwarn("`gausschebyshev(n, 1)` is deprecated and will be removed in the next breaking release. Please use `gausschebyshev1(n)` instead.", :gausschebyshev)
        return gausschebyshev1(n)
    elseif kind == 2
        # Gauss-ChebyshevU quadrature, i.e., w(x) = sqrt(1-x^2)
        Base.depwarn("`gausschebyshev(n, 2)` is deprecated and will be removed in the next breaking release. Please use `gausschebyshev2(n)` instead.", :gausschebyshev)
        return gausschebyshev2(n)
    elseif kind == 3
        # Gauss-ChebyshevV quadrature, i.e., w(x) = sqrt((1+x)/(1-x))
        Base.depwarn("`gausschebyshev(n, 3)` is deprecated and will be removed in the next breaking release. Please use `gausschebyshev3(n)` instead.", :gausschebyshev)
        return gausschebyshev3(n)
    elseif kind == 4
        # Gauss-ChebyshevW quadrature, i.e., w(x) = sqrt((1-x)/(1+x))
        Base.depwarn("`gausschebyshev(n, 4)` is deprecated and will be removed in the next breaking release. Please use `gausschebyshev4(n)` instead.", :gausschebyshev)
        return gausschebyshev4(n)
    else
        throw(ArgumentError("Chebyshev kind should be 1, 2, 3, or 4"))
    end
end
