@doc raw"""
    gausschebyshev(n::Integer) -> Tuple{Vector{Float64},Vector{Float64}}
    gausschebyshev(n::Integer, 1) -> Tuple{Vector{Float64},Vector{Float64}}

Return nodes and weights of [Gauss-Chebyshev quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) of the 1st kind.

```math
\int_{-1}^{1} \frac{f(x)}{\sqrt{1-x^2}} dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gausschebyshev(3);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ 3π/8
true
```

---

    gausschebyshev(n::Integer, 2) -> Tuple{Vector{Float64},Vector{Float64}}

Return nodes and weights of [Gauss-Chebyshev quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) of the 2nd kind.

```math
\int_{-1}^{1} f(x)\sqrt{1-x^2} dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gausschebyshev(3, 2);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ π/16
true
```

---

    gausschebyshev(n::Integer, 3) -> Tuple{Vector{Float64},Vector{Float64}}

Return nodes and weights of [Gauss-Chebyshev quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) of the 3rd kind.

```math
\int_{-1}^{1} f(x)\sqrt{1-x^2} dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gausschebyshev(3, 3);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ 3π/8
true
```

---

    gausschebyshev(n::Integer, 4) -> Tuple{Vector{Float64},Vector{Float64}}

Return nodes and weights of [Gauss-Chebyshev quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) of the 4th kind.

```math
\int_{-1}^{1} f(x)\sqrt{1-x^2} dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gausschebyshev(3, 4);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ 3π/8
true
```
"""
function gausschebyshev(n::Integer, kind::Integer=1)
    # GAUSS-CHEBYSHEV NODES AND WEIGHTS.

    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    end

    # Use known explicit formulas. Complexity O(n).
    if kind == 1
        # Gauss-ChebyshevT quadrature, i.e., w(x) = 1/sqrt(1-x^2)
        return ([cos((2 * k - 1) * π / (2 * n)) for k = n:-1:1], fill(π / n, n))
    elseif kind == 2
        # Gauss-ChebyshevU quadrature, i.e., w(x) = sqrt(1-x^2)
        return ([cos(k * π / (n + 1)) for k = n:-1:1],
                [π/(n + 1) * sin(k / (n + 1) * π)^2 for k = n:-1:1])
    elseif kind == 3
        # Gauss-ChebyshevV quadrature, i.e., w(x) = sqrt((1+x)/(1-x))
        return ([cos((k - .5) * π / (n + .5)) for k = n:-1:1],
                [2π / (n + .5) * cos((k - .5) * π / (2 * (n + .5)))^2 for k = n:-1:1])
    elseif kind == 4
        # Gauss-ChebyshevW quadrature, i.e., w(x) = sqrt((1-x)/(1+x))
        return ([cos(k * π / (n + .5)) for k = n:-1:1],
                [2π / (n + .5) * sin(k * π / (2 * (n + .5)))^2 for k = n:-1:1])
    else
        throw(ArgumentError("Chebyshev kind should be 1, 2, 3, or 4"))
    end
end
