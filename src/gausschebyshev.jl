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
function gausschebyshevt(::Type{T}, n::Integer) where {T}
    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    elseif n == 0
        x, w = T[], T[]
        validate_rule(x, w)
        return x, w
    elseif precision(T) < precision(Float32)
        x, w = gausschebyshevt(Float32, n)
        x, w = T.(x), T.(w)
        validate_rule(x, w)
        return x, w
    end
    # For i = 1,...,n, θ_i = (2i - 1)π/(2n):
    # x_i = -cos(θ_i), w_i = π/n, with nodes in ascending order.
    # Use x_{n+1-i} = -x_i to fill both halves with one trig evaluation.
    x = Vector{T}(undef, n)
    w = fill(T(π) / n, n)
    # Switch from cos(θ_i) to sin(π/2 - θ_i) near θ_i = π/4.
    @inbounds for i in 1:div(n, 4)
        x[n + 1 - i] = cospi(T(2i - 1) / T(2n))
        x[i] = -x[n + 1 - i]
    end
    @inbounds for i in (div(n, 4) + 1):div(n, 2)
        # φ_i = π/2 - θ_i = (n + 1 - 2i)π/(2n), so cos(θ_i) = sin(φ_i).
        # Form the small complementary phase directly near the centre, avoiding
        # subtraction of rounded angles close to π/2.
        x[n + 1 - i] = sinpi(T(n + 1 - 2i) / T(2n))
        x[i] = -x[n + 1 - i]
    end
    # For odd n, the central angle is exactly π/2 and its node is zero.
    isodd(n) && (x[div(n, 2) + 1] = zero(T))
    validate_rule(x, w)
    return x, w
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
function gausschebyshevu(::Type{T}, n::Integer) where {T}
    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    elseif n == 0
        x, w = T[], T[]
        validate_rule(x, w)
        return x, w
    elseif precision(T) < precision(Float32)
        x, w = gausschebyshevu(Float32, n)
        x, w = T.(x), T.(w)
        validate_rule(x, w)
        return x, w
    end
    # For i = 1,...,n, θ_i = iπ/(n + 1):
    # x_i = -cos(θ_i), w_i = π/(n + 1) * sin(θ_i)^2.
    # Symmetry gives x_{n+1-i} = -x_i and w_{n+1-i} = w_i.
    x = Vector{T}(undef, n)
    w = Vector{T}(undef, n)
    factor = T(π) / T(n + 1)
    # Use θ_i directly for θ_i ≤ π/4, and its complement nearer the centre.
    @inbounds for i in 1:div(n + 1, 4)
        si, ci = sincospi(T(i) / T(n + 1))
        x[i], x[n + 1 - i] = -ci, ci
        w[i] = w[n + 1 - i] = factor * si^2
    end
    # Near the centre, φ_i = π/2 - θ_i = (n + 1 - 2i)π/(2(n + 1)).
    # Thus x_i = -sin(φ_i), w_i = factor * cos(φ_i)^2. Forming φ_i
    # directly preserves small nodes that would lose accuracy in cos(θ_i).
    @inbounds for i in (div(n + 1, 4) + 1):div(n, 2)
        si, ci = sincospi(T(n + 1 - 2i) / (2 * T(n + 1)))
        x[i], x[n + 1 - i] = -si, si
        w[i] = w[n + 1 - i] = factor * ci^2
    end
    if isodd(n)
        # θ_i = π/2: x_i = 0 and sin(θ_i)^2 = 1.
        x[div(n, 2) + 1] = zero(T)
        w[div(n, 2) + 1] = factor
    end
    validate_rule(x, w)
    return x, w
end
gausschebyshevu(n::Integer) = gausschebyshevu(Float64, n)

@doc raw"""
    gausschebyshevv([T=Float64,] n::Integer) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Chebyshev quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) of the 3rd kind and type `T`.

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
function gausschebyshevv(::Type{T}, n::Integer) where {T}
    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    elseif n == 0
        x, w = T[], T[]
        validate_rule(x, w)
        return x, w
    elseif precision(T) < precision(Float32)
        x, w = gausschebyshevv(Float32, n)
        x, w = T.(x), T.(w)
        validate_rule(x, w)
        return x, w
    end
    # The phase below can be negative even when n is an unsigned integer.
    n = Int(n)
    # For i = 1,...,n, θ_i = 2iπ/(2n + 1):
    # x_i = -cos(θ_i), w_i = 2π/(2n + 1) * (1 + x_i).
    x = Vector{T}(undef, n)
    w = Vector{T}(undef, n)
    # Use half angles for θ_i ≤ π/4, where 1 + x_i is small.
    @inbounds for i in 1:n
        if i ≤ div(2n + 1, 8)
            # Near x = -1, let s² = sin(θ_i/2)^2. Then x_i = 2s² - 1
            # and w_i = 4π/(2n + 1) * s², avoiding cancellation in 1 + x_i.
            s² = sinpi(T(i) / T(2n + 1))^2
            x[i] = muladd(T(2), s², -one(T))
            w[i] = 2 * (2T(π) / T(2n + 1)) * s²
        else
            # x_i = sin(θ_i - π/2), where
            # θ_i - π/2 = (4i - 2n - 1)π/(2(2n + 1)). Form this phase
            # directly for accuracy near zero; here 1 + x_i is safe to use.
            x[i] = sinpi(T(4i - 2n - 1) / (2 * T(2n + 1)))
            w[i] = (2T(π) / T(2n + 1)) * (1 + x[i])
        end
    end
    validate_rule(x, w)
    return x, w
end
gausschebyshevv(n::Integer) = gausschebyshevv(Float64, n)

@doc raw"""
    gausschebyshevw([T=Float64,] n::Integer) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Chebyshev quadrature](https://en.wikipedia.org/wiki/Chebyshev%E2%80%93Gauss_quadrature) of the 4th kind and type `T`.

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
function gausschebyshevw(::Type{T}, n::Integer) where {T}
    if n > 0 && precision(T) < precision(Float32)
        # Reflect in the working precision before converting the output.
        x, w = gausschebyshevw(Float32, n)
        x, w = T.(x), T.(w)
        validate_rule(x, w)
        return x, w
    end
    # For i = 1,...,n, θ_i = (2i - 1)π/(2n + 1):
    # x_i = -cos(θ_i), w_i = 2π/(2n + 1) * (1 - x_i).
    # Reflect the third-kind rule and reverse it to keep ascending nodes:
    # x_i^(IV) = -x_{n+1-i}^(III), w_i^(IV) = w_{n+1-i}^(III).
    x, w = gausschebyshevv(T, n)
    reverse!(x)
    reverse!(w)
    x .= .-x
    validate_rule(x, w)
    return x, w
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
