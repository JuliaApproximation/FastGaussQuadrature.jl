@doc raw"""
    gaussradau([T=Float64,] n::Integer) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Radau quadrature](https://mathworld.wolfram.com/RadauQuadrature.html) and type `T`.
The order `n` must be positive. Use `gaussradau(n, α, β)` for a Jacobi weight.

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
function gaussradau(::Type{T}, n::Integer) where {T}
    n > 0 || throw(DomainError(n, "Input N must be a positive integer"))
    if n == 1
        x, w = T[-1], T[2]
    elseif precision(T) < precision(Float32)
        x, w = gaussradau(Float32, n)
        x, w = T.(x), T.(w)
    else
        a = b = zero(T)
        # RADAUPTS   Gauss-Legendre-Radau Quadrature Nodes and Weights
        if n == 2
            x, w = [-1, one(T) / 3], [one(T) / 2, convert(T, 3) / 2]
        elseif (T <: AbstractFloat && !(T <: Union{Float32, Float64})) || (T !== Float64 && n > 101)
            x, w = gaussradau_jacobi(n, a, b)
        else
            # Compute via GaussJacobi:
            if T == Float64 && n > 101
                x, w = jacobi_asy(n - 1, a, b + 1, Val(:radau))
            else
                x, w = gaussjacobi(n - 1, a, b + 1)
                @inbounds for i in 1:length(w)
                    w[i] = w[i] / (1 + x[i])
                end
            end
            pushfirst!(x, -1)
            pushfirst!(w, convert(T, 2) / n^2)
        end
    end
    validate_rule(x, w)
    return x, w
end
gaussradau(n::Integer) = gaussradau(Float64, n)
@deprecate gaussradau(n::Integer, T::Type) gaussradau(T, n)

@doc raw"""
    gaussradau(n::Integer, α, β) -> x, w  # nodes, weights

Return Gauss–Radau nodes and weights for the Jacobi weight, with the first
node fixed at ``-1``:

```math
\int_{-1}^{1} f(x) (1-x)^\alpha(1+x)^\beta dx \approx \sum_{i=1}^{n} w_i f(x_i).
```

The parameters must satisfy `n ≥ 1` and `α, β > -1`. Both output vectors use
the common type of `promote(float(α), float(β))`; types with fewer than 24
significand bits compute in `Float32` before conversion. The rule integrates
polynomials of degree up to `2n-2` exactly in exact arithmetic.

# Examples
```jldoctest
julia> x, w = gaussradau(4, 0.0f0, 1.0f0);

julia> eltype(x) === eltype(w) === Float32 && x[1] == -1
true

julia> dot(w, x.^2) ≈ 2/3
true
```
"""
@inline function gaussradau(n::Integer, α, β)
    if n ≤ 0
        throw(DomainError(n, "Input N must be a positive integer"))
    end
    α, β = promote(float(α), float(β))
    α > -1 && β > -1 || throw(DomainError((α, β), "The Jacobi parameters correspond to a nonintegrable weight function"))
    T = typeof(α)
    if precision(T) < precision(Float32)
        x, w = gaussradau(n, Float32(α), Float32(β))
        x, w = T.(x), T.(w)
    elseif n == 1
        x, w = [-one(α)], [iszero(α) && iszero(β) ? 2one(α) : jacobimoment(α, β)]
    else
        x, w = gaussradau_jacobi(n, α, β)
    end
    validate_rule(x, w)
    return x, w
end

function gaussradau_jacobi(n::Integer, α, β)
    T = typeof(α)
    if T === Float64 && iszero(α) && iszero(β) && (n ≤ 101 || n > 500)
        return gaussradau(n)
    end
    m = n - 1
    if n == 2
        # Positive shifted parameters retain the small mass on either side
        # as a parameter approaches -1. Do not subtract a weight from the mass.
        a, b = α + 1, β + 1
        s = a + b
        mass = jacobi_mass(α, β)
        isfinite(mass) && mass > 0 || error("The quadrature normalization is not representable in $T.")
        x = T[-1, (b - α) / (s + 1)]
        w = [mass * (a / s) / (b + 1), mass * (b / (b + 1)) * (1 + 1 / s)]
        validate_rule(x, w)
        return x, w
    elseif T === Float64 && m ≥ 8 && (max(α, β + 1) < 5 || jacobi_gst_regime(m, α, β + 1))
        # Gil--Segura--Temme (2026), §5.1: the interior is a Jacobi rule
        # for (1-x)^α(1+x)^(β+1). Transform weights before rounding endpoint
        # distances when the underlying kernel retains them. At intermediate
        # degrees GST avoids the setup cost of the asymptotic coefficients;
        # at larger fixed-parameter degrees the asymptotic weights are preferred.
        # Near α=-1 the Bessel boundary coefficients lose accuracy. Match
        # the lower-parameter limit of the public Jacobi asymptotic route.
        if α < -0.99 || jacobi_gst_regime(m, α, β + 1) || (m ≥ 20 && n ≤ 500)
            x, w = jacobi_gst(m, α, β + 1, Val(:radau))
        elseif m > 100
            x, w = jacobi_asy(m, α, β + 1, Val(:radau))
        else
            x, w = gaussjacobi(m, α, β + 1)
            @inbounds for i in eachindex(w)
                w[i] /= 1 + x[i]
            end
        end
    else
        x, w = jacobi_gst(m, α, β + 1, Val(:radau))
    end
    pushfirst!(x, -one(T))
    pushfirst!(w, radau_endpoint_weight(n, α, β))
    validate_rule(x, w)
    return x, w
end

function radau_endpoint_weight(n, α::T, β::T) where {T <: AbstractFloat}
    iszero(α) && iszero(β) && return T(2) / n^2
    weight = jacobi_mass(α, β)
    isfinite(weight) && weight > 0 || error("The quadrature normalization is not representable in $T.")
    # Apply each decreasing factor to the mass separately. Combining the two
    # factors first can underflow even when the endpoint weight is representable.
    # The k-th update is k/(k+β+1) * (k+α)/(k+α+β+1).
    @inbounds for k in 1:(n - 1)
        f1 = T(k) / (T(k) + (β + 1))
        f2 = (T(k - 1) + (α + 1)) / (T(k - 1) + (α + 1) + (β + 1))
        if (iszero(f1) || iszero(f2)) && !iszero(weight)
            error("The Radau endpoint weight cannot be resolved in $T.")
        end
        weight *= f1
        weight *= f2
    end
    return weight
end

function radau_endpoint_weight(n, α::Float64, β::Float64)
    iszero(α) && iszero(β) && return 2.0 / n^2
    if α == β && n ≤ α + 1
        # Avoid cancelling logarithms of order α when the symmetric mass
        # has a simple beta-function expression. This product decreases.
        weight = beta(0.5, α + 1)
        @inbounds for k in 1:(n - 1)
            weight *= (k / (k + β + 1)) * ((k + α) / (k + α + β + 1))
        end
        return weight
    end
    # Equivalent to the gamma formula, without separately overflowing its
    # numerator or underflowing beta factors before their powers of two.
    return exp(
        (α + β + 1) * log(2.0) + log(β + 1) +
            logbeta(Float64(n), β + 1) + logbeta(n + α, β + 1)
    )
end
