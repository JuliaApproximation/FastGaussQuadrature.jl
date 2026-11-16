@doc raw"""
    gausshermite([T=Float64,] n::Integer; normalize = false) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Hermite quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Hermite_quadrature) with type `T`.

```math
\int_{-\infty}^{+\infty} f(x) \exp(-x^2) dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

The option `normalize=true` instead computes
```math
\int_{-\infty}^{+\infty} f(x) \phi(x) dx \approx \sum_{i=1}^{n} w_i f(x_i),
```
where $\phi$ is the standard normal density function.

# Examples
```jldoctest
julia> x, w = gausshermite(3);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ 3(√π)/4
true
```
"""
function gausshermite(::Type{T}, n::Integer; normalize = false) where {T}
    if precision(T) < precision(Float32)
        x, w = gausshermite(Float32, n; normalize)
        x, w = T.(x), T.(w)
    else
        x, w = _hermite_rule(T, n, true)
        w .*= exp.(-x .^ 2)
        if normalize
            x *= sqrt(T(2))
            w /= sqrt(T(π))
        end
    end
    validate_rule(x, w)
    return x, w
end
gausshermite(n::Integer; normalize = false) = gausshermite(Float64, n; normalize = normalize)

unweightedgausshermite(::Type{T}, n::Integer) where {T <: Real} = _hermite_rule(T, n, false)

function _hermite_rule(::Type{T}, n::Integer, physical::Bool) where {T <: Real}
    # In scaled mode, return s_i = exp(x_i²) w_i for the unweighted rule.
    # Physical mode prepares weights for multiplication by exp(-x_i²) in
    # gausshermite, including the GST residual correction before normalization.
    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    elseif n == 0
        nodes, weights = T[], T[]
        validate_rule(nodes, weights)
        return nodes, weights
    elseif precision(T) < precision(Float32)
        x, w = _hermite_rule(Float32, n, physical)
        nodes, weights = T.(x), T.(w)
        validate_rule(nodes, weights)
        return nodes, weights
    elseif n == 1
        nodes, weights = T[0], T[sqrt(T(π))]
        validate_rule(nodes, weights)
        return nodes, weights
    elseif n ≤ 20
        x = hermite_rec(T, n)
    else
        x = hermite_gst(T, n, physical)
    end

    # Fold out directly into the final arrays, preserving exact symmetry.
    _x, _w = Vector{T}(undef, n), Vector{T}(undef, n)
    m = length(x[1])
    @inbounds for i in 1:m
        _x[i] = -x[1][m + 1 - i]
        _x[n + 1 - i] = x[1][m + 1 - i]
        _w[i] = _w[n + 1 - i] = x[2][m + 1 - i]
    end
    isodd(n) && (_x[m] = zero(T))
    # Normalize the scaled weights so Σ exp(-x_i²) * scaled_w_i = √π.
    _w .*= sqrt(T(π)) / sum(i -> exp(-_x[i]^2) * _w[i], eachindex(_x))

    validate_rule(_x, _w)
    return _x, _w
end
unweightedgausshermite(n::Integer) = unweightedgausshermite(Float64, n)

# At small order, a classical recurrence is cheaper than high-precision
# Taylor propagation. GST supplies seeds, which are then refined at precision T.
function hermite_rec(::Type{T}, n::Integer) where {T <: AbstractFloat}
    # Seeds never use wider arithmetic than T; Newton recovers higher precision.
    W = precision(T) < precision(Float64) ? Float32 : Float64
    nodes = T.(hermite_gst(W, n)[1])
    weights = similar(nodes)
    function polynomial_derivative(x)
        previous, p = one(T), 2x
        for k in 1:(n - 1)
            previous, p = p, 2x * p - 2k * previous
        end
        return p, 2n * previous
    end
    @inbounds for i in eachindex(nodes)
        x = nodes[i]
        for iteration in 1:max(10, precision(T))
            p, derivative = polynomial_derivative(x)
            step = p / derivative
            x -= step
            abs(step) ≤ eps(T) * max(one(T), abs(x)) && break
            iteration == max(10, precision(T)) && error("Hermite recurrence did not converge")
        end
        _, derivative = polynomial_derivative(x)
        nodes[i] = x
        # The common normalization constant is supplied by the zeroth-moment
        # calculation in the caller.
        weights[i] = exp(x^2) / derivative^2
    end
    validate_rule(nodes, weights)
    return nodes, weights
end

function hermpoly_rec(n::Integer, x0)
    # HERMPOLY_rec evaluation of scaled Hermite poly using recurrence
    n < 0 && throw(DomainError(n, "Input n must be a non-negative integer"))
    # evaluate:
    w = exp(-x0^2 / (4 * n))
    wc = 0 # 0 times we've applied wc
    Hold = one(x0)
    # n == 0 && return (Hold, 0)
    H = x0
    for k in 1:(n - 1)
        Hold, H = H, (x0 * H / sqrt(oftype(x0, k + 1)) - Hold / sqrt(oftype(x0, k + 1) / k))
        while abs(H) ≥ 100 && wc < n  # regularise
            H *= w
            Hold *= w
            wc += 1
        end
        k += 1
    end
    for _ in (wc + 1):n
        H *= w
        Hold *= w
    end

    return H, -x0 * H + sqrt(oftype(x0, n)) * Hold
end

function hermpoly_rec(r::Base.OneTo, x0)
    isempty(r) && return [1.0]
    n = maximum(r)
    # HERMPOLY_rec evaluation of scaled Hermite poly using recurrence
    n < 0 && throw(DomainError(n, "Input n must be a non-negative integer"))
    n == 0 && return [exp(-x0^2 / 4)]
    p = max(1, floor(Int, x0^2 / 100))
    w = exp(-x0^2 / 4p)
    wc = 0 # 0 times we've applied wc
    ret = Vector{Float64}()
    Hold = one(x0)
    push!(ret, Hold)
    H = x0
    push!(ret, H)
    for k in 1:(n - 1)
        Hold, H = H, (x0 * H / sqrt(k + 1) - Hold / sqrt(1 + 1 / k))
        while abs(H) ≥ 100 && wc < p  # regularise
            ret .*= w
            H *= w
            Hold *= w
            wc += 1
        end
        push!(ret, H)
        k += 1
    end
    ret .*= w^(p - wc)

    return ret
end

hermpoly_rec(r::AbstractRange, x0) = hermpoly_rec(Base.OneTo(maximum(r)), x0)[r .+ 1]

# Gil--Segura--Temme, Numer. Math. 143 (2019), §3, and
# SIAM J. Sci. Comput. 48 (2026), §4: march outwards from exact central data
# for y = exp(-x^2/2)H_n(x), y'' + (2n+1-x^2)y = 0. The arbitrary initial
# normalization is removed by unweightedgausshermite's zeroth-moment scaling.
function hermite_gst(::Type{T}, n::Integer, physical::Bool = false) where {T <: AbstractFloat}
    λ = T(2) * n + 1
    m = (n + 1) ÷ 2
    nodes, weights = Vector{T}(undef, m), Vector{T}(undef, m)
    x = zero(T)
    y = iseven(n) ? one(T) : zero(T)
    derivative = iseven(n) ? zero(T) : one(T)
    if isodd(n)
        nodes[1], weights[1] = zero(T), one(T)
    end
    # The leading relative error after a fourth-order step is -step^4/6.
    @inbounds for i in (isodd(n) ? 2 : 1):m
        frequency = sqrt(muladd(-x, x, λ))
        step = iszero(derivative) ? T(π) / (2frequency) :
            (T(π) - atan(frequency * y / derivative)) / frequency
        next = x + step
        # Use the actual difference of the rounded iterates. Carrying the
        # residual y between roots prevents accumulated phase error (§3.1.1).
        y, derivative = hermite_gst_taylor(λ, x, next - x, y, derivative)
        x = next
        for iteration in 1:max(10, precision(T))
            frequency = sqrt(muladd(-x, x, λ))
            step = -atan(frequency * y / derivative) / frequency
            next = x + step
            y, derivative = hermite_gst_taylor(λ, x, next - x, y, derivative)
            x = next
            abs(step) ≤ sqrt(sqrt(eps(T))) / 2 && break
            iteration == max(10, precision(T)) && error("Hermite iteration did not converge")
        end
        nodes[i] = x
        # y'' vanishes at every root, so this derivative-based scaled weight
        # is stationary there and avoids exponential polynomial growth.
        weight = inv(derivative)^2
        if physical
            # Transport exp(-x^2) from the rounded iterate to the root.
            # y/y' is the remaining root displacement; no coordinate pair is kept.
            weight *= muladd(2x, y / derivative, one(T))
        end
        weights[i] = weight
    end
    validate_rule(nodes, weights)
    return nodes, weights
end

function hermite_gst_taylor(λ, x, h, y::T, derivative, compensated = true) where {T <: AbstractFloat}
    iszero(h) && return y, derivative
    # t_k = y^(k)(x) h^k/k!. Differentiating the ODE gives
    # (k+1)(k+2)t_{k+2} = (x^2-λ)h^2 t_k + 2xh^3 t_{k-1} + h^4 t_{k-2}.
    h2 = h * h
    # Cache the recurrence coefficients shared by every Taylor term.
    a, b, c = muladd(x, x, -λ) * h2, 2x * h2 * h, h2 * h2
    t0, t1 = y, derivative * h
    tm1, tm2 = zero(T), zero(T)
    value, scaled_derivative = y + t1, t1
    value_error, derivative_error = zero(T), zero(T)
    small_terms = 0
    @inbounds for k in 0:max(100, 2precision(T))
        t2 = muladd(a, t0, muladd(b, tm1, c * tm2)) / ((k + 1) * (k + 2))
        if compensated
            # Compensated sums reduce cancellation on a whole root-to-root step.
            increment = t2 - value_error
            updated = value + increment
            value_error = (updated - value) - increment
            value = updated
            increment = (k + 2) * t2 - derivative_error
            updated = scaled_derivative + increment
            derivative_error = (updated - scaled_derivative) - increment
            scaled_derivative = updated
        else
            value += t2
            scaled_derivative += (k + 2) * t2
        end
        small_terms = abs(t2) * (k + 2) < eps(T) * max(abs(value), abs(scaled_derivative)) / 8 ? small_terms + 1 : 0
        # Four successive terms cover all lags, including parity zeros at x=0.
        small_terms ≥ 4 && return value, scaled_derivative / h
        tm2, tm1, t0, t1 = tm1, t0, t1, t2
    end
    error("Hermite Taylor series did not converge")
end
