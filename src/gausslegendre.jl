@doc raw"""
    gausslegendre(n::Integer) -> x, w

Return nodes `x` and weights `w` of [Gauss-Legendre quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature).

```math
\int_{-1}^{1} f(x) dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gausslegendre(3);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ 2/5
true
```
"""
@inline function gausslegendre(n::Integer)
    # GAUSSLEGENDRE(n) COMPUTE THE GAUSS-LEGENDRE NODES AND WEIGHTS IN O(n) time.

    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    elseif n == 0
        return Float64[], Float64[]
    elseif n == 1
        return [0.0], [2.0]
    elseif n == 2
        return [-1 / sqrt(3), 1 / sqrt(3)], [1.0, 1.0]
    elseif n == 3
        return [-sqrt(3 / 5), 0.0, sqrt(3 / 5)], [5 / 9, 8 / 9, 5 / 9]
    elseif n == 4
        a = 2 / 7 * sqrt(6 / 5)
        return ([-sqrt(3 / 7 + a), -sqrt(3/7-a), sqrt(3/7-a), sqrt(3/7+a)],
         [(18 - sqrt(30)) / 36, (18 + sqrt(30)) / 36,
          (18 + sqrt(30)) / 36, (18 - sqrt(30)) / 36])
    elseif n == 5
        b = 2 * sqrt(10 / 7)
        return ([-sqrt(5 + b) / 3, -sqrt(5 - b) / 3, 0.0,
          sqrt(5 - b) / 3, sqrt(5 + b) / 3],
         [(322 - 13 * sqrt(70)) / 900, (322 + 13 * sqrt(70)) / 900, 128 / 225,
          (322 + 13 * sqrt(70)) / 900, (322 - 13 * sqrt(70)) / 900])
    elseif n ≤ 60
        # NEWTON'S METHOD WITH THREE-TERM RECURRENCE:
        return rec(n)
    else
        # USE ASYMPTOTIC EXPANSIONS:
        return asy(n)
    end
end

function asy(n)
    # COMPUTE GAUSS-LEGENDRE NODES AND WEIGHTS USING ASYMPTOTIC EXPANSIONS.
    # COMPLEXITY O(n).

    # Nodes and weights:
    m = (n + 1) >> 1
    a = besselZeroRoots(m)
    rmul!(a, 1 / (n + 0.5))
    x = legpts_nodes(n, a)
    w = legpts_weights(n, m, a)
    # Use symmetry to get the others:
    resize!(x, n)
    resize!(w, n)
    @inbounds for i in 1:m
        x[n + 1 - i] = x[i]
        w[n + 1 - i] = w[i]
    end
    @inbounds for i in 1:m
        x[i] = -x[i]
    end
    @inbounds isodd(n) && (x[m] = 0.0)

    return x, w
end

function legpts_nodes(n, a)
    # ASYMPTOTIC EXPANSION FOR THE GAUSS-LEGENDRE NODES.
    vn = 1 / (n + 0.5)
    m = length(a)
    nodes = cot.(a)
    vn² = vn * vn
    vn⁴ = vn² * vn²
    @inbounds if n ≤ 255
        vn⁶ = vn⁴ * vn²
        for i in 1:m
            u = nodes[i]
            u² = u^2
            ai = a[i]
            ai² = ai * ai
            ai³ = ai² * ai
            ai⁵ = ai² * ai³
            node = ai + (u - 1 / ai) / 8 * vn²
            v1 = (6 * (1 + u²) / ai + 25 / ai³ - u * muladd(31, u², 33)) / 384
            v2 = u * @evalpoly(u², 2595 / 15360, 6350 / 15360, 3779 / 15360)
            v3 = (1 + u²) * (-muladd(31 / 1024, u², 11 / 1024) / ai +
                             u / 512 / ai² + -25 / 3072 / ai³)
            v4 = (v2 - 1073 / 5120 / ai⁵ + v3)
            node = muladd(v1, vn⁴, node)
            node = muladd(v4, vn⁶, node)
            nodes[i] = node
        end
    elseif n ≤ 3950
        for i in 1:m
            u = nodes[i]
            u² = u^2
            ai = a[i]
            ai² = ai * ai
            ai³ = ai² * ai
            node = ai + (u - 1 / ai) / 8 * vn²
            v1 = (6 * (1 + u²) / ai + 25 / ai³ - u * muladd(31, u², 33)) / 384
            node = muladd(v1, vn⁴, node)
            nodes[i] = node
        end
    else
        for i in 1:m
            u = nodes[i]
            ai = a[i]
            node = ai + (u - 1 / ai) / 8 * vn²
            nodes[i] = node
        end
    end
    @inbounds for jj = 1:m
        nodes[jj] = cos(nodes[jj])
    end

    return nodes
end

function legpts_weights(n, m, a)
    # ASYMPTOTIC EXPANSION FOR THE GAUSS-LEGENDRE WEIGHTS.
    vn = 1 / (n + 0.5)
    vn² = vn^2
    weights = Array{Float64}(undef, m)
    if n ≤ 850000
        @inbounds for i in 1:m
            weights[i] = cot(a[i])
        end
    end
    # Split out the part that can be vectorized by llvm
    @inbounds if n ≤ 170
        for i in 1:m
            u = weights[i]
            u² = u^2
            ai = a[i]
            ai⁻¹ = 1 / ai
            ai² = ai^2
            ai⁻² = 1 / ai²
            ua = u * ai
            W1 = muladd(ua - 1, ai⁻², 1.0) / 8
            W2 = @evalpoly(ai⁻², @evalpoly(u², -27.0, -84.0, -56.0),
                           muladd(-3.0, muladd(u², -2.0, 1.0), 6 * ua),
                           muladd(ua, -31.0, 81.0)) / 384
            W3 = @evalpoly(ai⁻¹, @evalpoly(u², 153 / 1024, 295 / 256, 187 / 96,
                                           151 / 160),
                           @evalpoly(u², -65 / 1024, -119 / 768, -35 / 384) * u,
                           @evalpoly(u², 5 / 512, 15 / 512, 7 / 384),
                           muladd(u², 1 / 512, -13 / 1536) * u,
                           muladd(u², -7 / 384, + 53 / 3072),
                           3749 / 15360 * u, -1125 / 1024)
            weights[i] = @evalpoly(vn², 1 / vn² + W1, W2, W3)
        end
    elseif n ≤ 1500
        for i in 1:m
            u = weights[i]
            u² = u^2
            ai = a[i]
            ai² = ai^2
            ai⁻² = 1 / ai²
            ua = u * ai
            W1 = muladd(ua - 1, ai⁻², 1.0) / 8
            W2 = @evalpoly(ai⁻², @evalpoly(u², -27.0, -84.0, -56.0),
                           muladd(-3.0, muladd(u², -2.0, 1.0), 6 * ua),
                           muladd(ua, -31.0, 81.0)) / 384
            weights[i] = muladd(vn², W2, 1 / vn² + W1)
        end
    elseif n ≤ 850000
        for i in 1:m
            u = weights[i]
            u² = u^2
            ai = a[i]
            ai² = ai^2
            ai⁻² = 1 / ai²
            ua = u * ai
            W1 = muladd(ua - 1, ai⁻², 1.0) / 8
            weights[i] = 1 / vn² + W1
        end
    else
        for i in 1:m
            weights[i] = 1 / vn²
        end
    end
    bJ1 = besselJ1(m)
    @inbounds for i in 1:m
        weights[i] = 2 / (bJ1[i] * (a[i] / sin(a[i])) * weights[i])
    end

    return weights
end

function rec(n)
    # COMPUTE GAUSS-LEGENDRE NODES AND WEIGHTS USING NEWTON'S METHOD.
    # THREE-TERM RECURENCE IS USED FOR EVALUATION. COMPLEXITY O(n^2).
    # Initial guesses:
    x0 = asy(n)[1]
    x = x0[1:n ÷ 2 + 1]
    # Perform Newton to find zeros of Legendre polynomial:
    PP1, PP2 = innerRec(n, x)
    @inbounds @simd for i in 1:length(x)
        x[i] -= PP1[i] / PP2[i]
    end
        # One more Newton for derivatives:
    PP1, PP2 = innerRec(n, x)
    @inbounds @simd for i in 1:length(x)
        x[i] -= PP1[i] / PP2[i]
    end

    # Use symmetry to get the other Legendre nodes and weights:
    m = length(x)
    resize!(x, n)
    resize!(PP2, n)
    @inbounds for i in 1:m-1
        x[n + 1 - i] = -x[i]
        PP2[n + 1 - i] = -PP2[i]
    end
    @inbounds for i in 1:n
        PP2[i] = 2 / ((1 - x[i]^2) * PP2[i]^2)
    end

    return x, PP2
end

function innerRec(n, x)
    # EVALUATE LEGENDRE AND ITS DERIVATIVE USING THREE-TERM RECURRENCE RELATION.
    N = size(x, 1)
    myPm1 = Array{Float64}(undef, N)
    myPPm1 = Array{Float64}(undef, N)
    @inbounds for j = 1:N
        xj = x[j]
        Pm2 = 1.0
        Pm1 = xj
        PPm1 = 1.0
        PPm2 = 0.0
        for k = 1:(n - 1)
            Pm2, Pm1 = Pm1, muladd((2 * k + 1) * Pm1, xj, - k * Pm2) / (k + 1)
            PPm2, PPm1 = PPm1, ((2 * k + 1) * muladd(xj, PPm1, Pm2) -
                                k * PPm2) / (k + 1)
        end
        myPm1[j] = Pm1
        myPPm1[j] = PPm1
    end
    return myPm1, myPPm1
end
