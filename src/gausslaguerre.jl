@doc raw"""
    gausslaguerre(n::Integer) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of [Gauss-Laguerre quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature).

```math
\int_{0}^{+\infty} f(x) e^{-x} dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gausslaguerre(3);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ 24
true
```
"""
function gausslaguerre(n::Integer)
    return gausslaguerre(n, 0.0)
end


@doc raw"""
    gausslaguerre(n::Integer, α::Real) -> x, w  # nodes, weights

Return nodes `x` and weights `w` of generalized [Gauss-Laguerre quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature).

```math
\int_{0}^{+\infty} f(x) x^\alpha e^{-x} dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gausslaguerre(3, 1.0);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ 120
true
```

Optionally, a reduced quadrature rule can be computed. In that case, only those
points and weights are computed for which the weight does not underflow in the
floating point precision type. Supply the optional argument `reduced = true`.

Though the code is generic, heuristical choices on the choice of the algorithm
are based on achieving machine precision accuracy only for `Float64` type. In
case the default choice of algorithm does not achieve the desired accuracy, the
user can manually invoke the following routines:
- `gausslaguerre_GW`: computation based on Golub-Welsch
- `gausslaguerre_rec`: computation based on Newton iterations applied to evaluation
   using the recurrence relation
- `gausslaguerre_asy`: the asymptotic expansions
"""
function gausslaguerre(n::Integer, α::Real; reduced = false)
    if α ≤ -1
        throw(DomainError(α, "The Laguerre parameter α ≤ -1 corresponds to a nonintegrable weight function"))
    end
    if n < 0
        throw(DomainError(n, "gausslaguerre($n,$α) not defined: n must be positive."))
    end

    # Guess the numerical type from the supplied type of α
    # Although the code is generic, the heuristics are derived for Float64 precision
    T = typeof(float(α))
    if n == 0
        T[],T[]
    elseif n == 1
        [1+α], [gamma(1+α)]
    elseif n == 2
        [α + 2-sqrt(α+2),α+2+sqrt(α+2)],
        [((α-sqrt(α+2)+2)*gamma(α+2))/(2*(α+2)*(sqrt(α+2)-1)^2),
         ((α+sqrt(α+2)+2)*gamma(α+2))/(2*(α+2)*(sqrt(α+2)+1)^2)]
    elseif n < 15
        # Use Golub-Welsch for small n
        gausslaguerre_GW(n, α)
    elseif n < 128
        # Use the recurrence relation for moderate n
        gausslaguerre_rec(n, α)
    else
        # Use explicit asymptotic expansions for larger n
        # The restriction to α comes from the restriction on ν in besselroots
        if α < 5
            gausslaguerre_asy(n, α, reduced=reduced, T=-1, recompute=true)
        else
            gausslaguerre_rec(n, α)
        end
    end
end

# Our threshold for deciding on underflow
underflow_threshold(x) = underflow_threshold(typeof(x))
underflow_threshold(::Type{T}) where {T <: AbstractFloat} = 10floatmin(T)


"""
Compute the Gauss-Laguerre rule using explicit asymptotic expansions for the nodes and weights.
Optional parameters are:
- `reduced`: compute a reduced quadrature rule, discarding all points and weights as soon as the weights underflow
- `T`: the order of the expansion. Set `T=-1` to determine the order adaptively depending on the size of the terms in the expansion
- `recompute`: if a crude measure of the error is larger than a tolerance, the point and weight are recomputed using the (slower) recursion+newton approach, yielding more reliable accurate results.
"""
function gausslaguerre_asy(n::Integer, α;
    reduced = false,
    T = max(1, ceil(Int, 50/log(n))),  # Heuristic for number of terms
    recompute = false)

    if α^2/n > 1
        @warn "A large value of α may lead to inaccurate results."
    end

    ELT = typeof(float(α))

    n_alloc = reduced ? 0 : n
    x = zeros(ELT, n_alloc)
    w = zeros(ELT, n_alloc)

    # The expansions are given in powers of 1/(4n+2α+2)
    d = one(ELT)/(4n+2α+2)

    # Heuristical indices for Bessel and Airy regions
    k_bessel = max(ceil(Int, sqrt(n) ), 7)
    k_airy = floor(Int, 0.9*n)

    # The Bessel region
    # First compute the roots of the Bessel function of order α
    jak_vector = approx_besselroots(α, k_bessel)

    bessel_wins = true
    k = 0
    while bessel_wins && k < n
        k += 1
        # We iterate until the estimated error of the bulk expansion is smaller
        # than the one of the Bessel expansion
        jak = (k < k_bessel) ? jak_vector[k] : McMahon(α, k)

        xk, wk, δ_bessel = gausslaguerre_asy_bessel(n, α, jak, d, T)
        xkb, wkb, δ_bulk = gausslaguerre_asy_bulk(n, α, k, d, T)
        if δ_bulk < δ_bessel
            bessel_wins = false
            xk = xkb
            wk = wkb
        end
        if recompute
            δ = min(δ_bessel,δ_bulk)
            if δ > 1e-13
                xk_rec, wk_rec = gl_rec_newton(xk, n, α)
                if abs(xk_rec-xk) < 100δ
                    xk = xk_rec
                    wk = wk_rec
                end
            end
        end
        if reduced
            if abs(wk) < underflow_threshold(ELT)
                return x, w
            else
                push!(x, xk); push!(w, wk)
            end
        else
            x[k] = xk; w[k] = wk
        end
    end

    # The bulk region
    # - First we go from where we left of to our heuristic
    while k < k_airy-1
        k += 1
        xk, wk, δ_bulk = gausslaguerre_asy_bulk(n, α, k, d, T)
        if recompute
            if δ_bulk > 1e-13
                xk_rec, wk_rec = gl_rec_newton(xk, n, α)
                if abs(xk_rec-xk) < 100δ_bulk
                    xk = xk_rec
                    wk = wk_rec
                end
            end
        end
        if reduced
            if abs(wk) < underflow_threshold(ELT)
                return x, w
            else
                push!(x, xk); push!(w, wk)
            end
        else
            x[k] = xk; w[k] = wk
        end
    end

    # - Then we compare to Airy until it wins, and then we switch to just Airy
    bulk_wins = true
    while bulk_wins && k < n
        k += 1
        xk, wk, δ_bulk = gausslaguerre_asy_bulk(n, α, k, d, T)
        xka, wka, δ_airy = gausslaguerre_asy_airy(n, α, k, d, T)
        if δ_airy < δ_bulk
            bulk_wins = false
            xk = xka
            wk = wka
        end
        if recompute
            δ = min(δ_airy,δ_bulk)
            if δ > 1e-13
                xk_rec, wk_rec = gl_rec_newton(xk, n, α)
                if abs(xk_rec-xk) < 100δ
                    xk = xk_rec
                    wk = wk_rec
                end
            end
        end
        if reduced
            if abs(wk) < underflow_threshold(ELT)
                return x, w
            else
                push!(x, xk); push!(w, wk)
            end
        else
            x[k] = xk; w[k] = wk
        end
    end

    # The Airy region
    while k < n
        k += 1
        xk, wk, δ_airy = gausslaguerre_asy_airy(n, α, k, d, T)
        if recompute
            if δ_airy > 1e-13
                xk_rec, wk_rec = gl_rec_newton(xk, n, α)
                if abs(xk_rec-xk) < 100δ_airy
                    xk = xk_rec
                    wk = wk_rec
                end
            end
        end
        if reduced
            if abs(wk) < underflow_threshold(ELT)
                return x, w
            else
                push!(x, xk); push!(w, wk)
            end
        else
            x[k] = xk; w[k] = wk
        end
    end

    # Sanity check
    if ( minimum(x) < 0.0 ) || ( maximum(x) > 4*n + 2*α + 2 ) ||  ( minimum(diff(x)) ≤ 0.0 ) || (minimum(w) < 0.0)
        @warn "Unexpected inconsistency in the computation of nodes and weights"
    end

    return x, w
end

## Expansion coefficients
# These are explicit formulas of the coefficients, up to a simple postprocessing
# that is common to all factors and not included here (see below).
#
# General expressions are given in terms of α, more specific expressions
# follow for the special case α = 0.

## The bulk

# Note: there is always one division by an integer, placed such that it preserves the type of `d`
function gl_bulk(t, d, α)
    α² = α * α
    tinv = inv(1-t)
    _t = (1-t) / t
    _t² = _t * _t
    _t³ = _t² * _t
    d² = d * d
    d⁴ = d² * d²
    d⁶ = d⁴ * d²

    c1 = evalpoly(α², (-4, 12))
    x3 = -evalpoly(tinv, (c1, -4, 5)) * d / 12
    c1, c2 = evalpoly(α², (224, -960, 480)), evalpoly(α², (7, -30, 15))
    x5 = d² * d * _t / 720 * evalpoly(tinv, (c1, -48*c2, -16, -576, 2814, -3815, 1600))
    c1, c2, c3 = evalpoly(α², (-285696, 1354752, -967680, 193536)), evalpoly(α², (-31, 147, -105, 21)), evalpoly(α², (-1346, 6405, -4620, 945))
    c4, c5 = evalpoly(α², (43, -126, 63)), evalpoly(α², (-221, -630, 315))
    x7 = -d⁴ * d / 181440 * _t² * evalpoly(tinv, (c1, -23040*c2, 384*c3, -320*c4, 80*c5, -1727136, 16131880, -48469876, 66424575, -43122800, 10797500))
    c1, c2, c3 = evalpoly(α², (210677760, -1028505600, 812851200, -232243200, 24883200)), evalpoly(α², (127, -620, 490, -140, 15)), evalpoly(α², (1193053, -5826660, 4613070, -1324260, 143325))
    c4, c5, c6 = evalpoly(α², (555239, -2716980, 2163630, -631260, 70875)), evalpoly(α², (-641, 2960, -2155, 450)), evalpoly(α², (-1598, 17685, -13905, 3375))
    c7, c8, c9 = evalpoly(α², (-7823, -9042, 4521)), evalpoly(α², (15948182, -206850, 103425)), evalpoly(α², (64957561, -24000, 12000))
    x9 = d⁶ * d / 10886400 * _t³ * evalpoly(tinv, (c1, -5806080*c2, 768*c3, -768*c4, 16128*c5, -1792*c6, 3360*c7, -192*c8, 672*c9, -212307298152, 518401904799, -714465642135, 566519158800, -241928673000, 43222750000))

    w3 = d² / 6 * (2*t + 3) / (t-1)^3
    c1 = evalpoly(α², (7, -30, 15))
    w5 = d⁴ * _t² / 720  * evalpoly(tinv, (0, 0, 16*c1, 32, 1712, -12408, 27517, -24860, 8000))
    c1, c2, c3 = evalpoly(α², (-31, 147, -105, 21)), evalpoly(α², (-416, 1995, -1470, 315)), evalpoly(α², (43, -126, 63))
    c4, c5 = evalpoly(α², (-89, -378, 189)), evalpoly(α², (53752, -630, 315))
    w7 = -d⁶ * _t³ / 90720 * evalpoly(tinv, (0, 0, 2304*c1, -384*c2, 480*c3, -320*c4, 80*c5, -50986344, 201908326, -386872990, 393326325, -204917300, 43190000))
    return (x3, x5, x7, x9), (w3, w5, w7)
end

function gl_bulk(t, d)
    tinv = inv(1-t)
    _t = (1-t) / t
    _t² = _t * _t
    _t³ = _t² * _t
    d² = d * d
    d⁴ = d² * d²
    d⁶ = d⁴ * d²

    x3 = -d / 12 * evalpoly(tinv, (-4, -4, 5))
    x5 = d² * d * _t / 720 * evalpoly(tinv, (224, -336, -16, -576, 2814, -3815, 1600))
    x7 = -d⁴ * d / 181440 * _t² * evalpoly(tinv, (-285696, 714240, -516864, -13760, -17680, -1727136, 16131880, -48469876, 66424575, -43122800, 10797500))
    x9 = d⁶ * d / 10886400 * _t³ * evalpoly(tinv, (210677760, -737372160, 916264704, -426423552, -10338048, 2863616, -26285280, -3062050944, 43651480992, -212307298152, 518401904799, -714465642135, 566519158800, -241928673000, 43222750000))

    w3 = d² / 6 * (2*t + 3) / (t-1)^3
    w5 = d⁴ / 720 * _t² * evalpoly(tinv, (0, 0, 112, 32, 1712, -12408, 27517, -24860, 8000))
    w7 = -d⁶ * _t³ / 90720 * evalpoly(tinv, (0, 0, -71424, 159744, 20640, 28480, 4300160, -50986344, 201908326, -386872990, 393326325, -204917300, 43190000))
    return (x3, x5, x7, x9), (w3, w5, w7)
end

## The hard edge (Bessel region)
function gl_bessel(jak, d, α)
    jak² = jak * jak
    d² = d * d
    d⁴ = d² * d²
    d⁸ = d⁴ * d⁴
    α² = α * α

    c1 = evalpoly(α², (-2, 2))
    x3 = d² * evalpoly(jak², (c1, 1)) / 3
    w3 = d² * evalpoly(jak², (c1, 2)) / 3
    c1, c2 = evalpoly(α², (94, -140, 46)), evalpoly(α², (-19, 11)) 
    x5 = d⁴ * evalpoly(jak², (c1, 3*c2, 11)) / 45
    w5 = d⁴ * evalpoly(jak², (c1, 6*c2, 33)) / 45
    c1, c2, c3 = evalpoly(α², (-12077, 19887, -9303, 1493)), evalpoly(α², (14051, -10750, 2459)), evalpoly(α², (-181, 73))
    x7 = d⁴ * d² * evalpoly(jak², (4*c1, 2*c2, 36*c3, 657)) / 2835
    w7 = 4 * d⁴ * d² * evalpoly(jak², (c1, c2, 27*c3, 657)) / 2835
    c1, c2, = evalpoly(α², (6029959, -10087180, 5095482, -1146220, 107959)), evalpoly(α², (-2201939, 1678761, -507801, 63299))
    c3, c4 = evalpoly(α², (1456807, -729422, 125671)), evalpoly(α², (-2879, 887))
    x9 = d⁸ * evalpoly(jak², (2*c1, 3*c2, c3, 60*c4, 10644)) / 42525
    w9 = d⁸ * evalpoly(jak², (2*c1, 6*c2, 3*c3, 240*c4, 53220)) / 42525
    return (x3, x5, x7, x9), (w3, w5, w7, w9)
end

# And for α = 0:
function gl_bessel(jak, d)
    jak² = jak * jak
    d² = d * d
    d⁴ = d² * d²
    d⁸ = d⁴ * d⁴

    x3 = d² * evalpoly(jak², (-2, 1)) / 3
    x5 = d⁴ * evalpoly(jak², (94, -57, 11)) / 45
    x7 = d⁴ * d² * evalpoly(jak², (-48308, 28102, -6516, 657)) / 2835
    x9 = d⁸ * evalpoly(jak², (12059918, -6605817, 1456807, -172740, 10644)) / 42525
    x11 = d⁸ * d² * evalpoly(jak², (-11427291076, 6028914206, -1248722004, 138902061, -9908262, 410649)) / 1403325

    w3 = 2 * d² * evalpoly(jak², (-1, 1)) / 3
    w5 = d⁴ * evalpoly(jak², (94, -114, 33)) / 45
    w7 = 4 * d⁴ * d² * evalpoly(jak², (-12077, 14051, -4887, 657)) / 2835
    w9 = d⁸ * evalpoly(jak², (12059918, -13211634, 4370421, -690960, 53220)) / 42525
    w11 = 2 * d⁸ * d² * evalpoly(jak², (-5713645538, 6028914206, -1873083006, 277804122, -24770655, 1231947)) / 1403325
    return (x3, x5, x7, x9, x11), (w3, w5, w7, w9, w11)
end

## The soft edge (Airy region)

gl_airy_x1(ak, d, α) = 1/d + ak*(d/4)^(-1/3)
gl_airy_x3(ak, d, α) = ak^2*(d*16)^(1/3)/5 + (11/35-α^2-12/175*ak^3)*d + (16/1575*ak+92/7875*ak^4)*2^(2/3)*d^(5/3)
gl_airy_x5(ak, d, α) = -(15152/3031875*ak^5+1088/121275*ak^2)*2^(1/3)*d^(7/3)

gl_airy_x1(ak, d) = 1/d + ak*(d/4)^(-1/3)
gl_airy_x3(ak, d) = ak^2*(d*16)^(1/3)/5 + (11/35-12/175*ak^3)*d + (16/1575*ak+92/7875*ak^4)*2^(2/3)*d^(5/3)
gl_airy_x5(ak, d) = -(15152/3031875*ak^5+1088/121275*ak^2)*2^(1/3)*d^(7/3)


# Sum the first Q elements of vals, and return the absolute value of the next
# value in the list (or of the last value in the list)
function sum_explicit(vals, Q)
    T = eltype(vals[1])
    z = zero(T)
    for q = min(Q,length(vals)):-1:1
        z += vals[q]
    end
    if Q < length(vals)
        delta = abs(vals[Q+1])
    else
        delta = abs(vals[end])
    end
    z, delta
end

function sum_adaptive(vals)
    z = vals[1]
    i = 1
    while (i < length(vals)) && (abs(vals[i+1]) < abs(vals[i]))
        i += 1
        z += vals[i]
    end
    delta = abs(vals[min(i+1,length(vals))])
    z, delta
end

function gl_bulk_solve_t(n, k, d)
    T = typeof(d)
    pt = (4n-4k+3)*d
    t = T(π)^2/16*(pt-1)^2
    diff = 100
    iter = 0
    maxiter = 20
    while (abs(diff) > 100eps(T)) && (iter < maxiter)
        iter += 1
        diff = (pt*π +2*sqrt(t-t^2) -acos(2*t-1) )*sqrt(t/(1-t))/2
        t -= diff
    end
    if iter == maxiter
        @warn "Maximal number of iterations reached in the computation of t for the bulk"
    end
    t
end

function gausslaguerre_asy_bulk(n, α, k, d, T)
    if α == 0
        return gausslaguerre_asy0_bulk(n, k, d, T)
    end

    t = gl_bulk_solve_t(n, k, d)
    xs, ws = gl_bulk(t, d, α)

    xk, xdelta = (T > 0) ? sum_explicit(xs, (T-1)>>1) : sum_adaptive(xs)
    wk, wdelta = (T > 0) ? sum_explicit(ws, (T-1)>>1) : sum_adaptive(ws)

    xk += t/d

    wfactor = xk^α * exp(-xk) * 2π * sqrt(t/(1-t))
    wk = wfactor * (1+wk)
    wdelta *= wfactor

    xk, wk, max(xdelta,wdelta)
end


function gausslaguerre_asy0_bulk(n, k, d, T)
    t = gl_bulk_solve_t(n, k, d)
    xs, ws = gl_bulk(t, d)


    xk, xdelta = (T > 0) ? sum_explicit(xs, (T-1)>>1) : sum_adaptive(xs)
    wk, wdelta = (T > 0) ? sum_explicit(ws, (T-1)>>1) : sum_adaptive(ws)

    xk += t/d

    wfactor = exp(-xk) * 2π * sqrt(t/(1-t))
    wk = wfactor * (1+wk)
    wdelta *= wfactor

    xk, wk, max(xdelta,wdelta)
end


function gausslaguerre_asy_bessel(n, α, jak, d, T)
    if α == 0
        return gausslaguerre_asy0_bessel(n, jak, d, T)
    end

    xs, ws = gl_bessel(jak, d, α)

    xk, xdelta = (T > 0) ? sum_explicit(xs, (T-1)>>1) : sum_adaptive(xs)
    wk, wdelta = (T > 0) ? sum_explicit(ws, (T-1)>>1) : sum_adaptive(ws)

    xfactor = jak^2 * d
    xk = xfactor * (1 + xk)
    xdelta *= xfactor

    # Invoking the besselj function below is the cause of memory
    # allocation of this routine
    wfactor = 4d * xk^α * exp(-xk) / besselj(α-1, jak)^2
    wk = wfactor * (1 + wk)
    wdelta *= wfactor

    return xk, wk, max(xdelta,wdelta)
end

function gausslaguerre_asy0_bessel(n, jak, d, T)
    xs, ws = gl_bessel(jak, d)
  
    xk, xdelta = (T > 0) ? sum_explicit(xs, (T-1)>>1) : sum_adaptive(xs)
    wk, wdelta = (T > 0) ? sum_explicit(ws, (T-1)>>1) : sum_adaptive(ws)

    xfactor = jak^2 * d
    xk = xfactor * (1 + xk)
    xdelta *= xfactor

    wfactor = 4d * exp(-xk) / besselj(-1, jak)^2
    wk = wfactor * (1 + wk)
    wdelta *= wfactor

    return xk, wk, max(xdelta,wdelta)
end

function compute_airy_root(n, k)
    index = n-k+1
    if index ≤ 11
        ak = AIRY_ROOTS[index]
    else
        t = 3 * π/2 * (index-0.25)
        ak = -t^(2/3) * evalpoly(inv(t*t), (6967296, 725760, -967680, 6478500, -10856875)) / 6967296
    end
    ak
end

function gausslaguerre_asy_airy(n, α, k, d, T)
    if α == 0
        return gausslaguerre_asy0_airy(n, k, d, T)
    end

    ak = compute_airy_root(n, k)
    x1 = gl_airy_x1(ak, d, α)
    x3 = gl_airy_x3(ak, d, α)
    x5 = gl_airy_x5(ak, d, α)

    xs = (x1, x3, x5)

    xk, xdelta = (T > 0) ? sum_explicit(xs, (T+1)>>1) : sum_adaptive(xs)

    wk = 4^(1/3)*xk^(α+1/3)*exp(-xk)/(airyaiprime(ak))^2
    wdelta = abs(wk)

    return xk, wk, max(xdelta,wdelta)
end

function gausslaguerre_asy0_airy(n, k, d, T)
    ak = compute_airy_root(n, k)

    x1 = gl_airy_x1(ak, d)
    x3 = gl_airy_x3(ak, d)
    x5 = gl_airy_x5(ak, d)

    xs = (x1, x3, x5)

    xk, xdelta = (T > 0) ? sum_explicit(xs, (T+1)>>1) : sum_adaptive(xs)

    wk = 4^(1/3) * xk^(1/3) * exp(-xk) / (airyaiprime(ak))^2
    wdelta = abs(wk)

    return xk, wk, max(xdelta,wdelta)
end


"""
Calculate Gauss-Laguerre nodes and weights from the eigenvalue decomposition of
the Jacobi matrix.
"""
function gausslaguerre_GW(n, α)
    _α = 2*(1:n) .+ (α-1)  # 3-term recurrence coeffs a and b
    β = sqrt.( (1:n-1).*(α .+ (1:n-1)) )
    T = SymTridiagonal(promote(collect(_α), β)...)
    x, V = eigen(T)  # eigenvalue decomposition
    w = gamma(α+1)*V[1,:].^2  # Quadrature weights
    return x, vec(w)
end


########################## Routines for the forward recurrence ##########################

function gl_rec_newton(x0, n, α; maxiter = 20, computeweight = true)
    T = eltype(x0)
    step = x0
    iter = 0
    xk = x0

    xk_prev = xk
    pn_prev = floatmax(T)
    pn_deriv = zero(T)
    while (abs(step) > 40eps(T)*xk) && (iter < maxiter)
        iter += 1
        pn, pn_deriv = evalLaguerreRec(n, α, xk)
        if abs(pn) ≥ abs(pn_prev)*(1-50eps(T))
            # The function values do not decrease enough any more due to roundoff errors.
            xk = xk_prev # Set to the previous value and quit.
            break
        end
        step = pn / pn_deriv
        xk_prev = xk
        xk -= step
        pn_prev = pn
    end
    if ( xk < 0 ) || ( xk > 4n + 2α + 2 ) || ( iter == maxiter )
        @warn "Newton method may not have converged in gausslaguerre_rec($n,$α)."
    end
    wk = oftype(xk, 0)
    if computeweight
        pn_min1, ~ = evalLaguerreRec(n-1, α, xk)
        wk = (n^2 +α*n)^(-1/2)/pn_min1/pn_deriv
    end

    return xk, wk
end

"Compute Gauss-Laguerre rule based on the recurrence relation, using Newton iterations on an initial guess."
function gausslaguerre_rec(n, α; reduced = false)
    T = typeof(float(α))

    n_alloc = reduced ? 0 : n
    w = zeros(T, n_alloc)
    x = zeros(T, n_alloc)

    # We compute up to 7 starting values for the Newton iterations
    n_pre = min(n, 7)

    ν = 4n + 2α + 2
    x_pre = T.(approx_besselroots(α, n_pre)).^2 / ν # this is a lower bound by [DLMF 18.16.10]

    noUnderflow = true  # this flag turns false once the weights start to underflow
    for k in 1:n
        # Use sextic extrapolation for a new initial guess
        xk = (k ≤ n_pre) ? x_pre[k] : 7*x[k-1] -21*x[k-2] +35*x[k-3] -35*x[k-4] +21*x[k-5] -7*x[k-6] +x[k-7]

        xk, wk = gl_rec_newton(xk, n, α, maxiter = 20, computeweight = noUnderflow)
        if noUnderflow && abs(wk) < underflow_threshold(T)
            noUnderflow = false
        end

        if reduced
            if !noUnderflow
                return x, w
            else
                push!(x, xk); push!(w, wk)
            end
        else
            x[k] = xk
            w[k] = wk
        end
    end

    return x, w
end


"""
Evaluate the orthonormal associated Laguerre polynomial with positive leading coefficient,
as well as its derivative, in the point x using the recurrence relation.
"""
function evalLaguerreRec(n, α, x)
    T = typeof(α)
    pnprev = zero(T)
    pn = 1/sqrt(gamma(α+1))
    pndprev = zero(T)
    pnd = zero(T)
    for k in 1:n
        pnold = pn
        pn = (x -2*k -α+1)/sqrt(k*(α+k))*pnold-sqrt((k-1+α)*(k-1)/k/(k+α))*pnprev
        pnprev = pnold
        pndold = pnd
        pnd = (pnold+(x-2*k-α+1)*pndold)/sqrt(k*(α+k)) -sqrt((k-1+α)*(k-1)/k/(α+k))*pndprev
        pndprev = pndold
    end

    return pn, pnd
end
