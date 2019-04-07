"""
(x,w) = gausslaguerre(n) returns n Gauss-Laguerre nodes and weights.
(x,w) = gausslaguerre(n, alpha) allows generalized Gauss-Laguerre quadrature.

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
function gausslaguerre(n::Integer, alpha = 0.0; reduced = false)
    if alpha <= -1
        error("The Laguerre parameter α <= -1 corresponds to a nonintegrable weight function")
    end
    if n < 0
        error("gausslaguerre($n,$alpha) not defined: n must be positive.")
    end

    # Guess the numerical type from the supplied type of alpha
    # Although the code is generic, the heuristics are derived for Float64 precision
    T = typeof(float(alpha))
    if n == 0
        T[],T[]
    elseif n == 1
        [1+alpha], [gamma(1+alpha)]
    elseif n == 2
        [alpha + 2-sqrt(alpha+2),alpha+2+sqrt(alpha+2)],
        [((alpha-sqrt(alpha+2)+2)*gamma(alpha+2))/(2*(alpha+2)*(sqrt(alpha+2)-1)^2),
         ((alpha+sqrt(alpha+2)+2)*gamma(alpha+2))/(2*(alpha+2)*(sqrt(alpha+2)+1)^2)]
    elseif n < 15
        # Use Golub-Welsch for small n
        gausslaguerre_GW(n, alpha)
    elseif n < 128
        # Use the recurrence relation for moderate n
        gausslaguerre_rec(n, alpha)
    else
        # Use explicit asymptotic expansions for larger n
        # The restriction to alpha comes from the restriction on nu in besselroots
        if alpha < 5
            gausslaguerre_asy(n, alpha, reduced=reduced, T=-1, recompute=true)
        else
            gausslaguerre_rec(n, alpha)
        end
    end
end

# Our threshold for deciding on underflow
underflow_threshold(x) = underflow_threshold(typeof(x))
underflow_threshold(::Type{T}) where {T <: AbstractFloat} = 10floatmin(T)

# We explicitly store the first 11 roots of the Airy function in double precision
const airy_roots = [-2.338107410459767, -4.08794944413097, -5.520559828095551,
    -6.786708090071759, -7.944133587120853, -9.02265085340981, -10.04017434155809,
    -11.00852430373326, -11.93601556323626, -12.828776752865757, -13.69148903521072]



"""
Compute the Gauss-Laguerre rule using explicit asymptotic expansions for the nodes
and weights.
Optional parameters are:
- `reduced`: compute a reduced quadrature rule, discarding all points and weights
as soon as the weights underflow
- `T`: the order of the expansion. Set `T=-1` to determine the order adaptively
depending on the size of the terms in the expansion
- `recompute`: if a crude measure of the error is larger than a tolerance,
the point and weight are recomputed using the (slower) recursion+newton approach,
yielding more reliable accurate results.
"""
function gausslaguerre_asy(n::Integer, alpha;
    reduced = false,
    T = max(1, ceil(Int, 50/log(n))),       # Heuristic for number of terms
    recompute = false)

    if alpha^2/n > 1
        Compat.@warn "A large value of alpha may lead to inaccurate results."
    end

    ELT = typeof(float(alpha))

    n_alloc = reduced ? 0 : n
    x = zeros(ELT, n_alloc)
    w = zeros(ELT, n_alloc)

    # The expansions are given in powers of 1/(4n+2α+2)
    d = one(ELT)/(4n+2alpha+2)

    # Heuristical indices for Bessel and Airy regions
    k_bessel = max(ceil(Int, sqrt(n) ), 7)
    k_airy = floor(Int, 0.9*n)

    # The Bessel region
    # First compute the roots of the Bessel function of order alpha
    jak_vector = besselroots(alpha, k_bessel)

    bessel_wins = true
    k = 0
    while bessel_wins && k < n
        k += 1
        # We iterate until the estimated error of the bulk expansion is smaller
        # than the one of the Bessel expansion
        jak = (k < k_bessel) ? jak_vector[k] : jak = FastGaussQuadrature.McMahon(alpha, k)

        xk, wk, δ_bessel = gausslaguerre_asy_bessel(n, alpha, jak, d, T)
        xkb, wkb, δ_bulk = gausslaguerre_asy_bulk(n, alpha, k, d, T)
        if δ_bulk < δ_bessel
            bessel_wins = false
            xk = xkb
            wk = wkb
        end
        if recompute
            δ = min(δ_bessel,δ_bulk)
            if δ > 1e-13
                xk_rec, wk_rec = gl_rec_newton(xk, n, alpha)
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
        xk, wk, δ_bulk = gausslaguerre_asy_bulk(n, alpha, k, d, T)
        if recompute
            if δ_bulk > 1e-13
                xk_rec, wk_rec = gl_rec_newton(xk, n, alpha)
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
        xk, wk, δ_bulk = gausslaguerre_asy_bulk(n, alpha, k, d, T)
        xka, wka, δ_airy = gausslaguerre_asy_airy(n, alpha, k, d, T)
        if δ_airy < δ_bulk
            bulk_wins = false
            xk = xka
            wk = wka
        end
        if recompute
            δ = min(δ_airy,δ_bulk)
            if δ > 1e-13
                xk_rec, wk_rec = gl_rec_newton(xk, n, alpha)
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
        xk, wk, δ_airy = gausslaguerre_asy_airy(n, alpha, k, d, T)
        if recompute
            if δ_airy > 1e-13
                xk_rec, wk_rec = gl_rec_newton(xk, n, alpha)
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
    if ( minimum(x) < 0.0 ) || ( maximum(x) > 4*n + 2*alpha + 2 ) ||  ( minimum(diff(x)) <= 0.0 ) || (minimum(w) < 0.0)
        Compat.@warn "Unexpected inconsistency in the computation of nodes and weights"
    end

    x, w
end

## Expansion coefficients
# These are explicit formulas of the coefficients, up to a simple postprocessing
# that is common to all factors and not included here (see below).
#
# General expressions are given in terms of alpha, more specific expressions
# follow for the special case alpha = 0.

## The bulk

# Note: there is always one division by an integer, placed such that it preserves the type of `d`
gl_bulk_x3(t, d, alpha) = -(12*alpha^2 + 5*(1-t)^(-2) - 4*(1-t)^(-1) - 4) * d / 12
gl_bulk_x5(t, d, alpha) = d^3*(1-t)/t/720*(1600*(1-t)^(-6) - 3815*(1-t)^(-5)
    + 480*alpha^4 +2814*(1-t)^(-4) - 576*(1-t)^(-3) - 960*alpha^2 - 48*(15*alpha^4 - 30*alpha^2 + 7)*(1-t)^(-1)
    -16*(1-t)^(-2) + 224)
gl_bulk_x7(t, d, alpha) = -d^5/181440*(1-t)^2/t^2*(10797500*(1-t)^(-10) - 43122800*(1-t)^(-9)
    + 66424575*(1-t)^(-8) -48469876*(1-t)^(-7) + 193536*alpha^6 + 16131880*(1-t)^(-6)
    + 80*(315*alpha^4 - 630*alpha^2 -221)*(1-t)^(-4) - 1727136*(1-t)^(-5)
    - 967680*alpha^4 - 320*(63*alpha^4 - 126*alpha^2 +43)*(1-t)^(-3)
    + 384*(945*alpha^6 - 4620*alpha^4 + 6405*alpha^2 - 1346)*(1-t)^(-2)
    + 1354752*alpha^2   - 23040*(21*alpha^6 - 105*alpha^4 + 147*alpha^2 - 31)*(1-t)^(-1) -285696)
gl_bulk_x9(t, d, alpha) = d^7/10886400*(1-t)^3/t^3*(43222750000*(1-t)^(-14) - 241928673000*(1-t)^(-13)
    + 566519158800*(1-t)^(-12) -714465642135*(1-t)^(-11) + 518401904799*(1-t)^(-10)
    + 672*(12000*alpha^4 - 24000*alpha^2 +64957561)*(1-t)^(-8) - 212307298152*(1-t)^(-9)
    + 24883200*alpha^8 - 192*(103425*alpha^4 -206850*alpha^2 + 15948182)*(1-t)^(-7)
    + 3360*(4521*alpha^4 - 9042*alpha^2 - 7823)*(1-t)^(-6) -232243200*alpha^6
    - 1792*(3375*alpha^6 - 13905*alpha^4  + 17685*alpha^2 - 1598)*(1-t)^(-5)
    + 16128*(450*alpha^6 - 2155*alpha^4 + 2960*alpha^2 - 641)*(1-t)^(-4)
    + 812851200*alpha^4 -768*(70875*alpha^8 - 631260*alpha^6 + 2163630*alpha^4
    - 2716980*alpha^2 +555239)*(1-t)^(-3)  + 768*(143325*alpha^8 - 1324260*alpha^6
    + 4613070*alpha^4 -5826660*alpha^2 + 1193053)*(1-t)^(-2) - 1028505600*alpha^2
    - 5806080*(15*alpha^8 -140*alpha^6 + 490*alpha^4 - 620*alpha^2 + 127)*(1-t)^(-1) + 210677760)

gl_bulk_w3(t, d, alpha) = d^2/6*(2*t + 3)/(t-1)^3
gl_bulk_w5(t, d, alpha) = (1-t)^2/720/t^2*d^4*(8000*(1-t)^(-8) - 24860*(1-t)^(-7) + 27517*(1-t)^(-6)
    - 12408*(1-t)^(-5) + 1712*(1-t)^(-4) +16*(15*alpha^4 - 30*alpha^2 + 7)*(1-t)^(-2) + 32*(1-t)^(-3))
gl_bulk_w7(t, d, alpha) = -(1-t)^3/90720/t^3*d^6*(43190000*(1-t)^(-12) -204917300*(1-t)^(-11)
    + 393326325*(1-t)^(-10) - 386872990*(1-t)^(-9) + 201908326*(1-t)^(-8)
    + 80*(315*alpha^4 - 630*alpha^2 + 53752)*(1-t)^(-6)  - 50986344*(1-t)^(-7)
    - 320*(189*alpha^4 -378*alpha^2 - 89)*(1-t)^(-5) + 480*(63*alpha^4 - 126*alpha^2
    + 43)*(1-t)^(-4)  -384*(315*alpha^6 - 1470*alpha^4 + 1995*alpha^2 - 416)*(1-t)^(-3)
    + 2304*(21*alpha^6 -105*alpha^4 + 147*alpha^2 - 31)*(1-t)^(-2) )

# And for alpha = 0
gl_bulk_x3(t, d) = -d/12*(5*(1-t)^(-2) - 4*(1-t)^(-1) - 4)
gl_bulk_x5(t, d) = d^3*(1-t)/t/720*(1600*(1-t)^(-6) - 3815*(1-t)^(-5) + 2814*(1-t)^(-4)
    - 576*(1-t)^(-3) - 48*7*(1-t)^(-1) -16*(1-t)^(-2) + 224)
gl_bulk_x7(t, d) = -d^5/181440*(1-t)^2/t^2*(10797500*(1-t)^(-10)
    - 43122800*(1-t)^(-9) + 66424575*(1-t)^(-8) -48469876*(1-t)^(-7)
    + 16131880*(1-t)^(-6) - 80*221*(1-t)^(-4) - 1727136*(1-t)^(-5) - 320*43*(1-t)^(-3)
    - 384*1346*(1-t)^(-2) + 23040*31*(1-t)^(-1) -285696)
gl_bulk_x9(t, d) = d^7/10886400*(1-t)^3/t^3*(43222750000*(1-t)^(-14)
    - 241928673000*(1-t)^(-13) + 566519158800*(1-t)^(-12) -714465642135*(1-t)^(-11)
    + 518401904799*(1-t)^(-10) + 672*64957561*(1-t)^(-8)   - 212307298152*(1-t)^(-9)
    - 192*15948182*(1-t)^(-7)  - 3360*7823*(1-t)^(-6) + 1792*1598*(1-t)^(-5)
    + 16128*(- 641)*(1-t)^(-4)  -768*555239*(1-t)^(-3)  + 768*1193053*(1-t)^(-2)
    - 5806080*127*(1-t)^(-1) + 210677760)

gl_bulk_w3(t, d) = d^2/6*(2*t + 3)/(t-1)^3
gl_bulk_w5(t, d) = (1-t)^2/720/t^2*d^4*(8000*(1-t)^(-8) - 24860*(1-t)^(-7)
    + 27517*(1-t)^(-6) - 12408*(1-t)^(-5) + 1712*(1-t)^(-4) +16*7*(1-t)^(-2)
    + 32*(1-t)^(-3))
gl_bulk_w7(t, d) = -(1-t)^3/90720/t^3*d^6*(43190000*(1-t)^(-12)
    - 204917300*(1-t)^(-11) + 393326325*(1-t)^(-10) - 386872990*(1-t)^(-9)
    + 201908326*(1-t)^(-8) +80*53752*(1-t)^(-6)
    - 50986344*(1-t)^(-7) + 320*89*(1-t)^(-5)
    + 480*43*(1-t)^(-4) + 384*416*(1-t)^(-3) - 2304*31*(1-t)^(-2) )


## The hard edge (Bessel region)

gl_bessel_x3(jak, d, alpha) = (jak^2 + 2*alpha^2 - 2)*d^2 / 3
gl_bessel_x5(jak, d, alpha) = (11*jak^4 +3*jak^2*(11*alpha^2-19) +46*alpha^4 -140*alpha^2 +94)*d^4 / 45
gl_bessel_x7(jak, d, alpha) = (657*jak^6 +36*jak^4*(73*alpha^2-181) +2*jak^2*(2459*alpha^4 -10750*alpha^2 +14051)
    + 4*(1493*alpha^6 -9303*alpha^4 +19887*alpha^2 - 12077) )*d^6 / 2835
gl_bessel_x9(jak, d, alpha) = (10644*jak^8 + 60*(887*alpha^2 - 2879)*jak^6 + (125671*alpha^4 -729422*alpha^2 + 1456807)*jak^4
    + 3*(63299*alpha^6 - 507801*alpha^4 + 1678761*alpha^2 - 2201939)*jak^2 + 2*(107959*alpha^8
    - 1146220*alpha^6 + 5095482*alpha^4 -10087180*alpha^2 + 6029959) )*d^8 / 42525

gl_bessel_w3(jak, d, alpha) = (alpha^2 + jak^2 -1)*2*d^2 / 3
gl_bessel_w5(jak, d, alpha) = (46*alpha^4 + 33*jak^4 +6*jak^2*(11*alpha^2 -19) -140*alpha^2 +94)*d^4 / 45
gl_bessel_w7(jak, d, alpha) = (1493*alpha^6 + 657*jak^6 + 27*(73*alpha^2 - 181)*jak^4 - 9303*alpha^4
    + (2459*alpha^4 -10750*alpha^2 + 14051)*jak^2 + 19887*alpha^2 - 12077)*4*d^6 / 2835
gl_bessel_w9(jak, d, alpha) = (215918*alpha^8 + 53220*jak^8 + 240*(887*alpha^2 - 2879)*jak^6 -2292440*alpha^6 +
    3*(125671*alpha^4 - 729422*alpha^2 + 1456807)*jak^4 + 10190964*alpha^4  +
    6*(63299*alpha^6 - 507801*alpha^4 + 1678761*alpha^2 -2201939)*jak^2 -
    20174360*alpha^2 + 12059918)*d^8 / 42525

# And for alpha = 0:
gl_bessel_x3(jak, d) = (jak^2 - 2)*d^2 / 3
gl_bessel_x5(jak, d) = (11*jak^4 - 57*jak^2 + 94)d^4 / 45
gl_bessel_x7(jak, d) = (657*jak^6 - 6516*jak^4 + 28102*jak^2 - 48308)*d^6 / 2835
gl_bessel_x9(jak, d) = (10644*jak^8 - 172740*jak^6 + 1456807*jak^4 -  6605817*jak^2 + 12059918)*d^8 / 42525
gl_bessel_x11(jak, d) = (410649*jak^10 -  9908262*jak^8 + 138902061*jak^6 - 1248722004*jak^4 + 6028914206*jak^2 - 11427291076)*d^10 / 1403325

gl_bessel_w3(jak, d) = (jak^2 - 1)*2*d^2 / 3
gl_bessel_w5(jak, d) = (33*jak^4 -114*jak^2 + 94)*d^4 / 45
gl_bessel_w7(jak, d) = (657*jak^6 - 4887*jak^4 + 14051*jak^2 - 12077)*4*d^6 / 2835
gl_bessel_w9(jak, d) = (53220*jak^8 - 690960*jak^6 + 4370421*jak^4 - 13211634*jak^2 + 12059918)*d^8 / 42525
gl_bessel_w11(jak, d) = (1231947*jak^10 - 24770655*jak^8 + 277804122*jak^6 - 1873083006*jak^4 + 6028914206*jak^2 - 5713645538)*2*d^10 / 1403325


## The soft edge (Airy region)

gl_airy_x1(ak, d, alpha) = 1/d + ak*(d/4)^(-1/3)
gl_airy_x3(ak, d, alpha) = ak^2*(d*16)^(1/3)/5 + (11/35-alpha^2-12/175*ak^3)*d + (16/1575*ak+92/7875*ak^4)*2^(2/3)*d^(5/3)
gl_airy_x5(ak, d, alpha) = -(15152/3031875*ak^5+1088/121275*ak^2)*2^(1/3)*d^(7/3)

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
    t = T(pi)^2/16*(pt-1)^2
    diff = 100
    iter = 0
    maxiter = 20
    while (abs(diff) > 100eps(T)) && (iter < maxiter)
        iter += 1
        diff = (pt*pi +2*sqrt(t-t^2) -acos(2*t-1) )*sqrt(t/(1-t))/2
        t -= diff
    end
    if iter == maxiter
        Compat.@warn "Maximal number of iterations reached in the computation of t for the bulk"
    end
    t
end

function gausslaguerre_asy_bulk(n, alpha, k, d, T)
    if alpha == 0
        return gausslaguerre_asy0_bulk(n, k, d, T)
    end

    t = gl_bulk_solve_t(n, k, d)
    x3 = gl_bulk_x3(t, d, alpha)
    x5 = gl_bulk_x5(t, d, alpha)
    x7 = gl_bulk_x7(t, d, alpha)
    x9 = gl_bulk_x9(t, d, alpha)
    w3 = gl_bulk_w3(t, d, alpha)
    w5 = gl_bulk_w5(t, d, alpha)
    w7 = gl_bulk_w7(t, d, alpha)

    xs = (x3, x5, x7, x9)
    ws = (w3, w5, w7)

    xk, xdelta = (T > 0) ? sum_explicit(xs, (T-1)>>1) : sum_adaptive(xs)
    wk, wdelta = (T > 0) ? sum_explicit(ws, (T-1)>>1) : sum_adaptive(ws)

    xk += t/d

    wfactor = xk^alpha * exp(-xk) * 2pi * sqrt(t/(1-t))
    wk = wfactor * (1+wk)
    wdelta *= wfactor

    xk, wk, max(xdelta,wdelta)
end


function gausslaguerre_asy0_bulk(n, k, d, T)
    t = gl_bulk_solve_t(n, k, d)
    x3 = gl_bulk_x3(t, d)
    x5 = gl_bulk_x5(t, d)
    x7 = gl_bulk_x7(t, d)
    x9 = gl_bulk_x9(t, d)
    w3 = gl_bulk_w3(t, d)
    w5 = gl_bulk_w5(t, d)
    w7 = gl_bulk_w7(t, d)

    xs = (x3, x5, x7, x9)
    ws = (w3, w5, w7)

    xk, xdelta = (T > 0) ? sum_explicit(xs, (T-1)>>1) : sum_adaptive(xs)
    wk, wdelta = (T > 0) ? sum_explicit(ws, (T-1)>>1) : sum_adaptive(ws)

    xk += t/d

    wfactor = exp(-xk) * 2pi * sqrt(t/(1-t))
    wk = wfactor * (1+wk)
    wdelta *= wfactor

    xk, wk, max(xdelta,wdelta)
end


function gausslaguerre_asy_bessel(n, alpha, jak, d, T)
    if alpha == 0
        return gausslaguerre_asy0_bessel(n, jak, d, T)
    end
    x3 = gl_bessel_x3(jak, d, alpha)
    x5 = gl_bessel_x5(jak, d, alpha)
    x7 = gl_bessel_x7(jak, d, alpha)
    x9 = gl_bessel_x9(jak, d, alpha)
    w3 = gl_bessel_w3(jak, d, alpha)
    w5 = gl_bessel_w5(jak, d, alpha)
    w7 = gl_bessel_w7(jak, d, alpha)
    w9 = gl_bessel_w9(jak, d, alpha)

    xs = (x3, x5, x7, x9)
    ws = (w3, w5, w7, w9)

    xk, xdelta = (T > 0) ? sum_explicit(xs, (T-1)>>1) : sum_adaptive(xs)
    wk, wdelta = (T > 0) ? sum_explicit(ws, (T-1)>>1) : sum_adaptive(ws)

    xfactor = jak^2 * d
    xk = xfactor * (1 + xk)
    xdelta *= xfactor

    # Invoking the besselj function below is the cause of memory
    # allocation of this routine
    wfactor = 4d * xk^alpha * exp(-xk) / besselj(alpha-1, jak)^2
    wk = wfactor * (1 + wk)
    wdelta *= wfactor

    xk, wk, max(xdelta,wdelta)
end

function gausslaguerre_asy0_bessel(n, jak, d, T)
    x3 = gl_bessel_x3(jak, d)
    x5 = gl_bessel_x5(jak, d)
    x7 = gl_bessel_x7(jak, d)
    x9 = gl_bessel_x9(jak, d)
    x11 = gl_bessel_x11(jak, d)
    w3 = gl_bessel_w3(jak, d)
    w5 = gl_bessel_w5(jak, d)
    w7 = gl_bessel_w7(jak, d)
    w9 = gl_bessel_w9(jak, d)
    w11 = gl_bessel_w11(jak, d)

    xs = (x3, x5, x7, x9, x11)
    ws = (w3, w5, w7, w9, w11)

    xk, xdelta = (T > 0) ? sum_explicit(xs, (T-1)>>1) : sum_adaptive(xs)
    wk, wdelta = (T > 0) ? sum_explicit(ws, (T-1)>>1) : sum_adaptive(ws)

    xfactor = jak^2 * d
    xk = xfactor * (1 + xk)
    xdelta *= xfactor

    wfactor = 4d * exp(-xk) / besselj(-1, jak)^2
    wk = wfactor * (1 + wk)
    wdelta *= wfactor

    xk, wk, max(xdelta,wdelta)
end

function compute_airy_root(n, k)
    index = n-k+1
    if index <= 11
        ak = airy_roots[index]
    else
        t = 3 * pi/2 * (index-0.25)
        ak = -t^(2/3)*(1 + 5/48/t^2 - 5/36/t^4 + 77125/82944/t^6 -10856875/6967296/t^8)
    end
    ak
end

function gausslaguerre_asy_airy(n, alpha, k, d, T)
    if alpha == 0
        return gausslaguerre_asy0_airy(n, k, d, T)
    end

    ak = compute_airy_root(n, k)
    x1 = gl_airy_x1(ak, d, alpha)
    x3 = gl_airy_x3(ak, d, alpha)
    x5 = gl_airy_x5(ak, d, alpha)

    xs = (x1, x3, x5)

    xk, xdelta = (T > 0) ? sum_explicit(xs, (T+1)>>1) : sum_adaptive(xs)

    wk = 4^(1/3)*xk^(alpha+1/3)*exp(-xk)/(airyaiprime(ak))^2
    wdelta = abs(wk)

    xk, wk, max(xdelta,wdelta)
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

    xk, wk, max(xdelta,wdelta)
end




"""
Calculate Gauss-Laguerre nodes and weights from the eigenvalue decomposition of
the Jacobi matrix.
"""
function gausslaguerre_GW(n, alpha)
    alph = 2*(1:n) .+ (alpha-1)     # 3-term recurrence coeffs a and b
    beta = sqrt.( (1:n-1).*(alpha .+ (1:n-1)) )
    T = SymTridiagonal(Vector(alph), beta)  # Jacobi matrix
    x, V = eigen(T)                 # eigenvalue decomposition
    w = gamma(alpha+1)*V[1,:].^2    # Quadrature weights
    x, vec(w)
end


########################## Routines for the forward recurrence ##########################

function gl_rec_newton(x0, n, alpha; maxiter = 20, computeweight = true)
    T = eltype(x0)
    step = x0
    iter = 0
    xk = x0

    xk_prev = xk
    pn_prev = floatmax(T)
    pn_deriv = zero(T)
    while (abs(step) > 40eps(T)*xk) && (iter < maxiter)
        iter += 1
        pn, pn_deriv = evalLaguerreRec(n, alpha, xk)
        if abs(pn) >= abs(pn_prev)*(1-50eps(T))
            # The function values do not decrease enough any more due to roundoff errors.
            xk = xk_prev # Set to the previous value and quit.
            break
        end
        step = pn / pn_deriv
        xk_prev = xk
        xk -= step
        pn_prev = pn
    end
    if ( xk < 0 ) || ( xk > 4n + 2alpha + 2 ) || ( iter == maxiter )
        Compat.@warn "Newton method may not have converged in gausslaguerre_rec($n,$alpha)."
    end
    wk = oftype(xk, 0)
    if computeweight
        pn_min1, ~ = evalLaguerreRec(n-1, alpha, xk)
        wk = (n^2 +alpha*n)^(-1/2)/pn_min1/pn_deriv
    end
    xk, wk
end

"Compute Gauss-Laguerre rule based on the recurrence relation, using Newton iterations on an initial guess."
function gausslaguerre_rec(n, alpha; reduced = false)
    T = typeof(float(alpha))

    n_alloc = reduced ? 0 : n
    w = zeros(T, n_alloc)
    x = zeros(T, n_alloc)

    # We compute up to 7 starting values for the Newton iterations
    n_pre = min(n, 7)

    nu = 4n + 2alpha + 2
    x_pre = T.(besselroots(alpha, n_pre)).^2 / nu # this is a lower bound by [DLMF 18.16.10]

    noUnderflow = true      # this flag turns false once the weights start to underflow
    for k in 1:n
        local pn_deriv

        # Use sextic extrapolation for a new initial guess
        xk = (k <= n_pre) ? x_pre[k] : 7*x[k-1] -21*x[k-2] +35*x[k-3] -35*x[k-4] +21*x[k-5] -7*x[k-6] +x[k-7]

        xk, wk = gl_rec_newton(xk, n, alpha, maxiter = 20, computeweight = noUnderflow)
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
            x[k] = xk; w[k] = wk
        end
    end
    x, w
end


"""
Evaluate the orthonormal associated Laguerre polynomial with positive leading coefficient,
as well as its derivative, in the point x using the recurrence relation.
"""
function evalLaguerreRec(n, alpha, x)
    T = typeof(alpha)
    pnprev = zero(T)
    pn = 1/sqrt(gamma(alpha+1))
    pndprev = zero(T)
    pnd = zero(T)
    for k in 1:n
        pnold = pn
        pn = (x -2*k -alpha+1)/sqrt(k*(alpha+k))*pnold-sqrt((k-1+alpha)*(k-1)/k/(k+alpha))*pnprev
        pnprev = pnold
        pndold = pnd
        pnd = (pnold+(x-2*k-alpha+1)*pndold)/sqrt(k*(alpha+k)) -sqrt((k-1+alpha)*(k-1)/k/(alpha+k))*pndprev
        pndprev = pndold
    end
    pn, pnd
end
