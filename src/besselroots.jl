@doc raw"""
    approx_besselroots(ν::Real, n::Integer) -> Vector{Float64}

Return the first ``n`` roots of [Bessel function](https://en.wikipedia.org/wiki/Bessel_function).
Note that this function is only 12-digits of precision.

```math
J_{\nu}(x) = \sum_{m=0}^{\infty}\frac{(-1)^j}{\Gamma(\nu+j+1)j!} \left(\frac{x}{2}\right)^{2j+\nu}
```

# Examples
```jldoctest
julia> ν = 0.3;

julia> roots = approx_besselroots(ν, 10);

julia> zeros = (x -> besselj(ν, x)).(roots);

julia> all(zeros .< 1e-12)
true
```
"""
function approx_besselroots(ν::Real, n::Integer)
    # FIXME (related issue #22 and #80)
    return approx_besselroots(Float64(ν), n)
end

function approx_besselroots(ν::Float64, n::Integer)
# DEVELOPERS NOTES:
#   ν = 0 --> Full Float64 precision for n ≤ 20 (Wolfram Alpha), and very
#     accurate approximations for n > 20 (McMahon's expansion)
#   -1 ≤ ν ≤ 5 : ν ~= 0 --> 12 decimal figures for the 6 first zeros
#     (Piessens's Chebyshev series approximations), and very accurate
#     approximations for the others (McMahon's expansion)
#   ν > 5 --> moderately accurate for the 6 first zeros and good
#     approximations for the others (McMahon's expansion)

# This code was originally written by L. L. Peixoto in MATLAB.
# Later modified by A. Townsend to work in Julia

    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    end

    x = zeros(n)
    if ν == 0
        for k in 1:min(n,20)
            x[k] = BESSELJ0_ROOTS[k]
        end
        for k in min(n,20)+1:n
            x[k] = McMahon(ν, k)
        end
    elseif -1 ≤ ν ≤ 5
        correctFirstFew = piessens(ν)
        for k in 1:min(n,6)
            x[k] = correctFirstFew[k]
        end
        for k in min(n,6)+1:n
            x[k] = McMahon(ν, k)
        end
    elseif 5 < ν
        for k in 1:n
            x[k] = McMahon(ν, k)
        end
    end
    return x
end

function McMahon(ν::Real, k::Integer)
    # FIXME (related issue #22 and #80)
    return McMahon(Float64(ν), k)
end

function McMahon(ν::Float64, k::Integer)
    # McMahon's expansion. This expansion gives very accurate approximation
    # for the sth zero (s ≥ 7) in the whole region ν ≥ -1, and moderate
    # approximation in other cases.
    μ = 4ν^2
    a1 = 1 / 8
    a3 = (7μ-31) / 384
    a5 = 4*(3779+μ*(-982+83μ)) / 61440 # Evaluate via Horner's method.
    a7 = 6*(-6277237+μ*(1585743+μ*(-153855+6949μ))) / 20643840
    a9 = 144*(2092163573+μ*(-512062548+μ*(48010494+μ*(-2479316+70197μ)))) / 11890851840
    a11 = 720 *(-8249725736393+μ*(1982611456181+μ*(-179289628602+μ*(8903961290 +
          μ*(-287149133 + 5592657μ))))) / 10463949619200
    a13 = 576 *(423748443625564327 + μ*(-100847472093088506 + μ*(8929489333108377 +
        μ*(-426353946885548+μ*(13172003634537+μ*(-291245357370 + 4148944183μ)))))) / 13059009124761600
    b = 0.25 * (2ν+4k-1)*π
    # Evaluate using Horner's scheme:
    x = b - (μ-1)*( ((((((a13/b^2 + a11)/b^2 + a9)/b^2 + a7)/b^2 + a5)/b^2 + a3)/b^2 + a1)/b)
    return x
end

function _piessens_chebyshev30(ν::Float64)
    # Piessens's Chebyshev series approximations (1984). Calculates the 6 first
    # zeros to at least 12 decimal figures in region -1 ≤ ν ≤ 5:
    pt = (ν-2)/3

    T1 = 1.0
    T2 = pt
    T3 = 2pt*T2 - T1
    T4 = 2pt*T3 - T2
    T5 = 2pt*T4 - T3
    T6 = 2pt*T5 - T4
    T7 = 2pt*T6 - T5
    T8 = 2pt*T7 - T6
    T9 = 2pt*T8 - T7
    T10 = 2pt*T9 - T8
    T11 = 2pt*T10 - T9
    T12 = 2pt*T11 - T10
    T13 = 2pt*T12 - T11
    T14 = 2pt*T13 - T12
    T15 = 2pt*T14 - T13
    T16 = 2pt*T15 - T14
    T17 = 2pt*T16 - T15
    T18 = 2pt*T17 - T16
    T19 = 2pt*T18 - T17
    T20 = 2pt*T19 - T18
    T21 = 2pt*T20 - T19
    T22 = 2pt*T21 - T20
    T23 = 2pt*T22 - T21
    T24 = 2pt*T23 - T22
    T25 = 2pt*T24 - T23
    T26 = 2pt*T25 - T24
    T27 = 2pt*T26 - T25
    T28 = 2pt*T27 - T26
    T29 = 2pt*T28 - T27
    T30 = 2pt*T29 - T28

    T = SVector(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17,T18,T19,T20,T21,T22,T23,T24,T25,T26,T27,T28,T29,T30)
    return T
end

function piessens(ν::Float64)
    # Piessens's Chebyshev series approximations (1984). Calculates the 6 first
    # zeros to at least 12 decimal figures in region -1 ≤ ν ≤ 5:
    C = PIESSENS_C
    T = _piessens_chebyshev30(ν)
    y = C'*T
    _y = collect(y)
    _y[1] *= sqrt(ν+1)  # Scale the first root.
    return _y
end


function besselZeroRoots(m)
    # BESSEL0ROOTS ROOTS OF BESSELJ(0,x). USE ASYMPTOTICS.
    # Use McMahon's expansion for the remainder (NIST, 10.21.19):
    jk = Array{Float64}(undef, m)
    p = (1071187749376 / 315, 0.0, -401743168 / 105, 0.0, 120928 / 15,
         0.0, -124 / 3, 0.0, 1.0, 0.0)
    # First 20 are precomputed:
    @inbounds for jj = 1:min(m, 20)
        jk[jj] = BESSELJ0_ROOTS[jj]
    end
    @inbounds for jj = 21:min(m, 47)
        ak = π * (jj - .25)
        ak82 = (.125 / ak)^2
        jk[jj] = ak + .125 / ak * @evalpoly(ak82, 1.0, p[7], p[5], p[3])
    end
    @inbounds for jj = 48:min(m, 344)
        ak = π * (jj - .25)
        ak82 = (.125 / ak)^2
        jk[jj] = ak + .125 / ak * @evalpoly(ak82, 1.0, p[7], p[5])
    end
    @inbounds for jj = 345:min(m,13191)
        ak = π * (jj - .25)
        ak82 = (.125 / ak)^2
        jk[jj] = ak + .125 / ak * muladd(ak82, p[7], 1.0)
    end
    @inbounds for jj = 13192:m
        ak = π * (jj - .25)
        jk[jj] = ak + .125 / ak
    end
    return jk
end

function besselJ1(m)
    # BESSELJ1 EVALUATE BESSELJ(1,x)^2 AT ROOTS OF BESSELJ(0,x).
    # USE ASYMPTOTICS. Use Taylor series of (NIST, 10.17.3) and McMahon's
    # expansion (NIST, 10.21.19):
    Jk2 = Array{Float64}(undef, m)
    c = (-171497088497 / 15206400, 461797 / 1152, -172913 / 8064,
         151 / 80, -7 / 24, 0.0, 2.0)
    # First 10 are precomputed:
    @inbounds for jj = 1:min(m, 10)
        Jk2[jj] = BESSELJ1_ON_BESSELJ0_ROOTS[jj]
    end
    @inbounds for jj = 11:min(m, 15)
        ak = π * (jj - .25)
        ak2 = (1 / ak)^2
        Jk2[jj] = 1 / (π * ak) * muladd(@evalpoly(ak2, c[5], c[4], c[3],
                                                  c[2], c[1]), ak2^2, c[7])
    end
    @inbounds for jj = 16:min(m, 21)
        ak = π * (jj - .25)
        ak2 = (1 / ak)^2
        Jk2[jj] = 1 / (π * ak) * muladd(@evalpoly(ak2, c[5], c[4], c[3], c[2]),
                                        ak2^2, c[7])
    end
    @inbounds for jj = 22:min(m,55)
        ak = π * (jj - .25)
        ak2 = (1 / ak)^2
        Jk2[jj] = 1 / (π * ak) * muladd(@evalpoly(ak2, c[5], c[4], c[3]),
                                        ak2^2, c[7])
    end
    @inbounds for jj = 56:min(m,279)
        ak = π * (jj - .25)
        ak2 = (1 / ak)^2
        Jk2[jj] = 1 / (π * ak) * muladd(muladd(ak2, c[4], c[5]), ak2^2, c[7])
    end
    @inbounds for jj = 280:min(m,2279)
        ak = π * (jj - .25)
        ak2 = (1 / ak)^2
        Jk2[jj] = 1 / (π * ak) * muladd(ak2^2, c[5], c[7])
    end
    @inbounds for jj = 2280:m
        ak = π * (jj - .25)
        Jk2[jj] = 1 / (π * ak) * c[7]
    end
    return Jk2
end

function besselroots(ν::Real, n::Integer)
    @warn "`besselroots` has been renamed to `approx_besselroots`, and `besselroots` will be removed in the next major release."
    return approx_besselroots(ν, n)
end
