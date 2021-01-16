@doc raw"""
    gausslegendre(ν::Real, n::Integer) -> Vector{Float64}

Return the first ``n`` roots of [Bessel function](https://en.wikipedia.org/wiki/Bessel_function).

```math
J_{\nu}(x) = \sum_{m=0}^{\infty}\frac{(-1)^j}{\Gamma(\nu+j+1)j!} \left(\frac{x}{2}\right)^{2j+\nu}
```

# Examples
```jldoctest; setup = :(using FastGaussQuadrature, SpecialFunctions)
julia> ν = 0.3;

julia> roots = besselroots(ν, 10);

julia> zeros = (x -> besselj(ν, x)).(roots);

julia> all(zeros .< 1e-12)
true
```
"""
function besselroots(ν::Real, n::Integer)
    # FIXME (related issue #22 and #80)
    return besselroots(Float64(ν), n)
end

function besselroots(ν::Float64, n::Integer)
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
            x[k] = J0_roots[k]
        end
        for k in min(n,20)+1:n
            x[k] = McMahon(ν, k)
        end
    elseif -1 ≤ ν ≤ 5
        correctFirstFew = Piessens(ν)
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

# Roots of Bessel funcion ``J_0`` in Float64.
# https://mathworld.wolfram.com/BesselFunctionZeros.html
const J0_roots = @SVector [
    2.4048255576957728,
    5.5200781102863106,
    8.6537279129110122,
    11.791534439014281,
    14.930917708487785,
    18.071063967910922,
    21.211636629879258,
    24.352471530749302,
    27.493479132040254,
    30.634606468431975,
    33.775820213573568,
    36.917098353664044,
    40.058425764628239,
    43.199791713176730,
    46.341188371661814,
    49.482609897397817,
    52.624051841114996,
    55.765510755019979,
    58.906983926080942,
    62.048469190227170,
]


const Piessens_C = @SMatrix [
       2.883975316228  8.263194332307 11.493871452173 14.689036505931 17.866882871378 21.034784308088
       0.767665211539  4.209200330779  4.317988625384  4.387437455306  4.435717974422  4.471319438161
      -0.086538804759 -0.164644722483 -0.130667664397 -0.109469595763 -0.094492317231 -0.083234240394
       0.020433979038  0.039764618826  0.023009510531  0.015359574754  0.011070071951  0.008388073020
      -0.006103761347 -0.011799527177 -0.004987164201 -0.002655024938 -0.001598668225 -0.001042443435
       0.002046841322  0.003893555229  0.001204453026  0.000511852711  0.000257620149  0.000144611721
      -0.000734476579 -0.001369989689 -0.000310786051 -0.000105522473 -0.000044416219 -0.000021469973
       0.000275336751  0.000503054700  0.000083834770  0.000022761626  0.000008016197  0.000003337753
      -0.000106375704 -0.000190381770 -0.000023343325 -0.000005071979 -0.000001495224 -0.000000536428
       0.000042003336  0.000073681222  0.000006655551  0.000001158094  0.000000285903  0.000000088402
      -0.000016858623 -0.000029010830 -0.000001932603 -0.000000269480 -0.000000055734 -0.000000014856
       0.000006852440  0.000011579131  0.000000569367  0.000000063657  0.000000011033  0.000000002536
      -0.000002813300 -0.000004672877 -0.000000169722 -0.000000015222 -0.000000002212 -0.000000000438
       0.000001164419  0.000001903082  0.000000051084  0.000000003677  0.000000000448  0.000000000077
      -0.000000485189 -0.000000781030 -0.000000015501 -0.000000000896 -0.000000000092 -0.000000000014
       0.000000203309  0.000000322648  0.000000004736  0.000000000220  0.000000000019  0.000000000002
      -0.000000085602 -0.000000134047 -0.000000001456 -0.000000000054 -0.000000000004               0
       0.000000036192  0.000000055969  0.000000000450  0.000000000013               0               0
      -0.000000015357 -0.000000023472 -0.000000000140 -0.000000000003               0               0
       0.000000006537  0.000000009882  0.000000000043  0.000000000001               0               0
      -0.000000002791 -0.000000004175 -0.000000000014               0               0               0
       0.000000001194  0.000000001770  0.000000000004               0               0               0
      -0.000000000512 -0.000000000752               0               0               0               0
       0.000000000220  0.000000000321               0               0               0               0
      -0.000000000095 -0.000000000137               0               0               0               0
       0.000000000041  0.000000000059               0               0               0               0
      -0.000000000018 -0.000000000025               0               0               0               0
       0.000000000008  0.000000000011               0               0               0               0
      -0.000000000003 -0.000000000005               0               0               0               0
       0.000000000001  0.000000000002               0               0               0               0
]


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

function Piessens(ν::Float64)
    # Piessens's Chebyshev series approximations (1984). Calculates the 6 first
    # zeros to at least 12 decimal figures in region -1 ≤ ν ≤ 5:
    C = Piessens_C
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
        jk[jj] = J0_roots[jj]
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

const besselJ1_10 = @SVector [
    0.2695141239419169,
    0.1157801385822037,
    0.07368635113640822,
    0.05403757319811628,
    0.04266142901724309,
    0.03524210349099610,
    0.03002107010305467,
    0.02614739149530809,
    0.02315912182469139,
    0.02078382912226786,
]

function besselJ1(m)
    # BESSELJ1 EVALUATE BESSELJ(1,x)^2 AT ROOTS OF BESSELJ(0,x).
    # USE ASYMPTOTICS. Use Taylor series of (NIST, 10.17.3) and McMahon's
    # expansion (NIST, 10.21.19):
    Jk2 = Array{Float64}(undef, m)
    c = (-171497088497 / 15206400, 461797 / 1152, -172913 / 8064,
         151 / 80, -7 / 24, 0.0, 2.0)
    # First 10 are precomputed:
    @inbounds for jj = 1:min(m, 10)
        Jk2[jj] = besselJ1_10[jj]
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
