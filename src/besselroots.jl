function besselroots(nu::Real, n::Integer)
    # FIXME (related issue #22 and #80)
    return besselroots(Float64(nu), n)
end

function besselroots(nu::Float64, n::Integer)
#BESSELROOTS    The first N roots of the function J_v(x)

# DEVELOPERS NOTES:
#   V = 0 --> Full Float64 precision for N <= 20 (Wolfram Alpha), and very
#     accurate approximations for N > 20 (McMahon's expansion)
#   -1 <= V <= 5 : V ~= 0 --> 12 decimal figures for the 6 first zeros
#     (Piessens's Chebyshev series approximations), and very accurate
#     approximations for the others (McMahon's expansion)
#   V > 5 --> moderately accurate for the 6 first zeros and good
#     approximations for the others (McMahon's expansion)

# This code was originally written by L. L. Peixoto in MATLAB.
# Later modified by A. Townsend to work in Julia

    if n < 0
        throw(DomainError(n, "Input N must be a non-negative integer"))
    end

    x = zeros(n)
    if n > 0 && nu == 0
        for k in 1:min(n,20)
            x[k] = J0_roots[k]
        end
        for k in min(n,20)+1:n
            x[k] = McMahon(nu, k)
        end
    elseif n > 0 && nu >= -1 && nu <= 5
        correctFirstFew = Piessens( nu )
        for k in 1:min(n,6)
            x[k] = correctFirstFew[k]
        end
        for k in min(n,6)+1:n
            x[k] = McMahon(nu, k)
        end
    elseif nu > 5
        for k in 1:n
            x[k] = McMahon(nu, k)
        end
    end
    return x
end

function McMahon(nu::Real, k::Integer)
    # FIXME (related issue #22 and #80)
    return McMahon(Float64(nu), k)
end

function McMahon(nu::Float64, k::Integer)
    # McMahon's expansion. This expansion gives very accurate approximation
    # for the sth zero (s >= 7) in the whole region V >=- 1, and moderate
    # approximation in other cases.
    mu = 4nu^2
    a1 = 1 / 8
    a3 = (7mu-31) / 384
    a5 = 4*(3779+mu*(-982+83mu)) / 61440 # Evaluate via Horner's method.
    a7 = 6*(-6277237+mu*(1585743+mu*(-153855+6949mu))) / 20643840
    a9 = 144*(2092163573+mu*(-512062548+mu*(48010494+mu*(-2479316+70197mu)))) / 11890851840
    a11 = 720 *(-8249725736393+mu*(1982611456181+mu*(-179289628602+mu*(8903961290 +
          mu*(-287149133 + 5592657mu))))) / 10463949619200
    a13 = 576 *(423748443625564327 + mu*(-100847472093088506 + mu*(8929489333108377 +
        mu*(-426353946885548+mu*(13172003634537+mu*(-291245357370 + 4148944183mu)))))) / 13059009124761600
    b = 0.25 * (2nu+4k-1)*pi
    # Evaluate using Horner's scheme:
    x = b - (mu-1)*( ((((((a13/b^2 + a11)/b^2 + a9)/b^2 + a7)/b^2 + a5)/b^2 + a3)/b^2 + a1)/b)
    return x
end

# Roots of Bessel funcion ``J_0`` in Float64.
# https://mathworld.wolfram.com/BesselFunctionZeros.html
const J0_roots =
    [   2.4048255576957728
        5.5200781102863106
        8.6537279129110122
        11.791534439014281
        14.930917708487785
        18.071063967910922
        21.211636629879258
        24.352471530749302
        27.493479132040254
        30.634606468431975
        33.775820213573568
        36.917098353664044
        40.058425764628239
        43.199791713176730
        46.341188371661814
        49.482609897397817
        52.624051841114996
        55.765510755019979
        58.906983926080942
        62.048469190227170  ]


const Piessens_C = [
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
       0.000000000001  0.000000000002               0               0               0               0]


function Piessens(nu::Float64)
    # Piessens's Chebyshev series approximations (1984). Calculates the 6 first
    # zeros to at least 12 decimal figures in region -1 <= V <= 5:
    C = Piessens_C
    T = Array{Float64}(undef,size(C,1))
    pt = (nu-2)/3
    T[1], T[2] = 1., pt
    for k = 2:size(C,1)-1
        T[k+1] = 2*pt*T[k] - T[k-1]
    end
    y = C'*T
    y[1] *= sqrt(nu+1)  # Scale the first root.
    return y
end
