@doc raw"""
    approx_besselroots(ν::Real, n::Integer) -> Vector{Float64}

Approximate the first ``n`` positive roots of the Bessel function ``J_\nu``,
for `n ≥ 0` and `ν > -1`. The result always uses `Float64`. The approximations
target about 12 decimal digits for the first few roots when `-1 < ν ≤ 5`;
accuracy depends on the parameter and root index, and can be lower for `ν > 5`.

```math
J_{\nu}(x) = \sum_{j=0}^{\infty}\frac{(-1)^j}{\Gamma(\nu+j+1)j!} \left(\frac{x}{2}\right)^{2j+\nu}
```

# Examples
```jldoctest
julia> ν = 0.3;

julia> roots = approx_besselroots(ν, 10);

julia> zeros = (x -> besselj(ν, x)).(roots);

julia> all(abs.(zeros) .< 1e-12)
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
        for k in 1:min(n, 20)
            x[k] = BESSELJ0_ROOTS[k]
        end
        for k in (min(n, 20) + 1):n
            x[k] = McMahon(ν, k)
        end
    elseif -1 ≤ ν ≤ 5
        correctFirstFew = piessens(ν)
        for k in 1:min(n, 6)
            x[k] = correctFirstFew[k]
        end
        for k in (min(n, 6) + 1):n
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
    a3 = (7μ - 31) / 384
    a5 = 4 * (3779 + μ * (-982 + 83μ)) / 61440 # Evaluate via Horner's method.
    a7 = 6 * (-6277237 + μ * (1585743 + μ * (-153855 + 6949μ))) / 20643840
    a9 = 144 * (2092163573 + μ * (-512062548 + μ * (48010494 + μ * (-2479316 + 70197μ)))) / 11890851840
    a11 = 720 * (
        -8249725736393 + μ * (
            1982611456181 + μ * (
                -179289628602 + μ * (
                    8903961290 + μ * (
                        -287149133 + μ * (
                            5592657
                        )
                    )
                )
            )
        )
    ) / 10463949619200
    a13 = 576 * (
        423748443625564327 + μ * (
            -100847472093088506 + μ * (
                8929489333108377 + μ * (
                    -426353946885548 + μ * (
                        13172003634537 + μ * (
                            -291245357370 + μ * (
                                4148944183
                            )
                        )
                    )
                )
            )
        )
    ) / 13059009124761600
    # β=(k+ν/2-1/4)π; the corrections use odd inverse powers of β.
    b = 0.25 * (2ν + 4k - 1) * π
    # Evaluate using Horner's scheme:
    return b - (μ - 1) * (((((((a13 / b^2 + a11) / b^2 + a9) / b^2 + a7) / b^2 + a5) / b^2 + a3) / b^2 + 1 / 8) / b)
end

function _piessens_chebyshev30(ν::Float64)
    # Piessens's Chebyshev series approximations (1984). Calculates the 6 first
    # zeros to at least 12 decimal figures in region -1 ≤ ν ≤ 5:
    pt = (ν - 2) / 3

    # T₀(pt), …, T₂₉(pt), using T_{j+1}=2pt*T_j-T_{j-1}.
    T1 = 1.0
    T2 = pt
    T3 = 2pt * T2 - T1
    T4 = 2pt * T3 - T2
    T5 = 2pt * T4 - T3
    T6 = 2pt * T5 - T4
    T7 = 2pt * T6 - T5
    T8 = 2pt * T7 - T6
    T9 = 2pt * T8 - T7
    T10 = 2pt * T9 - T8
    T11 = 2pt * T10 - T9
    T12 = 2pt * T11 - T10
    T13 = 2pt * T12 - T11
    T14 = 2pt * T13 - T12
    T15 = 2pt * T14 - T13
    T16 = 2pt * T15 - T14
    T17 = 2pt * T16 - T15
    T18 = 2pt * T17 - T16
    T19 = 2pt * T18 - T17
    T20 = 2pt * T19 - T18
    T21 = 2pt * T20 - T19
    T22 = 2pt * T21 - T20
    T23 = 2pt * T22 - T21
    T24 = 2pt * T23 - T22
    T25 = 2pt * T24 - T23
    T26 = 2pt * T25 - T24
    T27 = 2pt * T26 - T25
    T28 = 2pt * T27 - T26
    T29 = 2pt * T28 - T27
    T30 = 2pt * T29 - T28

    return SVector(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25, T26, T27, T28, T29, T30)
end

function piessens(ν::Float64)
    # Piessens's Chebyshev series approximations (1984). Calculates the 6 first
    # zeros to at least 12 decimal figures in region -1 ≤ ν ≤ 5:
    y = PIESSENS_C' * _piessens_chebyshev30(ν)
    return setindex(y, y[1] * sqrt(ν + 1), 1)
end


function besselZeroRoots(m, ::Type{T} = Float64) where {T}
    # Bessel J₀ roots for small-order Legendre seeds (m ≤ 31, at most 53 bits).
    # Use McMahon's expansion for the remainder (NIST, 10.21.19):
    jk = Vector{T}(undef, m)
    # First 20 are precomputed:
    @inbounds for jj in 1:min(m, 20)
        jk[jj] = BESSELJ0_ROOTS[jj]
    end
    c = T(0.125)
    coefficients = map(T, (1.0, -124 / 3, 120928 / 15, -401743168 / 105))
    @inbounds for jj in 21:m
        ak = T(π) * (T(jj) - T(0.25))
        jk[jj] = ak + c / ak * evalpoly((c / ak)^2, coefficients)
    end
    return jk
end

Base.@deprecate besselroots(ν::Real, n::Integer) approx_besselroots(ν, n)
