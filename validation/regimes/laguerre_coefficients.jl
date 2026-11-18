# Reference-only copy of the MIT package coefficient formulas.
# Opsomer & Huybrechs2023 / arXiv2208.11423; preserved independently of outputs.
# Coefficients are evaluated in high reference precision; no numerical errors enter.
function gl_bulk(t, d, α)
    α² = α * α
    tinv = inv(1 - t)
    _t = (1 - t) / t
    _t² = _t * _t
    d² = d * d
    d⁴ = d² * d²

    c1 = evalpoly(α², (-4, 12))
    x3 = -evalpoly(tinv, (c1, -4, 5)) * d / 12
    c1, c2 = evalpoly(α², (224, -960, 480)), evalpoly(α², (7, -30, 15))
    x5 = d² * d * _t / 720 * evalpoly(tinv, (c1, -48 * c2, -16, -576, 2814, -3815, 1600))
    c1, c2, c3 = evalpoly(α², (-285696, 1354752, -967680, 193536)), evalpoly(α², (-31, 147, -105, 21)), evalpoly(α², (-1346, 6405, -4620, 945))
    c4, c5 = evalpoly(α², (43, -126, 63)), evalpoly(α², (-221, -630, 315))
    x7 = -d⁴ * d / 181440 * _t² * evalpoly(tinv, (c1, -23040 * c2, 384 * c3, -320 * c4, 80 * c5, -1727136, 16131880, -48469876, 66424575, -43122800, 10797500))
    c1, c2, c3 = evalpoly(α², (210677760, -1028505600, 812851200, -232243200, 24883200)), evalpoly(α², (127, -620, 490, -140, 15)), evalpoly(α², (1193053, -5826660, 4613070, -1324260, 143325))
    c4, c5, c6 = evalpoly(α², (555239, -2716980, 2163630, -631260, 70875)), evalpoly(α², (-641, 2960, -2155, 450)), evalpoly(α², (-1598, 17685, -13905, 3375))
    c7, c8, c9 = evalpoly(α², (-7823, -9042, 4521)), evalpoly(α², (15948182, -206850, 103425)), evalpoly(α², (64957561, -24000, 12000))
    x9 = (d⁴ * d²) * d / 10886400 * (_t² * _t) * evalpoly(tinv, (c1, -5806080 * c2, 768 * c3, -768 * c4, 16128 * c5, -1792 * c6, 3360 * c7, -192 * c8, 672 * c9, -212307298152, 518401904799, -714465642135, 566519158800, -241928673000, 43222750000))

    w3 = d² / 6 * (2 * t + 3) / (t - 1)^3
    c1 = evalpoly(α², (7, -30, 15))
    w5 = d⁴ * _t² / 720 * evalpoly(tinv, (0, 0, 16 * c1, 32, 1712, -12408, 27517, -24860, 8000))
    c1, c2, c3 = evalpoly(α², (-31, 147, -105, 21)), evalpoly(α², (-416, 1995, -1470, 315)), evalpoly(α², (43, -126, 63))
    c4, c5 = evalpoly(α², (-89, -378, 189)), evalpoly(α², (53752, -630, 315))
    w7 = -(d⁴ * d²) * (_t² * _t) / 90720 * evalpoly(tinv, (0, 0, 2304 * c1, -384 * c2, 480 * c3, -320 * c4, 80 * c5, -50986344, 201908326, -386872990, 393326325, -204917300, 43190000))
    return (x3, x5, x7, x9), (w3, w5, w7)
end

function gl_bulk(t, d)
    tinv = inv(1 - t)
    _t = (1 - t) / t
    _t² = _t * _t
    d² = d * d
    d⁴ = d² * d²

    x3 = -d / 12 * evalpoly(tinv, (-4, -4, 5))
    x5 = d² * d * _t / 720 * evalpoly(tinv, (224, -336, -16, -576, 2814, -3815, 1600))
    x7 = -d⁴ * d / 181440 * _t² * evalpoly(tinv, (-285696, 714240, -516864, -13760, -17680, -1727136, 16131880, -48469876, 66424575, -43122800, 10797500))
    x9 = (d⁴ * d²) * d / 10886400 * (_t² * _t) * evalpoly(tinv, (210677760, -737372160, 916264704, -426423552, -10338048, 2863616, -26285280, -3062050944, 43651480992, -212307298152, 518401904799, -714465642135, 566519158800, -241928673000, 43222750000))

    w3 = d² / 6 * (2 * t + 3) / (t - 1)^3
    w5 = d⁴ / 720 * _t² * evalpoly(tinv, (0, 0, 112, 32, 1712, -12408, 27517, -24860, 8000))
    w7 = -(d⁴ * d²) * (_t² * _t) / 90720 * evalpoly(tinv, (0, 0, -71424, 159744, 20640, 28480, 4300160, -50986344, 201908326, -386872990, 393326325, -204917300, 43190000))
    return (x3, x5, x7, x9), (w3, w5, w7)
end

## The hard edge (Bessel region)
function gl_bessel(jak, d, α)
    jak² = jak * jak
    d² = d * d
    d⁴ = d² * d²
    α² = α * α

    c1 = evalpoly(α², (-2, 2))
    x3 = d² * evalpoly(jak², (c1, 1)) / 3
    w3 = d² * evalpoly(jak², (c1, 2)) / 3
    c1, c2 = evalpoly(α², (94, -140, 46)), evalpoly(α², (-19, 11))
    x5 = d⁴ * evalpoly(jak², (c1, 3 * c2, 11)) / 45
    w5 = d⁴ * evalpoly(jak², (c1, 6 * c2, 33)) / 45
    c1, c2, c3 = evalpoly(α², (-12077, 19887, -9303, 1493)), evalpoly(α², (14051, -10750, 2459)), evalpoly(α², (-181, 73))
    x7 = d⁴ * d² * evalpoly(jak², (4 * c1, 2 * c2, 36 * c3, 657)) / 2835
    w7 = 4 * d⁴ * d² * evalpoly(jak², (c1, c2, 27 * c3, 657)) / 2835
    c1, c2, = evalpoly(α², (6029959, -10087180, 5095482, -1146220, 107959)), evalpoly(α², (-2201939, 1678761, -507801, 63299))
    c3, c4 = evalpoly(α², (1456807, -729422, 125671)), evalpoly(α², (-2879, 887))
    x9 = (d⁴ * d⁴) * evalpoly(jak², (2 * c1, 3 * c2, c3, 60 * c4, 10644)) / 42525
    w9 = (d⁴ * d⁴) * evalpoly(jak², (2 * c1, 6 * c2, 3 * c3, 240 * c4, 53220)) / 42525
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
function gl_airy(ak, d, α)
    T = eltype(α)
    cbrt_d = cbrt(d)
    ak² = ak * ak

    # Airy corrections in powers of d^(1/3), starting from x≈d⁻¹+a_k(4/d)^(1/3).
    x1 = 1 / d + ak * cbrt(T(4)) / cbrt_d
    x3 = ak² * cbrt(T(16)) / 5 * cbrt_d + (T(11) / 35 - α^2 - T(12) / 175 * (ak² * ak)) * d + (T(16) / 1575 * ak + T(92) / 7875 * (ak² * ak²)) * cbrt(T(4)) * d * (cbrt_d * cbrt_d)
    x5 = -(T(15152) / 3031875 * (ak² * ak²) * ak + T(1088) / 121275 * ak²) * cbrt(T(2)) * d * d * cbrt_d
    return x1, x3, x5
end
