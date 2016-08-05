function gausslegendre(n::Int)
    # GAUSSLEGENDRE(n) COMPUTE THE GAUSS-LEGENDRE NODES AND WEIGHTS IN O(n) time.

  if n <= 0
        Float64[], Float64[]
    elseif n == 1
        [0.0], [2.0]
    elseif n == 2
        [-1 / sqrt(3), 1 / sqrt(3)], [1.0, 1.0]
    elseif n == 3
        [-sqrt(3 / 5), 0.0, sqrt(3 / 5)], [5 / 9, 8 / 9, 5 / 9]
    elseif n == 4
        a = 2 / 7 * sqrt(6 / 5)
        ([-sqrt(3 / 7 + a), -sqrt(3/7-a), sqrt(3/7-a), sqrt(3/7+a)],
         [(18 - sqrt(30)) / 36, (18 + sqrt(30)) / 36,
          (18 + sqrt(30)) / 36, (18 - sqrt(30)) / 36])
    elseif n == 5
        b = 2 * sqrt(10 / 7)
        ([-sqrt(5 + b) / 3, -sqrt(5 - b) / 3, 0.0,
          sqrt(5 - b) / 3, sqrt(5 + b) / 3],
         [(322 - 13 * sqrt(70)) / 900, (322 + 13 * sqrt(70)) / 900, 128 / 225,
          (322 + 13 * sqrt(70)) / 900, (322 - 13 * sqrt(70)) / 900])
    elseif n <= 60
        # NEWTON'S METHOD WITH THREE-TERM RECURRENCE:
        rec(n)
    else
        # USE ASYMPTOTIC EXPANSIONS:
        asy(n)
    end
end

function asy(n)
    # COMPUTE GAUSS-LEGENDRE NODES AND WEIGHTS USING ASYMPTOTIC EXPANSIONS.
    # COMPLEXITY O(n).

    # Nodes and weights:
    m = (n + 1) >> 1
    a = besselZeroRoots(m)
    scale!(a, 1 / (n + 0.5))
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
    @inbounds mod(n, 2) == 1 && (x[m] = 0.0)
    x, w
end

function legpts_nodes(n, a)
    # ASYMPTOTIC EXPANSION FOR THE GAUSS-LEGENDRE NODES.
    vn = 1 / (n + 0.5)
    m = length(a)
    nodes = cot(a)
    vn² = vn * vn
    vn⁴ = vn² * vn²
    @inbounds if n <= 255
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
    elseif n <= 3950
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
    nodes
end

function legpts_weights(n, m, a)
    # ASYMPTOTIC EXPANSION FOR THE GAUSS-LEGENDRE WEIGHTS.
    vn = 1 / (n + 0.5)
    vn² = vn^2
    weights = Array(Float64, m)
    if n <= 850000
        @inbounds for i in 1:m
            weights[i] = cot(a[i])
        end
    end
    # Split out the part that can be vectorized by llvm
    @inbounds if n <= 170
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
    elseif n <= 1500
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
    elseif n <= 850000
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
    x,PP2
end

function innerRec(n, x)
    # EVALUATE LEGENDRE AND ITS DERIVATIVE USING THREE-TERM RECURRENCE RELATION.
    N = size(x, 1)
    myPm1 = Array(Float64, N)
    myPPm1 = Array(Float64, N)
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

const besselZeros_20 = [2.4048255576957728, 5.5200781102863106,
                        8.6537279129110122, 11.791534439014281,
                        14.930917708487785, 18.071063967910922,
                        21.211636629879258, 24.352471530749302,
                        27.493479132040254, 30.634606468431975,
                        33.775820213573568, 36.917098353664044,
                        40.058425764628239, 43.199791713176730,
                        46.341188371661814, 49.482609897397817,
                        52.624051841114996, 55.765510755019979,
                        58.906983926080942, 62.048469190227170]

function besselZeroRoots(m)
    # BESSEL0ROOTS ROOTS OF BESSELJ(0,x). USE ASYMPTOTICS.
    # Use McMahon's expansion for the remainder (NIST, 10.21.19):
    jk = Array(Float64, m)
    p = (1071187749376 / 315, 0.0, -401743168 / 105, 0.0, 120928 / 15,
         0.0, -124 / 3, 0.0, 1.0, 0.0)
    # First 20 are precomputed:
    @inbounds for jj = 1:min(m, 20)
        jk[jj] = besselZeros_20[jj]
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

const besselJ1_10 = [0.2695141239419169, 0.1157801385822037,
                     0.07368635113640822, 0.05403757319811628,
                     0.04266142901724309, 0.03524210349099610,
                     0.03002107010305467, 0.02614739149530809,
                     0.02315912182469139, 0.02078382912226786]

function besselJ1(m)
    # BESSELJ1 EVALUATE BESSELJ(1,x)^2 AT ROOTS OF BESSELJ(0,x).
    # USE ASYMPTOTICS. Use Taylor series of (NIST, 10.17.3) and McMahon's
    # expansion (NIST, 10.21.19):
    Jk2 = Array(Float64, m)
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
