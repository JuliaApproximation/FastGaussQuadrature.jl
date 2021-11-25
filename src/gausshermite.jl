@doc raw"""
    gausshermite(n::Integer) -> Tuple{Vector{Float64},Vector{Float64}}

Return nodes and weights of [Gauss-Hermite quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Hermite_quadrature).

```math
\int_{-\infty}^{+\infty} f(x) \exp(-x^2) dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest
julia> x, w = gausshermite(3);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ 3(√π)/4
true
```
"""
function gausshermite(n::Integer)
    x,w = unweightedgausshermite(n)
    w .*= exp.(-x.^2)
    x, w
end

function unweightedgausshermite(n::Integer)
    # GAUSSHERMITE(n) COMPUTE THE GAUSS-HERMITE NODES AND WEIGHTS IN O(n) time.
    if n < 0
        throw(DomainError(n, "Input n must be a non-negative integer"))
    elseif n == 0
        return Float64[],Float64[]
    elseif n == 1
        return [0.0],[sqrt(π)]
    elseif n ≤ 20
       # GW algorithm
       x = hermite_gw(n)
    elseif n ≤ 200
       # REC algorithm
       x = hermite_rec(n)
    else
       # ASY algorithm
       x = hermite_asy(n)
    end

    # fold out
    if isodd(n)
        _w = [reverse(x[2][:]); x[2][2:end]]
        _x = [-reverse(x[1]) ; x[1][2:end]]
    else
        _w = [reverse(x[2][:]); x[2][:]]
        _x = [-reverse(x[1]) ; x[1]]
    end
    _w .*= sqrt(π)/sum(exp.(-_x.^2).*_w)

    return _x, _w
end

function hermite_asy(n::Integer)
    # Compute Hermite nodes and weights using asymptotic formula

    x0 = hermite_initialguess(n)  # get initial guesses
    t0 = x0./sqrt(2n+1)
    θ = acos.(t0)  # convert to θ-variable
    val = x0
    for _ in 1:20
        val = hermpoly_asy_airy(n, θ)
        dθ = -val[1]./(sqrt(2).*sqrt(2n+1).*val[2].*sin.(θ))
        θ .-= dθ  # Newton update
        if norm(dθ,Inf) < sqrt(eps(Float64))/10
           break
        end
    end
    t0 = cos.(θ)
    x = sqrt(2n+1)*t0  #back to x-variable
    w = x.*val[1] .+ sqrt(2).*val[2]
    w .= 1 ./ w.^2  # quadrature weights

    return x, w
end

function hermite_rec(n::Integer)
    # Compute Hermite nodes and weights using recurrence relation.

    x0 = hermite_initialguess(n)
    x0 .*= sqrt(2)
    val = x0
    for _ in 1:10
        val = hermpoly_rec.(n, x0)
        dx = first.(val)./last.(val)
        dx[ isnan.( dx ) ] .= 0
        x0 .= x0 .- dx
        if norm(dx, Inf) < sqrt(eps(Float64))
            break
        end
    end
    x0 ./= sqrt(2)
    w = 1 ./ last.(val).^2  # quadrature weights

    return x0, w
end

function hermpoly_rec(n::Integer, x0)
    # HERMPOLY_rec evaluation of scaled Hermite poly using recurrence
    n < 0 && throw(DomainError(n, "Input n must be a non-negative integer"))
    # evaluate:
    w = exp(-x0^2 / (4*n))
    wc = 0 # 0 times we've applied wc
    Hold = one(x0)
    # n == 0 && return (Hold, 0)
    H = x0
    for k = 1:n-1
        Hold, H = H, (x0*H/sqrt(k+1) - Hold/sqrt(1+1/k))
        while abs(H) ≥ 100 && wc < n  # regularise
            H *= w
            Hold *= w
            wc += 1
        end
        k += 1
    end
    for _ = wc+1:n
        H *= w
        Hold *= w
    end

    return H, -x0*H + sqrt(n)*Hold
end

function hermpoly_rec(r::Base.OneTo, x0)
    isempty(r) && return [1.0]
    n = maximum(r)
    # HERMPOLY_rec evaluation of scaled Hermite poly using recurrence
    n < 0 && throw(DomainError(n, "Input n must be a non-negative integer"))
    n == 0 && return [exp(-x0^2 / 4)]
    p = max(1,floor(Int,x0^2/100))
    w = exp(-x0^2 / 4p)
    wc = 0 # 0 times we've applied wc
    ret = Vector{Float64}()
    Hold = one(x0)
    push!(ret, Hold)
    H = x0
    push!(ret, H)
    for k = 1:n-1
        Hold, H = H, (x0*H/sqrt(k+1) - Hold/sqrt(1+1/k))
        while abs(H) ≥ 100 && wc < p  # regularise
            ret .*= w
            H *= w
            Hold *= w
            wc += 1
        end
        push!(ret, H)
        k += 1
    end
    ret .*= w^(p-wc)

    return ret
end

hermpoly_rec(r::AbstractRange, x0) = hermpoly_rec(Base.OneTo(maximum(r)), x0)[r.+1]

function hermpoly_asy_airy(n::Integer, θ::AbstractVector)
    # HERMPOLY_ASY evaluation hermite poly using Airy asymptotic formula in θ-space.

    musq = 2n+1
    cosT = cos.(θ)
    sinT = sin.(θ)
    sin2T = 2 .* cosT.*sinT
    η = 0.5 .* θ .- 0.25 .* sin2T
    χ = -(3*η/2).^(2/3)
    φ = (-χ./sinT.^2).^(1/4)
    C = 2*sqrt(π)*musq^(1/6)*φ
    Airy0 = real.(airyai.(musq.^(2/3).*χ))
    Airy1 = real.(airyaiprime.(musq.^(2/3).*χ))

    # Terms in (12.10.43):
    a0 = 1
    b0 = 1
    a1 = 15/144
    b1 = -7/5*a1
    a2 = 5*7*9*11/2/144^2
    b2 = -13/11*a2
    a3 = 7*9*11*13*15*17/6/144^3
    b3 = -19/17*a3

    # u polynomials in (12.10.9)
    u0 = 1
    u1 = (cosT.^3-6*cosT)/24
    u2 = (-9*cosT.^4 + 249*cosT.^2 .+ 145)/1152
    u3 = (-4042*cosT.^9+18189*cosT.^7-28287*cosT.^5-151995*cosT.^3-259290*cosT)/414720

    # first term
    A0 = 1
    val = A0*Airy0

    # second term
    B0 = -(a0*u1.*φ.^6 .+ a1*u0) ./ χ.^2
    val .+=  B0.*Airy1./musq.^(4/3)

    # third term
    A1 = (b0*u2.*φ.^12 + b1*u1.*φ.^6 .+ b2*u0) ./ χ.^3
    val .+= A1.*Airy0/musq.^2

    # fourth term
    B1 = -(u3.*φ.^18 + a1*u2.*φ.^12 + a2*u1.*φ.^6 .+ a3*u0) ./ χ.^5
    val .+= B1.*Airy1./musq.^(4/3+2)

    val .= C.*val

    ## Derivative
    η = .5*θ - .25*sin2T
    χ = -(3*η/2).^(2/3)
    φ = (-χ./sinT.^2).^(1/4)
    C = sqrt(2*π)*musq^(1/3)./φ

    # v polynomials in (12.10.10)
    v0 = 1
    v1 = (cosT.^3+6cosT)/24
    v2 = (15*cosT.^4 - 327*cosT.^2 .- 143)/1152
    v3 = (259290*cosT + 238425*cosT.^3 - 36387*cosT.^5 + 18189*cosT.^7 - 4042*cosT.^9)/414720

    # first term
    C0 = -(b0*φ.^6 .* v1 .+ b1.*v0)./χ
    dval = C0.*Airy0/musq.^(2/3)

    # second term
    D0 =  a0*v0
    dval = dval + D0*Airy1

    # third term
    C1 = -(v3.*φ.^18 + b1*v2.*φ.^12 + b2*v1.*φ.^6 .+ b3*v0) ./ χ.^4
    dval = dval + C1.*Airy0/musq.^(2/3+2)

    # fourth term
    D1 = (a0*v2.*φ.^12 + a1*v1.*φ.^6 .+ a2*v0) ./ χ.^3
    dval = dval + D1.*Airy1/musq.^2

    dval = C.*dval

    return val, dval
end

function hermite_initialguess(n::Integer)
    # HERMITEINTITIALGUESSES(N), Initial guesses for Hermite zeros.
    #
    # [1] L. Gatteschi, Asymptotics and bounds for the zeros of Laguerre
    # polynomials: a survey, J. Comput. Appl. Math., 144 (2002), pp. 7-27.
    #
    # [2] F. G. Tricomi, Sugli zeri delle funzioni di cui si conosce una
    # rappresentazione asintotica, Ann. Mat. Pura Appl. 26 (1947), pp. 283-300.

    # Error if n < 20 because initial guesses are based on asymptotic expansions:
    @assert n ≥ 20

    # Gatteschi formula involving airy roots [1].
    # These initial guess are good near x = sqrt(n+1/2);
    if isodd(n)
        m = (n-1)>>1
        # bess = (1:m)*π
        a = .5
    else
        m = n>>1
        # bess = ((0:m-1) .+ 0.5)*π
        a = -.5
    end
    ν = 4m + 2a + 2

    T(t) = t^(2/3)*(1+5/48*t^(-2)-5/36*t^(-4)+(77125/82944)*t^(-6) -108056875/6967296*t^(-8)+162375596875/334430208*t^(-10))
    airyrts = -T.(3π/8*(4*(1:m) .- 1))

    airyrts[1:10] .= AIRY_ROOTS[1:10]  # correct first 10.

    x_init = sqrt.(abs.(ν .+ (2^(2/3)).*airyrts.*ν^(1/3) .+ (1/5*2^(4/3)).*airyrts.^2 .* ν^(-1/3) .+
        (11/35-a^2-12/175).*airyrts.^3 ./ ν .+ ((16/1575).*airyrts.+(92/7875).*airyrts.^4).*2^(2/3).*ν^(-5/3) .-
        ((15152/3031875).*airyrts.^5 .+ (1088/121275).*airyrts.^2).*2^(1/3).*ν^(-7/3)))
    x_init_airy = reverse(x_init)

    # Tricomi initial guesses. Equation (2.1) in [1]. Originally in [2].
    # These initial guesses are good near x = 0 . Note: zeros of besselj(+/-.5,x)
    # are integer and half-integer multiples of π.
    # x_init_bess =  bess/sqrt(ν).*sqrt((1+ (bess.^2+2*(a^2-1))/3/ν^2) );
    Tnk0 = fill(π/2,m)
    rhs = ((4m+3) .- 4*(1:m))/ν*π

    for k = 1:7
        val = Tnk0 .- sin.(Tnk0) .- rhs
        dval = 1 .- cos.(Tnk0)
        dTnk0 = val./dval
        Tnk0 = Tnk0 .- dTnk0
    end

    tnk = cos.(Tnk0/2).^2
    x_init_sin = sqrt.(ν*tnk - ((tnk.+1/4)./(tnk.-1).^2 .+ (3a^2-1))/3ν)

    # Patch together
    p = 0.4985+eps(Float64)
    x_init = [x_init_sin[1:convert(Int,floor(p*n))]
    x_init_airy[convert(Int,ceil(p*n)):end]]

    if isodd(n)
        x_init = [0 ; x_init]
        x_init = x_init[1:m+1]
    else
        x_init = x_init[1:m]
    end

    return x_init
end

function hermite_gw(n::Integer)
    # Golub--Welsch algorithm. Used here for n ≤ 20.

    beta = sqrt.((1:n-1)/2)  # 3-term recurrence coeffs
    T = SymTridiagonal(zeros(n), beta)  # Jacobi matrix
    D, V = eigen(T)  # Eigenvalue decomposition
    ind = sortperm(D)  # Hermite points
    x = D[ind]  # nodes
    w = sqrt(π)*V[1,ind].^2  # weights

    # Enforce symmetry:
    i = floor(Int, n/2)+1:n
    x = x[i]
    w = w[i]
    return x, exp.(x.^2).*w
end
