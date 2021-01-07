@doc raw"""
    gaussjacobi(n::Integer) -> Tuple{Vector{Float64},Vector{Float64}}

Return nodes and weights of [Gauss-Jacobi quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Jacobi_quadrature).

```math
\int_{-1}^{1} f(x) (1-x)^\alpha(1+x)^\beta dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

# Examples
```jldoctest; setup = :(using FastGaussQuadrature, LinearAlgebra)
julia> x, w = gaussjacobi(3, 1/3, -1/3);

julia> f(x) = x^4;

julia> I = dot(w, f.(x));

julia> I ≈ 268π/729(√3)
true
```
"""
function gaussjacobi(n::Integer, a::Real, b::Real)
    #GAUSS-JACOBI QUADRATURE NODES AND WEIGHTS
    if n < 0
        throw(DomainError(n, "gaussjacobi($n,$a,$b) not defined: n must be non-negative."))
    elseif a == 0. && b == 0.
        gausslegendre(n)
    elseif a == -0.5 && b == -0.5
        gausschebyshev(n, 1)
    elseif a == 0.5 && b == 0.5
        gausschebyshev(n, 2)
    elseif a == -0.5 && b == 0.5
        gausschebyshev(n, 3)
    elseif a == 0.5 && b == -0.5
        gausschebyshev(n, 4)
    elseif n == 0
        Float64[], Float64[]
    elseif n == 1
        [(b - a) / (a + b + 2)], [2^(a + b + 1) * beta(a + 1, b + 1)]
    elseif min(a,b) <= -1.
        throw(DomainError((a,b), "The Jacobi parameters correspond to a nonintegrable weight function"))
    elseif n <= 100 && max(a,b) < 5.
        jacobi_rec(n, a, b)
    elseif n > 100 && max(a,b) < 5.
        jacobi_asy(n, a, b)
    elseif n <= 4000 && max(a,b)>=5.
        jacobi_gw(n, a, b)
    else
        error("gaussjacobi($n,$a,$b) is not implemented: n must be ≤ 4000 for max(a,b)≥5.")
    end
end

# Convenience function: convert any kind of numbers a and b to a joint floating point type
jacobi_rec(n::Integer, a::Real, b::Real) = jacobi_rec(n, promote(float(a), float(b))...)

function jacobi_rec(n::Integer, a::T, b::T) where {T <: AbstractFloat}
    #Compute nodes and weights using recurrrence relation.
    x11, x12 = HalfRec(n, a, b, 1)
    x21, x22 = HalfRec(n, b, a, 0)

    x = Array{T}(undef,n)
    w = Array{T}(undef,n)
    m1 = length(x11)
    m2 = length(x21)
    sum_w = zero(T)
    @inbounds for i in 1:m2
        idx = m2 + 1 - i
        xi = -x21[i]
        der = x22[i]
        wi = 1 / ((1 - xi^2) * der^2)
        w[idx] = wi
        x[idx] = xi
        sum_w += wi
    end
    @inbounds for i in 1:m1
        idx = m2 + i
        xi = x11[i]
        der = x12[i]
        wi = 1 / ((1 - xi^2) * der^2)
        w[idx] = wi
        x[idx] = xi
        sum_w += wi
    end
    c = (2^(a + b + 1) * gamma(2 + a) *
         gamma(2 + b) / (gamma(2 + a + b) * (a + 1) * (b + 1)))
    rmul!(w, c / sum_w)
    return x, w
end

function HalfRec(n::Integer, a::T, b::T, flag) where {T <: AbstractFloat}
    # HALFREC  Jacobi polynomial recurrence relation.
    # Asymptotic formula - only valid for positive x.
    r = (flag == 1) ? (ceil(n / 2):-1:1) : (floor(n / 2):-1:1)
    m = length(r)
    c1 = 1 / (2 * n + a + b + 1)
    a1 = one(T)/4 - a^2
    b1 = one(T)/4 - b^2
    c1² = c1^2
    x = Array{T}(undef,m)
    @inbounds for i in 1:m
        C = muladd(2*one(T), r[i], a - one(T)/2) * (T(π) * c1)
        C_2 = C/2
        x[i] = cos(muladd(c1², muladd(-b1, tan(C_2), a1 * cot(C_2)), C))
    end

    P1 = Array{T}(undef,m)
    P2 = Array{T}(undef,m)
    # Loop until convergence:
    for _ in 1:10
        innerjacobi_rec!(n, x, a, b, P1, P2)
        dx2 = 0.0
        @inbounds for i in 1:m
            dx = P1[i] / P2[i]
            _dx2 = abs2(dx)
            dx2 = ifelse(_dx2 > dx2, _dx2, dx2)
            x[i] = x[i] - dx
        end
        dx2 > eps(T) / 1e6 || break
    end
    # Once more for derivatives:
    innerjacobi_rec!(n, x, a, b, P1, P2)
    return x, P2
end

function innerjacobi_rec!(n, x, a::T, b::T, P, PP) where {T <: AbstractFloat}
    # EVALUATE JACOBI POLYNOMIALS AND ITS DERIVATIVE USING THREE-TERM RECURRENCE.
    N = length(x)
    @inbounds for j = 1:N
        xj = x[j]

        Pj = (a - b + (a + b + 2) * xj)/2
        Pm1 = one(T)
        PPj = (a + b + 2)/2
        PPm1 = zero(T)
        for k in 1:n-1
            k0 = muladd(2*one(T), k, a + b)
            k1 = k0 + 1
            k2 = k0 + 2
            A = 2 * (k + 1) * (k + (a + b + 1)) * k0
            B = k1 * (a^2 - b^2)
            C = k0 * k1 * k2
            D = 2 * (k + a) * (k + b) * k2
            c1 = muladd(C, xj, B)
            Pm1, Pj = Pj, muladd(-D, Pm1, c1 * Pj) / A
            PPm1, PPj = PPj, muladd(c1, PPj, muladd(-D, PPm1, C * Pm1)) / A
        end
        P[j] = Pj
        PP[j] = PPj
    end
    nothing
end

function innerjacobi_rec(n, x, a::T, b::T) where {T <: AbstractFloat}
    # EVALUATE JACOBI POLYNOMIALS AND ITS DERIVATIVE USING THREE-TERM RECURRENCE.
    N = length(x)
    P = Array{T}(undef,N)
    PP = Array{T}(undef,N)
    innerjacobi_rec!(n, x, a, b, P, PP)
    return P, PP
end

function weightsConstant(n, a, b)
    # Compute the constant for weights:
    M = min(20, n - 1)
    C = 1.0
    p = -a * b / n
    for m = 1:M
        C += p
        p *= -(m + a) * (m + b) / (m + 1) / (n - m)
        abs(p / C) < eps(Float64) / 100 && break
    end
    return 2^(a + b + 1) * C
end

function jacobi_asy(n, a, b)
    # ASY  Compute nodes and weights using asymptotic formulae.

    # Determine switch between interior and boundary regions:
    nbdy = 10
    bdyidx1 = n - (nbdy - 1):n
    bdyidx2 = nbdy:-1:1

    # Interior formula:
    x, w = asy1(n, a, b, nbdy)

    # Boundary formula (right):
    xbdy = boundary(n, a, b, nbdy)
    x[bdyidx1], w[bdyidx1] = xbdy

    # Boundary formula (left):
    if a != b
        xbdy = boundary(n, b, a, nbdy)
    end
    x[bdyidx2] = -xbdy[1]
    w[bdyidx2] = xbdy[2]

    rmul!(w, weightsConstant(n, a, b))
    return x, w
end

function asy1(n::Integer, a::Float64, b::Float64, nbdy::Integer)
    # Algorithm for computing nodes and weights in the interior.

    # Approximate roots via asymptotic formula: (Gatteschi and Pittaluga, 1985)
    K = (2(n:-1:1).+a.-0.5).*pi./(2n+a+b+1)
    tt = K .+ (1/(2n+a+b+1)^2).*((0.25-a^2).*cot.(K/2).-(0.25-b^2).*tan.(K/2))

    # First half (x > 0):
    t = tt[tt .<= pi/2]
    mint = t[end-nbdy+1]
    idx = 1:max(findfirst(t .< mint)-1, 1)

    # Newton iteration
    for _ in 1:10
        vals, ders = feval_asy1(n, a, b, t, idx)         # Evaluate
        dt = vals./ders
        t += dt                                     # Next iterate
        if norm(dt[idx],Inf) > sqrt(eps(Float64))/100
            break
        end
    end
    vals, ders = feval_asy1(n, a, b, t, idx)      # Once more for luck
    t += vals./ders

    # Store
    x_right = cos.(t)
    w_right = 1 ./ ders.^2

    # Second half (x < 0):
    a, b = b, a
    t = pi .- tt[1:(n-length(x_right))]
    mint = t[nbdy]
    idx = max(findfirst(t .> mint), 1):length(t)

    # Newton iteration
    for _ in 1:10
        vals, ders = feval_asy1(n, a, b, t, idx)  # Evaluate.
        dt = vals./ders                                # Newton update.
        t += dt
        if norm(dt[idx],Inf) > sqrt(eps(Float64))/100
            break
        end
    end
    vals, ders = feval_asy1(n, a, b, t, idx)     # Once more for luck.
    t += vals./ders                                 # Newton update.

    # Store
    x_left = cos.(t)
    w_left = 1 ./ ders.^2

    return vcat(-x_left, x_right), vcat(w_left, w_right)
end

"""
Evaluate the interior asymptotic formula at x = cos(t).
Assumption:
* length(t) == n ÷ 2
"""
function feval_asy1(n::Integer, a::Float64, b::Float64, t::AbstractVector, idx)
    # Number of terms in the expansion:
    M = 20

    # Number of elements in t:
    N = length(t)

    # Some often used vectors/matrices:
    onesM = ones(M)
    MM = collect(1:M)

    # The sine and cosine terms:
    alpha = repeat((2n+a+b).+MM,1,N).*repeat(t',M)/2 .- (a+0.5)*pi/2  # M × N matrix
    cosA = cos.(alpha)
    sinA = sin.(alpha)

    sinT = repeat(sin.(t)',M)
    cosT = repeat(cos.(t)',M)
    cosA2 = cosA.*cosT .+ sinA.*sinT
    sinA2 = sinA.*cosT .- cosA.*sinT

    sinT = hcat(ones(N), cumprod(repeat((csc.(t/2)/2),1,M-1),2))  # M × N matrix
    cosT = sec.(t/2)/2

    j = 0:M-2
    _vec = @. (0.5+a+j)*(0.5-a+j)/(j+1)/(2n+a+b+j+2)
    P1 = [1;cumprod(_vec,1)]
    P1[3:4:end] = -P1[3:4:end]
    P1[4:4:end] = -P1[4:4:end]
    P2 = Matrix(1.0I, M, M)
    for l = 0:M-1
        j = 0:M-l-2
        _vec = @. (0.5+b+j)*(0.5-b+j)/(j+1)/(2n+a+b+j+l+2)
        P2[l+1,(l+1).+(1:length(j))] = cumprod(_vec,1)
    end
    PHI = repeat(P1,1,M).*P2

    j = 0:M-2
    _vec = @. (0.5+a+j)*(0.5-a+j)/(j+1)/(2*(n-1)+a+b+j+2)
    P1 = [1;cumprod(_vec,1)]
    P1[3:4:end] = -P1[3:4:end]
    P1[4:4:end] = -P1[4:4:end]
    P2 = Matrix(1.0I, M, M)
    for l = 0:M-1
        j = 0:M-l-2
        _vec = @. (0.5+b+j)*(0.5-b+j)/(j+1)/(2*(n-1)+a+b+j+l+2)
        P2[l+1,(l+1).+(1:length(j))] = cumprod(_vec,1)
    end
    PHI2 = repeat(P1,1,M).*P2

    S = S2 = zeros(N)
    SC = sinT
    for m = 1:M
        l = 1:2:m
        phi = PHI[l, m]
        dS1 = (SC[:, l]*phi) .* cosA[m, :]
        phi2 = PHI2[l, m]
        dS12 = (SC[:, l]*phi2) .* cosA2[m, :]
        l = 2:2:m
        phi = PHI[l, m]
        dS2 = (SC[:, l]*phi) .* sinA[m, :]
        phi2 = PHI2[l, m]
        dS22 = (SC[:, l]*phi2) .* sinA2[m, :]
        if m - 1 > 10 && norm(dS1[idx] + dS2[idx], Inf) < eps(Float64) / 100
            break
        end
        S = S .+ dS1 .+ dS2
        S2 = S2 .+ dS12 .+ dS22
        SC[:,1:m] = SC[:,1:m].*cosT
    end

    # Constant out the front:
    dsa = a^2/2n
    dsb = b^2/2n
    dsab = (a+b)^2/4n
    ds = dsa + dsb - dsab
    s = ds
    j = 1
    dsold = ds # to fix a = -b bug.
    while (abs(ds/s) + dsold) > eps(Float64)/10
        dsold = abs(ds/s)
        j += 1
        tmp = -(j-1)/(j+1)/n
        dsa = tmp*dsa*a
        dsb = tmp*dsb*b
        dsab = tmp*dsab*(a+b)/2
        ds = dsa + dsb - dsab
        s = s + ds
    end
    p2 = exp(s)*sqrt(2*pi)*sqrt((n+a)*(n+b)/(2n+a+b))/(2n+a+b+1)
    # g is a coefficients in
    # ``\Gamma(z) = \frac{z^{z-1/2}}{\exp(z)}\sqrt{2\pi} \left(\sum_{i} B_i z^{-i}\right)``, where B_{i-1} = g[i].
    # https://math.stackexchange.com/questions/1714423/what-is-the-pattern-of-the-stirling-series
    g = [1, 1/12, 1/288, -139/51840, -571/2488320, 163879/209018880,
         5246819/75246796800, -534703531/902961561600,
         -4483131259/86684309913600, 432261921612371/514904800886784000]
    f(g,z) = sum(g.*[1;cumprod(ones(9)./z)])
    C = p2*(f(g,n+a)*f(g,n+b)/f(g,2n+a+b))*2/pi
    C2 = C*(a+b+2n).*(a+b+1+2n)./(4*(a+n).*(b+n))

    _vals = vec(C*S')
    _S2 = vec(C2*S2')

    # Use relation for derivative:
    ders = (n.*(a.-b.-(2n+a+b).*cos.(t)).*_vals .+ 2 .*(n+a).*(n+b).*_S2)./(2n+a+b)./sin.(t)
    t_abs = abs.(t)
    denom = 1 ./ (sin.(t_abs/2).^(a+0.5).*cos.(t_abs/2).^(b+0.5))
    vals = _vals.*denom
    ders = ders.*denom

    return vals, ders
end

function boundary(n::Integer, a::Float64, b::Float64, npts)
# Algorithm for computing nodes and weights near the boundary.

    # Use Newton iterations to find the first few Bessel roots:
    smallK = min(30, npts)
    jk = besselroots(a, min(npts, smallK))

    # Approximate roots via asymptotic formula: (see Olver 1974)
    phik = jk/(n + .5*(a + b + 1))
    x = cos.( phik .+ ((a^2-0.25).*(1 .-phik.*cot.(phik))./(8*phik) .- 0.25.*(a^2-b^2).*tan.(0.5.*phik))./(n + 0.5*(a + b + 1))^2 )

    # Newton iteration:
    for _ in 1:10
        vals, ders = innerjacobi_rec(n, x, a, b)   # Evaluate via asymptotic formula.
        dx = -vals./ders                   # Newton update.
        x += dx                             # Next iterate.
        if norm(dx,Inf) > sqrt(eps(Float64))/200
            break
        end
    end
    vals, ders = innerjacobi_rec(n, x, a, b);     # Evaluate via asymptotic formula.
    dx = -vals./ders
    x += dx

    # flip:
    x = x[npts:-1:1]
    ders = ders[npts:-1:1]

    # Revert to x-space:
    w = @. 1 /((1-x^2)*ders^2)
    return x, w
end

function jacobi_jacobimatrix(n, a, b)
    ab = a + b
    ii = 2:n-1
    abi = 2*ii .+ ab
    aa = [(b - a)/(2 + ab);
          (b^2 - a^2) ./ ((abi .- 2).*abi);
          (b^2 - a^2) ./ ((2n - 2+ab).*(2n+ab))]
    bb = [2*sqrt( (1 + a)*(1 + b)/(ab + 3))/(ab + 2) ;
          2 .*sqrt.(ii.*(ii .+ a).*(ii .+ b).*(ii .+ ab)./(abi.^2 .- 1))./abi]
    return SymTridiagonal(aa, bb)
end


function jacobimoment(a,b)
    ab = a + b
    T = float(typeof(ab))
    # Same as 2^(a+b+1) * beta(a+1,b+1)
    return exp((ab+1)*log(convert(T,2)) + loggamma(a+1)+loggamma(b+1)-loggamma(2+ab))
end

function jacobi_gw(n::Integer, a, b)
    # Golub-Welsh for Gauss--Jacobi quadrature. This is used when max(a,b)>5.
    ab = a + b
    x, V = eigen(jacobi_jacobimatrix(n,a,b))  # Eigenvalue decomposition.
    # Quadrature weights:
    w = V[1,:].^2 .* jacobimoment(a,b)
    return x, w
end

