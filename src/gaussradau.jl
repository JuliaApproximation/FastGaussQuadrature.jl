"""
   gaussradau(n)

Creates the n-point Gauss-Radau quadrature rule, with the first node fixed at -1.
"""
function gaussradau(n::Integer, T::Type=Float64)
    a = b = zero(T)
    # RADAUPTS   Gauss-Legendre-Radau Quadrature Nodes and Weights
    if n == 1
        return T[-1], T[2]
    elseif n == 2
        return [-1, one(T)/3], [one(T)/2, convert(T,3)/2]
    else
        # Compute via GaussJacobi:
        x, w = gaussjacobi(n - 1, a, b+1)
        @inbounds for i in 1:length(w)
            w[i] = w[i] / (1 + x[i])
        end
        pushfirst!(x, -1)
        pushfirst!(w, convert(T, 2) / n^2)
        return x, w
    end
end

function gaussradau(n::Integer, a, b)
    if n ≤ 0
        throw(DomainError(n, "Input N must be a positive integer"))
    end
    m = n - 1
    ab = a + b
    T = float(eltype(ab))
    μ = jacobimoment(a, b)
    n == 0 && return T[], T[]
    n == 1 && return [-one(T)], [μ]
    J = jacobi_jacobimatrix(n, a, b)
    aᴿ = -1 + 2m*convert(T,m+a)/((2m+ab)*(2m+ab+1))
    J.dv[end] = aᴿ
    x, V = eigen(J)
    w = V[1,:].^2 .* μ
    x[1] = -1 # fix rounding from eigen computation
    return x, w
end
