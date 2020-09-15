function gaussradau(n::Integer, a=0.0, b=0.0)
    T = float(promote_type(eltype(a), eltype(b)))
    # RADAUPTS   Gauss-Legendre-Radau Quadrature Nodes and Weights
    if n == 1
        T[-1], T[2]
    elseif n == 2
        [-1, one(T)/3], [one(T)/2, convert(T,3)/2]
    else
        # Compute via GaussJacobi:
        x, w = gaussjacobi(n - 1, a, b+1)
        @inbounds for i in 1:length(w)
            w[i] = w[i] / (1 + x[i])
        end
        pushfirst!(x, -1)
        pushfirst!(w, convert(T, 2) / n^2)
        x, w
    end
end