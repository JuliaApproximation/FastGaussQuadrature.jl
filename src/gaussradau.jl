function gaussradau(n::Int)
    # RADAUPTS   Gauss-Legendre-Radau Quadrature Nodes and Weights
    if n == 1
        [-1.0], [2.0]
    elseif n == 2
        [-1.0, 1/3], [.5, 1.5]
    else
        # Compute via GaussJacobi:
        x, w = gaussjacobi(n - 1, 0.0, 1.0)
        @inbounds for i in 1:length(w)
            w[i] = w[i] / (1.0 + x[i])
        end
        unshift!(x, -1.0)
        unshift!(w, 2.0 / n^2)
        x, w
    end
end
