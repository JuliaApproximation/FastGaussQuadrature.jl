function gausslobatto(n::Integer)
    # Gauss-Legendre-Lobatto Quadrature Nodes and Weights
    if n ≤ 1
        throw(DomainError(n, "Lobatto undefined for n ≤ 1."))
    elseif n == 2
        return [-1.0, 1.0], [1.0, 1.0]
    elseif n == 3
        return [-1.0, 0.0, 1.0], [1.0 / 3, 4.0 / 3, 1.0 / 3]
    else
        # Compute via GaussJacobi:
        x, w = gaussjacobi(n - 2, 1.0, 1.0)
        @inbounds for i in 1:length(x)
            w[i] = w[i] / (1 - x[i]^2)
        end
        pushfirst!(x, -1.0)
        push!(x, 1.0)
        pushfirst!(w, 2 / (n * (n - 1)))
        push!(w, 2 / (n * (n - 1)))
        return x, w
    end
end
