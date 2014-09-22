function gaussradau( n::Int64 )
#RADAUPTS   Gauss-Legendre-Radau Quadrature Nodes and Weights

    if ( n == 1 )
        x = ([-1.0], [2.0])
    elseif ( n == 2 )
        x = ([-1.0, 1/3], [.5, 1.5])
    else
    # Compute via GaussJacobi:
        x = gaussjacobi( n - 1, 0.0, 1.0 )
        w = x[2]; x = x[1];
        x = [-1.0, x]
        w = [2.0/n^2, w./(1.0 + x[2:end])]
        x = (x, w)
    end
    return x
end
