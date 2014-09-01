
n testGaussLegendre()
    # Test for GaussLegendre quadrature. 
    
    tol = 1e-14
    ntests = 15
    pass = Array(Float64, ntests)
    # Test a small n (using REC)
    n = 42
    x = GaussLegendre(n);
    pass[1] = all(length(x[1]) == n);
    x = GaussLegendre(n)
    pass[2] = all(length(x[1]) == n) && all(length(x[2]) == n);
    pass[3] = abs(x[1][37] - 0.910959724904127) < tol
    pass[4] = abs(x[2][37] - 0.030479240699603) < tol

    # Test a larger n (using ASY)
    n = 251
    x = GaussLegendre(n)
    w = x[2]; x = x[1]
    pass[5] = all(length(x) == n) && all(length(w) == n)
    pass[6] = abs(x[37] + 0.896467746955729) < tol
    pass[7] = abs(w[37] - 0.005535005742012) < tol

    x = GaussLegendre(1013); 
    pass[8] = norm( x[1][2] - -0.999985167586110, Inf ) < tol
    pass[9] = norm( x[1][13] - -0.999218995240887, Inf ) < tol
    pass[10] = norm( x[2][2] - 1.681691163200592e-05, Inf ) < tol
    pass[11] = norm( x[2][13] - 1.224755309137936e-04, Inf ) < tol
    
    x = GaussLegendre(10013); 
    pass[12] = norm( x[1][2] - -0.999999848054223, Inf ) < tol
    pass[13] = norm( x[1][13] - -0.999991998242661, Inf ) < tol
    pass[14] = norm( x[2][2] - 1.722757320118474e-07, Inf ) < tol
    pass[15] = norm( x[2][13] - 1.254980540032470e-06, Inf ) < tol
    
    return ( sum(pass) == ntests )
end
