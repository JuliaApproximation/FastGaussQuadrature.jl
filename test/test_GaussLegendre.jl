# Test for GaussLegendre(). 
    
tol = 1e-14
    pass = Array(Float64, ntests)
    n = 42
    x = GaussLegendre(n);
    pass[1] = all(length(x[1]) == n);
    x = GaussLegendre(n)
    @test length(x[1]) == n && length(x[2]) == n
    @test abs(x[1][37] - 0.910959724904127) < tol
    @test abs(x[2][37] - 0.030479240699603) < tol

    # Test a larger n (using ASY)
    n = 251
    x = GaussLegendre(n)
    w = x[2]; x = x[1]
    @test all(length(x) == n) && all(length(w) == n)
    @test abs(x[37] + 0.896467746955729) < tol
    @test abs(w[37] - 0.005535005742012) < tol

    x = GaussLegendre(1013); 
    @test norm( x[1][2] - -0.999985167586110, Inf ) < tol
    @test norm( x[1][13] - -0.999218995240887, Inf ) < tol
    @test norm( x[2][2] - 1.681691163200592e-05, Inf ) < tol
    @test norm( x[2][13] - 1.224755309137936e-04, Inf ) < tol
    
    x = GaussLegendre(10013); 
    @test norm( x[1][2] - -0.999999848054223, Inf ) < tol
    @test norm( x[1][13] - -0.999991998242661, Inf ) < tol
    @test norm( x[2][2] - 1.722757320118474e-07, Inf ) < tol
    @test norm( x[2][13] - 1.254980540032470e-06, Inf ) < tol
