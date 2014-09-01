    # Test for GaussRadau() 
    tol = 1e-14
    n = 42
    ntests = 4
    x = GaussRadau(n)
    w = x[2]; x = x[1]
    @test ( length(x) == n) && ( length(w) == n )   
    @test abs(x[37] - 0.908847278001044) < tol
    @test abs(w[37] - 0.031190846817016) < tol
    @test ( x[1] == -1 )
