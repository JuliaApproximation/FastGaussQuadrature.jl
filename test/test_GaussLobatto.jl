    # Test for GaussLobatto() 
    tol = 1e-14
    n = 42
    x = GaussLobatto(n)
    w = x[2]; x = x[1]
    @test ( length(x) == n) && ( length(w) == n )   
    @test abs(x[37] - 0.922259214258616) < tol
    @test abs(w[37] - 0.029306411216166) < tol
    @test ( x[1] == -1 && x[n] == 1 )
