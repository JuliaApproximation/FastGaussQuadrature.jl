function test_GaussLobatto()
    # Test for GaussLobatto() 
    tol = 1e-14
    n = 42
    ntests = 4
    pass = Array(Float64,ntests)
    x = GaussLobatto(n)
    w = x[2]; x = x[1]
    pass[1] = ( length(x) == n) && ( length(w) == n )   
    pass[2] = abs(x[37] - 0.922259214258616) < tol
    pass[3] = abs(w[37] - 0.029306411216166) < tol
    pass[4] = ( x[1] == -1 && x[n] == 1 )
    return ( sum( pass ) == ntests )
end
