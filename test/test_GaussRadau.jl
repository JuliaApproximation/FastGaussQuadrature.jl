function testGaussRadau()
    # Test for GaussRadau() 
    tol = 1e-14
    n = 42
    ntests = 4
    pass = Array(Float64,ntests)
    x = GaussRadau(n)
    w = x[2]; x = x[1]
    pass[1] = ( length(x) == n) && ( length(w) == n )   
    pass[2] = abs(x[37] - 0.908847278001044) < tol
    pass[3] = abs(w[37] - 0.031190846817016) < tol
    pass[4] = ( x[1] == -1 )
    return ( sum( pass ) == ntests )
end
