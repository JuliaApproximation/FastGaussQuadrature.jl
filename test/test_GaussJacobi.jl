function testGaussJacobi()

    tol = 1e-14
    ntests = 16
    pass = Array(Float64,ntests)
    # Test a small n:
    n = 42; a = -.1; b = .3
    x = GaussJacobi(n, -.1, .3)
    w = x[2]; x = x[1]
    pass[1] = all(length(x) == n) && all(length(w) == n)
    pass[2] = abs(x[37] - 0.912883347814032) < tol;
    pass[3] = abs(w[37] - 0.046661910947553) < tol;

    # Test a larger n (using ASY)
    a = -.7; b = 1.3;
    n = 251;
    x = GaussJacobi(n,-.7,1.3);
    w = x[2]; x = x[1]
    pass[4]= all(length(x) == n) && all(length(w) == n) 
    pass[5] = abs(x[37] + 0.893103435898983) < tol
    pass[6] = abs(w[37] - 1.962534523788093e-04) < tol

    # Test n = 1: 
    a = 1.0; b = 2.0; 
    x = GaussJacobi( 1, a, b ); 
    pass[7] = abs( x[1] - (b-a)/(a+b+2) ) < tol 
    pass[8] = abs( x[2] - 2^(a+b+1)*beta(a+1, b+1) ) < tol
    
    x = GaussJacobi( 1013, .9, -.1 )
    pass[9] = abs( x[1][2]  + 0.999986012231899 ) < tol 
    pass[10] = abs( x[1][13]  + 0.999225722939832 ) < tol 
    pass[11] = abs( x[2][2] - 9.314674169892358e-05 ) < tol 
    pass[12] = abs( x[2][13] - 4.654651764553262e-04 ) < tol 

    x = GaussJacobi( 10013, .9, -.1 )
    pass[13] = abs( x[1][2]  +0.999999856605293 ) < tol 
    pass[14] = abs( x[1][13]  + 0.999992061552711 ) < tol 
    pass[15] = abs( x[2][2] - 1.509654630405615e-06 ) < tol 
    pass[16] = abs( x[2][13] - 7.548275262993863e-06 ) < tol 

    return ( sum( pass ) == ntests )
end
