# Test for GaussJacobi(). 
    tol = 1e-14
    # Test a small n:
    n = 42; a = -.1; b = .3
    x = GaussJacobi(n, -.1, .3)
    w = x[2]; x = x[1]
    @test all(length(x) == n) && all(length(w) == n)
    @test abs(x[37] - 0.912883347814032) < tol
    @test abs(w[37] - 0.046661910947553) < tol

    # Test a larger n (using ASY)
    a = -.7; b = 1.3;
    n = 251;
    x = GaussJacobi(n,-.7,1.3);
    w = x[2]; x = x[1]
    @test all(length(x) == n) && all(length(w) == n) 
    @test abs(x[37] + 0.893103435898983) < tol
    @test abs(w[37] - 1.962534523788093e-04) < tol

    # Test n = 1: 
    a = 1.0; b = 2.0; 
    x = GaussJacobi( 1, a, b ); 
    @test abs( x[1][1] - (b-a)/(a+b+2) ) < tol 
    @test abs( x[2][1] - 2^(a+b+1)*beta(a+1, b+1) ) < tol
    
    x = GaussJacobi( 1013, .9, -.1 )
    @test abs( x[1][2]  + 0.999986012231899 ) < tol 
    @test abs( x[1][13]  + 0.999225722939832 ) < tol 
    @test abs( x[2][2] - 9.314674169892358e-05 ) < tol 
    @test abs( x[2][13] - 4.654651764553262e-04 ) < tol 

    x = GaussJacobi( 10013, .9, -.1 )
    @test abs( x[1][2]  +0.999999856605293 ) < tol 
    @test abs( x[1][13]  + 0.999992061552711 ) < tol 
    @test abs( x[2][2] - 1.509654630405615e-06 ) < tol 
    @test abs( x[2][13] - 7.548275262993863e-06 ) < tol 
