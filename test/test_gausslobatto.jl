# Test for GaussLobatto()
n = 2
x,w = gausslobatto(n)
@test_approx_eq x[1] -1.
@test_approx_eq x[2] 1.
@test_approx_eq w[1] 1.
@test_approx_eq w[2] 1.

n = 3
x,w = gausslobatto(n)
@test_approx_eq x[1] -1.
@test abs(x[2])<1e-15
@test_approx_eq x[3] 1.
@test_approx_eq w[1] 1./3
@test_approx_eq w[2] 4./3
@test_approx_eq w[3] 1./3

tol = 1e-14
n = 42
x,w = gausslobatto(n)
@test ( length(x) == n) && ( length(w) == n )
@test abs(x[37] - 0.922259214258616) < tol
@test abs(w[37] - 0.029306411216166) < tol
@test ( x[1] == -1 && x[n] == 1 )

@test_approx_eq dot( w,(x.^2)) 2/3
@test_approx_eq dot( w,exp(x)) exp(1)-exp(-1)
