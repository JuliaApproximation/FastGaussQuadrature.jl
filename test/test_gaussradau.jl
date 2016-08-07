# Test for GaussRadau()

n = 1 
x, w = gaussradau(n)
@test_approx_eq x[1] -1.
@test_approx_eq w[1] 2.

n = 2 
x, w = gaussradau(n)
@test_approx_eq x[1] -1.
@test_approx_eq x[2] 1./3.
@test_approx_eq w[1] .5
@test_approx_eq w[2] 1.5

tol = 1e-14
n = 42
ntests = 4
x, w = gaussradau(n)
@test length(x) == n && length(w) == n
@test abs(x[37] - 0.908847278001044) < tol
@test abs(w[37] - 0.031190846817016) < tol
@test x[1] == -1

@test_approx_eq dot( w,exp(x)) exp(1)-exp(-1)