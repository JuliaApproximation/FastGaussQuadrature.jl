# Test for GaussRadau()
tol = 1e-14
n = 42
ntests = 4
x, w = gaussradau(n)
@test length(x) == n && length(w) == n
@test abs(x[37] - 0.908847278001044) < tol
@test abs(w[37] - 0.031190846817016) < tol
@test x[1] == -1

@test_approx_eq dot( w,exp(x)) exp(1)-exp(-1)
