# Test gausslaguerre.jl 

tol = 3.e-13

# Test a small n (using GQ)
n = 42
x, w = gausslaguerre( n )
@test abs(dot(x,w) - 1) <= tol
@test abs(dot(x.^2,w) - 2) <= tol
@test abs(x[37] - 98.388267163326702) < tol
@test abs(w[7] - 0.055372813167092) < tol

# Test a larger n (using GLR)
n = 251
x, w = gausslaguerre( n )

@test abs(dot(x,w) - 1) < 100*tol
@test abs(dot(x.^2,w) - 2) < 100*tol
@test abs(x[37] - 13.309000189442097) < tol
@test abs(w[3] - 0.050091759039996) < tol

# Test an even larger n (using RH)
n = 5000
x, w = gausslaguerre( n )

@test abs(dot(x,w) - 1) < 100*tol
@test abs(dot(x.^2,w) - 2) < 100*tol
