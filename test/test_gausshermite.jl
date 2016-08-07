# Test for gausshermite()

tol = 1e-14

n = 18; 
x,w = gausshermite( n )
@test (length(x) == n && length(w) == n)
@test (dot(w,x) < tol && abs(dot(w,x.^2) - sqrt(pi)/2) < tol)

# Test a small n:
n = 42
x,w = gausshermite( n )
@test (length(x) == n && length(w) == n)
@test (dot(w,x) < tol && abs(dot(w,x.^2) - sqrt(pi)/2) < tol)
@test abs(x[37] - 5.660357581283058) < tol
@test abs(w[17] - 0.032202101288908) < tol

# Test a larger n:
n = 251
x,w = gausshermite( n )
@test (length(x) == n && length(w) == n)
@test (dot(w,x) < tol && abs(dot(w,x.^2) - sqrt(pi)/2) < 300*tol)
@test abs(x[37] - -13.292221459334638) < tol
@test abs(w[123] - 0.117419270715955) < 2*tol