# Test for gausschebyshev().

n = 10 

# Test x.^2:
x, w = gausschebyshev(n,1)
@test_approx_eq dot(x.^2,w) pi/2
x, w = gausschebyshev(n,2)
@test_approx_eq dot(x.^2,w) pi/8
x, w = gausschebyshev(n,3)
@test_approx_eq dot(x.^2,w) pi/2
x, w = gausschebyshev(n,4)
@test_approx_eq dot(x.^2,w) pi/2

# Test x^3:
x, w = gausschebyshev(n,1)
@test abs(dot(x.^3,w))<1e-15
x, w = gausschebyshev(n,2)
@test abs(dot(x.^3,w))<1e-15
x, w = gausschebyshev(n,3)
@test_approx_eq dot(x.^3,w) 3*pi/8
x, w = gausschebyshev(n,4)
@test_approx_eq dot(x.^3,w) -3*pi/8