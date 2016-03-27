# Test for gaussjacobi().
tol = 1e-14

# Check if besselj(nu, besselroots(nu, n) ) is small: 

for nu = 0.:0.1:10 
n = 10; 
@test norm( besselj(nu, besselroots(nu, n) ) ) < tol; 

end