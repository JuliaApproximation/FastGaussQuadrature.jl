# Test gausslaguerre.jl

# Test integration
tol = 4.0e-10

# Evaluate the exact integral of x^alpha * p(x) *exp(-x) on the positive halfline,
# where p(x) = sum( coef[i+1]*x^i, i=0..degree(p)) is a polynomial given by its
# monomial coefficients
function exact_integral_poly(alpha, coef)
    # do calculations in BigFloat precision
    z = big(0.0)
    for i in 1:length(coef)
        z += gamma(big(alpha)+big(i)) * big(coef[i])
    end
    Float64(z)
end

function polyval(x, coef)
    z = zeros(size(x))
    d = length(coef)-1
    for i in 1:length(x)
        z[i] = sum(coef .* x[i].^(0:d))
    end
    z
end

# Evaluate the exact integral of cos(a*x) * x^alpha * exp(-x) on the positive halfline
# - For alpha = 0
exact_integral_cos1(a) = 1/(a^2+1)
# - For alpha = 1/2
exact_integral_cos2(a) = sqrt(pi)*cos(3/2*atan(a))/2 / (a^2+1)^(3/4)

if VERSION < v"0.7"
    srand(0)
else
    Random.seed!(0)
end

##########
# Test the special cases
##########

x, w = gausslaguerre(0)
@test x == Float64[]
@test w == Float64[]

x, w = gausslaguerre(1, 0.4)
@test x == [1.4]
@test w == [gamma(1.4)]

x, w = gausslaguerre(2, -0.37)
@test x ≈ [0.35328546651962944; 2.90671453348037]
@test w ≈ [1.269857265034167; 0.15433992954023468]


##########
n = 42
##########
alpha = 0.0
coef = rand(17)
Z = exact_integral_poly(alpha, coef)

x, w = gausslaguerre(n, alpha)
@test isa(x, Vector{Float64})
@test isa(w, Vector{Float64})
@test abs(sum(w) - 1) < tol
@test abs(dot(w, x) - 1) < tol
@test abs(x[37] - 98.388267163326702) < tol
@test abs(w[7] - 0.055372813167092) < tol
@test abs(dot(w, polyval(x, coef)) - Z)/abs(Z) < tol

x_gw, w_gw = FastGaussQuadrature.gausslaguerre_GW(n, alpha)
@test abs(x[37] - 98.388267163326702) < tol
@test abs(w[7] - 0.055372813167092) < tol
@test abs(dot(w, polyval(x, coef)) - Z)/abs(Z) < tol

x_rec, w_rec = FastGaussQuadrature.gausslaguerre_rec(n, alpha)
@test abs(x[37] - 98.388267163326702) < tol
@test abs(w[7] - 0.055372813167092) < tol
@test abs(dot(w, polyval(x, coef)) - Z)/abs(Z) < tol


alpha = 0.5
coef = rand(17)
Z = exact_integral_poly(alpha, coef)
x, w = gausslaguerre(n, alpha)
@test abs(dot(w, polyval(x, coef)) - Z)/abs(Z) < tol



##########
n = 251
##########
alpha = 0.0
x, w = gausslaguerre(n, alpha)
a = 4
Z = exact_integral_cos1(a)
@test isa(x, Vector{Float64})
@test isa(w, Vector{Float64})
@test abs(x[37] - 13.309000189442097) < tol
@test abs(w[3] - 0.050091759039996) < tol
@test abs(dot(w, cos.(a*x)) - Z) < tol

alpha = 0.5
x, w = gausslaguerre(n, alpha)
@test abs(x[44] -19.095577327730616 ) < tol
@test abs(w[18] - 0.026245779174690266) < tol
a = 4
Z = exact_integral_cos2(a)
@test abs(dot(w, cos.(a*x)) - Z) < tol


############
n = 350000
############
alpha = 0.0
a = 50
Z = exact_integral_cos1(a)
x, w = gausslaguerre(n, alpha)
@test isa(x, Vector{Float64})
@test isa(w, Vector{Float64})
@test abs(dot(w, cos.(a*x)) - Z) < tol

alpha = 0.44
x, w = gausslaguerre(n, alpha; reduced = true)
@test w[end] < 100*floatmin(Float64)
