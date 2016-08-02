function gausshermite( n::Int )
# GAUSSHERMITE(n) COMPUTE THE GAUSS-HERMITE NODES AND WEIGHTS IN O(n) time.


if n < 0
    x = (Float64[],Float64[])
    return x
elseif n == 0
    x = (Float64[],Float64[])
    return x
elseif n == 1
    x = ([0.0],[sqrt(pi)])
    return x
elseif n <= 20
   # GW algorithm
   x = hermpts_gw( n )
elseif n <= 200
   # REC algorithm
   x = hermpts_rec( n )
else
   # ASY algorithm
   x = hermpts_asy( n )
end

if mod(n,2) == 1                              # fold out
    w = [flipdim(x[2][:],1); x[2][2:end]]
    w = (sqrt(pi)/sum(w))*w
    x = ([-flipdim(x[1],1) ; x[1][2:end]], w)
else
    w = [flipdim(x[2][:],1); x[2][:]]
    w = (sqrt(pi)/sum(w))*w
    x = ([-flipdim(x[1],1) ; x[1]], w)
end

end

function hermpts_asy( n::Int )
# Compute Hermite nodes and weights using asymptotic formula

x0 = HermiteInitialGuesses( n ) # get initial guesses
t0 = x0./sqrt(2n+1)
theta0 = acos(t0)               # convert to theta-variable
val = x0;
for k = 1:20
    val = hermpoly_asy_airy(n, theta0);
    dt = -val[1]./(sqrt(2)*sqrt(2n+1)*val[2].*sin(theta0))
    theta0 = theta0 - dt;                        # Newton update
    if norm(dt,Inf) < sqrt(eps(Float64))/10
       break
    end
end
t0 = cos(theta0)
x = sqrt(2n+1)*t0                          #back to x-variable
ders = x.*val[1] + sqrt(2)*val[2]
w = (exp(-x.^2)./ders.^2)';            # quadrature weights

x = (x, w)
end

function hermpts_rec( n::Int )
# Compute Hermite nodes and weights using recurrence relation.

x0 = HermiteInitialGuesses( n )
x0 = x0.*sqrt(2)
val = x0
for kk = 1:10
    val = hermpoly_rec(n, x0)
    dx = val[1]./val[2]
    dx[ isnan( dx ) ] = 0
    x0 = x0 - dx
    if norm(dx, Inf)<sqrt(eps(Float64))
        break
    end
end
x = x0/sqrt(2)
w = exp(-x.^2)./val[2].^2           # quadrature weights

x = (x, w)
end

function hermpoly_rec( n::Int, x0)
# HERMPOLY_rec evaluation of scaled Hermite poly using recurrence

# evaluate:
Hold = exp(-x0.^2/4)
H = x0.*exp(-x0.^2/4)
for k = 1:n-1
    Hold, H = H, (x0.*H./sqrt(k+1) - Hold./sqrt(1+1/k))
end
# return (value, derivative):
val = (H, (-x0.*H + sqrt(n)*Hold))
end


function hermpoly_asy_airy(n::Int, theta)
# HERMPOLY_ASY evaluation hermite poly using Airy asymptotic formula in
# theta-space.

musq = 2n+1;
cosT = cos(theta)
sinT = sin(theta)
sin2T = 2*cosT.*sinT
eta = .5*theta - .25*sin2T
chi = -(3*eta/2).^(2/3)
phi = (-chi./sinT.^2).^(1/4)
C = 2*sqrt(pi)*musq^(1/6)*phi
Airy0 = real(airy(musq.^(2/3)*chi))
Airy1 = real(airy(1,musq.^(2/3)*chi))

# Terms in (12.10.43):
a0 = 1; b0 = 1
a1 = 15/144; b1 = -7/5*a1
a2 = 5*7*9*11/2/144^2; b2 = -13/11*a2
a3 = 7*9*11*13*15*17/6/144^3
b3 = -19/17*a3

# u polynomials in (12.10.9)
u0 = 1; u1 = (cosT.^3-6*cosT)/24
u2 = (-9*cosT.^4 + 249*cosT.^2 + 145)/1152
u3 = (-4042*cosT.^9+18189*cosT.^7-28287*cosT.^5-151995*cosT.^3-259290*cosT)/414720

#first term
A0 = 1
val = A0*Airy0

#second term
B0 = -(a0*phi.^6.*u1+a1*u0)./chi.^2
val = val + B0.*Airy1./musq.^(4/3)

# third term
A1 = (b0*phi.^12.*u2 + b1*phi.^6.*u1 + b2*u0)./chi.^3
val = val + A1.*Airy0/musq.^2

# fourth term
B1 = -(phi.^18.*u3 + a1*phi.^12.*u2 + a2*phi.^6.*u1 + a3*u0)./chi.^5
val = val + B1.*Airy1./musq.^(4/3+2)

val = C.*val

## Derivative

eta = .5*theta - .25*sin2T
chi = -(3*eta/2).^(2/3)
phi = (-chi./sinT.^2).^(1/4)
C = sqrt(2*pi)*musq^(1/3)./phi

# v polynomials in (12.10.10)
v0 = 1; v1 = (cosT.^3+6*cosT)/24
v2 = (15*cosT.^4-327*cosT.^2-143)/1152
v3 = (259290*cosT + 238425*cosT.^3 - 36387*cosT.^5 + 18189*cosT.^7 - 4042*cosT.^9)/414720

# first term
C0 = -(b0*phi.^6.*v1 + b1.*v0)./chi
dval = C0.*Airy0/musq.^(2/3)

# second term
D0 =  a0*v0
dval = dval + D0*Airy1

# third term
C1 = -(phi.^18.*v3 + b1*phi.^12.*v2 + b2*phi.^6.*v1 + b3*v0)./chi.^4
dval = dval + C1.*Airy0/musq.^(2/3+2)

#fourth term
D1 = (a0*phi.^12.*v2 + a1*phi.^6.*v1 + a2*v0)./chi.^3
dval = dval + D1.*Airy1/musq.^2

dval = C.*dval

val = (val, dval)
end

function HermiteInitialGuesses( n::Int )
#HERMITEINTITIALGUESSES(N), Initial guesses for Hermite zeros.
#
# [1] L. Gatteschi, Asymptotics and bounds for the zeros of Laguerre
# polynomials: a survey, J. Comput. Appl. Math., 144 (2002), pp. 7-27.
#
# [2] F. G. Tricomi, Sugli zeri delle funzioni di cui si conosce una
# rappresentazione asintotica, Ann. Mat. Pura Appl. 26 (1947), pp. 283-300.

# Error if n < 20 because initial guesses are based on asymptotic expansions: 
@assert n>=20    
    
# Gatteschi formula involving airy roots [1].
# These initial guess are good near x = sqrt(n+1/2);
if mod(n,2) == 1
    m = (n-1)>>1
    bess = collect(1:m)*pi
    a = .5
else
    m = n>>1
    bess = (collect(0:m-1)+.5)*pi
    a = -.5
end
nu = 4*m + 2*a + 2
T(t) = t.^(2/3).*(1+5/48*t.^(-2)-5/36*t.^(-4)+(77125/82944)*t.^(-6) -108056875/6967296*t.^(-8)+162375596875/334430208*t.^(-10));
airyrts = -T(3/8*pi*(4*collect(1:m)-1))

airyrts_exact = [-2.338107410459762           # Exact Airy roots.
    -4.087949444130970
    -5.520559828095555
    -6.786708090071765
    -7.944133587120863
    -9.022650853340979
    -10.040174341558084
    -11.008524303733260
    -11.936015563236262
    -12.828776752865757]
airyrts[1:10] = airyrts_exact  # correct first 10.

x_init = sqrt(abs(nu + 2^(2/3)*airyrts*nu^(1/3) + 1/5*2^(4/3)*airyrts.^2*nu^(-1/3) +
    (11/35-a^2-12/175*airyrts.^3)/nu + (16/1575*airyrts+92/7875*airyrts.^4)*2^(2/3)*nu^(-5/3) -
    (15152/3031875*airyrts.^5+1088/121275*airyrts.^2)*2^(1/3)*nu^(-7/3)))
x_init_airy = real( flipdim(x_init,1) )

# Tricomi initial guesses. Equation (2.1) in [1]. Originally in [2].
# These initial guesses are good near x = 0 . Note: zeros of besselj(+/-.5,x)
# are integer and half-integer multiples of pi.
# x_init_bess =  bess/sqrt(nu).*sqrt((1+ (bess.^2+2*(a^2-1))/3/nu^2) );
Tnk0 = pi/2*ones(m)
nu = (4*m+2*a+2)
rhs = (4*m-4*collect(1:m)+3)./nu*pi

for k = 1:7
    val = Tnk0 - sin(Tnk0) - rhs
    dval = 1 - cos(Tnk0)
    dTnk0 = val./dval
    Tnk0 = Tnk0 - dTnk0
end

tnk = cos(Tnk0/2).^2
x_init_sin = sqrt(nu*tnk - (5./(4*(1-tnk).^2) - 1./(1-tnk)-1+3*a^2)/3/nu)

# Patch together
p = 0.4985+eps(Float64)
x_init = [x_init_sin[1:convert(Int,floor(p*n))] ;
x_init_airy[convert(Int,ceil(p*n)):end]]

if mod(n, 2) == 1
    x_init = [0 ; x_init]
    x_init = x_init[1:m+1]
else
    x_init = x_init[1:m]
end

return x_init
end


function hermpts_gw( n::Int )
# Golub--Welsch algorithm. Used here for n<=20.

    beta = sqrt(.5*(1:n-1))              # 3-term recurrence coeffs
    T = diagm(beta, 1) + diagm(beta, -1)   # Jacobi matrix
    (D, V) = eig(T)                      # Eigenvalue decomposition
    indx = sortperm(D)                  # Hermite points
    x = D[indx]
    w = sqrt(pi)*V[1,indx].^2            # weights

    # Enforce symmetry:
    ii = floor(Int, n/2)+1:n
    x = x[ii]
    w = w[ii]
    return (x,w)
end