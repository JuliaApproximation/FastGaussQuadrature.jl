function gausslegendre( n::Int64 )
# GAUSSLEGENDRE(n)  COMPUTE THE GAUSS-LEGENDRE NODES AND WEIGHTS IN O(n) time.

if n < 0
    x = (Float64[],Float64[])
elseif n == 0 
    x = (Float64[],Float64[])
elseif n == 1
    x = ([0.0],[2.0])
elseif n == 2
    x = Array(Float64,2)
    w = Array(Float64,2)
    x[1] = -1/sqrt(3)
    x[2] = 1/sqrt(3)
    w[1] = 1.0 
    w[2] =  1.0
    x = (x,w)
elseif n == 3
    x = Array(Float64,3)
    w = Array(Float64,3)
    x[1] = -sqrt(3/5) 
    x[2] = 0.0
    x[3] = sqrt(3/5) 
    w[1] = 5/9
    w[2] = 8/9
    w[3] = 5/9
    x = (x,w)
elseif n == 4 
    x = Array(Float64,4)
    w = Array(Float64,4)
    const a::Float64 = 2/7*sqrt(6/5)
    x[1] = -sqrt(3/7+a)
    x[2] = -sqrt(3/7-a)
    x[3] =  sqrt(3/7-a)
    x[4] =  sqrt(3/7+a)
    w[1] = (18 - sqrt(30))/36;
    w[2] = (18 + sqrt(30))/36;
    w[3] = (18 + sqrt(30))/36;
    w[4] = (18 - sqrt(30))/36;
    x = (x, w)
elseif n == 5
    x = Array(Float64, 5)
    w = Array(Float64, 5)
    const b::Float64 = 2sqrt(10/7)
    x[1] = -sqrt(5+b)/3 
    x[2] = -sqrt(5-b)/3
    x[3] = 0.0
    x[4] = sqrt(5-b)/3 
    x[5] = sqrt(5+b)/3
    w[1] = (322-13sqrt(70))/900
    w[2] = (322+13sqrt(70))/900
    w[3] = 128/225
    w[4] = (322+13sqrt(70))/900
    w[5] = (322-13sqrt(70))/900
    x = (x, w);
elseif n <= 60 
# NEWTON'S METHOD WITH THREE-TERM RECURRENCE:   
    x = rec( n::Int64 )
else
# USE ASYMPTOTIC EXPANSIONS:
    x = asy( n::Int64 )
end
    return x   # nodes and weights in tuple. 
end

function asy( n::Int64 )
# COMPUTE GAUSS-LEGENDRE NODES AND WEIGHTS USING ASYMPTOTIC EXPANSIONS. COMPLEXITY O(n). 
    # Nodes and weights:
    m = (mod(n,2)==0) ? n>>1 : (n+1)>>1
    a = besselZeroRoots(m)
    scale!(a,1/(n + 0.5))
    x = legpts_nodes(n, a); 
    w = legpts_weights(n, m, a)
    # Use symmetry to get the others:
    if ( mod(n, 2) == 1 )
        x = [-x[1:end-1] ; 0.0 ; x[end-1:-1:1] ]
        w = [w ; w[end-1:-1:1] ]
    else
        x = [-x ; x[end:-1:1] ]
        w = [w ; w[end:-1:1] ]
    end
    return x, w
end

function legpts_nodes(n::Int64, a::Array{Float64})
# ASYMPTOTIC EXPANSION FOR THE GAUSS-LEGENDRE NODES.
    vn = 1/(n + 0.5)
    m = length(a)
    u = Array(Float64,m); nodes = Array(Float64,m)
    for jj=1:m u[jj] = cot(a[jj]) end
    for jj=1:m nodes[jj] = a[jj] + ((u[jj]-1/a[jj])/8*vn^2) end
    nodes += ((n <= 3950) ? (6(1+u.^2)./a + 25./a.^3 - u.*(31u.^2+33))/384 : 0.0)*vn^4 +
            ((n <= 255) ? u.*(2595 + 6350*u.^2 + 3779*u.^4)/15360 -1073/5120./a.^5 +
            (1+u.^2).*(-(31*u.^2 + 11)/1024./a + u/512./a.^2 + -25/3072./a.^3) : 0.0)*vn^6
    for jj=1:m nodes[jj] = cos(nodes[jj]) end
    return nodes
end

function legpts_weights(n::Int64, m::Int64, a::Array{Float64})
# ASYMPTOTIC EXPANSION FOR THE GAUSS-LEGENDRE WEIGHTS.
    vn = 1/(n + 0.5);
    u = Array(Float64,m); ua = Array(Float64,m); weights = Array(Float64,m)
    for jj=1:m u[jj] = cot(a[jj]) end
    for jj=1:m ua[jj] = u[jj]*a[jj] end
     # Assemble coefficients:
    W1 = Array(Float64,m); W2 = Array(Float64,m); W3 = Array(Float64,m)
    for jj=1:m W1[jj] = (n <= 850000) ? (ua[jj] + a[jj]^2 - 1)/a[jj]^2/8 : 0.0 end
    for jj=1:m W2[jj] = (n <= 1500) ? ( 81 - 31*ua[jj] - 3*(1-2*u[jj]^2)*a[jj]^2 + 6*u[jj]*a[jj]^3 -
                                        (27 + 84*u[jj]^2 + 56*u[jj]^4)*a[jj]^4 )/a[jj]^4/384 : 0.0 end
    for jj=1:m W3[jj] = ( n <= 170 ) ? (187/96*u[jj]^4 + 295/256*u[jj]^2 + 151/160*u[jj]^6 + 153/1024) +
                                    (-119/768*u[jj]^2 -35/384*u[jj]^4 - 65/1024)*u[jj]/a[jj] +
                                    (5/512 + 7/384*u[jj]^4 + 15/512*u[jj]^2)/a[jj]^2 +
                                    (u[jj]^3/512 - 13/1536*u[jj])/a[jj]^3 +
                                    (-7/384*u[jj]^2 + 53/3072)/a[jj]^4 +
                                    (3749/15360*u[jj])/a[jj]^5 -1125/1024/a[jj]^6 : 0.0 end
    bJ1 = besselJ1(m)
    for jj=1:m weights[jj] = 2/(bJ1[jj]*(a[jj]/sin(a[jj]))*(1/vn^2 + W1[jj] + W2[jj]*vn^2 + W3[jj]*vn^4)) end
    return weights
end

function rec( n::Int64 ) 
# COMPUTE GAUSS-LEGENDRE NODES AND WEIGHTS USING NEWTON'S METHOD. THREE-TERM RECURENCE 
# IS USED FOR EVALUATION. COMPLEXITY O(n^2). 
    hN = mod(n,2)
    # Initial guesses: 
    x0 = asy(n)[1]; x = x0[(n-hN)/2+1:n]   
    # Perform Newton to find zeros of Legendre polynomial:
    PP = innerRec( n, x ); dx = -PP[1]./PP[2]; x += dx
    # One more Newton for derivatives: 
    PP = innerRec( n, x ); dx = -PP[1]./PP[2]; x += dx
    #PP = (n<45)? innerRec( n, x ) : PP
    # Use symmetry to get the other Legendre nodes and weights: 
    nodes = vcat( -x[(n+hN)/2:-1:hN+1] , x )
    weights = 2./((1-nodes.^2).*vcat( PP[2][(n+hN)/2:-1:1+hN] , PP[2] ).^2)
    return nodes, weights
end

function innerRec( n::Int64, x )
# EVALUATE LEGENDRE AND ITS DERIVATIVE USING THREE-TERM RECURRENCE RELATION.
    N = size(x,1) 
    myPm1 = Array(Float64,N); myPPm1 = Array(Float64,N) 
    for j = 1:N
        xj = x[j]; Pm2 = 1.0; Pm1 = xj; PPm1 = 1.0; PPm2 = 0.0
        for k = 1:n-1
            Pm2, Pm1 = Pm1, ((2k+1)*Pm1*xj - k*Pm2)/(k+1)
            PPm2, PPm1 = PPm1, ((2k+1)*(Pm2 + xj*PPm1)-k*PPm2)/(k+1)
        end
        myPm1[j] = Pm1
        myPPm1[j] = PPm1 
    end
    return myPm1, myPPm1
end


function besselZeroRoots(m::Int64)
#BESSEL0ROOTS ROOTS OF BESSELJ(0,x). USE ASYMPTOTICS.
# Use McMahon's expansion for the remainder (NIST, 10.21.19):
    jk = Array(Float64,m); ak = Array(Float64,m); ak82 = Array(Float64,m)
    p = [1071187749376/315 0.0 -401743168/105 0.0 120928/15 0.0 -124/3 0.0 1.0 0.0]
    for jj=1:m ak[jj] = pi*( jj -.25 ) end
    for jj=1:m ak82[jj] = (.125/ak[jj])^2 end
    # First 20 are precomputed:
    twenty = [2.4048255576957728, 5.5200781102863106, 8.6537279129110122,
            11.791534439014281, 14.930917708487785, 18.071063967910922,
            21.211636629879258, 24.352471530749302, 27.493479132040254,
            30.634606468431975, 33.775820213573568, 36.917098353664044,
            40.058425764628239, 43.199791713176730, 46.341188371661814,
            49.482609897397817, 52.624051841114996, 55.765510755019979,
            58.906983926080942, 62.048469190227170]
    for jj=1:min(m,20) jk[jj] = twenty[jj] end
    for jj=21:min(m,47) jk[jj] = ak[jj] + .125/ak[jj]*(1 + ak82[jj]*(p[7] + ak82[jj]*(p[5] + ak82[jj]*p[3]))) end
    for jj=48:min(m,344) jk[jj] = ak[jj] + .125/ak[jj]*(1 + ak82[jj]*(p[7]+ ak82[jj]*p[5])) end
    for jj=345:min(m,13191) jk[jj] = ak[jj] + .125/ak[jj]*(1 + ak82[jj]*p[7]) end
    for jj=13192:m jk[jj] = ak[jj] + .125/ak[jj] end
    return jk
end


function besselJ1( m::Int64 )
# BESSELJ1 EVALUATE BESSELJ(1,x)^2 AT ROOTS OF BESSELJ(0,x). USE ASYMPTOTICS.
# Use Taylor series of (NIST, 10.17.3) and McMahon's expansion (NIST, 10.21.19):
    Jk2 = Array(Float64,m); ak = Array(Float64,m); ak2 = Array(Float64,m)
    c = [-171497088497/15206400, 461797/1152, -172913/8064, 151/80, -7/24, 0.0, 2.0]
    for jj=1:m ak[jj] = pi*(jj-.25) end
    for jj=1:m ak2[jj] = (1/ak[jj])^2 end
    # First 10 are precomputed:
    ten = [ 0.2695141239419169, 0.1157801385822037, 0.07368635113640822,
            0.05403757319811628, 0.04266142901724309, 0.03524210349099610,
            0.03002107010305467, 0.02614739149530809, 0.02315912182469139, 0.02078382912226786]
    for jj=1:min(m,10) Jk2[jj] = ten[jj] end
    for jj=11:min(m,15) Jk2[jj] = 1/(pi*ak[jj])*(c[7] + ak2[jj]^2*(c[5] + ak2[jj]*(c[4] + 
                ak2[jj]*(c[3] + ak2[jj]*(c[2]+ak2[jj]*c[1]))))) end
    for jj=16:min(m,21) Jk2[jj] = 1/(pi*ak[jj])*(c[7] + ak2[jj]^2*(c[5] + ak2[jj]*(c[4] + ak2[jj]*(c[3] + ak2[jj]*c[2])))) end
    for jj=22:min(m,55) Jk2[jj] =  1/(pi*ak[jj])*(c[7] + ak2[jj]^2*(c[5] + ak2[jj]*(c[4] + ak2[jj]*c[3]))) end
    for jj=56:min(m,279) Jk2[jj] = 1/(pi*ak[jj])*(c[7] + ak2[jj]^2*(c[5] + ak2[jj]*c[4])) end
    for jj=280:min(m,2279) Jk2[jj] = 1/(pi*ak[jj])*(c[7] + ak2[jj]^2*c[5]) end
    for jj=2280:m Jk2[jj] = 1/(pi*ak[jj])*c[7] end
    return Jk2
end
