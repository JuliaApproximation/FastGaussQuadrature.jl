function gausslaguerre( n::Int64 )
    
# TODO: This code needs optimizing. 
if  n == 0
    Float64[], Float64[]
elseif n == 1
    [1.0], [1.0]
elseif n == 2
    [2.-sqrt(2.) 2.+sqrt(2.)], [(2.+sqrt(2.))/4 (2.-sqrt(2.))/4]
elseif n <= 128 
    laguerreGW( n )   # Use Golub-Welsch
else  
    laguerreGLR( n )  # Use GLR 
end
    
end

function laguerreGW( n::Int64 )
# Calculate Gauss-Laguerre nodes and weights based on Golub-Welsch       
    
    alpha = 2*(1:n)-1       # 3-term recurrence coeffs
    beta = 1:n-1       
    T = diagm(beta,1) + diagm(alpha) + diagm(beta,-1)  # Jacobi matrix
    x, V = eig( T )                  # eigenvalue decomposition
    w = V[1,:].^2                     # Quadrature weights
    x, vec(w)
    
end

function laguerreGLR( n::Int64 )
    
    x, ders = ScaledLaguerreGLR( n )  # Nodes and L_n'(x)
    w = exp(-x)./(x.*ders.^2)         # Quadrature weights
    x, w
end

function  ScaledLaguerreGLR( n::Int64 )
    # Calculate the nodes and L_n'(x). 
    
    ders = Array(Float64, n )
    x = Array(Float64,n)
    xs = 1/(2*n+1)
    n1 = 20
    n1 = min(n1, n)
    for k = 1:n1
        xs, dxs = BoundaryLaguerreGLR(n, xs)
        ders[k] = dxs
        x[k] = xs
        xs *= 1.1
    end
    x, ders = InteriorLaguerreGLR(n, x, ders, n1)   
end

function InteriorLaguerreGLR(n::Int64, roots, ders, n1::Int64)
    
    m = 30
    hh1 = ones(m+1) 
    zz = Array(Float64,m) 
    u = Array(Float64,m+1) 
    up = Array(Float64,m+1)
    x = roots[ n1 ] 
    for j = n1 : (n - 1)
        # initial approx
        h = rk2_Lag(pi/2, -pi/2, x, n)
        h = h - x

        # scaling:
        M = 1/h 
        M2 = M^2 
        M3 = M^3 
        M4 = M^4

        # recurrence relation for Laguerre polynomials
        r = x*(n + .5 - .25*x) 
        p = x^2
        u[1] = 0.
        u[2] = ders[j]/M
        u[3] = -.5*u[2]/(M*x) - (n + .5 - .25*x)*u[1]/(x*M2)
        u[4] = -u[3]/(M*x) + ( -(1+r)*u[2]/6/M2 - (n+.5-.5*x)*u[1]/M3 ) / p
        up[1] = u[2]; up[2] = 2*u[3]*M; up[3] = 3*u[4]*M

        for k = 2:(m - 2)
              u[k+3] = ( -x*(2*k+1)*(k+1)*u[k+2]/M - (k*k+r)*u[k+1]/M2 - 
                (n+.5-.5*x)*u[k]/M3 + .25*u[k-1]/M4 ) / (p*(k+2)*(k+1))
              up[k+2] = (k+2)*u[k+3]*M
        end
        up[m+1] = 0

        # Flip for more accuracy in inner product calculation.
        u = u[m+1:-1:1]  
        up = up[m+1:-1:1]

        # Newton iteration
        hh = hh1 
        hh[m+1] = M    
        step = Inf  
        l = 0
        if M == 1
            Mhzz = M*h + zz
            hh = [M ; cumprod( Mhzz )]
            hh = hh[end:-1:1]
        end
        while ( (abs(step) > eps(Float64)) && (l < 10) )
            l = l + 1
            step = dot(u,hh)/dot(up,hh)
            h = h - step
            Mhzz = (M*h) + zz
            # Powers of h (This is the fastest way!)
            hh = [M ; cumprod(Mhzz)]     
            # Flip for more accuracy in inner product
            hh = hh[end:-1:1]          
        end

        # Update
        x = x + h
        roots[j+1] = x
        ders[j+1] = dot(up,hh)
    end
    roots, ders
end

function BoundaryLaguerreGLR(n::Int64, xs)
    
    u, up = eval_Lag(xs, n)
    theta = atan(sqrt(xs/(n + .5 - .25*xs))*up/u)
    x1 = rk2_Lag(theta, -pi/2, xs, n)

    # Newton iteration
    step = Inf  
    l = 0
    while ( (abs(step) > eps(Float64) || abs(u) > eps(Float64)) && (l < 200) )
        l = l + 1
        u, up = eval_Lag(x1, n)
        step = u/up
        x1 = x1 - step
    end

    ignored, d1 = eval_Lag(x1, n)
    
    x1, d1
end

function eval_Lag(x, n::Int64)
    # Evauate Laguerre polynomial via recurrence

    L = 0.
    Lp = 0.
    Lm2 = 0.
    Lm1 = exp(-x/2) 
    Lpm2 = 0.
    Lpm1 = 0.
    for k = 0:n-1
        L = ( (2*k+1-x).*Lm1 - k*Lm2 ) / (k + 1)
        Lp = ( (2*k+1-x).*Lpm1 - Lm1 - k*Lpm2 ) / (k + 1)
        Lm2 = Lm1 
        Lm1 = L
        Lpm2 = Lpm1 
        Lpm1 = Lp
    end
    L, Lp
end

function rk2_Lag(t, tn, x, n::Int64)
    # Runge-Kutta for Laguerre Equation

    m = 10
    h = (tn - t)/m
    for j = 1:m
        f1 = (n + .5 - .25*x)
        k1 = -h/( sqrt(abs(f1/x)) + .25*(1/x-.25/f1)*sin(2*t) )
        t = t + h
        x = x + k1   
        f1 = (n + .5 - .25*x)
        k2 = -h/( sqrt(abs(f1/x)) + .25*(1/x-.25/f1)*sin(2*t) )
        x = x + .5*(k2 - k1)
    end
    x
end
