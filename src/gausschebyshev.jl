function gausschebyshev( n::Int64, kind::Int64=1 )
# GAUSS-CHEBYSHEV NODES AND WEIGHTS. 

x = (Array(Float64,n), Array(Float64,n))
# Use known explicit formulas. Complexity O(n).
if kind == 1 
    # Gauss-ChebyshevT quadrature, i.e., w(x) = 1/sqrt(1-x^2)
    x = (cos((2*[n:-1:1]-1)*pi/2n), pi./n*ones(n))
elseif kind == 2 
    # Gauss-ChebyshevU quadrature, i.e., w(x) = sqrt(1-x^2)
    x = (cos([n:-1:1]*pi./(n+1)), pi/(n+1)*sin([n:-1:1]./(n+1)*pi).^2 )
elseif kind == 3 
    # Gauss-ChebyshevV quadrature, i.e., w(x) = sqrt((1+x)/(1-x))
    x = (cos(([n:-1:1]-.5)*pi/(n+.5)), 2*pi/(n+.5)*cos(([n:-1:1]-.5)*pi/(2(n+.5))).^2)
elseif kind == 4 
    # Gauss-ChebyshevW quadrature, i.e., w(x) = sqrt((1-x)/(1+x))
    x = (cos([n:-1:1]*pi/(n+.5)), 2*pi/(n+.5)*sin([n:-1:1]*pi/(2(n+.5))).^2)
else 
   throw(ArgumentError("Chebyshev kind should be 1, 2, 3, or 4")) 
end
return x 
end
