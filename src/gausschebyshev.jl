function gausschebyshev( n::Int64, kind::Int64=1 )
    # GAUSS-CHEBYSHEV NODES AND WEIGHTS. 
    
    # Use known explicit formulas. Complexity O(n).
    if kind == 1 
        # Gauss-ChebyshevT quadrature, i.e., w(x) = 1/sqrt(1-x^2)
        ([cos((2*k-1)*pi/2n) for k=n:-1:1], fill(pi./n,n))
    elseif kind == 2 
        # Gauss-ChebyshevU quadrature, i.e., w(x) = sqrt(1-x^2)
        ([cos(k*pi./(n+1)) for k=n:-1:1], [pi/(n+1)*sin(k./(n+1)*pi).^2 for k=n:-1:1])
    elseif kind == 3 
        # Gauss-ChebyshevV quadrature, i.e., w(x) = sqrt((1+x)/(1-x))
        ([cos((k-.5)*pi/(n+.5)) for k=n:-1:1], [2*pi/(n+.5)*cos((k-.5)*pi/(2(n+.5))).^2 for k=n:-1:1])
    elseif kind == 4 
        # Gauss-ChebyshevW quadrature, i.e., w(x) = sqrt((1-x)/(1+x))
        ([cos(k*pi/(n+.5)) for k=n:-1:1], [2*pi/(n+.5)*sin(k*pi/(2(n+.5))).^2 for k=n:-1:1])
    else 
       throw(ArgumentError("Chebyshev kind should be 1, 2, 3, or 4")) 
    end
end

