function gausschebyshev(n::Int, kind::Int=1)
    # GAUSS-CHEBYSHEV NODES AND WEIGHTS.

    # Use known explicit formulas. Complexity O(n).
    if kind == 1
        # Gauss-ChebyshevT quadrature, i.e., w(x) = 1/sqrt(1-x^2)
        [cos((2 * k - 1) * π / (2 * n)) for k = n:-1:1], fill(π / n, n)
    elseif kind == 2
        # Gauss-ChebyshevU quadrature, i.e., w(x) = sqrt(1-x^2)
        ([cos(k * π / (n + 1)) for k = n:-1:1],
         [π/(n + 1) * sin(k / (n + 1) * π)^2 for k = n:-1:1])
    elseif kind == 3
        # Gauss-ChebyshevV quadrature, i.e., w(x) = sqrt((1+x)/(1-x))
        ([cos((k - .5) * π / (n + .5)) for k = n:-1:1],
         [2π / (n + .5) * cos((k - .5) * π / (2 * (n + .5)))^2 for k = n:-1:1])
    elseif kind == 4
        # Gauss-ChebyshevW quadrature, i.e., w(x) = sqrt((1-x)/(1+x))
        ([cos(k * π / (n + .5)) for k = n:-1:1],
         [2π / (n + .5) * sin(k * π / (2 * (n + .5)))^2 for k = n:-1:1])
    else
        throw(ArgumentError("Chebyshev kind should be 1, 2, 3, or 4"))
    end
end
