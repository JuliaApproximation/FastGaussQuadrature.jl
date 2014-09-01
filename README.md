FastGauss
=========
A Julia package to compute 16-digit accurate Gauss nodes and weights. At the moment we have Gauss-Legendre and Gauss-Chebyshev implemented. 

## Our Aim 
To be the fastest Julia code for computing 16-digit accurate Gauss nodes and weights (without tabulation).

## Examples 
```
tic(), GaussChebyshev( 100000 ); toc()
elapsed time: 0.007636825 seconds

tic(), GaussLegendre( 100000 ); toc() 
elapsed time: 0.092068965 seconds

tic(), GaussJacobi( 100000,.9,-.1 ); toc() 
elapsed time: 6.078796565 seconds

tic(), GaussRadau( 100000 ); toc() 
elapsed time: 6.1875638 seconds

tic(), GaussLobatto( 100000 ); toc() 
elapsed time: 4.901654062 second
```

## The algorithm for Gauss-Chebyshev
There are four kinds of Gauss-Chebyshev quadrature rules: 
1. 1st kind, weight function w(x) = 1/sqrt(1-x^2)
2. 2nd kind, weight function w(x) = sqrt(1-x^2) 
3. 3rd kind, weight function  w(x) = sqrt((1+x)/(1-x))
4. 4th kind, weight function  w(x) = sqrt((1-x)/(1+x)) 

They are all have explicit simple formulas for the nodes and weights. 
## The algorithm for Gauss-Legendre
 For `n<=5`: Use an analytic expression.
 
 For `n<=60`: Use Newton's method to solve Pn(x)=0. Evaluate Pn and Pn' by 3-term recurrence. Weights are related to Pn'. 
 
 For `n>60`: Use asymptotic expansions for the Legendre nodes and weights, see [1].  

## The algorithm for Gauss-Jacobi
  For `n<=100`: Use Newton's method to solve Pn(x)=0. Evaluate Pn and Pn' by 3-term recurrence.

  For `n>100`: Use Newton's method to solve Pn(x)=0. Evaluate Pn and Pn' by asymptotics expansion (in the interior of [-1,1]) and 3-term very near the boundary.

## The algorithm for Gauss-Radau

## The algorithm for Gauss-Lobatto


## References:
1. I. Bogaert, "Iteration-free computation of Gauss-Legendre quadrature
       nodes and weights", SIAM J. Sci. Comput., 36(3), A1008-A1026, 2014.
       (Describes the algorithm implemented here for n>60. A few modifications were made 
         to optimize for Julia. )
2. N. Hale and A. Townsend, "Fast computation of Gauss-Jacobi quadrature 
       nodes and weights", SIAM J. Sci. Comput., 2012.
       (Describes the algorithm implemented here for 5<n<=60 and an O(n) algorithm for 
        computing Gauss-Jacobi nodes and weights.) 

For a review and comparison of the recent methods for computing Gauss nodes, see [2].  

