FastGauss
=========
A Julia package to compute 16-digit accurate Gauss nodes and weights. At the moment we have Gauss-Legendre and Gauss-Chebyshev implemented. 

## Our Aim 
To be the fastest Julia code for computing 16-digit accurate Gauss nodes and weights (without tabulation).

## Example 
```
tic(), GaussLegendre( 1000 ); toc()
elapsed time: 0.000829828 seconds

tic(), GaussLegendre( 1000000 ); toc()
elapsed time: 0.310792745 seconds

tic(), GaussChebyshev( 1000000 ); toc() 
elapsed time: 0.04614016 seconds
```

## The underlying algorithm for Gauss-Legendre
 For n<=5: Use an analytic expression.
 
 For n<=60: Use Newton's method to solve Pn(x)=0. Weights are related to P'n(x), see [2].  
 
 For n>60: Use asymptotic expansions for the Legendre nodes and weights, see [1].  

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

