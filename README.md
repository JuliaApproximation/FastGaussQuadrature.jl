FastGauss
=========
A Julia package to compute 16-digit accurate Gauss nodes and weights. At the moment we only have Gauss-Legendre implemented. 

## Our Aim 
To be the fastest Julia code for computing 16-digit accurate Gauss nodes and weights (without tabulation).

<<<<<<< HEAD
## Example 
```
tic(), GaussLegendre( 1000 ); toc()
elapsed time: 0.000829828 seconds
=======
  GAUSS-LEGENDRE EXAMPLE: 
  
          tic(), GaussLegendre( 1000 ); toc()
          elapsed time: 0.000829828 seconds
>>>>>>> 9d9759a9c38bfba2cc4713c210fcccb989a6b29d

tic(), GaussLegendre( 1000000 ); toc()
elapsed time: 0.310792745 seconds
```

## The underlying algorithm
 For n<=5: Use an analytic expression. <return>
 For n<=60: Use Newton's method to solve P_n(x)=0. Weights are related to P'_n(x), see [2].<return> 
 For n>60: Use asymptotic expansions for the Legendre nodes and weights, see [1].<return>

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

