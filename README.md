FastGaussQuadrature
=========
A Julia package to compute `n`-point Gauss quadrature nodes and weights to 16-digit accuracy and in `O(n)` time. So far the package includes `gausschebyshev()`, `gausslegendre()`, `gaussjacobi()`, `gaussradau()`, and `gausslobatto()`. This package is heavily influenced by <a href="http://www.chebfun.org">Chebfun</a>. 

An introduction to Gauss quadrature can be found <a href="http://en.wikipedia.org/wiki/Gaussian_quadrature">here</a>.

## Our Aims 

* The fastest Julia code for Gauss quadrature nodes and weights (without tabulation). 

* Change the perception that Gauss quadrature rules are expensive to compute. 

## Examples 
Here we compute `100000` nodes and weights of the Gauss rules. Try a million or ten million. 

```
tic(), gausschebyshev( 100000 ); toc()
elapsed time: 0.007636825 seconds

tic(), gausslegendre( 100000 ); toc() 
elapsed time: 0.017749388 seconds

tic(), gaussjacobi( 100000, .9, -.1 ); toc() 
elapsed time: 4.670444327 seconds

tic(), gaussradau( 100000 ); toc() 
elapsed time: 4.51240011 seconds

tic(), gausslobatto( 100000 ); toc() 
elapsed time: 3.989099163 seconds
```

Here is (probably) a world record: The largest Gauss-Legendre quadrature rule ever computed: 
```
tic(), gausslegendre( 100000001 ); toc()
elapsed time: 98.554081035
```
(A little larger and the nodes coalesce in 16-digits of precision.)

## The algorithm for Gauss-Chebyshev
There are four kinds of Gauss-Chebyshev quadrature rules, corresponding to four weight functions: 

1. 1st kind, weight function `w(x) = 1/sqrt(1-x^2)`

2. 2nd kind, weight function `w(x) = sqrt(1-x^2)` 

3. 3rd kind, weight function `w(x) = sqrt((1+x)/(1-x))`

4. 4th kind, weight function `w(x) = sqrt((1-x)/(1+x))` 

They are all have explicit simple formulas for the nodes and weights <a href="http://books.google.com/books?id=8FHf0P3to0UC&lpg=PP1&pg=PA180#v=onepage&q&f=false">[3]</a>. 
## The algorithm for Gauss-Legendre
Gauss quadrature for the weight function `w(x) = 1`. 

* For `n<=5`: Use an analytic expression.
 
* For `n<=60`: Use Newton's method to solve `Pn(x)=0`. Evaluate `Pn` and `Pn'` by 3-term recurrence. Weights are related to `Pn'`. 
 
* For `n>60`: Use asymptotic expansions for the Legendre nodes and weights <a href="http://epubs.siam.org/doi/abs/10.1137/140954969">[1]</a>.  

## The algorithm for Gauss-Jacobi
Gauss quadrature for the weight functions `w(x) = (1-x)^a(1+x)^b`, `a,b>-1`.

*  For `n<=100`: Use Newton's method to solve `Pn(x)=0`. Evaluate `Pn` and `Pn'` by three-term recurrence.

*  For `n>100`: Use Newton's method to solve `Pn(x)=0`. Evaluate `Pn` and `Pn'` by an asymptotic expansion (in the interior of `[-1,1]`) and the three-term recurrence `O(n^-2)` close to the endpoints. (This is a small modification to the algorithm described in <a href="http://epubs.siam.org/doi/abs/10.1137/120889873">[2]</a>.) 

*Warning:* `a` and `b` need to be relatively small `(-1<a,b<10)`. 

## The algorithm for Gauss-Radau
Gauss quadrature for the weight function `w(x)=1`, except the endpoint `-1` is included as a quadrature node. 

The Gauss-Radau nodes and weights can be computed via the `(0,1)` Gauss-Jacobi nodes and weights <a href="http://epubs.siam.org/doi/abs/10.1137/120889873">[2]</a>. 
 
## The algorithm for Gauss-Lobatto
Gauss quadrature for the weight function `w(x)=1`, except the endpoints `-1` and `1` are included as nodes. 

The Gauss-Lobatto nodes and weights can be computed via the `(1,1)` Gauss-Jacobi nodes and weights <a href="http://epubs.siam.org/doi/abs/10.1137/120889873">[2]</a>. 

## References:
[1] I. Bogaert, <a href="http://epubs.siam.org/doi/abs/10.1137/140954969">"Iteration-free computation of Gauss-Legendre quadrature
       nodes and weights"</a>, SIAM J. Sci. Comput., 36(3), A1008-A1026, 2014.

[2] N. Hale and A. Townsend, <a href="http://epubs.siam.org/doi/abs/10.1137/120889873">"Fast and accurate computation of Gauss-Legendre and Gauss-Jacobi quadrature 
       nodes and weights"</a>, SIAM J. Sci. Comput., 2012.

[3] J. C. Mason and D. C. Handscomb, <a href="http://books.google.com/books?id=8FHf0P3to0UC&lpg=PP1&dq=Mason%20and%20Handscomb&pg=PP1#v=onepage&q=Mason%20and%20Handscomb&f=false">"Chebyshev Polynomials"</a>, CRC Press, 2002.
