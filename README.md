FastGaussQuadrature.jl
=========
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaApproximation.github.io/FastGaussQuadrature.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaApproximation.github.io/FastGaussQuadrature.jl/dev)
[![Build Status](https://travis-ci.com/JuliaApproximation/FastGaussQuadrature.jl.svg?branch=master)](https://travis-ci.com/JuliaApproximation/FastGaussQuadrature.jl)
[![codecov](https://codecov.io/gh/JuliaApproximation/FastGaussQuadrature.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaApproximation/FastGaussQuadrature.jl)

A Julia package to compute `n`-point Gauss quadrature nodes and weights to 16-digit accuracy and in `O(n)` time.
So far the package includes `gausschebyshev()`, `gausslegendre()`, `gaussjacobi()`, `gaussradau()`, `gausslobatto()`, `gausslaguerre()`, and `gausshermite()`.
This package is heavily influenced by [Chebfun](http://www.chebfun.org).

An introduction to Gauss quadrature can be found [here](http://en.wikipedia.org/wiki/Gaussian_quadrature).
For a quirky account on the history of computing Gauss-Legendre quadrature, see [[6]](http://pi.math.cornell.edu/~ajt/papers/QuadratureEssay.pdf).

## Our Aims

* The fastest Julia code for Gauss quadrature nodes and weights (without tabulation).
* Change the perception that Gauss quadrature rules are expensive to compute.

## Examples
Here we compute `100000` nodes and weights of the Gauss rules.
Try a million or ten million.

```julia
julia> @time gausschebyshev( 100000 );
  0.001788 seconds (4 allocations: 1.526 MiB)

julia> @time gausslegendre( 100000 );
  0.002976 seconds (10 allocations: 2.289 MiB)

julia> @time gaussjacobi( 100000, .9, -.1 );
  0.894373 seconds (3.59 k allocations: 1.255 GiB, 36.38% gc time)

julia> @time gaussradau( 100000 );
  0.684122 seconds (3.59 k allocations: 1.256 GiB, 21.71% gc time)

julia> @time gausslobatto( 100000 );
  0.748166 seconds (3.57 k allocations: 1.256 GiB, 27.78% gc time)

julia> @time gausslaguerre( 100000 );
  0.156867 seconds (7 allocations: 2.292 MiB)

julia> @time gausshermite( 100000 );
  0.175055 seconds (386 allocations: 67.916 MiB, 9.18% gc time)
```

The paper [[1]](http://epubs.siam.org/doi/abs/10.1137/140954969) computed a billion Gauss-Legendre nodes.
So here we will do a billion + 1.
```julia
julia> @time gausslegendre( 1000000001 );
 24.441304 seconds (10 allocations: 22.352 GiB, 2.08% gc time)
```
(The nodes near the endpoints coalesce in 16-digits of precision.)

## The algorithm for Gauss-Chebyshev
There are four kinds of Gauss-Chebyshev quadrature rules, corresponding to four weight functions:

1. 1st kind, weight function `w(x) = 1/sqrt(1-x^2)`

2. 2nd kind, weight function `w(x) = sqrt(1-x^2)`

3. 3rd kind, weight function `w(x) = sqrt((1+x)/(1-x))`

4. 4th kind, weight function `w(x) = sqrt((1-x)/(1+x))`

They are all have explicit simple formulas for the nodes and weights [[4]](https://books.google.co.jp/books?id=8FHf0P3to0UC).

## The algorithm for Gauss-Legendre
Gauss quadrature for the weight function `w(x) = 1`.

* For `n ≤ 5`: Use an analytic expression.
* For `n ≤ 60`: Use Newton's method to solve `Pₙ(x)=0`. Evaluate Legendre polynomials `Pₙ` and their derivatives `Pₙ'` by 3-term recurrence. Weights are related to `Pₙ'`.
* For `n > 60`: Use asymptotic expansions for the Legendre nodes and weights [[1]](http://epubs.siam.org/doi/abs/10.1137/140954969).

## The algorithm for Gauss-Jacobi
Gauss quadrature for the weight functions `w(x) = (1-x)^a(1+x)^b`, `a,b > -1`.

* For `n ≤ 100`: Use Newton's method to solve `Pₙ(x)=0`. Evaluate `Pₙ` and `Pₙ'` by three-term recurrence.
* For `n > 100`: Use Newton's method to solve `Pₙ(x)=0`. Evaluate `Pₙ` and `Pₙ'` by an asymptotic expansion (in the interior of `[-1,1]`) and the three-term recurrence `O(n^-2)` close to the endpoints. (This is a small modification to the algorithm described in [[3]](http://epubs.siam.org/doi/abs/10.1137/120889873).)
* For `max(a,b) > 5`: Use the Golub-Welsch algorithm requiring `O(n^2)` operations.

## The algorithm for Gauss-Radau
Gauss quadrature for the weight function `w(x)=1`, except the endpoint `-1` is included as a quadrature node.

The Gauss-Radau nodes and weights can be computed via the `(0,1)` Gauss-Jacobi nodes and weights[[3]](http://epubs.siam.org/doi/abs/10.1137/120889873).

## The algorithm for Gauss-Lobatto
Gauss quadrature for the weight function `w(x)=1`, except the endpoints `-1` and `1` are included as nodes.

The Gauss-Lobatto nodes and weights can be computed via the `(1,1)` Gauss-Jacobi nodes and weights[[3]](http://epubs.siam.org/doi/abs/10.1137/120889873).

## The algorithm for Gauss-Laguerre
Gauss quadrature for the weight function `w(x) = exp(-x)` on `[0,Inf)`

* For `n < 128`: Use the Golub-Welsch algorithm.
* For `method=GLR`: Use the Glaser-Lui-Rohklin algorithm. Evaluate Laguerre polynomials `Lₙ` and their derivatives `Lₙ'` by using Taylor series expansions near roots generated by solving the second-order differential equation that `Lₙ` satisfies, see [[2]](http://epubs.siam.org/doi/pdf/10.1137/06067016X).
* For `n ≥ 128`: Use a Newton procedure on Riemann-Hilbert asymptotics of Laguerre polynomials, see [5], based on [8]. There are some heuristics to decide which expression to use, it allows a general weight `w(x) = x^alpha exp(-q_m x^m)` and this is O(sqrt(n)) when allowed to stop when the weights are below the smallest positive floating point number.

## The algorithm for Gauss-Hermite
Gauss quadrature for the weight function `w(x) = exp(-x^2)` on the real line.

* For `n < 200`: Use Newton's method to solve `Hₙ(x)=0`. Evaluate Hermite polynomials `Hₙ` and their derivatives `Hₙ'` by three-term recurrence.
* For `n ≥ 200`: Use Newton's method to solve `Hₙ(x)=0`. Evaluate `Hₙ` and `Hₙ'` by a uniform asymptotic expansion, see [[7]](http://arxiv.org/abs/1410.5286).

The paper [[7]](http://arxiv.org/abs/1410.5286) also derives an `O(n)` algorithm for generalized Gauss-Hermite nodes and weights associated to weight functions of the form `exp(-V(x))`, where `V(x)` is a real polynomial.

## Example usage
```julia
julia> @time nodes, weights = gausslegendre( 100000 );
  0.002192 seconds (10 allocations: 2.289 MiB)

# integrates f(x) = x^2 from -1 to 1
julia> @time dot( weights, nodes.^2 )
  0.000184 seconds (7 allocations: 781.422 KiB)
0.6666666666666665
```

## References:
[1] I. Bogaert, ["Iteration-free computation of Gauss-Legendre quadrature nodes and weights"](http://epubs.siam.org/doi/abs/10.1137/140954969), SIAM J. Sci. Comput., 36(3), A1008-A1026, 2014.

[2] A. Glaser, X. Liu, and V. Rokhlin. ["A fast algorithm for the calculation of the roots of special functions."](http://epubs.siam.org/doi/pdf/10.1137/06067016X) SIAM J. Sci. Comput., 29 (2007), 1420-1438.

[3] N. Hale and A. Townsend, ["Fast and accurate computation of Gauss-Legendre and Gauss-Jacobi quadrature nodes and weights"](http://epubs.siam.org/doi/abs/10.1137/120889873), SIAM J. Sci. Comput., 2012.

[4] J. C. Mason and D. C. Handscomb, ["Chebyshev Polynomials"](https://books.google.co.jp/books?id=8FHf0P3to0UC), CRC Press, 2002.

[5] P. Opsomer, (in preparation).

[6] A. Townsend, [The race for high order Gauss-Legendre quadrature](http://pi.math.cornell.edu/~ajt/papers/QuadratureEssay.pdf), in SIAM News, March 2015.

[7] A. Townsend, T. Trogdon, and S. Olver, ["Fast computation of Gauss quadrature nodes and weights on the whole real line"](http://arxiv.org/abs/1410.5286), to appear in IMA Numer. Anal., 2014.

[8] M. Vanlessen, "Strong asymptotics of Laguerre-Type orthogonal polynomials and applications in Random Matrix Theory", Constr. Approx., 25:125-175, 2007.
