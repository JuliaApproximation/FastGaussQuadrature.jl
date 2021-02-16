# FastGaussQuadrature.jl

## Abstract

FastGaussQuadrature.jl is a Julia package to compute `n`-point Gauss quadrature nodes and weights to 16-digit accuracy and in `O(n)` time.
So far the package includes `gausschebyshev()`, `gausslegendre()`, `gaussjacobi()`, `gaussradau()`, `gausslobatto()`, `gausslaguerre()`, and `gausshermite()`.
This package is heavily influenced by [Chebfun](http://www.chebfun.org).

An introduction to Gauss quadrature can be found [here](http://en.wikipedia.org/wiki/Gaussian_quadrature).
For a quirky account on the history of computing Gauss-Legendre quadrature, see [[6]](http://pi.math.cornell.edu/~ajt/papers/QuadratureEssay.pdf).

## Our Aims

* The fastest Julia code for Gauss quadrature nodes and weights (without tabulation).
* Change the perception that Gauss quadrature rules are expensive to compute.

## First example
To check an integral
```math
\int_{-1}^{1} x^4 dx = \frac{2}{5}
```
by numerically, try following code.

```@repl
using FastGaussQuadrature, LinearAlgebra
x, w = gausslegendre(3)
f(x) = x^4
I = dot(w, f.(x))
I ≈ 2/5
```
