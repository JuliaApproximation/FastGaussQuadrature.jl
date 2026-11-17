# FastGaussQuadrature.jl

## Abstract

FastGaussQuadrature.jl computes `n`-point Gauss quadrature nodes and weights. Large-order Float64 rules use linear-work asymptotic or fourth-order iterative methods where applicable. Accuracy and computational cost depend on the rule, parameters and numeric type; see [Gaussian Quadrature](@ref). The package includes `gausslegendre`, `gaussjacobi`, `gaussradau`, `gausslobatto`, `gausslaguerre`, `gausshermite`, and the four Chebyshev rules `gausschebyshevt`, `gausschebyshevu`, `gausschebyshevv`, and `gausschebyshevw`. This package is heavily influenced by [Chebfun](http://www.chebfun.org).

An introduction to Gauss quadrature can be found [here](http://en.wikipedia.org/wiki/Gaussian_quadrature). For a quirky account on the history of computing Gauss-Legendre quadrature, see [[6]](https://www.siam.org/publications/siam-news/articles/the-race-to-compute-high-order-gauss-legendre-quadrature/).

## Our Aims

* The fastest Julia code for Gauss quadrature nodes and weights (without tabulation).
* Change the perception that Gauss quadrature rules are expensive to compute.

## First example

To check the integral
```math
\int_{-1}^{1} x^4 dx = \frac{2}{5}
```
numerically, try the following code.

```@repl
using FastGaussQuadrature, LinearAlgebra
x, w = gausslegendre(3)
f(x) = x^4
I = dot(w, f.(x))
I ≈ 2/5
```

The rules also support `BigFloat`, low-precision types, and compatible custom floating-point types. See [Floating-point types](@ref) for examples and the required arithmetic and auxiliary functions.
