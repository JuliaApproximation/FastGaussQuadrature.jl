FastGaussQuadrature.jl
=========
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaApproximation.github.io/FastGaussQuadrature.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaApproximation.github.io/FastGaussQuadrature.jl/dev)
[![Build Status](https://github.com/JuliaApproximation/FastGaussQuadrature.jl/workflows/CI/badge.svg)](https://github.com/JuliaApproximation/FastGaussQuadrature.jl/actions)
[![codecov](https://codecov.io/gh/JuliaApproximation/FastGaussQuadrature.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaApproximation/FastGaussQuadrature.jl)

A Julia package to compute `n`-point Gauss quadrature nodes and weights to 16-digit accuracy and in `O(n)` time.
So far the package includes `gausschebyshev()`, `gausslegendre()`, `gaussjacobi()`, `gaussradau()`, `gausslobatto()`, `gausslaguerre()`, and `gausshermite()`.
This package is heavily influenced by [Chebfun](http://www.chebfun.org).

An introduction to Gauss quadrature can be found [here](http://en.wikipedia.org/wiki/Gaussian_quadrature).
For a quirky account on the history of computing Gauss-Legendre quadrature, see [[6]](http://pi.math.cornell.edu/~ajt/papers/QuadratureEssay.pdf).

## Our Aims

* The fastest Julia code for Gauss quadrature nodes and weights (without tabulation).
* Change the perception that Gauss quadrature rules are expensive to compute.

## Example usage
```julia
julia> @time nodes, weights = gausslegendre( 100000 );
  0.002192 seconds (10 allocations: 2.289 MiB)

# integrates f(x) = x^2 from -1 to 1
julia> @time dot( weights, nodes.^2 )
  0.000184 seconds (7 allocations: 781.422 KiB)
0.6666666666666665
```
