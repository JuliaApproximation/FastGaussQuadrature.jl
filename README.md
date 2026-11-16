FastGaussQuadrature.jl
=========
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaApproximation.github.io/FastGaussQuadrature.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaApproximation.github.io/FastGaussQuadrature.jl/dev)
[![Build Status](https://github.com/JuliaApproximation/FastGaussQuadrature.jl/workflows/CI/badge.svg)](https://github.com/JuliaApproximation/FastGaussQuadrature.jl/actions)
[![codecov](https://codecov.io/gh/JuliaApproximation/FastGaussQuadrature.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaApproximation/FastGaussQuadrature.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl)

A Julia package for fast computation of `n`-point Gauss quadrature nodes and weights. The default large-order rules use linear-work sweeps at fixed working precision. Small rules use explicit formulas, recurrence, or native-precision iteration. Accuracy depends on the rule, order, parameters, and numeric type; see the [algorithm documentation](docs/src/gaussquadrature.md).

The package includes `gausslegendre`, `gaussjacobi`, `gaussradau`, `gausslobatto`, `gausslaguerre`, `gausshermite`, and the four Chebyshev rules `gausschebyshevt`, `gausschebyshevu`, `gausschebyshevv`, and `gausschebyshevw`.

This package is heavily influenced by [Chebfun](http://www.chebfun.org).

An introduction to Gauss quadrature can be found [here](http://en.wikipedia.org/wiki/Gaussian_quadrature). For a quirky account on the history of computing Gauss-Legendre quadrature, see [[6]](https://www.siam.org/publications/siam-news/articles/the-race-to-compute-high-order-gauss-legendre-quadrature/).

## Our Aims

* The fastest Julia code for Gauss quadrature nodes and weights (without tabulation).
* Change the perception that Gauss quadrature rules are expensive to compute.

## Example usage

```julia
julia> using FastGaussQuadrature, LinearAlgebra

julia> nodes, weights = gausslegendre(100000);

# integrates f(x) = x^2 from -1 to 1
julia> dot(weights, nodes.^2) ≈ 2/3
true
```

## Floating-point types

The default rules support `BigFloat` and other compatible floating-point types:

```julia
using FastGaussQuadrature

nodes, weights = setprecision(BigFloat, 256) do
    gausslaguerre(BigFloat, 32)
end
```

Support for other `AbstractFloat` types depends on the arithmetic and special functions they provide; see the [type requirements](docs/src/gaussquadrature.md#type-requirements). For parameterized rules, pass parameters of the desired type, for example `gaussjacobi(32, big"0.1", big"0.2")` or `gausslaguerre(32, 0.5f0)`. The minimum working precision is `Float32` (24 significand bits). Requests for lower precision, including low-bit `BigFloat`, generally compute in `Float32` and convert the completed nodes and weights to the requested type. For example, `gausslaguerre(Float16, 32)` returns `Float16` vectors computed in `Float32`. Legendre rules with `n > 60` and 24–53 bits of precision evaluate their expansion in the requested working precision; lower precisions compute a completed `Float32` rule before conversion.

The requested type specifies the output format, not a uniform relative-error bound. Node and weight errors can grow differently with order, parameters, and position within the rule. Accurate moments alone do not establish the accuracy of individual nodes and weights; see the [accuracy discussion](docs/src/gaussquadrature.md).
