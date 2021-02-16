# Benchmark

Here we compute `100000` nodes and weights of the Gauss rules.
Try a million or ten million.

```@repl
using LinearAlgebra, BenchmarkTools, FastGaussQuadrature
@btime gausschebyshev(100000);
@btime gausslegendre(100000);
@btime gaussjacobi(100000, 0.9, -0.1);
@btime gaussradau(100000);
@btime gausslobatto(100000);
@btime gausslaguerre(100000);
@btime gausshermite(100000);
```

The paper [[1]](http://epubs.siam.org/doi/abs/10.1137/140954969) computed a billion Gauss-Legendre nodes.
So here we will do a billion + 1.

```julia
julia> @btime gausslegendre(1000000001);
  23.363 s (10 allocations: 22.35 GiB)
```

(The nodes near the endpoints coalesce in 16-digits of precision.)
