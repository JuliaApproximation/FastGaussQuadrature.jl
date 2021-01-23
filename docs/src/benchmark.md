# Benchmark

Here we compute `100000` nodes and weights of the Gauss rules.
Try a million or ten million.

```@repl
using LinearAlgebra, FastGaussQuadrature
@time gausschebyshev(100000);
@time gausslegendre(100000);
@time gaussjacobi(100000, 0.9, -0.1);
@time gaussradau(100000);
@time gausslobatto(100000);
@time gausslaguerre(100000);
@time gausshermite(100000);
```

The paper [[1]](http://epubs.siam.org/doi/abs/10.1137/140954969) computed a billion Gauss-Legendre nodes.
So here we will do a billion + 1.

```@repl
using LinearAlgebra, FastGaussQuadrature
@time gausslegendre(1000000001);
```

(The nodes near the endpoints coalesce in 16-digits of precision.)
