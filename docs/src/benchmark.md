# Benchmark

The following examples use BenchmarkTools to measure `Float64` rules with `100000` nodes. Timings and allocations are generated when this documentation is built; they depend on Julia, dependencies, hardware, order, and parameters. Run the examples locally to measure your workload.

```@repl
using BenchmarkTools, FastGaussQuadrature
@btime gausschebyshevt(100000);
@btime gausslegendre(100000);
@btime gaussjacobi(100000, 0.9, -0.1);
@btime gaussradau(100000);
@btime gausslobatto(100000);
@btime gausslaguerre(100000);
@btime gausshermite(100000);
```

Reduced Laguerre rules omit the tail below a weight cutoff in supported `Float64` regimes; see [Gauss-Laguerre quadrature](@ref). They return fewer nodes and weights, so compare both output size and accuracy when measuring performance. The algorithm and its cost also depend on the Laguerre parameter.

```@repl
using BenchmarkTools, FastGaussQuadrature
@btime gausslaguerre(10000, 0.5);
@btime gausslaguerre(10000, 0.5; reduced=true);
```

At fixed precision and bounded per-node iteration, the default large-order algorithms use linear work and storage. Two full `Float64` output vectors alone require approximately `16n` bytes, in addition to any working storage. Very large rules can also lose distinct endpoint nodes when rounded to the requested type.
