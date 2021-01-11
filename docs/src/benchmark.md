# Benchmark

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
