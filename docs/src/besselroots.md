# Roots of Bessel function

Since SpecialFunctions.jl doesn't have a method to calculate roots of [Bessel function](https://en.wikipedia.org/wiki/Bessel_function), we implemented `besselroots`.

```@docs
besselroots(nu::Float64, n::Integer)
```

This method `besselroots` is used in `gaussjacobi` and `gausslaguerre`.
