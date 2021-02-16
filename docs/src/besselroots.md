# Roots of Bessel function

Since [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl) doesn't have a method to calculate roots of [Bessel function](https://en.wikipedia.org/wiki/Bessel_function), we implemented `approx_besselroots`.

```@docs
approx_besselroots(ν::Real, n::Integer)
```

This method `approx_besselroots` is used to calculate `gaussjacobi` and `gausslaguerre`.
