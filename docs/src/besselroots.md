# Roots of the Bessel function

`approx_besselroots` supplies approximate positive zeros of the Bessel function ``J_\nu`` for the asymptotic quadrature routines. It returns `Float64` values, even when the input parameter has higher precision. Accuracy depends on the parameter and root index; this helper is not an arbitrary-precision root solver.

```@docs
approx_besselroots(ν::Real, n::Integer)
```

Some Jacobi and Laguerre algorithms use these approximations as starting values. The generic default quadrature algorithms refine roots at their working precision; the precision of this helper does not set their output precision.
