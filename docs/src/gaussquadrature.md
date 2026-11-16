# Gaussian Quadrature

The rules return ascending nodes and corresponding weights. Accuracy depends on
order, parameters, numeric type, position within the rule, and whether the nodes
and weights can be represented. Selecting an output type does not guarantee
relative errors bounded by a fixed multiple of `eps(T)` across these regimes.
Node and weight errors can grow differently; accurate moments alone do not
establish the accuracy of each node and weight.

Accuracy comparisons should use independent higher-precision references for the
exactly represented input parameters. Absolute errors are useful near zero,
relative errors for normal nonzero weights, and errors relative to the output
spacing for subnormal weights. Underflow and indistinguishable rounded nodes
are separate limitations from error in the numerical algorithm.

Large-order asymptotic methods are designed for `Float64`;
Legendre also uses its expansion for lower-precision output types. Higher
precisions use native Julia recurrence and iteration. See
[Floating-point types](@ref) for type selection and precision requirements.

## Gauss-Legendre quadrature

Gauss quadrature for ``w(x)=1`` on ``[-1,1]``. Rules through ``n=5`` use analytic
expressions. For ``6\le n\le60``, Newton iteration uses a three-term recurrence
in the working type. Initial guesses use an asymptotic expansion based on Bessel
zeros, evaluated in the working type below 53 significand bits. At 53 bits and
above, `Float64` initial guesses are converted to the working type before Newton
refinement.

For ``n>60`` and
`precision(T) <= precision(Float64)` (53 bits), rules use the ``O(n)`` expansions
of [Bogaert (2014)](https://doi.org/10.1137/140954969), evaluated in the requested
working precision at 24–53 bits. The number of retained terms depends on both
order and precision, with an omission bound that accounts for small central
nodes. The stored coefficients limit this expansion to 53 bits. Lower
precisions, including `Float16` and
low-precision `BigFloat`, compute the completed rule in `Float32` before conversion.
Endpoint weights retain small corrections separately to avoid rounding them away
near one. Above 53 bits, large rules use fourth-order Jacobi iteration at
the requested precision. Conversion is followed by validation; an error is
thrown if distinct nodes round to the same value in `T`.

```@docs
gausslegendre(n::Integer)
```

## Gauss-Hermite quadrature

Gauss quadrature for ``w(x)=\exp(-x^2)`` on the real line. Small rules use
explicit formulas or recurrence through ``n=20`` for all working types.
The cutoff keeps the polynomial recurrence and
squared derivative within the exponent range of `Float32`. Initial guesses use
GST iteration in `Float32` below 53 significand bits and `Float64` otherwise.
The guesses are converted to the working type before native Newton refinement;
the final derivatives and weights are also computed in that type.
Other rules use the fourth-order iteration of
[Gil, Segura and Temme](https://arxiv.org/html/2509.16716v2#S4),
with the final root residual correcting the exponential factor in physical
weights. This correction uses the same working precision and leaves the scaled
weights unchanged. Returned physical weights may underflow to zero.
With `normalize=true`, nodes and weights correspond to the
standard normal density.

```@docs
gausshermite(n::Integer)
```

## Gauss-Laguerre quadrature

Gauss quadrature for ``w(x)=x^\alpha\exp(-x)`` on ``[0,\infty)``, with
``\alpha>-1``. The default is ``\alpha=0``. Orders one and two use explicit
formulas; larger rules use the fourth-order iteration of
[Gil, Segura and Temme](https://arxiv.org/html/2509.16716v2#S3).
For large `Float64` rules with ``\alpha=\pm1/2``, Bessel expansions based on
[Opsomer and Huybrechs](https://doi.org/10.1016/j.cam.2023.115317) improve the
first few nodes and weights. The bounded recurrence and full asymptotic
implementation remain available as specialized manual alternatives.

When the working type is `Float64`, with ``n ≥ 128`` and ``α < 5``,
`reduced=true` returns a prefix of the rule, stopping when weights fall below
`10floatmin(Float64)`; other default paths return the full rule.
Specify the parameter explicitly to use this keyword, for example
`gausslaguerre(1000, 0.0; reduced=true)`.
The iterative path avoids computing the far tail, using positivity and exact
polynomial moments to bound the
omitted weights before stopping. This cutoff can exclude representable normal
and subnormal weights. It is not an error tolerance for the integral: the
effect of truncation also depends on the integrand.

Tiny weights in full rules may underflow to zero. Large parameters can also be
limited by overflow in the normalization or intermediate computations.

```@docs
gausslaguerre(n::Integer)
gausslaguerre(n::Integer, α::Real)
```

## Gauss-Chebyshev quadrature

There are four weight functions on ``[-1,1]``:

* First kind: ``1/\sqrt{1-x^2}``.
* Second kind: ``\sqrt{1-x^2}``.
* Third kind: ``\sqrt{(1+x)/(1-x)}``.
* Fourth kind: ``\sqrt{(1-x)/(1+x)}``.

All four use explicit trigonometric formulas in ``O(n)`` work; see
[DLMF §18.5(i)](https://dlmf.nist.gov/18.5#i).

```@docs
gausschebyshevt(n::Integer)
gausschebyshevu(n::Integer)
gausschebyshevv(n::Integer)
gausschebyshevw(n::Integer)
```

## Gauss-Jacobi quadrature

Gauss quadrature for ``w(x)=(1-x)^\alpha(1+x)^\beta`` on ``[-1,1]``, with
``\alpha,\beta>-1``. Legendre and Chebyshev parameter pairs use their dedicated
rules. Large `Float64` rules use the fixed-parameter asymptotic expansions of
[Hale and Townsend (2013)](https://doi.org/10.1137/120889873) for ``n\ge768``,
``\min(\alpha,\beta)\ge -0.99`` and ``\max(\alpha,\beta)<5``.
The iteration handles parameters closer to ``-1``, where the asymptotic
endpoint weights lose relative accuracy. Recurrence handles parameters
``-0.5\le\alpha,\beta\le2`` through ``n=24`` for `Float64`, and through
``n=20`` for other working types, including `Float32` and `BigFloat`. Other cases use
fourth-order iteration of
[Gil, Segura and Temme](https://arxiv.org/html/2509.16716v2#S2.SS3).
The cutoffs retain recurrence for small rules and avoid the setup cost of
asymptotic coefficients at intermediate orders. The recurrence parameter range
also limits loss of relative accuracy in endpoint weights computed from rounded
nodes. Low-precision requests use the routing of their `Float32` working type.

The iteration evaluates its ODE coefficients, local Taylor series, and root
updates in the requested working type at 24 significand bits and above.
Lower-precision rules compute in `Float32` before conversion. A single local
coordinate measured from ``0`` or ``\pm1`` retains endpoint distances for weight
evaluation without subtracting rounded output nodes. Roundoff can accumulate
during Taylor propagation, so node and weight accuracy depends on order,
parameters, and node location.

The iterative recurrence and normalization form the positive quantities
``\alpha+1`` and ``\beta+1`` before adding them, avoiding cancellation when
both parameters approach ``-1``. Large parameters and parameters close to
``-1`` can challenge floating-point accuracy, and endpoint nodes need not remain
distinct when rounded.

```@docs
gaussjacobi(n::Integer, α::Real, β::Real)
```

## Gauss-Radau quadrature

The Legendre–Radau rule has weight ``w(x)=1`` and includes the prescribed node
``-1``. The method `gaussradau(n, α, β)` also supports the Jacobi weight. The
rules use explicit small-order formulas or a Jacobi reduction, without an
eigensolver. Interior weights are transformed using retained endpoint distances
or angles where the Jacobi kernel supplies them; see
[Gil, Segura and Temme, §5.1](https://arxiv.org/html/2509.16716v2#S5.SS1).

```@docs
gaussradau
```

## Gauss-Lobatto quadrature

The Legendre–Lobatto rule has weight ``w(x)=1`` and includes both endpoints
``-1`` and ``1``. It requires ``n≥2``. Interior nodes come from the ``(1,1)``
Jacobi rule, and both endpoint weights are ``2/(n(n-1))``.
For larger orders, the interior weights are transformed using the Jacobi
iteration's endpoint distances or the asymptotic expansion's angles. This
preserves their accuracy near the endpoints when the returned nodes have
already lost some relative accuracy in ``1-x^2``.

```@docs
gausslobatto(n::Integer)
```

## Removed Golub–Welsch helpers

Golub–Welsch implementations and automatic eigensolver fallbacks have been
removed. The public quadrature functions keep their calling conventions and
do not require an external eigensolver.

Code calling `FastGaussQuadrature.gausslaguerre_GW` or `jacobi_gw` must use
`gausslaguerre` or `gaussjacobi`, respectively. Use `gausshermite` for Hermite
quadrature: unlike the removed `hermite_gw` helper, it returns the full rule
with physical weights, rather than the nonnegative half with weights scaled
by `exp(x^2)`.
The internal `radau_gw` helper is also removed; use `gaussradau(n, α, β)`.
The shared `golubwelsch`, `has_eigen`, and `jacobi_jacobimatrix` helpers are no
longer available. There is no replacement helper for constructing arbitrary
rules from Jacobi matrix coefficients. This is a compatibility change for
callers of these non-exported helpers.

## Floating-point types

Select a type with a leading type argument for unparameterized rules. For
parameterized rules, supply parameters of the desired type:

| Rule | Type selection |
|:-----|:---------------|
| Legendre, Hermite, Lobatto, Legendre–Radau | `gausslegendre(T, n)`, `gausshermite(T, n; normalize=false)`, `gausslobatto(T, n)`, `gaussradau(T, n)` |
| Chebyshev | `gausschebyshevt(T, n)` and similarly for `u`, `v`, and `w` |
| Laguerre | `gausslaguerre(T, n)` for ``α=0``; `gausslaguerre(n, T(α))` otherwise |
| Jacobi, Jacobi–Radau | `gaussjacobi(n, T(α), T(β))`, `gaussradau(n, T(α), T(β))` |

Omitting the type in an unparameterized call selects `Float64`. Jacobi and
generalized Radau apply `float` to both parameters and promote them to a common
type. To request `Float32` throughout, for example, pass both parameters as
`Float32`; mixing `Float32` and `Float64` parameters selects `Float64`.
The type of a floating-point Laguerre parameter determines both output types.

### Higher precision

The default rules support `BigFloat` and other compatible `AbstractFloat` types
without an external eigensolver:

```jldoctest
julia> using FastGaussQuadrature

julia> x, w = setprecision(BigFloat, 256) do
           gaussjacobi(32, BigFloat("0.1"), BigFloat("0.2"))
       end;

julia> eltype(x) === eltype(w) === BigFloat
true

julia> all(t -> precision(t) == 256, x) && all(t -> precision(t) == 256, w)
true
```

At fixed working precision and bounded per-node iteration, the default
large-order paths use ``O(n)`` arithmetic work and storage. The cost of
arbitrary-precision arithmetic also depends on the precision chosen by the
caller. Generic iterations for `BigFloat` and other higher-precision types use
the requested precision. The caller may need to increase precision
for difficult orders or parameters. Bounded small-order paths may use
recurrence. The explicit `FastGaussQuadrature.gausslaguerre_rec` helper accepts
orders through 127; larger orders throw `DomainError`. The optional
`recompute=true` mode of `FastGaussQuadrature.gausslaguerre_asy` accepts orders
through 1000. These limits keep quadratic recurrence work bounded. Use
`gausslaguerre` for larger orders and generic precision. Both explicit helpers
require parameters whose floating-point type is `Float64`; other formats throw
`ArgumentError` rather than computing with fixed-precision seeds or coefficients.
The explicit recurrence can lose relative node accuracy, particularly as
``\alpha`` approaches ``-1``. The fixed asymptotic expansion can also lose tail
node accuracy at moderate orders, even when the associated weights underflow.
Automatic selection uses GST to avoid those losses; the manual helpers retain
their specialized accuracy limitations.

Normalization also limits the supported exponent range. In particular, the
Laguerre GST path requires the zeroth moment ``\Gamma(\alpha+1)`` to be finite
in the working type. An overflowing moment can cause rejection even when some
or all individual weights would fit; it is not by itself proof that the output
rule is unrepresentable. Use a type with a wider exponent range in that regime.

### Low precision

The minimum working precision is `Float32` (24 significand bits). For an output
type `T` with fewer than 24 bits, the default rules generally compute in `Float32`,
including any endpoint-weight transformation, then convert the completed nodes
and weights to `T`. Parameterized rules first promote their parameters to a
common output type, then convert them to `Float32` for computation when that
type has fewer than 24 bits. For example, `Float16` requests return
`Vector{Float16}` computed in `Float32`, including Legendre rules with ``n>60``.

Empty rules and exact special cases use direct formulas, including one-point
Legendre and Legendre–Radau, two-point Lobatto, and normalized one-point Hermite.

No arbitrary order cap is imposed for `Float16`. The completed working-precision rule
is converted to the requested output type and validated before returning.
All rules require equal numbers of nodes and weights, finite strictly increasing
nodes, and finite nonnegative weights. A nonempty rule must have at least one
positive weight. Validation throws an error if rounding merges nodes or a weight
overflows. Individual weights may underflow to zero. For example:

```jldoctest
julia> using FastGaussQuadrature

julia> x, w = gausslaguerre(Float16, 32);

julia> eltype(x) === eltype(w) === Float16
true

julia> length(x) == length(w) == 32 && all(diff(x) .> 0)
true

julia> any(iszero, w) && isapprox(sum(Float64, w), 1; atol=eps(Float16))
true
```

The same policy applies to `BigFloat` with fewer than 24 bits: computation uses
`Float32`, and conversion uses the requested `BigFloat` precision. The package
does not change `BigFloat` precision internally. Construct parameters inside
the `setprecision` block so that their values also have the intended precision.
These computations have the precision and exponent range of `Float32`, even
when the output type has a wider exponent range.

### Type requirements

The working `AbstractFloat` type must provide ordered arithmetic and promotion,
conversions from integers and `π`, `eps`, `precision` in bits, `floatmin`,
`floatmax`, `isfinite`, `copysign`, and `muladd`. Depending on the requested rule,
the algorithms use `sqrt`, `exp`, `log`, `log1p`, `sin`, `cos`, `sinpi`, `cospi`,
`sincospi`, `atan`, `atanh`, and `tanh`. The Hermite, Laguerre, Jacobi, and
endpoint-rule algorithms also use special functions such as `gamma`, `loggamma`,
`beta`, and `logbeta` where required by their normalizations or iterations.
These operations must preserve the working type's precision and range.

A custom output type with fewer than 24 significand bits must provide
`precision`, conversions to and from `Float32`, `isfinite`, and the comparisons
with other output values and zero used to validate the rule. Parameterized rules
also require `float`, promotion, and the comparisons used for domain checks.
`Float32` supplies the arithmetic, elementary functions, and special functions
for these computations.

Providing these operations enables the generic methods, but does not guarantee
correctly rounded results or uniform accuracy at every order and parameter.
Conditioning, exponent range, indistinguishable rounded nodes, and iteration
convergence remain limitations.
