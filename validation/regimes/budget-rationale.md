# Finite formula and recurrence criteria

The coefficients below are derived before new broad candidate measurements.
They are engineering allowances for the stated finite regimes, not universal
forward-error theorems. Existing narrower regressions and previously frozen
models remain unchanged. Reference arguments are exact represented inputs;
`e=2^(1-p)`, `u=e/2` and all model arithmetic uses independent high precision.

For Chebyshev I–IV let `q=max(abs(x),1/n)`. Use absolute node allowance `6e*q`,
with an exact zero at the odd centers of I/II. A rational phase formed from
integer numerator/denominator has relative error at most about `3u` including
integer conversion. The code forms small complementary phases directly near
zero, so no subtraction of rounded half-pi angles occurs. Differentiating
sinpi/cospi gives the phase sensitivity; on each selected half interval this
and a `2e` relative elementary-function allowance fit within `6e*q`. The
third-kind half-angle branch has `x=2s²-1`; its square/FMA propagation fits the
same absolute allowance. IV's reflection adds no rounding.

Normal-weight relative allowances are `4e` for I, `12e` for II and the
third-kind half-angle branch. These reserve pi/division, the phase-conditioned
trigonometric value, its square and final multiplication. Elsewhere III uses
`e*(8+6q/(1+x))`, explicitly propagating its returned-node error through the
safe `1+x` factor. The branch guarantees `1+x` is bounded away from zero.
IV uses the reflected III criterion. These small constants are justified by
stable closed forms, not extended to iterative quadratures.

Legendre's explicit orders2:5 use `8e*max(abs(x),1/n)` nodes and `16e` weights;
exact zeros and order1/2 integer weights remain exact. The only subtractive
radicals are bounded fixed constants: the order4 smaller squared node has
subtraction condition below7, and order5 weight subtraction below3. Propagating
unit-roundoff through square roots/divisions fits the stated reserve. Hermite
order1 uses an exact zero and `4e` for sqrt(pi), with exact normalized weight1.
Standard Radau1/2 and Lobatto2/3 use exact endpoint/zero/integer values and
rational division allowances; the existing Lobatto4e endpoint screen remains.

Jacobi order1 uses `4e*abs(x)` for the mean ratio; its positive shifted
denominator avoids cancellation and equal parameters produce exact zero.
The mass uses the unchanged earlier coefficient8 times one plus the absolute
log-gamma/power terms. Laguerre1/2 use `2e*abs(x)`/`8e*abs(x)` nodes. Their
gamma allowance is `8e*(1+abs(loggamma(a+1))+(a+1)*abs(digamma(a+1)))`, including
rounded argument sensitivity. The positive rationalized two-node factors add
16e. Generalized Radau's existing positive-factor model is retained verbatim.
An overflowed mass is classified independently, not passed by a relative bound.

For bounded Jacobi REC, retain the Lobatto audit's `X=8max(abs(x),1/nu)` and
derivative allowance `A=8(n+1)`, where `nu=n+(a+b+1)/2`. The Jacobi ODE gives
raw-weight logarithmic sensitivity
`K=2abs(a-b+(a+b+1)x)/(1-x²)`. Thus
`L=2A+KX+2/(1-x²)+8`; the subtraction in the actual `1-x²` formula is explicit.
Positive-sum normalization adds a global envelope for L and `gamma_(n-1)`.
The endpoints of that envelope follow from
`x>=-1+2(b+1)/(n(n+a+b+1))` and the reflected upper bound, obtained from the
sum of reciprocal endpoint distances. The actual direct gamma-ratio mass
formula contributes its log magnitudes and digamma argument sensitivities.
This extends the prior small-order mechanism to parameters in[-1/2,2]; it
does not enlarge or replace the frozen Lobatto(1,1) consumer criterion.

For the Float64 Laguerre REC route, the implemented Newton trigger is
`40e*x`; the node allowance `48e*x` reserves that trigger and an8e finite
residual/subtraction margin. The orthonormal derivative recurrence has an
operational relative allowance `8(n+1)(1+abs(a))*e`, doubled in the inverse
square. The Laguerre ODE gives weight sensitivity
`abs(2-(2a+1)/x)`, which multiplies the node allowance. Gamma initialization
and final divisions are separate. This is restricted to the actual bounded
automatic REC region through127, plus independently classified ASY recurrence
recomputation through1000. No unbounded recurrence is certified. Failures or
inconclusive allowances require diagnosis, not coefficient inflation.

For every criterion, required normal weights, independently indexed roots,
finite ordered outputs, symmetry/exact centers where applicable, requested
type/precision and restored caller precision are separate gates. Below24bits,
test the completed Float32 calculation and conversion separately. Relative
allowance at least1/8 is inconclusive. Sparse large-order references are labeled
sampled. ASY needs its own published-tail and conditioning context; it is not
silently assigned one of these recurrence budgets.
