# Jacobi fixed-depth ASY finite model

This is a reference-only engineering model for the existing Float64 ASY route,
not a new production algorithm or a general asymptotic remainder theorem.
The public domain is degree>=768, represented parameters in[-0.99,5); the
standard Radau consumer also uses degree>=101 with(a,b)=(0,1). Generalized
Radau uses degree>=500 with -1<a<5 and 0<b<5, where b is the mathematical
beta+1 shift; this consumer does not inherit Jacobi's -0.99 lower restriction.
Lobatto uses degree>=768,(1,1). No other precision is assigned this fixed-depth model.
The first/last ten Jacobi roots use the boundary chart; all remaining roots
use the interior chart. The index passed to the model excludes fixed consumer
endpoints. Its weight argument includes the consumer transform.

The primary source is Hale–Townsend2013, section3.3, equations(3.22)–(3.24),
the derivative identity immediately following(3.24), and section3.3.1:
https://appliedmaths.sun.ac.za/~nhale/publications/HaleTownsend2013a.pdf.
Only published formulas are used; no external implementation code is copied.
The source retains20 Hahn terms and B0/B1/B2,A1/A2/A3 boundary terms.

## Truncation before observing errors

For the interior series, take the absolute coefficient envelope of each Hahn
degree m, including every binomial cross-term. Its two factors are the
Pochhammer products(.5±a)_l and(.5±b)_(m-l), divided by the corresponding
factorials, powers of2sin(theta/2),2cos(theta/2), and(2nu+1)_m.
The implementation can stop before adding degree11 when a signed term is
small. Therefore this conservative model covers all potentially omitted
degrees11:20, then uses the degree21 envelope and the22/21 ratio. It also
covers incomplete cross-terms in production degrees18/19. The ratio must be
below1/2. The factor2 is consistent with the known Hahn bound in its restricted
parameter range; outside that range it is explicitly a finite next-term
engineering assumption, to be challenged, not an extension of that theorem.
The source's stopping uses a vector maximum, so this pointwise allowance may
include terms actually retained. No candidate error selects a truncation depth.

For the boundary expansion, derive the next B3 and A4 coefficient functions
independently from the published recursion. Write q=theta² and f=g'/2, where
g=(1/4-a²)(cot(theta/2)-2/theta)-(1/4-b²)tan(theta/2).
Its exact Bernoulli series is analytic inside|theta|<pi. For coefficient k,

```
B_m,k = -(k+1)[1+(1+2a)/(2k+1)] A_m,k+1
        +(f*A_m)_k/[2(2k+1)]
A_m+1,k = (k-a)B_m,k -(f*B_m)_k-1/(4k), k>=1
A_m+1,0 = 0, A_0 = 1.
```

The helper uses33 exact rational Bernoulli factors, retaining29 coefficients
after four recurrences. It evaluates absolute power-series norms, with a
factor2 analytic-tail reserve on theta<=1/2. The first omitted terms are
theta*B3/nu^7 and A4/nu^8. The odd term cannot be discarded: it has inverse
eighth-order scale only because the selected boundary roots have theta=O(1/nu).
Eight times their absolute norm, a finite Bessel-amplitude reserve1+|a|,
and the geometric reserve with ratio(1+|a|+|b|)^2/nu² form the boundary
truncation estimate. This is a finite next-term estimate, not a proved uniform
remainder bound. The exact-series construction avoids merely guessing a
parameter polynomial D^4 for the unknown coefficient.

The production boundary coefficient functions use degree-nine Chebyshev
collocation. The finite rounding allowance below also reserves its numerical
differentiation/integration/interpolation error on theta<=1/2; the power-series
tail alone does not certify that separate coefficient construction. Smooth
analytic coefficients, constant interpolation degree, and the nu^-3 and
higher multipliers make this a bounded engineering assumption. The reserved
tests must include first/last roots and both sides of the ten-root splice.
Any failure requires diagnosis rather than adjusting these coefficients.

## Independent sensitivity transfer

Let K be the exact gamma-ratio weight constant. Reconstruct the base Jacobi
weight by multiplying the consumer weight by1,1+x,or1-x². Then
|P_theta'|=sqrt(K/wJacobi). Multiply the interior or boundary series tail by
its exact high-precision prefactor to obtain absolute polynomial errors for
degree n and n-1. The published derivative identity propagates both errors
to an absolute derivative allowance; division by|P_theta'| gives a relative
derivative contribution. Node displacement is the polynomial error divided
by|P_theta'|. Thus derivative cancellation is exposed rather than treating
a polynomial-value tail as a relative weight bound.

The native phase arithmetic contributes an angle reserve
e[8theta+8(1+D+chi)/nu], D=1+a²+b²,
chi=D/(nu*sin(theta)^2) in the interior and D/(nu²*sin(theta)^2)
in the boundary chart. Direct cancellation in g' contributes O(e/theta²),
then its A1/nu² factor contributes O(e/(nu*theta)²); the correlated Bessel
shift avoids the additional interior derivative-subtraction loss. The coherent phase perturbation moves a root;
the chi term describes the remaining angle-derivative cancellation.
Amplitude rounding adds8e(1+D+chi+L), where L is the sum of absolute
logarithmic power-prefactor terms, parameter magnitudes and log(nu).
Eight reserves the few operations in each stable scalar factor and bounded
coefficient construction; this is an operational envelope, not a claim of
correct rounding of the elementary/special functions. It is not a universal
C*n allowance: its growth follows the actual chart and local conditioning.

The Jacobi ODE gives |d log|P_theta'||/dtheta| at a root as
|a-b+(a+b+1)x|/sin(theta). Multiply this by the angle allowance. The Radau
and Lobatto divisors add sensitivities sin(theta)/(1+x) and2|cot(theta)|,
plus4e native trigonometric operations. Returned node error includes cosine
rounding and the finite second-order angle term; exact symmetric centers
remain exact. The weight constant gets a separate native gamma-ratio/sum
reserve and the first omitted degree21 hypergeometric term. A finite inverse
square/divisor quotient combines these errors. Require all tail ratios<1/2,
derivative and divisor errors<1/8, nu*angle<1/8, and final relative allowance
<1/8; otherwise the result is inconclusive. Structural/representation gates
remain separate.

The earlier Lobatto8e absolute-node/32e normal-weight diagnostic is returned
as an additional independent condition. This model must not loosen it.

## Symbolic check and limitations

Wolfram's pure mathematical evaluator session z9l8BuBj independently returned
`{0,0,0}` for the A1 constant-at-zero identity and the two coefficient-factor
identities in the recurrence above. The Limit evaluation emitted removable
singularity warnings before returning0; these warnings are recorded rather
than called a clean numerical evaluation. No file/network/shell operation was
requested. This algebra check supplements numerical references and does not
establish the engineering constants or asymptotic tail assumptions.

The model deliberately states the remaining finite assumptions. It cannot be
used to certify untested parameters, arbitrary orders/precisions, interpolant
accuracy, or correctly rounded weights from source formulas alone. Frozen
calibration and reserved challenge evidence is required before reporting this
regime as assessed.

Pre-execution revision2: independent mathematical review removed an unnecessary
nu factor from boundary chi. The original frozen proposal is preserved in
`jacobi-asy-initial-model/`; no candidate errors had been assessed. This
correction changes neither the grid nor production. The extreme a+1≈eps
case can still be inconclusive and must be investigated independently.
