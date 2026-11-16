# Hermite GST finite engineering model v1

Scope: n>20 on the frozen finite grid, native precision p>=24, physical,
standard-normal and scaled output views. The n<=20 recurrence has its separate
unchanged frozen model. No candidate errors have been inspected for this model.
The constants below are operational rounding reserves, not a global theorem.

Write e=eps(T), u=e/2, lambda=2n+1, A_i=lambda-x_i^2, R=max|x_i|,
Amin=lambda-R^2 and k_i=ceil(|i-(n+1)/2|), the number of outward root steps.
All x_i and the extreme R are independently indexed reference roots. Set

```
X_i = 8|x_i| + 8k_i/sqrt(A_i), d_i=e X_i
Q_i = 2(A_i+2|x_i|d_i+d_i^2)d_i^2
L_i = 16(k_i+1)+8+2|x_i|X_i+x_i^2     (physical)
Ls_i = 16(k_i+1)+8                     (scaled)
```

Coordinate rounding and the fourth-order stop contribute a relative root
scale; accumulated phase errors are converted to coordinate error through the
local frequency sqrt(A_i). The 8 reserves several native operations per
transport/iteration and the stopped Taylor tail. Derivative amplitude is
propagated across k steps and squared inversely, giving 16(k+1). This is an
explicit accumulation model with location dependence, not a universal n*eps
claim. Its coefficients must be challenged without tuning to observed errors.

For y=exp(-x^2/2)H_n, y''+A y=0, so y'' vanishes at a root and the scaled
weight 1/y'^2 is stationary. Q_i reserves the second-order displacement term
in a small local phase tube. Physical weights add log sensitivity -2x and
the rounded exponent argument x^2. The odd center must be exactly zero.
The existing iteration's leading relative root error is -step^4/6; its
threshold eps^(1/4)/2 puts this leading term below eps/96. This accounts only
for local iteration truncation, not accumulated amplitude error.

The Taylor recurrence has four lags, stops after four consecutive sufficiently
small terms, and caps work by max(100,2p). This remains native arithmetic.
Its truncation/accumulation reserve is included above; validation includes
precision challenges and operation counts. A cap hit is a failure, not a
budget pass. Kahan summation is the unchanged local native sum algorithm,
not a two-part root coordinate or extended arithmetic type.

Normalize using physical probabilities p_i=w_i/sqrt(pi). Sturm comparison for
A<=lambda gives k_i<=1+sqrt(lambda)|x_i|/pi. Using the exact moments
E[x^2]=1/2 and E|x|<=1/sqrt(2), the weighted mean of physical L_i is bounded by

```
Kbar = 1+sqrt(lambda/2)/pi
Lbar = 16(Kbar+1)+8+8
       +16(sqrt(lambda)/(2pi)+1/sqrt(2))/sqrt(Amin)+1/2
```

For the stationary remainder use the conservative common bound
D=e(8R+8floor(n/2)/sqrt(Amin)) and
Qbar=2(lambda+2RD+D^2)D^2. Set rho=e Lbar+Qbar,
gamma=(n-1)u/(1-(n-1)u), eta=rho+gamma(1+rho). The normal-weight relative
allowance for either view is
`(1+e*L_view+Q_i)*(1+4e)/(1-eta)-1`.
This finite quotient covers local error, a positive normalization sum, mass
constant and final operations. Require eta<1/8, local phase<1/8 and allowance
<1/8. The node absolute allowance is d_i. Normalize=true uses the earlier
finite sqrt(2) node transfer with delta=2u+u^2 and adds2e for final weight
division. All arithmetic in this assessor is reference-only high precision.

For subnormal or truly rounded-zero reference weights, use the absolute
allowance `w_ref*relative_allowance + 2*minimum_subnormal_spacing(T)` and retain
the measured spacing error. No required normal weight may be lost. Full vector
ordering, root counts, symmetry, output/caller precision and finiteness are
independent acceptance conditions. Inconclusive models are not passes.

Primary basis: GST2019 section3, https://arxiv.org/abs/1906.05414, and
GST2025v2 section4.1, https://arxiv.org/html/2509.16716v2#S4.SS1.
These supply the normal-form iteration and Taylor propagation. The finite
rounding/error envelope above is derived here, not copied from authors' code.
The Sturm comparison and moment inequalities are mathematical bounds; the
8/16 operation reserves remain finite engineering assumptions.

Grid: fixed full indices through n128 and fixed endpoint/central/fraction
indices above. Additional indices at output representability transitions are
locator diagnostics, chosen without reference error magnitudes and logged.
They can detect a missing tail but do not turn sampled coverage into exhaustive
large-order certification. Calibration and reserved challenge are separate.
