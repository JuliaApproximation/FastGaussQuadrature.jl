# Native-precision accuracy models

These are finite engineering criteria, supported by the operation structure,
conditioning, independent references and reserved challenge cases. They are not
uniform forward-error theorems or correctly-rounded guarantees. The working
format has p significand bits and e=2^(1-p). Node allowances are absolute;
normal-weight allowances are relative. Their coefficients were frozen before
reserved validation, and were not enlarged to turn failures into passes.

The integrated assessor uses exactly represented input parameters. The committed
case definitions contain Float64 literals; each is converted to the requested
format before being lifted into reference arithmetic. Its classifier is scoped
to that matrix, not every possible public call. Separate regime studies cover
additional parameters, transitions and precision challenges.

| Regime | Model and numerical basis |
| --- | --- |
| Chebyshev I–IV and explicit small rules | [formulas_rec.jl](formulas_rec.jl), [rationale](budget-rationale.md): phase/radical/rational rounding, stable endpoint factors and mass sensitivity. Exact centers/endpoints retain exact checks. |
| Legendre recurrence, n=6:60 | Absolute node allowance `4e*max(abs(x),1/n)`; relative weight allowance `e*(8n+8+4κ*max(abs(x),1/n))`, `κ=2abs(x)/(1-x²)`. The degree term covers recurrence accumulation and the location term covers endpoint sensitivity. |
| Legendre expansion, 24–53 bits | [Frozen model](legendre-asy.json): absolute node allowance `4e*max(abs(x),1/n)`, relative weight allowance `8e`. The precision-dependent omission rule has a separate finite coefficient bound. The fixed table ceiling is retained; above 53 bits the Jacobi GST model applies. |
| Hermite recurrence, n=2:20 | [hermite_rec.jl](hermite_rec.jl), [frozen criterion](hermite-rec.json): recurrence growth, root sensitivity, physical/scaled weights and positive-sum normalization; normalized outputs add their coordinate/weight scaling. |
| Hermite GST, n>20 | [hermite_gst.jl](hermite_gst.jl), [rationale](hermite-gst.md): coordinate accumulation follows outward root count and the normal-form coefficient `2n+1-x²`; transformed-derivative error is stationary at a root. Physical exponential sensitivity and normalization are separate. The native residual correction improves physical weights without a larger allowance. |
| Jacobi GST | [jacobi_gst.jl](jacobi_gst.jl), [criterion](jacobi-gst.json): native coordinate error, accumulated derivative error, endpoint series, logarithmic factors, mass and positive-sum normalization. Radau and Lobatto preserve the accurate interior geometry when transforming weights. |
| Bounded Jacobi recurrence and endpoint consumers | [formulas_rec.jl](formulas_rec.jl) and [assess.jl](../assess.jl): recurrence/derivative conditioning, finite normalization and endpoint divisor sensitivity. Standard Radau's returned-node division has explicit finite-normal denominator/numerator guards. |
| Laguerre GST | [laguerre_gst.jl](laguerre_gst.jl), [rationale](laguerre-gst.md): square-root coordinates, derivative accumulation, physical/logarithmic factors and normalization. Near−1 uses independently evaluated terminating-series sensitivity and a second-root anchor. All automatic orders ≥3 now use this kernel, with the half-parameter exceptions below. |
| Jacobi ASY and Radau/Lobatto consumers | [jacobi_asy_budget.jl](jacobi_asy_budget.jl), [rationale](jacobi-asy-rationale.md): independent next boundary coefficients, Hahn tail envelopes, phase/amplitude rounding and derivative/divisor conditioning. The earlier Lobatto `8e` node/`32e` weight diagnostic is an additional condition, never replaced by the larger model. |
| Laguerre half-parameter Bessel replacements | [laguerre_asy_budget.jl](laguerre_asy_budget.jl), [rationale](laguerre-asy-rationale.md): independent constants, retained/next terms, scalar rounding and physical-weight sensitivity. Only the first 2 or 12 Float64 nodes use these formulas automatically. The trace records formula/depth and native seed metadata; measured errors never set allowances. |

The frozen `budget-rationale.md` predates the final Laguerre selector repair;
its phrase “automatic REC region through 127” describes the earlier selector.
That model now applies only to the manually requested bounded recurrence.

The manual Laguerre recurrence and complete asymptotic helper have separate
finite models and documented limits. The review found recurrence cancellation
and a fixed-depth Airy tail floor, so neither selects automatic full rules.
The stored asymptotic coefficients remain useful for half-parameter endpoints
and the explicitly requested manual method. No new kernel or wider working
arithmetic was introduced to hide those losses.

Allowances reaching 1/8 are generally inconclusive, rather than passes. A
relative tail estimate for a genuinely underflowed weight is instead evaluated
in absolute terms. Subnormal/rounded-zero weights use their reference weight
times the applicable relative allowance plus two minimum spacings. Structural
checks and loss of a required normal weight are independent failure gates.
Reduced-prefix boundary checks are recorded separately, including the first
omitted weight; a short returned vector alone does not establish correct
reduction semantics. The Float32 calculation and final conversion below 24 bits
are separate controls.

## Independent evidence and changes prompted by it

The final integrated matrix has 185 numeric passes and one independently proved
Float32 overflow case, with 4,815 separately assessed node/weight pairs. These
are 4,621 normal, 39 subnormal and 155 true-underflow reference weights. The
maximum node/weight fractions are 0.5032/0.2380; no assessed result is
inconclusive. Larger orders are sampled, not exhaustively certified.

Focused studies additionally cover the 20/21 Hermite boundary, 53/54-bit
Legendre transition, algorithm splices, parameters near−1, large parameters,
low-precision floors, actual 113-bit Float128, and BigFloat at multiple precisions.
Frozen calibration/challenge cases and post-repair replays are distinguished
in the preserved study reports. An examined case is not relabeled an unseen
holdout after a repair.

The review retained and repaired three kinds of conspicuous loss: small-order
Laguerre recurrence cancellation, moderate-order Laguerre asymptotic tail-node
error, and generalized Radau's near−1 upper weight. Native GST was an existing,
comparably simple alternative. The final selection changes retain the frozen
budgets and independent references. The retained Jacobi ASY scope has 21 passing
calibration/challenge rules; the old near−1 Radau failure is preserved separately,
with its repaired points and three new neighboring cases assessed by the
unchanged GST model. Hermite's 180 physical/normalized/scaled views pass its
unchanged model before and after the residual correction.

Precision and exponent-range qualifications still matter. In particular,
normalization may require a representable total mass even when individual
weights could fit. This is an algorithmic capability limit; only an independent
`mass/n > floatmax(T)` calculation proves unavoidable output overflow. The
reported Float32 Laguerre(39,44.8) exception satisfies that stronger condition.
