# Tests removed from master

The preservation baseline is master `66f752a`. Three routine assertions and
one internal-table docstring example were removed for deleted functionality.
Redundant development-only tests are reviewed separately and are not part of
this master-removal inventory.

## Routine assertions

| File/test context | Exact inputs and checks | Count | Justification |
| --- | --- | ---: | --- |
| `test_gausslaguerre.jl`, n=42 block | `gausslaguerre_GW(42,0.0)`: `abs(x_gw[37]-98.388267163326702)<tol`; `abs(w_gw[7]-0.055372813167092)<tol`; analytic degree-16 polynomial integral `abs(dot(w_gw,polyval(x_gw,coef))-Z)/abs(Z)<tol` | 3 | The documented GW helper is removed. Identical literals, polynomial integral and tolerance remain for public `gausslaguerre` and the explicit recurrence. |

No other original master assertion or inference line is absent. The two master
`jacobi_gw` expressions used only to choose polynomial coefficients were replaced
by `gausslegendre(BigFloat,...)`; the analytic integrals, assertions, degrees,
loops and tolerances are preserved. The master public Hermite testset whose old
label says “Golub-Welsch” remains because it calls public `gausshermite`, not a
removed helper. The BigFloat coefficient test inclusion moved to the end of
`runtests.jl`, so the retained no-GenericLinearAlgebra check also covers it.


## Removed docstring example

There is also **one removed master docstring example with one numerical
comparison**, separate from these routine assertion counts. It was already
absent from the initial uncommitted branch: `src/constants.jl` documented the
unused `BESSELJ1_ON_BESSELJ0_ROOTS` table with setup
`roots = approx_besselroots(0,10)` and the expected-true comparison
`(besselj1.(roots)).^2 ≈ FastGaussQuadrature.BESSELJ1_ON_BESSELJ0_ROOTS`.
Removing the unused `besselJ1` kernel and its table removes that table-specific
doctest. Public Bessel root residual tests remain. The public Bessel docstring's
`all(zeros .< 1e-12)` comparison became the stronger
`all(abs.(zeros) .< 1e-12)`; it was not removed.
