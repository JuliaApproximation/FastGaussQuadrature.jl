# Final validation and environment

Final production and test sources passed **7,121 routine assertions on both Julia 1.10.13 and Julia 1.13.3**, in fresh isolated Julia processes. All eight Aqua groups passed on both versions. The 701 recorded testset-count rows agree. Source/test SHA-256 hashes and exact resolved environment versions are in `quality-summary.json`.

| Check | Actual result |
|---|---|
| Runic 1.13.0 | All source/test files checked; final changed files rechecked, remaining files byte-identical to the checked snapshots. |
| ExplicitImports 1.15.0 | All four strict import checks pass. |
| JET 0.13.1 | 23 targeted calls, zero findings, including direct bounded REC/ASY helpers and the repaired Radau endpoint route. |
| CodeComplexity 0.2.0 | Complete production source: cognitive 835, cyclomatic 646. These are diagnostics, not correctness or complexity proofs. |
| Documenter 1.20.0 | Final documentation build and doctests pass; no deployment. |

The routine suites were executed by including `test/runtests.jl` through Julia MCP. `Pkg.test` enables `--depwarn=yes`, which activates 89 additional existing deprecation assertions. The final snapshot was not rerun in that execution mode, so its predicted total is not presented as a completed test result.

## Environment

Hardware: Apple M3 Ultra, arm64 macOS, 32 Julia threads. Quality checks used 24 BLAS threads; the separate timing protocol explicitly used one BLAS thread. The package remained the active project, with a separately resolved validation-tool environment on the effective load path. No personal tool environment is needed to reproduce the declared setup.

| Component | Julia 1.10.13 tests | Julia 1.13.3 tests / tools |
|---|---|---|
| SpecialFunctions | 2.9.2 | 2.9.2 |
| StaticArrays | 1.9.24 | 1.9.24 |
| StaticArraysCore | 1.4.5 | 1.4.5 |
| PrecompileTools | 1.2.1 | 1.4.0 |
| Aqua | 0.8.18 | 0.8.18 |

The validation environment pins Quadmath **1.1.0** (actual Float128, 113-bit significand), BenchmarkTools **1.8.0**, Runic **1.13.0**, JET **0.13.1**, ExplicitImports **1.15.0**, CodeComplexity **0.2.0**, Documenter **1.20.0**, Aqua **0.8.18**, Revise **3.18.0**, and TestEnv **1.103.7**. Its tools require Julia 1.12 or later; that is separate from the package's retained Julia 1.10 compatibility.

Runtime dependencies remain LinearAlgebra, SpecialFunctions and StaticArrays. Test extras, documentation dependencies and validation tools are declared separately. The Julia 1.10 environment uses its compatible PrecompileTools 1.2.1; the newer runtime environment was preserved. Exact transitive versions and manifest hashes are recorded in JSON.

## Tests removed from pinned master

Exactly three routine assertions were removed, all for the deliberately removed `gausslaguerre_GW(42,0.0)` helper:

1. Node 37 compared with `98.388267163326702` at the original `tol`.
2. Weight 7 compared with `0.055372813167092` at the original `tol`.
3. The analytic degree-16 polynomial integral compared at the original relative `tol`.

The same literals, integral, and tolerances remain for the public Laguerre rule and explicit recurrence. No other master assertion or inference line is removed. Two `jacobi_gw` expressions used only to choose polynomial coefficients were adapted to public BigFloat Legendre rules while retaining all independent integral checks.

One additional master **docstring example with one numerical comparison** was removed with the unused `besselJ1` helper and `BESSELJ1_ON_BESSELJ0_ROOTS` table: after `roots = approx_besselroots(0,10)`, it checked `(besselj1.(roots)).^2 ≈ FastGaussQuadrature.BESSELJ1_ON_BESSELJ0_ROOTS`. Public Bessel root residual tests remain.

## Logs

- `tests-julia-1.10.13.log` and `tests-julia-1.13.3.log`: final fresh-process suites and Aqua.
- `quality-julia-1.13.3.log`: imports, formatting, all 23 JET calls and complete per-function complexity diagnostics.
- `documentation.log`: final documentation build/doctest completion.

Published logs replace local directory prefixes with placeholders; results are unchanged. Both original and redacted log hashes are recorded. Numerical accuracy and timing acceptance are separate evidence, not consequences of these tooling checks.
