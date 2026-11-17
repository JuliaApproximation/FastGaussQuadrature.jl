# Reproducing the numerical review

These developer checks supplement the routine tests. They use independent classical-polynomial/derivative references, Sturm root counts, increasing reference precision, frozen engineering budgets, and warmed BenchmarkTools measurements. They are deliberately outside the routine CI suite and runtime dependencies. The finite cases support the documented regimes; they do not prove uniform accuracy or convergence for every representable input.

The baseline is master `66f752a`. The comparison uses the exact 60 cases in [cases.jl](cases.jl), including separate full/reduced Laguerre calls. Candidate and baseline run in the same Julia process and dependency environment. Invalid baseline outputs remain in the results, but their timing ratios are not acceptance passes. The nominal cap is candidate/baseline ≤2 for each valid agreed case. The user explicitly accepted the final 2.043554684× reduced Laguerre(10000,-0.99) result; the exact ratio and strict comparison remain in the report.

## Environment and commands

Use Julia 1.13.3 for the recorded tooling run. The package supports Julia 1.10; its routine tests also ran on 1.10.13. Tool versions in [Project.toml](Project.toml) are pinned independently of the runtime package, and some require Julia ≥ 1.12. Install these developer dependencies only in this environment, then keep the package project active with the validation environment on the load path:

```sh
julia --project=validation -e 'using Pkg; Pkg.instantiate()'
JULIA_LOAD_PATH="@:$PWD/validation:@stdlib" julia --threads=32 --project=.
```

In that Julia session:

```julia
using FastGaussQuadrature, Quadmath, LinearAlgebra
BLAS.set_num_threads(1)
include("validation/quality.jl")
QuadratureQuality.run()

# Use an existing checkout of exactly 66f752a; no checkout is downloaded here.
baseline = "/path/to/master-66f752a"
output = mktempdir()
include("validation/measurements.jl")
QuadratureMeasurements.run(baseline, output)

include("validation/accuracy.jl")
master = QuadratureMeasurements.load_snapshot(baseline)
accuracy = ValidationAccuracy.run(output;
    implementations=["master" => master, "candidate" => FastGaussQuadrature])

include("validation/assess.jl")
QuadratureAssessment.assess([accuracy.rulepath], [accuracy.samplepath],
    joinpath(output, "budget_points.csv"))
```

Timing runs need an otherwise idle machine. Float64 uses eight alternating batches of 30 samples, one evaluation per sample, with compilation excluded by warmup. Reported times are pooled medians in nanoseconds; the accompanying table converts to microseconds. Other formats and doubling-order scaling use three batches of 10 samples. Native-format results include actual 113-bit `Quadmath.Float128`, alongside Float32 and BigFloat128/256.

The accuracy driver checks every node through order 128 at ≤53 bits and through order 40 at higher precision; larger rules use fixed endpoint/center/fraction indices and both sides of each normal/subnormal/zero transition. The Float64 parameter literals in `cases.jl` are converted to the output type before lifting into reference arithmetic. References run at p+128 and p+192 bits and verify root indices. Separate focused regime studies cover additional transitions, near-singular parameters and precision challenges. Normal weight loss, invalid rules and unresolved roots are separate from numerical-error budgets.

Run `Pkg.test()` in the package environment for the routine suite. Its default `--depwarn=yes` adds 89 deprecation assertions to direct test-file execution; reported totals state which mode actually ran. Documentation uses the existing `docs/Project.toml` environment and `docs/make.jl`; deployment is not part of local validation.

## Evidence

[The final results report](results/README.md) records the final source hashes, environments, test/tool outcomes, accuracy/timing tables, budgets and limitations. Compressed raw CSVs can be unpacked into a separate directory with `gzip -dc` or Python's standard `gzip` module, preserving the original artifacts. `python3 validation/report.py` verifies the saved timing counts and medians and regenerates the joined tables. Reproduction outputs refuse to overwrite an earlier result. Keep baseline and candidate source identities with every run; branch HEAD alone does not identify an uncommitted source tree.

[Removed tests](removed-tests.md) lists every removed assertion against master and the deleted internal-table doctest, with their justifications.
