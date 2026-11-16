# Final numerical and performance evidence

Baseline: master `66f752a`. Final production source and test identities are in [source-sha256.json](source-sha256.json).

All 60 agreed Float64 cases were measured on the final production source. Of 59 valid baseline comparisons, 58 are at most 2×; reduced Laguerre(10000,−0.99) is **2.043554684×**, explicitly accepted by the user. The baseline Laguerre(39,44.8) output is invalid and is excluded from timing acceptance. Exact measurements remain visible.

The independent accuracy matrix has **185 passes and one proved Float32 representation limit**, covering **4,815 node/weight pairs**: 4,621 normal, 39 subnormal and 155 genuinely underflowed reference weights. All applicable frozen budgets, 28 additional Lobatto checks and eight reduced-prefix boundary checks pass. The largest node/weight budget fractions are 0.5032/0.2380. Larger orders are sampled.

Fresh suites pass **7,121 assertions on both Julia 1.10.13 and 1.13.3**, with all Aqua groups. Runic, strict imports, 23 targeted JET checks, and Documenter/doctests pass. [Quality details](quality-summary.md) include exact versions, source hashes and redacted logs. Developer tools are outside runtime dependencies.

## Selected Float64 comparisons

Node errors are maximum absolute errors; weight errors are maximum relative errors among checked normal weights. Every selected reference is independently indexed and stable at increased precision. Times are warmed pooled medians in microseconds. Full results, ulp diagnostics, allocations and statuses are in [the 60-case CSV](float64-accuracy-timing.csv).

| Case | Node error, master → candidate | Normal weight error, master → candidate | Master µs | Candidate µs | Candidate/master |
| --- | --- | --- | ---: | ---: | ---: |
| `legendre_6` | 3.21e-17 → 3.21e-17 | 5.8e-16 → 5.8e-16 | 0.209 | 0.209 | 1.000× |
| `legendre_10000` | 2.68e-16 → 9.98e-17 | 4.66e-16 → 2.4e-16 | 101.375 | 114.834 | 1.133× |
| `hermite_201` | 7.69e-15 → 9.52e-16 | 1.02e-11 → 3.24e-14 | 282.562 | 22.709 | 0.080× |
| `hermite_10000` | 6.74e-14 → 1.13e-14 | 4.64e-12 → 6.53e-14 | 7140.479 | 942.979 | 0.132× |
| `jacobi_25_quarter_zero` | 8.24e-17 → 1.19e-16 | 2.15e-14 → 5.09e-15 | 7.708 | 11.791 | 1.530× |
| `jacobi_4096_quarter_zero` | 1.45e-16 → 1.45e-16 | 8.56e-15 → 8.56e-15 | 11248.292 | 11291.521 | 1.004× |
| `laguerre_31_point2_full` | 6.28e-15 → 2.73e-14 | 3.02e-13 → 3.15e-14 | 28.041 | 10.166 | 0.363× |
| `laguerre_10000_zero_full` | 3.55e-12 → 1.14e-11 | 8.07e-13 → 1.48e-13 | 4956.021 | 2273.979 | 0.459× |
| `laguerre_10000_near_minus_one_reduced` | 1.01e-12 → 1.74e-13 | 1.41e-12 → 1.72e-13 | 225.292 | 460.396 | 2.044× |
| `radau_jacobi_100_large_parameters` | 6.6e-16 → 1.15e-16 | 1.39e-12 → 1.74e-14 | 282.646 | 41.417 | 0.147× |
| `lobatto_10000` | 1.75e-16 → 1.75e-16 | 8.97e-10 → 1.05e-15 | 14547.729 | 8026.292 | 0.552× |

## Other formats and scaling

[Native-format tables](native-accuracy-timing.csv) retain Float32, actual 113-bit Quadmath.Float128, and BigFloat at 128/256 bits, including allocations and all available matching accuracy rows. Native timings precede the final Float64-only selection repairs; their algorithms are unchanged by the reviewed Float64-only source changes. [Complete-output comparisons](source-output-bridge.csv) additionally verify the eight overlapping accuracy/timing inputs per format; they do not replay all 16 timed inputs per format. Missing baseline methods and invalid outputs are retained in the raw status files. A time without a matching reference row is not an accuracy pass.

[Doubling-order measurements](scaling.csv) cover 14 family/regime/semantics combinations at n=2048,4096,8192,16384. Full-rule timings support the reviewed linear-work implementation; individual doubling ratios still vary with allocation and runtime effects. Reduced Laguerre returns a growing prefix and is separate work. These measurements supplement bounded loops and call-chain review, not a proof of all-input convergence.

## Reproduction and limits

See [the reproducible drivers](../README.md), [frozen model explanations](../regimes/README.md), and [removed master tests](../removed-tests.md). `python3 validation/report.py` verifies raw timing medians/counts and regenerates these joined tables. Raw CSVs are gzip-compressed; no numerical values are redacted. Environment prose replaces local home/check-out paths with placeholders.

Float64 uses eight alternating batches ×30 samples with compilation excluded; other formats/scaling use three ×10. Hardware: Apple M3 Ultra, Julia 1.13.3, 32 Julia threads and one BLAS thread. Baseline and candidate share the same runtime/dependencies within each paired measurement. The native run was continued after a 300-second tool timeout, preserving completed cases; the final Float64/scaling refresh completed all observations before a timestamp-only logger error. Both events and raw completeness checks are retained in the environment records.

Accuracy means good stable native-precision behavior with justified order/parameter/location growth, not a uniform ulp or correctly-rounded promise. The Float32 floor, normalization exponent range, distinct-node resolution and reduced cutoff remain explicit limits. The Float32 Laguerre(39,44.8) exception has an independent mass/n overflow certificate. A normalization overflow alone is not such a certificate. Manual recurrence/asymptotic helpers retain documented limits; automatic Laguerre now uses GST to avoid their observed losses.

The exact four agreed reduced-rule boundaries are in [reduced-boundaries.csv](reduced-boundaries.csv): the last retained and first omitted references are on opposite sides of 10floatmin(Float64). Prefix nodes match exactly; prefix weights agree within the preserved 8eps allowance (one n100000 half-parameter case differs by 2.22e−16 relatively).

Earlier failing measurements are preserved in the local study archive, with selected repair records here. The final fixes use existing native formulas and unchanged budgets: Hermite physical-weight residual transport, generalized Radau near−1 routing, and removal of recurrence/full asymptotic selection from automatic Laguerre. No higher working precision or guard arithmetic was introduced.
