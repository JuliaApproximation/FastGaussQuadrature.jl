"""Regenerate review tables from the committed raw numerical/timing evidence."""
from pathlib import Path
import csv
import gzip
import json
import statistics

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "results"


def read(name):
    path = DATA / name
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as stream:
        return list(csv.DictReader(stream))


def write(name, rows):
    with (DATA / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def verify_timings(prefix):
    raw, summaries = read(prefix + "_raw.csv.gz"), read(prefix + "_summary.csv")
    groups = {}
    for row in raw:
        key = (row["id"], row["implementation"], row["type"], row["bits"])
        groups.setdefault(key, []).append(row)
    assert len(groups) == len(summaries)
    for summary in summaries:
        key = tuple(summary[k] for k in ("id", "implementation", "type", "bits"))
        rows = groups[key]
        times = [float(row["time_ns"]) for row in rows]
        assert len(times) == int(summary["samples"])
        assert statistics.median(times) == float(summary["median_ns"])
        assert min(times) == float(summary["minimum_ns"])
    return summaries


rules = read("accuracy_rules.csv.gz")
points = read("accuracy_samples.csv.gz")
budgets = read("budget_rules.csv")
assert len(budgets) == 186
assert sum(r["status"] == "pass" for r in budgets) == 185
assert sum(r["status"] == "proved_representation_limit" for r in budgets) == 1
assert len({(r["format"], r["case_id"]) for r in budgets}) == 186
budget_points = read("budget_points.csv.gz")
assert len(budget_points) == 4815 and all(r["status"] == "pass" for r in budget_points)
lookup = {(r["implementation"], r["format"], r["case_id"]): r for r in rules}
point_groups = {}
for row in points:
    point_groups.setdefault((row["implementation"], row["format"], row["case_id"]), []).append(row)


def metrics(implementation, format_name, case_id):
    rows = [r for r in point_groups.get((implementation, format_name, case_id), [])
            if r["sample_status"] == "ok" and r["reference_stable"] == "true" and r["root_index_ok"] == "true"]
    normal = [r for r in rows if r["reference_weight_class"] == "normal"]
    def maximum(rs, key):
        return max((float(r[key]) for r in rs), default=float("nan"))
    return {"node_absolute": maximum(rows, "node_abs"), "node_ulp": maximum(rows, "node_ulp"),
            "normal_weight_relative": maximum(normal, "weight_rel"), "normal_weight_ulp": maximum(normal, "weight_ulp")}


summary = verify_timings("float64")
times = {(r["id"], r["implementation"]): r for r in summary}
assert len(times) == 120
acceptance = json.loads((DATA / "acceptance.json").read_text())
rows = []
for candidate in (r for r in summary if r["implementation"] == "candidate"):
    name = candidate["id"]
    baseline = times[name, "master"]
    a, b = lookup["master", "Float64", name], lookup["candidate", "Float64", name]
    ratio = float(candidate["median_ns"]) / float(baseline["median_ns"])
    if a["status"] != "ok":
        status = "invalid_baseline"
    elif ratio <= 2:
        status = "within_nominal_target"
    else:
        assert name == acceptance["accepted_case"]["case_id"]
        assert ratio == acceptance["accepted_case"]["measured_ratio"]
        status = "explicitly_accepted_2.044x"
    row = {"case_id": name, "family": b["family"], "n": b["n"], "alpha": b["a"], "beta": b["b"],
           "reduced": b["reduced"], "normalized": b["normalize"], "baseline_status": a["status"],
           "candidate_status": b["status"], "baseline_microseconds": float(baseline["median_ns"]) / 1000,
           "candidate_microseconds": float(candidate["median_ns"]) / 1000, "candidate_over_baseline": ratio,
           "acceptance": status, "candidate_bytes": candidate["bytes"], "candidate_allocations": candidate["allocations"]}
    for implementation in ("master", "candidate"):
        row.update({implementation + "_" + k: v for k, v in metrics(implementation, "Float64", name).items()})
    rows.append(row)
write("float64-accuracy-timing.csv", rows)

# Native timings share eight cases/format with the independent accuracy grid;
# other timing rows remain available without inventing matching accuracy data.
native = []
for tag, format_name in (("float32", "Float32"), ("float128", "Float128"), ("bigfloat128", "BigFloat128"), ("bigfloat256", "BigFloat256")):
    summaries = verify_timings("native_" + tag)
    timed = {(r["id"], r["implementation"]): r for r in summaries}
    for row in summaries:
        family, n, a, b, reduced, normalized = row["id"].rsplit("_", 5)
        matching = [r for r in rules if r["implementation"] == row["implementation"] and r["format"] == format_name
                    and (r["family"], int(r["n"]), float(r["a"]), float(r["b"]), r["reduced"], r["normalize"])
                    == (family, int(n), float(a), float(b), reduced, normalized)]
        result = {"format": format_name, "case_id": row["id"], "implementation": row["implementation"],
                  "microseconds": float(row["median_ns"]) / 1000, "bytes": row["bytes"], "allocations": row["allocations"],
                  "accuracy_case": matching[0]["case_id"] if matching else "not_in_joint_accuracy_grid"}
        if matching:
            result.update(metrics(row["implementation"], format_name, matching[0]["case_id"]))
        else:
            result.update({k: "" for k in ("node_absolute", "node_ulp", "normal_weight_relative", "normal_weight_ulp")})
        native.append(result)
write("native-accuracy-timing.csv", native)

scaling = verify_timings("scaling")
scaling_rows = []
for r in scaling:
    family, n, a, b, reduced, normalized = r["id"].rsplit("_", 5)
    scaling_rows.append({"family": family, "n": int(n), "alpha": a, "beta": b, "reduced": reduced,
                         "nanoseconds": float(r["median_ns"]), "nanoseconds_per_requested_node": float(r["median_ns"]) / int(n),
                         "bytes": r["bytes"], "allocations": r["allocations"]})
write("scaling.csv", scaling_rows)

select = ["legendre_6", "legendre_10000", "hermite_201", "hermite_10000", "jacobi_25_quarter_zero",
          "jacobi_4096_quarter_zero", "laguerre_31_point2_full", "laguerre_10000_zero_full",
          "laguerre_10000_near_minus_one_reduced", "radau_jacobi_100_large_parameters", "lobatto_10000"]
byid = {r["case_id"]: r for r in rows}
lines = ["# Final numerical and performance evidence", "", "Baseline: master `66f752a`. Final production source and test identities are in [source-sha256.json](source-sha256.json).", "",
         "All 60 agreed Float64 cases were measured on the final production source. Of 59 valid baseline comparisons, 58 are at most 2×; reduced Laguerre(10000,−0.99) is **2.043554684×**, explicitly accepted by the user. The baseline Laguerre(39,44.8) output is invalid and is excluded from timing acceptance. Exact measurements remain visible.", "",
         "The independent accuracy matrix has **185 passes and one proved Float32 representation limit**, covering **4,815 node/weight pairs**: 4,621 normal, 39 subnormal and 155 genuinely underflowed reference weights. All applicable frozen budgets, 28 additional Lobatto checks and eight reduced-prefix boundary checks pass. The largest node/weight budget fractions are 0.5032/0.2380. Larger orders are sampled.", "",
         "Fresh suites pass **7,121 assertions on both Julia 1.10.13 and 1.13.3**, with all Aqua groups. Runic, strict imports, 23 targeted JET checks, and Documenter/doctests pass. [Quality details](quality-summary.md) include exact versions, source hashes and redacted logs. Developer tools are outside runtime dependencies.", "",
         "## Selected Float64 comparisons", "", "Node errors are maximum absolute errors; weight errors are maximum relative errors among checked normal weights. Every selected reference is independently indexed and stable at increased precision. Times are warmed pooled medians in microseconds. Full results, ulp diagnostics, allocations and statuses are in [the 60-case CSV](float64-accuracy-timing.csv).", "",
         "| Case | Node error, master → candidate | Normal weight error, master → candidate | Master µs | Candidate µs | Candidate/master |", "| --- | --- | --- | ---: | ---: | ---: |"]
for name in select:
    r = byid[name]
    lines.append(f"| `{name}` | {r['master_node_absolute']:.3g} → {r['candidate_node_absolute']:.3g} | {r['master_normal_weight_relative']:.3g} → {r['candidate_normal_weight_relative']:.3g} | {r['baseline_microseconds']:.3f} | {r['candidate_microseconds']:.3f} | {r['candidate_over_baseline']:.3f}× |")
lines += ["", "## Other formats and scaling", "", "[Native-format tables](native-accuracy-timing.csv) retain Float32, actual 113-bit Quadmath.Float128, and BigFloat at 128/256 bits, including allocations and all available matching accuracy rows. Native timings precede the final Float64-only selection repairs; their algorithms are unchanged by the reviewed Float64-only source changes. [Complete-output comparisons](source-output-bridge.csv) additionally verify the eight overlapping accuracy/timing inputs per format; they do not replay all 16 timed inputs per format. Missing baseline methods and invalid outputs are retained in the raw status files. A time without a matching reference row is not an accuracy pass.", "",
          "[Doubling-order measurements](scaling.csv) cover 14 family/regime/semantics combinations at n=2048,4096,8192,16384. Full-rule timings support the reviewed linear-work implementation; individual doubling ratios still vary with allocation and runtime effects. Reduced Laguerre returns a growing prefix and is separate work. These measurements supplement bounded loops and call-chain review, not a proof of all-input convergence.", "",
          "## Reproduction and limits", "", "See [the reproducible drivers](../README.md), [frozen model explanations](../regimes/README.md), and [removed master tests](../removed-tests.md). `python3 validation/report.py` verifies raw timing medians/counts and regenerates these joined tables. Raw CSVs are gzip-compressed; no numerical values are redacted. Environment prose replaces local home/check-out paths with placeholders.", "",
          "Float64 uses eight alternating batches ×30 samples with compilation excluded; other formats/scaling use three ×10. Hardware: Apple M3 Ultra, Julia 1.13.3, 32 Julia threads and one BLAS thread. Baseline and candidate share the same runtime/dependencies within each paired measurement. The native run was continued after a 300-second tool timeout, preserving completed cases; the final Float64/scaling refresh completed all observations before a timestamp-only logger error. Both events and raw completeness checks are retained in the environment records.", "",
          "Accuracy means good stable native-precision behavior with justified order/parameter/location growth, not a uniform ulp or correctly-rounded promise. The Float32 floor, normalization exponent range, distinct-node resolution and reduced cutoff remain explicit limits. The Float32 Laguerre(39,44.8) exception has an independent mass/n overflow certificate. A normalization overflow alone is not such a certificate. Manual recurrence/asymptotic helpers retain documented limits; automatic Laguerre now uses GST to avoid their observed losses.", "",
          "The exact four agreed reduced-rule boundaries are in [reduced-boundaries.csv](reduced-boundaries.csv): the last retained and first omitted references are on opposite sides of 10floatmin(Float64). Prefix nodes match exactly; prefix weights agree within the preserved 8eps allowance (one n100000 half-parameter case differs by 2.22e−16 relatively).", "",
          "Earlier failing measurements are preserved in the local study archive, with selected repair records here. The final fixes use existing native formulas and unchanged budgets: Hermite physical-weight residual transport, generalized Radau near−1 routing, and removal of recurrence/full asymptotic selection from automatic Laguerre. No higher working precision or guard arithmetic was introduced.", ""]
(DATA / "README.md").write_text("\n".join(lines))
print("Verified and generated 60 Float64 rows,", len(native), "native timing rows, and", len(scaling_rows), "scaling rows.")
