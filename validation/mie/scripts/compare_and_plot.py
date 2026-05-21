"""Compare TMatrix.jl and PyMieScatt outputs, generate report and plots.

Usage:
    python compare_and_plot.py [--rtol RTOL] [--atol ATOL]

Reads:
    - ../outputs/julia_results.csv
    - ../outputs/pymiescatt_results.csv

Writes:
    - ../outputs/comparison_report.txt
    - ../outputs/comparison_plot.png
"""

import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Paths
BASE_DIR = Path(__file__).parent.parent
JULIA_CSV = BASE_DIR / "outputs" / "julia_results.csv"
PYMIE_CSV = BASE_DIR / "outputs" / "pymiescatt_results.csv"
REPORT_TXT = BASE_DIR / "outputs" / "comparison_report.txt"
PLOT_PNG = BASE_DIR / "outputs" / "comparison_plot.png"

# Quantities to compare
QUANTITIES = ["Qext", "Qsca", "Qabs", "g", "Cext", "Csca", "Cabs"]


def read_csv(path: Path) -> dict:
    """Read CSV into dict: name -> {quantity: value}."""
    data = {}
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {q: float(row[q]) for q in QUANTITIES}
    return data


def approx(actual, expected, rtol=1e-6, atol=1e-10):
    """Check if actual ≈ expected, handling near-zero values."""
    if abs(expected) < atol or abs(actual) < atol:
        return abs(actual - expected) < atol
    return abs(actual - expected) <= rtol * abs(expected)


def compare(julia_data, pymie_data, rtol, atol):
    """Compare outputs and return results."""
    all_pass = True
    results = []

    for name in sorted(julia_data.keys()):
        if name not in pymie_data:
            print(f"  SKIP: {name} not found in PyMieScatt output")
            continue

        j = julia_data[name]
        p = pymie_data[name]
        case_pass = True
        case_results = {"name": name, "quantities": {}}

        for q in QUANTITIES:
            diff = abs(j[q] - p[q])
            rel_diff = diff / abs(p[q]) if p[q] != 0 else float("inf")
            passed = approx(j[q], p[q], rtol=rtol, atol=atol)
            case_pass = case_pass and passed
            case_results["quantities"][q] = {
                "julia": j[q],
                "pymie": p[q],
                "diff": diff,
                "rel_diff": rel_diff,
                "passed": passed,
            }

        case_results["passed"] = case_pass
        results.append(case_results)
        all_pass = all_pass and case_pass

    return all_pass, results


def generate_report(results, rtol, atol):
    """Generate text report."""
    lines = []
    lines.append("=" * 70)
    lines.append("TMatrix.jl vs PyMieScatt Comparison Report")
    lines.append("=" * 70)
    lines.append(f"Tolerance: rtol={rtol}, atol={atol}")
    lines.append("")

    total_passed = 0
    total_failed = 0

    for case in results:
        name = case["name"]
        passed = case["passed"]
        status = "PASS" if passed else "FAIL"
        lines.append(f"\n{name}: {status}")
        lines.append("-" * 40)

        for q in QUANTITIES:
            qres = case["quantities"][q]
            qp = "PASS" if qres["passed"] else "FAIL"
            if not qres["passed"]:
                total_failed += 1
            else:
                total_passed += 1
            lines.append(
                f"  {q:8s}: julia={qres['julia']:.6e}, "
                f"pymie={qres['pymie']:.6e}, "
                f"rel_diff={qres['rel_diff']:.6e} [{qp}]"
            )

    lines.append("")
    lines.append("=" * 70)
    lines.append(f"Summary: {total_passed} passed, {total_failed} failed")
    lines.append("=" * 70)

    return "\n".join(lines)


def generate_plot(results):
    """Generate comparison bar plot."""
    n_cases = len(results)
    n_quantities = len(QUANTITIES)

    fig, axes = plt.subplots(
        n_quantities, 1, figsize=(10, 3 * n_quantities), squeeze=False
    )

    for idx, q in enumerate(QUANTITIES):
        ax = axes[idx, 0]
        names = [r["name"] for r in results]
        julia_vals = [r["quantities"][q]["julia"] for r in results]
        pymie_vals = [r["quantities"][q]["pymie"] for r in results]

        x = np.arange(len(names))
        width = 0.35

        ax.bar(x - width / 2, julia_vals, width, label="TMatrix.jl", alpha=0.8)
        ax.bar(x + width / 2, pymie_vals, width, label="PyMieScatt", alpha=0.8)

        ax.set_ylabel(q)
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=15, ha="right")
        ax.legend()
        ax.set_yscale("log" if all(v > 0 for v in julia_vals + pymie_vals) else "linear")
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(PLOT_PNG, dpi=150, bbox_inches="tight")
    print(f"Comparison plot saved to {PLOT_PNG}")


def main():
    parser = argparse.ArgumentParser(description="Compare TMatrix.jl and PyMieScatt")
    parser.add_argument("--rtol", type=float, default=1e-6, help="Relative tolerance")
    parser.add_argument("--atol", type=float, default=1e-10, help="Absolute tolerance")
    args = parser.parse_args()

    if not JULIA_CSV.exists():
        print(f"Error: Julia output not found: {JULIA_CSV}")
        print("Run: julia --project=../.. generate_julia_output.jl")
        sys.exit(1)

    if not PYMIE_CSV.exists():
        print(f"Error: PyMieScatt output not found: {PYMIE_CSV}")
        print("Run: python generate_pymiescatt_output.py")
        sys.exit(1)

    julia_data = read_csv(JULIA_CSV)
    pymie_data = read_csv(PYMIE_CSV)

    print(f"Comparing {len(julia_data)} test cases...")
    all_pass, results = compare(julia_data, pymie_data, args.rtol, args.atol)

    report = generate_report(results, args.rtol, args.atol)
    print(report)

    with open(REPORT_TXT, "w") as f:
        f.write(report)
    print(f"\nReport saved to {REPORT_TXT}")

    generate_plot(results)

    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
