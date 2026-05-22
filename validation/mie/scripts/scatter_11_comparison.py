"""Scatter plot comparison (1:1 line) for parametric scan results.

Usage:
    python scatter_11_comparison.py [julia_csv] [pymie_csv] [output_png]

Defaults:
    julia_csv:  ../outputs/julia_parametric.csv
    pymie_csv:  ../outputs/pymiescatt_parametric.csv
    output_png: ../outputs/scatter_11_comparison.png
"""

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

QUANTITIES = ["Qext", "Qsca", "Qabs", "g", "Cext", "Csca", "Cabs"]
NEAR_ZERO_ATOL = 1e-6  # threshold for treating values as near-zero


def read_csv(path: Path) -> dict:
    data = {}
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            data[row["name"]] = {q: float(row[q]) for q in QUANTITIES}
    return data


def compute_rel_diff(j, p, atol=NEAR_ZERO_ATOL):
    """Compute relative difference, treating near-zero values as exact matches."""
    if abs(j) < atol and abs(p) < atol:
        return 0.0
    if abs(p) > 1e-12:
        return abs(j - p) / abs(p)
    return abs(j - p)


def main():
    base_dir = Path(__file__).parent.parent
    julia_csv = Path(sys.argv[1]) if len(sys.argv) > 1 else base_dir / "outputs" / "julia_parametric.csv"
    pymie_csv = Path(sys.argv[2]) if len(sys.argv) > 2 else base_dir / "outputs" / "pymiescatt_parametric.csv"
    output_png = Path(sys.argv[3]) if len(sys.argv) > 3 else base_dir / "outputs" / "scatter_11_comparison.png"

    if not julia_csv.exists():
        print(f"Error: Julia output not found: {julia_csv}")
        sys.exit(1)
    if not pymie_csv.exists():
        print(f"Error: PyMieScatt output not found: {pymie_csv}")
        sys.exit(1)

    julia_data = read_csv(julia_csv)
    pymie_data = read_csv(pymie_csv)

    names = sorted(set(julia_data.keys()) & set(pymie_data.keys()))
    print(f"Comparing {len(names)} cases...")

    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()

    max_rel_diff = 0.0
    failed_cases = []

    for idx, q in enumerate(QUANTITIES):
        ax = axes[idx]
        j_vals = np.array([julia_data[n][q] for n in names])
        p_vals = np.array([pymie_data[n][q] for n in names])

        # Compute relative differences with near-zero handling
        rel_diffs = []
        for n, j, p in zip(names, j_vals, p_vals):
            rd = compute_rel_diff(j, p)
            rel_diffs.append(rd)
            if rd > 1e-6:
                failed_cases.append((n, q, j, p, rd))

        local_max = max(rel_diffs) if rel_diffs else 0.0
        max_rel_diff = max(max_rel_diff, local_max)

        # Scatter plot: use log scale for quantities spanning many orders
        use_log = q in ("Cext", "Csca", "Cabs") and np.min(np.abs(j_vals) + 1e-30) > 0
        if use_log:
            # Filter out non-positive for log scale
            mask = (j_vals > 1e-30) & (p_vals > 1e-30)
            ax.scatter(p_vals[mask], j_vals[mask], s=15, alpha=0.5, edgecolors="none")
        else:
            ax.scatter(p_vals, j_vals, s=15, alpha=0.5, edgecolors="none")

        # 1:1 line
        vmin = min(np.min(j_vals), np.min(p_vals))
        vmax = max(np.max(j_vals), np.max(p_vals))
        margin = (vmax - vmin) * 0.05 if vmax != vmin else 1.0
        ax.plot([vmin - margin, vmax + margin], [vmin - margin, vmax + margin], "r--", lw=1.5, label="1:1")

        ax.set_xlabel(f"PyMieScatt {q}")
        ax.set_ylabel(f"TMatrix.jl {q}")
        ax.set_title(f"{q}  (max rel. diff: {local_max:.2e})")
        ax.legend(loc="upper left")
        ax.grid(True, alpha=0.3)

    # Hide unused subplot
    axes[-1].axis("off")

    # Summary statistics on unused subplot
    ax_summary = axes[-1]
    ax_summary.axis("on")
    ax_summary.set_xlim(0, 1)
    ax_summary.set_ylim(0, 1)
    ax_summary.axis("off")

    summary_text = (
        f"Total cases: {len(names)}\n"
        f"Quantities checked: {len(QUANTITIES)}\n"
        f"Failed checks (>{'1e-6'}): {len(failed_cases)}\n"
        f"Max rel. diff: {max_rel_diff:.2e}"
    )
    ax_summary.text(0.1, 0.5, summary_text, fontsize=14, family="monospace",
                    verticalalignment="center")

    fig.suptitle(
        f"TMatrix.jl vs PyMieScatt: Parametric Scan\n"
        f"Max relative difference (near-zero-aware): {max_rel_diff:.2e}",
        fontsize=12,
    )

    plt.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(output_png, dpi=150, bbox_inches="tight")
    print(f"Scatter plot saved to {output_png}")

    if max_rel_diff < 1e-6:
        print("\n✅ All points lie on the 1:1 line. Algorithms agree.")
    elif max_rel_diff < 1e-4:
        print(f"\n⚠️  Small deviations (max rel. diff = {max_rel_diff:.2e}).")
    else:
        print(f"\n❌ Significant deviations (max rel. diff = {max_rel_diff:.2e}).")
        if failed_cases:
            print(f"   First failure: {failed_cases[0]}")


if __name__ == "__main__":
    main()
