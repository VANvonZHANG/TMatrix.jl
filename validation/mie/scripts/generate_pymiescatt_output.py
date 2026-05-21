"""Generate PyMieScatt output from input CSV.

Usage:
    python generate_pymiescatt_output.py

Reads:  ../inputs/test_cases.csv
Writes: ../outputs/pymiescatt_results.csv
"""

import csv
import sys
from pathlib import Path

try:
    import PyMieScatt as ps
except ImportError:
    print("Error: PyMieScatt not installed. Run: pip install PyMieScatt")
    sys.exit(1)

# Paths
BASE_DIR = Path(__file__).parent.parent
INPUT_CSV = BASE_DIR / "inputs" / "test_cases.csv"
OUTPUT_CSV = BASE_DIR / "outputs" / "pymiescatt_results.csv"


def main():
    if not INPUT_CSV.exists():
        print(f"Error: Input file not found: {INPUT_CSV}")
        sys.exit(1)

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    results = []
    with open(INPUT_CSV, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            m = float(row["m_real"]) + float(row["m_imag"]) * 1j
            wavelength = float(row["wavelength_nm"])
            diameter = float(row["diameter_nm"])

            # Efficiencies
            Qext, Qsca, Qabs, g, Qpr, Qback, Qratio = ps.MieQ(
                m, wavelength, diameter
            )
            # Cross-sections
            Cext, Csca, Cabs, g_cs, Cpr, Cback, Cratio = ps.MieQ(
                m, wavelength, diameter, asCrossSection=True
            )

            results.append(
                {
                    "name": name,
                    "m_real": row["m_real"],
                    "m_imag": row["m_imag"],
                    "wavelength_nm": row["wavelength_nm"],
                    "diameter_nm": row["diameter_nm"],
                    "Qext": Qext,
                    "Qsca": Qsca,
                    "Qabs": Qabs,
                    "g": g,
                    "Cext": Cext,
                    "Csca": Csca,
                    "Cabs": Cabs,
                }
            )

    # Write CSV
    fieldnames = [
        "name",
        "m_real",
        "m_imag",
        "wavelength_nm",
        "diameter_nm",
        "Qext",
        "Qsca",
        "Qabs",
        "g",
        "Cext",
        "Csca",
        "Cabs",
    ]
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"PyMieScatt results written to {OUTPUT_CSV}")
    for r in results:
        print(f"  {r['name']}: Qext={r['Qext']:.6f}, Qsca={r['Qsca']:.6f}, g={r['g']:.6f}")


if __name__ == "__main__":
    main()
