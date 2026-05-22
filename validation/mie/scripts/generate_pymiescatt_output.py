"""Generate PyMieScatt output from input CSV.

Usage:
    python generate_pymiescatt_output.py [input_csv] [output_csv]

Defaults:
    input:  ../inputs/test_cases.csv
    output: ../outputs/pymiescatt_results.csv
"""

import csv
import sys
from pathlib import Path

try:
    import PyMieScatt as ps
except ImportError:
    print("Error: PyMieScatt not installed. Run: pip install PyMieScatt")
    sys.exit(1)


def main():
    base_dir = Path(__file__).parent.parent
    input_csv = Path(sys.argv[1]) if len(sys.argv) > 1 else base_dir / "inputs" / "test_cases.csv"
    output_csv = Path(sys.argv[2]) if len(sys.argv) > 2 else base_dir / "outputs" / "pymiescatt_results.csv"

    if not input_csv.exists():
        print(f"Error: Input file not found: {input_csv}")
        sys.exit(1)

    output_csv.parent.mkdir(parents=True, exist_ok=True)

    results = []
    with open(input_csv, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            m = float(row["m_real"]) + float(row["m_imag"]) * 1j
            wavelength = float(row["wavelength_nm"])
            diameter = float(row["diameter_nm"])

            Qext, Qsca, Qabs, g, _, _, _ = ps.MieQ(m, wavelength, diameter)
            Cext, Csca, Cabs, _, _, _, _ = ps.MieQ(
                m, wavelength, diameter, asCrossSection=True
            )

            results.append({
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
            })

    fieldnames = [
        "name", "m_real", "m_imag", "wavelength_nm", "diameter_nm",
        "Qext", "Qsca", "Qabs", "g", "Cext", "Csca", "Cabs",
    ]
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"PyMieScatt results written to {output_csv} ({len(results)} cases)")


if __name__ == "__main__":
    main()
