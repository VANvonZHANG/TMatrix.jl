"""Generate parametric scan input CSV for Mie validation.

Scans over refractive index (m_real, m_imag) and size parameter x,
keeping wavelength fixed at 532 nm.

Usage:
    python generate_parametric_scan.py [output_csv]

Default output: ../inputs/parametric_scan.csv
"""

import csv
import sys
from pathlib import Path

# Fixed wavelength
WAVELENGTH_NM = 532.0

# Scan ranges
M_REAL_VALUES = [1.1, 1.2, 1.33, 1.5, 1.7, 2.0, 2.5, 3.0]
M_IMAG_VALUES = [0.0, 0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0]
X_VALUES = [0.1, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0]


def main():
    base_dir = Path(__file__).parent.parent
    output_csv = Path(sys.argv[1]) if len(sys.argv) > 1 else base_dir / "inputs" / "parametric_scan.csv"

    output_csv.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    counter = 0
    for m_real in M_REAL_VALUES:
        for m_imag in M_IMAG_VALUES:
            for x in X_VALUES:
                counter += 1
                # Compute diameter from size parameter x = k*r = π*d/λ
                # d = x * λ / π
                diameter_nm = x * WAVELENGTH_NM / 3.141592653589793
                rows.append({
                    "name": f"scan_{counter:04d}",
                    "m_real": m_real,
                    "m_imag": m_imag,
                    "wavelength_nm": WAVELENGTH_NM,
                    "diameter_nm": round(diameter_nm, 3),
                })

    fieldnames = ["name", "m_real", "m_imag", "wavelength_nm", "diameter_nm"]
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Generated {len(rows)} scan cases: {output_csv}")
    print(f"  m_real: {M_REAL_VALUES}")
    print(f"  m_imag: {M_IMAG_VALUES}")
    print(f"  x: {X_VALUES}")


if __name__ == "__main__":
    main()
