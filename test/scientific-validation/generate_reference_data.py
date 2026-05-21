"""Generate Mie scattering benchmark data using PyMieScatt."""

import json
import PyMieScatt as ps

test_cases = [
    {"name": "typical_aerosol", "m": 1.5 + 0.01j, "wavelength": 532, "diameter": 100},
    {"name": "water_droplet", "m": 1.33 + 0.0j, "wavelength": 650, "diameter": 500},
    {"name": "soot", "m": 2.0 + 1.0j, "wavelength": 1000, "diameter": 200},
    {"name": "high_absorption", "m": 1.5 + 0.5j, "wavelength": 800, "diameter": 1000},
]

results = []
for case in test_cases:
    Qext, Qsca, Qabs, g, Qpr, Qback, Qratio = ps.MieQ(
        case["m"], case["wavelength"], case["diameter"]
    )
    Cext, Csca, Cabs, g_cs, Cpr, Cback, Cratio = ps.MieQ(
        case["m"], case["wavelength"], case["diameter"], asCrossSection=True
    )
    results.append({
        "name": case["name"],
        "m": {"real": case["m"].real, "imag": case["m"].imag},
        "wavelength_nm": case["wavelength"],
        "diameter_nm": case["diameter"],
        "Qext": Qext,
        "Qsca": Qsca,
        "Qabs": Qabs,
        "g": g,
        "Cext": Cext,
        "Csca": Csca,
        "Cabs": Cabs,
    })

with open("reference_data.json", "w") as f:
    json.dump(results, f, indent=2)

print("Reference data written to reference_data.json")
