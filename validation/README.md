# Scientific Validation

Validate TMatrix.jl's Mie solver against [PyMieScatt](https://pypi.org/project/PyMieScatt/) reference data.

## Directory Structure

```
validation/
├── inputs/
│   └── test_cases.csv          # Input parameters (shared by both solvers)
├── outputs/
│   ├── julia_results.csv       # TMatrix.jl output
│   ├── pymiescatt_results.csv  # PyMieScatt output
│   ├── comparison_report.txt   # Text comparison report
│   └── comparison_plot.png     # Bar chart comparison
├── scripts/
│   ├── generate_julia_output.jl      # Run TMatrix.jl solver
│   ├── generate_pymiescatt_output.py # Run PyMieScatt solver
│   └── compare_and_plot.py           # Compare outputs + generate plots
└── README.md
```

## Workflow

### 1. Define test cases

Edit `inputs/test_cases.csv`:

```csv
name,m_real,m_imag,wavelength_nm,diameter_nm
"typical_aerosol",1.5,0.01,532,100
```

Columns:
- `name`: test case identifier
- `m_real`, `m_imag`: real and imaginary parts of refractive index
- `wavelength_nm`: wavelength in nanometers
- `diameter_nm`: particle diameter in nanometers

### 2. Generate outputs

Run **both** solvers from the `scripts/` directory:

```bash
cd validation/scripts

# TMatrix.jl
julia --project=../.. generate_julia_output.jl

# PyMieScatt (requires: pip install PyMieScatt)
python generate_pymiescatt_output.py
```

Each script reads `../inputs/test_cases.csv` and writes to `../outputs/<solver>_results.csv`.

### 3. Compare and plot

```bash
python compare_and_plot.py
```

This produces:
- `../outputs/comparison_report.txt` — numerical comparison with pass/fail per quantity
- `../outputs/comparison_plot.png` — side-by-side bar charts

Adjust tolerances:
```bash
python compare_and_plot.py --rtol 1e-5 --atol 1e-8
```

## Output CSV Format

Both `julia_results.csv` and `pymiescatt_results.csv` share identical columns:

| Column | Description |
|--------|-------------|
| `name` | test case name |
| `m_real`, `m_imag` | refractive index components |
| `wavelength_nm` | wavelength (nm) |
| `diameter_nm` | diameter (nm) |
| `Qext` | extinction efficiency |
| `Qsca` | scattering efficiency |
| `Qabs` | absorption efficiency |
| `g` | asymmetry parameter |
| `Cext` | extinction cross-section (nm²) |
| `Csca` | scattering cross-section (nm²) |
| `Cabs` | absorption cross-section (nm²) |

## Notes

- These validation tests are **not** run in CI. They require both Julia + TMatrix.jl and Python + PyMieScatt.
- The comparison script uses smart near-zero tolerance: when the reference value is below `atol`, the test checks that the actual value is also near-zero rather than using relative tolerance.
