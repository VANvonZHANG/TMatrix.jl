# Scientific Validation

These tests validate T-Matrix.jl's Mie solver against PyMieScatt reference data.

## Running Validation

1. Install PyMieScatt:
   ```bash
   pip install PyMieScatt
   ```

2. Generate reference data:
   ```bash
   cd test/scientific-validation
   python generate_reference_data.py
   ```

3. Run Julia validation:
   ```bash
   cd test/scientific-validation
   julia --project=../.. -e 'using Pkg; Pkg.test()' test_pymiescatt_validation.jl
   ```

## Note

These tests are **not** run in CI. They require Python + PyMieScatt.
