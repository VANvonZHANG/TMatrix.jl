# Getting Started

This tutorial walks through your first T-matrix computation in Julia.

## Installation

Install TMatrix.jl from the GitHub repository:

```julia
using Pkg
Pkg.add("https://github.com/VanVonZhang/T-Matrix.jl")
```

## Your First Mie Calculation

Import the package:

```julia
using TMatrix
```

### Step 1: Define the particle

Create a spherical particle with radius 1.0 (in whatever length units you are using):

```julia
particle = Sphere(1.0)
```

### Step 2: Set physical parameters

```julia
wavelength = 0.5       # wavelength in the same units as radius
k = 2π / wavelength    # wave number in the surrounding medium
m = 1.5 + 0.01im       # refractive index of the sphere
```

### Step 3: Determine truncation order

The Mie series must be truncated at some maximum multipole order `N_max`.
Use the Wiscombe criterion:

```julia
x = k * 1.0            # size parameter
N_max = mie_nmax(x)    # automatic truncation
```

### Step 4: Solve the T-matrix

```julia
t_matrix = solve_tmatrix(particle, MieMethod(), k, m, N_max)
```

The result is a [`MieTMatrix`](@ref) struct containing the diagonal Mie coefficients.

### Step 5: Compute cross-sections

```julia
cs = calc_cross_sections(t_matrix)
```

This returns a [`CrossSections`](@ref) struct with:

| Field | Symbol | Description |
|-------|--------|-------------|
| `extinction` | `C_ext` | Extinction cross-section |
| `scattering` | `C_sca` | Scattering cross-section |
| `absorption` | `C_abs` | Absorption cross-section |
| `Q_ext` | — | Extinction efficiency |
| `Q_sca` | — | Scattering efficiency |
| `Q_abs` | — | Absorption efficiency |
| `asymmetry` | `g` | Asymmetry parameter |

Print the result:

```julia
println(cs)
```

## Complete Example

```julia
using TMatrix

particle = Sphere(1.0)
k = 2π / 0.5
m = 1.5 + 0.01im
N_max = mie_nmax(k * 1.0)

t_matrix = solve_tmatrix(particle, MieMethod(), k, m, N_max)
cs = calc_cross_sections(t_matrix)

println(cs)
```

## Next Steps

- Learn about [Mie Scattering](mie-scattering.md) with parametric scans and visualizations
- Browse the [API Reference](../api/mie-solver.md) for detailed function documentation
- Read the [Theory](../theory/mie-theory.md) section for the mathematical background
