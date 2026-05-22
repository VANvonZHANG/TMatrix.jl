# TMatrix.jl

A Julia library for T-matrix electromagnetic scattering computations, providing a unified framework for analytical, boundary-integral, and volume-differential solvers.

## Installation

```julia
using Pkg
Pkg.add("https://github.com/VanVonZhang/T-Matrix.jl")
```

## Quick Start

Compute Mie scattering from a homogeneous sphere:

```julia
using TMatrix

# Define a spherical particle
particle = Sphere(1.0)  # radius = 1.0 μm

# Solve for the T-matrix
k = 2π / 0.5  # wave number for λ = 0.5 μm
m = 1.5 + 0.01im  # refractive index
N_max = mie_nmax(k * 1.0)
t_matrix = solve_tmatrix(particle, MieMethod(), k, m, N_max)

# Compute cross-sections
cs = calc_cross_sections(t_matrix)
# CrossSections:
#   C_ext = ...
#   C_sca = ...
#   C_abs = ...
```

## Features

| Method | Status | Description |
|--------|--------|-------------|
| **Mie solver** | ✅ Implemented | Homogeneous spheres (exact diagonal T-matrix) |
| **SVM** | Planned | Spheroids via spheroidal wave functions |
| **EBCM** | Planned | Moderate aspect ratio non-spherical particles |
| **NFM-DS** | Planned | Extreme aspect ratio particles |
| **IITM** | Planned | Complex/inhomogeneous particles |
| **MSTM** | Planned | Multi-sphere aggregates |

See [ROADMAP.md](https://github.com/VanVonZhang/T-Matrix.jl/blob/main/ROADMAP.md) for the full development plan.

## Documentation Overview

- **Tutorials** — Step-by-step guides for getting started and worked examples
- **API Reference** — Complete function documentation auto-generated from docstrings
- **Theory** — Mathematical background following Sun et al. 2020 conventions
- **Developer Docs** — Code architecture, numerical stability notes, design rationale

## Citation

If you use TMatrix.jl in your research, please cite:

- Sun et al. (2020), *T-matrix Concept*, Chapter 3

## License

MIT License. See [LICENSE](https://github.com/VanVonZhang/T-Matrix.jl/blob/main/LICENSE).
