# TMatrix.jl

A Julia library for T-matrix electromagnetic scattering computations, providing a unified framework for analytical, boundary-integral, and volume-differential solvers.

## Installation

```julia
using Pkg
Pkg.add("https://github.com/VanVonZhang/T-Matrix.jl")
```

## Quick Start

```julia
using TMatrix

# Coming soon — see ROADMAP.md for planned features
```

## Features

See [ROADMAP.md](ROADMAP.md) for the full development plan.

- **Mie solver** — Homogeneous and multi-layered spheres
- **EBCM** — Extended Boundary Condition Method for moderate aspect ratios
- **IITM** — Invariant Imbedding T-matrix Method for complex/inhomogeneous particles
- **NFM-DS** — Null-Field Method with Discrete Sources for extreme aspect ratios
- **MSTM** — Superposition T-matrix for multi-sphere aggregates

## License

MIT License. See [LICENSE](LICENSE).
