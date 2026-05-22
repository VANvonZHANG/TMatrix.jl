# Solver Interface Design

This page documents the planned solver interface for Phase 2+ of TMatrix.jl.

## Abstract Type Hierarchy

```julia
abstract type AbstractTMatrixMethod end

struct MieMethod <: AbstractTMatrixMethod end
struct SVMMethod <: AbstractTMatrixMethod end
struct EBCMMethod <: AbstractTMatrixMethod end
struct NFMDSMethod <: AbstractTMatrixMethod end
struct IITMMethod <: AbstractTMatrixMethod end
struct MSTMMethod <: AbstractTMatrixMethod end
```

## Generic solve_tmatrix Interface

```julia
solve_tmatrix(particle, method::AbstractTMatrixMethod, k, refractive_index, N_max)
```

Each solver implements this method via multiple dispatch. The return type is
solver-specific (e.g., `MieTMatrix` for Mie, dense matrix for EBCM).

## Generic Post-Processing

Cross-section computation should work for any solver:

```julia
calc_cross_sections(t_matrix::MieTMatrix) -> CrossSections
# Future:
calc_cross_sections(t_matrix::EBCMTMatrix) -> CrossSections
```

This is achieved by defining `calc_cross_sections` for each concrete T-matrix type.

## Solver-Specific Options

Solvers may accept additional keyword arguments via a configuration struct:

```julia
struct EBCMOptions
    n_quadrature::Int      # Number of Gauss-Legendre points
    tolerance::Float64     # Convergence threshold
end

solve_tmatrix(particle, EBCMMethod(), k, m, N_max; options=EBCMOptions(100, 1e-10))
```

## Hybrid Solvers

The architecture supports solver fusion:

```julia
# SVM core + IITM rough shell
t_core = solve_tmatrix(spheroid_core, SVMMethod(), k, m, N_max)
t_total = solve_tmatrix(rough_shell, IITMMethod(), k, m, N_max; initial_t=t_core)
```

This is a key advantage of Julia's multiple dispatch system.

## Status

| Solver | Status | T-matrix Output |
|--------|--------|-----------------|
| Mie | Implemented | `MieTMatrix` (diagonal) |
| SVM | Planned | Dense (spheroidal basis) |
| EBCM | Planned | Dense |
| NFM-DS | Planned | Dense |
| IITM | Planned | Dense |
| MSTM | Planned | Dense (aggregate) |
