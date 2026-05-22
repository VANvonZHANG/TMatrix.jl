# Architecture

TMatrix.jl is organized in dependency layers. Each layer depends only on layers
below it, following a strict bottom-up design.

## Layer Structure

```
Layer 2+: Mie Solver (mie.jl)
    ↑ depends on
Layer 2: VSWFs & Geometry (vswf.jl, geometry.jl)
    ↑ depends on
Layer 1: Composed Functions (composed.jl)
    ↑ depends on
Layer 0: Direct Dependencies & Thin Compositions (special_functions.jl)
    ↑ depends on
External: SpecialFunctions.jl, LegendrePolynomials.jl, WignerSymbols.jl,
          FastGaussQuadrature.jl, StaticArrays.jl
```

## Layer 0: Mathematical Primitives

Thin compositions around mature Julia libraries. No wrappers — direct API calls
with only the functions that don't exist in any library.

Examples:
- `shankelh1` = composition of `sphericalbesselj` + `sphericalbessely`
- `dplmdtheta` = recurrence relation using `Plm`

## Layer 1: Composed Functions

Functions built from Layer 0 that have no direct library equivalent:

- `gaunt` — Wigner 3-j composition
- `P_vsh`, `B_vsh`, `C_vsh` — vector spherical harmonics

## Layer 2: VSWFs and Geometry

Complete physical objects assembled from lower layers:

- `vswf_M`, `vswf_N`, `vswf_L` — full wave functions
- `Sphere`, `Spheroid`, `ChebyshevParticle` — particle geometries

## Layer 2+: Solvers

The Mie solver is the first complete solver. Future solvers (EBCM, IITM, etc.)
will follow the same pattern: consume Layer 2 objects and produce a T-matrix.

## Solver Interface Design

The target interface for Phase 2+ is:

```julia
abstract type AbstractTMatrixMethod end
struct MieMethod <: AbstractTMatrixMethod end
struct EBCM <: AbstractTMatrixMethod end
# ... etc

function solve_tmatrix(particle::AbstractParticle, method::AbstractTMatrixMethod, args...)
    # Multiple dispatch routes to the correct solver
end
```

This allows post-processing code (`calc_cross_sections`, orientation averaging)
to operate generically on any solver output.

## Module Organization

`src/TMatrix.jl` controls load order via `include()`:

```julia
include("special_functions.jl")  # Layer 0
include("composed.jl")           # Layer 1
include("vswf.jl")               # Layer 2
include("geometry.jl")           # Layer 2
include("mie.jl")                # Layer 2+
```

Exports are grouped by layer for clarity.
