# Phase 1 Roadmap: Mathematical and Physical Primitives

> **Status:** Planning | **Reference Convention:** Sun et al. 2020, *T-matrix Concept*

This document provides the detailed technical roadmap for Phase 1 of T-Matrix.jl. Phase 1 is organized into three dependency layers: each layer depends only on the layers below it, and produces testable output before the next layer begins.

The guiding principle: **use existing Julia libraries directly, compose only what doesn't exist.**

---

## Dependency Overview

```
Layer 2: VSWF Module (M, N, L vectors + Geometry)
    ↑ depends on
Layer 1: Composed Functions (Gaunt coefficients, vector spherical harmonics P/B/C)
    ↑ depends on
Layer 0: Direct Dependencies (SpecialFunctions.jl, LegendrePolynomials.jl,
         WignerSymbols.jl, FastGaussQuadrature.jl)
```

---

## Layer 0: Direct Dependencies & Thin Compositions

Direct API calls to mature Julia libraries. No wrappers. Only thin compositions for functions not provided by any library.

### 0.1 Spherical Bessel & Hankel Functions

**Reference:** Sun et al. 2020, Eqs. 3.2.27–3.2.33

**Direct API (`SpecialFunctions.jl`):**

| Function | Julia API | Notes |
|----------|-----------|-------|
| j_n(x) | `sphericalbesselj(n, x)` | Spherical Bessel of the 1st kind |
| y_n(x) | `sphericalbessely(n, x)` | Spherical Bessel of the 2nd kind |

Both accept `Float64` and `ComplexF64` arguments natively (complex refractive index requires complex arguments).

**Thin compositions (must write):**

- `shankelh1(n, x) = sphericalbesselj(n, x) + im * sphericalbessely(n, x)` — outward-traveling wave, used for scattered fields (Sun Eq. 3.2.28c)
- `shankelh2(n, x) = sphericalbesselj(n, x) - im * sphericalbessely(n, x)` — inward-traveling wave (Sun Eq. 3.2.28d)
- Derivative functions via recurrence relations (Sun Eq. 3.2.29): `(2n+1)/x · z_n = z_{n-1} + z_{n+1}`, rather than numerical differentiation

**Known limitation:** `SpecialFunctions.jl` has known gaps in `BigFloat` support for spherical Bessel functions. For n > 100 with extreme size parameters, consider `Bessels.jl` as a fallback for real arguments, or direct `BigFloat` computation. This is not a blocker for Phase 1 but should be tracked.

**Validation:**

| Test | Expected | Tolerance |
|------|----------|-----------|
| `j_0(x) = sin(x)/x` (Sun Eq. 3.2.33a) | Exact closed form | `rtol=1e-14` |
| `y_0(x) = -cos(x)/x` (Sun Eq. 3.2.33b) | Exact closed form | `rtol=1e-14` |
| `h_0^{(1)}(x) = exp(ix)/(ix)` (Sun Eq. 3.2.33c) | Exact closed form | `rtol=1e-14` |
| Wronskian: `j_n y'_n - j'_n y_n = 1/x²` (Sun Eq. 3.2.32a) | Identity holds for all n | `rtol=1e-12` |
| Wronskian: `j_n h'_n - j'_n h_n = i/x²` (Sun Eq. 3.2.32b) | Identity holds for all n | `rtol=1e-12` |
| Recurrence: `(2n+1)/x · z_n = z_{n-1} + z_{n+1}` (Sun Eq. 3.2.29a) | Holds for n = 0..50 | `rtol=1e-12` |

---

### 0.2 Associated Legendre Polynomials

**Reference:** Sun et al. 2020, Eqs. 3.2.18–3.2.26

**Direct API (`LegendrePolynomials.jl`):**

| Function | Julia API | Notes |
|----------|-----------|-------|
| P_n^m(x) | `Plm(x, l, m)` | Includes Condon-Shortley phase (-1)^m |
| P_n(x) | `Pl(x, l)` | m=0 case |
| All P_n^m up to lmax | `collectPlm(x; m, lmax)` | Vectorized, efficient for batch evaluation |

Pure Julia implementation. Supports `BigFloat` for high-order calculations. The Condon-Shortley phase convention `(-1)^m` (Sun Eq. 3.2.18) is included by default, which matches the Sun et al. 2020 convention.

**Thin compositions (must write):**

- Derivative `dP_n^m/dθ` via recurrence relations (Sun Eq. 3.2.24), rather than numerical differentiation. This is critical for assembling the B and C vector harmonics in Layer 1.

**Overflow warning:** For n > 100, `Float64` overflow is possible in the associated Legendre polynomials. `collectPlm` with `BigFloat` input works correctly. The design should document this boundary and provide a graceful fallback path.

**Validation:**

| Test | Expected | Tolerance |
|------|----------|-----------|
| `P_n(1) = 1` for all n (Sun Eq. 3.2.15) | Exact | `rtol=1e-14` |
| `P_n(-1) = (-1)^n` (Sun Eq. 3.2.15) | Exact | `rtol=1e-14` |
| `P_2^2(cosθ) = 3 sin²θ` (Sun Eq. 3.2.22c) | Exact | `rtol=1e-14` |
| Orthogonality: `∫ P_n^m · P_{n'}^m dx` (Sun Eq. 3.2.23) | `(n+m)!/(n-m)! · 2/(2n+1) · δ_{nn'}` | `rtol=1e-10` |
| Recurrence (Sun Eq. 3.2.24a) for n = 0..50 | Identity holds | `rtol=1e-12` |

---

### 0.3 Wigner 3-j Symbols & Clebsch-Gordan Coefficients

**Reference:** Angular momentum coupling theory; used in surface integrals for mode coupling.

**Direct API (`WignerSymbols.jl`):**

| Function | Julia API | Notes |
|----------|-----------|-------|
| Wigner 3-j | `wigner3j(T, j1, j2, j3, m1, m2, m3)` | Exact rational via prime factorization |
| Clebsch-Gordan | `clebschgordan(T, j1, m1, j2, m2, j3, m3)` | Directly provided |
| 6-j, 9-j symbols | `wigner6j(...)`, not available for 9-j | Use `CGcoefficient.jl` for 9-j if needed |

Key advantage: **exact rational computation** — no floating-point overflow regardless of order. Results computed as `RationalRoot{BigInt}`, convertible to `BigFloat` with arbitrary precision. Thread-safe with LRU cache (default 10^6 entries).

**Note on 9-j symbols:** `WignerSymbols.jl` does not support 9-j symbols (listed as TODO). If Phase 3 (orientation averaging) requires them, use `CGcoefficient.jl` which supports 3-j, 6-j, 9-j and Moshinsky brackets with exact results via `SqrtRational` type.

**Validation:**

| Test | Expected | Tolerance |
|------|----------|-----------|
| Triangular selection rule: `|j1-j2| ≤ j3 ≤ j1+j2` | Zero when violated | Exact |
| m-sum rule: `m1+m2+m3 = 0` | Zero when violated | Exact |
| `{1/2, 1/2, 1; -1/2, 1/2, 0} = 1/√3` | Known closed form | `rtol=1e-14` |
| `{1, 1, 0; -1, 1, 0} = -1/√3` | Known closed form | `rtol=1e-14` |
| Permutation symmetry | `(-1)^{j1+j2+j3}` factor | Exact |

---

### 0.4 Gauss-Legendre Quadrature

**Reference:** Used for numerical integration over particle surfaces in EBCM.

**Direct API (`FastGaussQuadrature.jl`):**

| Function | Julia API | Notes |
|----------|-----------|-------|
| Gauss-Legendre | `gausslegendre(N)` | Returns `(nodes, weights)` on [-1, 1] |
| Gauss-Lobatto | `gausslobatto(N)` | Includes endpoints, useful for some formulations |

O(N) complexity. 100,000 points computed in ~2ms. Numerically stable even for very large N.

**Thin compositions (must write):**

- Interval mapping: `[0, π]` for θ integration, `[0, 2π]` for φ integration — trivial affine transformations, no separate module needed.

**Validation:**

| Test | Expected | Tolerance |
|------|----------|-----------|
| Sum of weights: `Σ w_i = 2` | Exact for [-1,1] | `rtol=1e-14` |
| Exact integration of `x^k` for k < 2N | Analytical value | `rtol=1e-12` |
| Node symmetry about origin | `x_i = -x_{N+1-i}` | Exact |

---

## Layer 1: Composed Functions

Functions that do not exist in any Julia library and must be built by combining Layer 0 primitives.

### 1.1 Gaunt Coefficients

**Reference:** Gaunt coefficient = integral of triple product of spherical harmonics over the full sphere.

The Gaunt coefficient `G(l1,l2,l3; m1,m2,m3)` decomposes into Wigner 3-j symbols via a known closed-form expression:

```
G(l1,l2,l3; m1,m2,m3) = √[(2l1+1)(2l2+1)(2l3+1)/(4π)] ×
    ( l1  l2  l3 )   ( l1  l2  l3 )
    ( 0   0   0  ) × ( m1  m2  m3 )
```

In T-matrix, Gaunt coefficients determine which modes couple in the surface integral — modes that are orthogonal integrate to zero.

**Dependencies:** `WignerSymbols.jl` (Layer 0.3)

**Implementation notes:**
- Selection rules allow early return of zero: triangular condition on l's, m-sum rule, parity check (l1+l2+l3 must be even for the "000" 3-j to be non-zero)
- Optional caching: Gaunt coefficients are constants for given (l1,l2,l3,m1,m2,m3), recomputation is wasteful in iterative solvers

**Validation:**

| Test | Expected | Tolerance |
|------|----------|-----------|
| `G(0,0,0; 0,0,0) = 1/√(4π)` | Known value | `rtol=1e-14` |
| Triangular selection rule | Zero when violated | Exact |
| m-selection rule: `m1+m2 = m3` | Zero when violated | Exact |
| Parity rule: l1+l2+l3 even | Zero when odd | Exact |

---

### 1.2 Vector Spherical Harmonics (P, B, C)

**Reference:** Sun et al. 2020, Eqs. 3.2.39–3.2.45

Three vector spherical harmonics form the angular building blocks of the VSWFs. They are defined as:

**Normalization constants (Sun Eq. 3.2.39):**

- γ'_{mn} = √[(2n+1)(n-m)! / (4π(n+m)!)]
- γ_{mn} = √[(2n+1)(n-m)! / (4πn(n+1)(n+m)!)]

**Scalar harmonic (Sun Eq. 3.2.40):**

- 𝕐_{mn}(θ,φ) = P_n^m(cosθ) · exp(imφ)

**Vector harmonics (Sun Eq. 3.2.41):**

| Symbol | Definition | Physical role |
|--------|-----------|---------------|
| **P_{mn}(θ,φ)** | r̂ · 𝕐_{mn}(θ,φ) | Radial component (Sun Eq. 3.2.41a) |
| **B_{mn}(θ,φ)** | [θ̂ dP_n^m/dθ + φ̂ im/sinθ · P_n^m] exp(imφ) | Tangential gradient type (Sun Eq. 3.2.41b) |
| **C_{mn}(θ,φ)** | [θ̂ im/sinθ · P_n^m - φ̂ dP_n^m/dθ] exp(imφ) | Curl type (Sun Eq. 3.2.41c) |

Note: B_{mn} = r̂ × C_{mn} and C_{mn} = B_{mn} × r̂ (Sun Eq. 3.2.41b-c).

**Dependencies:** `LegendrePolynomials.jl` + derivative recurrence (Layer 0.2)

**Implementation notes:**
- The derivative `dP_n^m/dθ` in B and C must use the stable recurrence from Sun Eq. 3.2.24, not numerical differentiation
- The `im/sinθ` factor requires careful handling near θ = 0 and θ = π (singularity); L'Hôpital's rule or limiting values must be applied
- Complex conjugate relation (Sun Eq. 3.2.41d): `V*_{mn} = (-1)^m (n+m)!/(n-m)! V_{m(-n)}`

**Validation:**

| Test | Expected | Reference |
|------|----------|-----------|
| Cross-orthogonality: ∫ B_{mn} · C*_{m'n'} dΩ = 0 | Zero | Sun Eq. 3.2.42a |
| Self-orthogonality: ∫ B_{mn} · B*_{m'n'} dΩ/(4π) = δ_{mm'}δ_{nn'}/γ²_{mn} | Kronecker delta | Sun Eq. 3.2.42b |
| Self-orthogonality of P: ∫ P_{mn} · P*_{m'n'} dΩ/(4π) = δ_{mm'}δ_{nn'}/(γ'_{mn})² | Kronecker delta | Sun Eq. 3.2.42c |
| Parity: P_{mn}(π-θ, π+φ) = (-1)^{n+1} P_{mn}(θ,φ) | Sign check | Sun Eq. 3.2.43a |
| Parity: C_{mn}(π-θ, π+φ) = (-1)^n C_{mn}(θ,φ) | Sign check | Sun Eq. 3.2.43c |

---

## Layer 2: Vector Spherical Wave Functions & Geometry

The final deliverable of Phase 1. Assembles complete VSWFs from Layer 0 + Layer 1, and provides geometry infrastructure for Phase 2 solvers.

### 2.1 Vector Spherical Wave Functions (M, N, L)

**Reference:** Sun et al. 2020, Eqs. 3.2.46–3.2.51

VSWFs are constructed by combining the vector spherical harmonics (P, B, C) with radial functions (spherical Bessel or Hankel):

**Regular (incident/internal field) VSWFs — use j_n(kr):**

- **RgM_{mn}(kr,θ,φ)** = γ_{mn} · j_n(kr) · C_{mn}(θ,φ) — regular magnetic multipole (TE)
- **RgN_{mn}(kr,θ,φ)** — regular electric multipole (TM); radial part involves `[kr·j_n(kr)]'/(kr)` combined with B_{mn}
- **RgL_{mn}(kr,θ,φ)** — longitudinal component; uses γ'_{mn} and P_{mn}

**Outgoing (scattered field) VSWFs — use h_n^{(1)}(kr):**

- **M_{mn}(kr,θ,φ)** = γ_{mn} · h_n^{(1)}(kr) · C_{mn}(θ,φ) — outgoing magnetic multipole
- **N_{mn}(kr,θ,φ)** — outgoing electric multipole; radial part involves `[kr·h_n^{(1)}(kr)]'/(kr)` combined with B_{mn}
- **L_{mn}(kr,θ,φ)** — longitudinal component

**Dependencies:** Layer 0.1 (Bessel/Hankel) + Layer 1.2 (P/B/C vector harmonics)

**Implementation notes:**
- The radial part of N involves the derivative of `[ρ·z_n(ρ)]` where ρ = kr and z_n is j_n or h_n^{(1)}. This can be computed via the recurrence `d/dρ[ρ·z_n(ρ)] = ρ·z_{n-1}(ρ) - n·z_n(ρ)` or equivalent, avoiding numerical differentiation.
- The "Rg" (regular) prefix denotes using j_n; absence of "Rg" denotes using h_n^{(1)}. This is the standard notation from Mishchenko's convention.

**Validation:**

| Test | Expected | Reference |
|------|----------|-----------|
| Far-field M: `M_{mn} → (-i)^{n+1} e^{iρ}/ρ · γ_{mn} · C_{mn}` | Asymptotic form | Sun Eq. 3.2.46b |
| Far-field N: `N_{mn} → (-i)^n e^{iρ}/ρ · γ_{mn} · B_{mn}` | Asymptotic form | Sun Eq. 3.2.46c |
| Divergence-free M: `∇ · M_{mn} = 0` | Zero | Definition property |
| Curl relation: `∇ × M_{mn} = k · N_{mn}` | Coupling identity | Definition property |
| Plane wave expansion coefficients (Sun Eq. 3.2.51) | Match Mie a_n, b_n for sphere | Cross-check with analytical Mie |
| Parity: `RgM_{mn}(-kr) = (-1)^n RgM_{mn}(kr)` | Sign check | Sun Eq. 3.2.44b |

---

### 2.2 Geometry & Surface Parameterization

Particle shape types for surface integrals in Phase 2 (EBCM, NFM-DS). Each type provides the geometric quantities needed for boundary condition matching.

**Abstract type hierarchy:**

```
AbstractParticle
├── Sphere
├── Spheroid
├── Cylinder
└── ChebyshevParticle
```

**Each concrete type must provide:**

| Function | Returns | Description |
|----------|---------|-------------|
| `radius(p, θ)` | `Float64` | Surface-to-origin distance r(θ) |
| `dradius_dtheta(p, θ)` | `Float64` | dr/dθ for surface normal computation |
| `surface_element(p, θ)` | `Float64` | Jacobian × sinθ for integration |
| `surface_normal(p, θ)` | `SVector{3}` | Outward unit normal vector |

**Shape definitions:**

| Type | Surface equation | Notes |
|------|-----------------|-------|
| Sphere | r(θ) = a | Trivial; benchmark case |
| Spheroid | r(θ) = ab / √(a²sin²θ + b²cos²θ) | a = equatorial, b = polar radius |
| Cylinder | Piecewise: top cap + side + bottom cap | Requires θ segmentation |
| ChebyshevParticle | r(θ) = a[1 + ε·T_n(cosθ)] | Deformed sphere; T_n = Chebyshev polynomial |

**Dependencies:** None beyond Julia Base. Analytic formulas only — no mesh libraries needed for Phase 1.

**Validation:**

| Test | Expected |
|------|----------|
| Sphere surface area = 4πa² | Analytical value |
| Sphere volume = 4πa³/3 | Analytical value |
| Spheroid surface area vs. known formula | Reference tables |
| Surface normal · r̂ > 0 everywhere | Outward-pointing check |
| Cylinder: continuity at segment boundaries | Smooth joining |

---

## Convention Summary

The following conventions are adopted throughout Phase 1, following Sun et al. 2020:

| Convention | Choice | Impact |
|-----------|--------|--------|
| Time-harmonic factor | exp(-iωt) | Standard in optics; determines h_n^{(1)} for outgoing waves |
| Condon-Shortley phase | Included in P_n^m | Matches LegendrePolynomials.jl default |
| Spherical harmonic normalization | Unnormalized 𝕐_{mn} = P_n^m exp(imφ) | Sun Eq. 3.2.40; standard Y_{lm} differs by a factor |
| VSWF normalization | γ_{mn}, γ'_{mn} per Sun Eq. 3.2.39 | Propagates into all cross-section formulas in Phase 3 |
| Mode ordering | Deferred to implementation | Will be specified per solver in Phase 2 |

---

## Existing Julia Ecosystem Note

**TransitionMatrices.jl** (v0.3.1, JuliaRemoteSensing organization) already implements EBCM T-matrix for spheroids, cylinders, and Chebyshev particles with arbitrary precision support. It uses `FastGaussQuadrature.jl`, `ForwardDiff.jl`, and `Wigxjpf` (C library for Wigner symbols). Key differences from our approach:

- Uses `Wigxjpf` (C FFI) rather than `WignerSymbols.jl` (pure Julia, exact rational)
- Does not appear to follow Sun et al. 2020 conventions explicitly
- Limited to axisymmetric shapes in EBCM; IITM for arbitrary shapes
- Still at v0.3.x maturity

Our library differentiates through: explicit reference to Sun et al. 2020 conventions, pure Julia Wigner symbols with exact arithmetic, cleaner Layer 0 dependency strategy (direct API calls, no wrappers), and a multi-solver architecture (Phase 2) including SVM, IITM, and MSTM beyond EBCM.

---

## Julia Dependency Summary

| Package | Version | Role | AD-compatible | BigFloat |
|---------|---------|------|---------------|----------|
| SpecialFunctions.jl | ≥ 2.7 | Spherical Bessel functions | No (C FFI) | Partial |
| LegendrePolynomials.jl | ≥ 0.4 | Associated Legendre polynomials | Likely (pure Julia) | Yes |
| WignerSymbols.jl | ≥ 2.0 | Wigner 3-j, CG coefficients | N/A (discrete inputs) | Yes (exact rational) |
| FastGaussQuadrature.jl | ≥ 1.0 | Gauss-Legendre quadrature | N/A (precomputed constants) | Limited |
