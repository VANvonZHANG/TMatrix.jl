# T-Matrix.jl Roadmap

> A comprehensive Julia library for T-matrix electromagnetic scattering computations, providing a unified framework for analytical, boundary-integral, and volume-differential solvers.

## Vision

Build a modern, extensible, and high-performance T-matrix scattering library in Julia that unifies multiple solution methods (Mie, SVM, EBCM, NFM-DS, IITM, MSTM) under a single coherent interface. The library targets applications in atmospheric science, remote sensing, and electromagnetic scattering research, offering both numerical stability and developer productivity that legacy Fortran implementations cannot match.

---

## Motivation: Why Julia

The choice of Julia addresses a fundamental tension in computational electromagnetics: the tradeoff between development efficiency (Python/MATLAB) and runtime performance (C++/Fortran). Julia resolves this two-language problem and provides specific advantages for T-matrix algorithms.

**World-class ODE ecosystem for IITM.** The Invariant Imbedding T-matrix Method (IITM) solves large-scale complex matrix Riccati ODEs along the radial direction. Julia's [DifferentialEquations.jl](https://diffeq.sciml.ai/) provides high-order Runge-Kutta methods, adaptive step control, stiff solvers, and GPU-accelerated ODE solving out of the box -- reducing IITM development from months to days.

**Multiple dispatch and arbitrary precision.** EBCM algorithms frequently encounter ill-conditioned matrix inversions for non-spherical particles. Julia's generic programming allows the same code to operate on `Float64`/`ComplexF64` for speed, and seamlessly switch to `BigFloat` when numerical instability is detected -- trading speed for robustness without code changes.

**Native complex number and linear algebra support.** T-matrix computations involve massive complex matrix multiplications, inversions, and SVDs. Julia has native complex types with zero overhead and seamless integration with OpenBLAS/MKL, achieving performance on par with hand-tuned C/Fortran.

**Automatic differentiation.** For future inverse problems (e.g., retrieving ice crystal shape from radar Mueller matrices), Julia's AD ecosystem ([Zygote.jl](https://github.com/FluxML/Zygote.jl), [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)) enables end-to-end differentiation of the T-matrix simulation pipeline -- something impractical in Fortran and difficult in Python/JAX for complex physical models.

---

## Development Phases

The library is built bottom-up in four phases. Each phase produces testable, usable output before the next begins.

### Phase 1: Foundation -- Mathematical and Physical Primitives

> **Status:** Not Started | **ETA:** 1-2 months

Building the numerical bedrock. Every downstream solver depends on these primitives being correct and numerically stable.

**Deliverables:**

- [ ] **Special Functions Module**
  - High-precision spherical Bessel and Hankel functions
  - Associated Legendre polynomials (using logarithmic derivatives or backward recurrence for high-order numerical stability)

- [ ] **Mode Coupling Coefficients**
  - Wigner 3-j symbols
  - Gaunt coefficients / Clebsch-Gordan coefficients

- [ ] **Vector Spherical Wave Functions (VSWFs)**
  - M_mn and N_mn mode functions composed from the above primitives

- [ ] **Geometry and Quadrature Infrastructure**
  - Gauss-Legendre quadrature nodes and weights
  - Particle shape types (`Sphere`, `Spheroid`, `Cylinder`, `CustomShape`) with surface parameterization

---

### Phase 2: Core Solvers -- Hierarchical Solver Family

> **Status:** Not Started | **ETA:** 3-4 months

Solvers are organized in a layered hierarchy. Higher-level solvers can consume results from lower-level ones or serve as fallbacks under extreme conditions.

#### L0: Analytical Solvers (Fastest, Narrowest Scope)

Benchmark standards for idealized geometries.

- [ ] **Mie Solver** -- Homogeneous and multi-layered concentric spheres; produces exact diagonal T-matrix
- [ ] **SVM (Separation of Variables) Solver** -- Perfect spheroids via vector spheroidal wave functions, with spheroidal-to-spherical wave translation interface

#### L1: Boundary Integral Solvers (Fast, Homogeneous/Layered Non-spherical)

The workhorses for meteorology and remote sensing, based on surface boundary matching.

- [ ] **EBCM (Extended Boundary Condition Method)** -- Moderate aspect ratios (< 3:1), smooth non-spherical particles (Chebyshev particles, short cylinders)
- [ ] **NFM-DS (Null-Field Method with Discrete Sources)** -- EBCM upgrade with multipole discrete source placement for extreme aspect ratios (e.g., 20:1 elongated ice crystals or dust)

#### L2: Volume Differential Solvers (Heavy, Complex/Inhomogeneous)

For complex materials and aerosols with 3D inhomogeneous properties.

- [ ] **IITM (Invariant Imbedding T-matrix Method)** -- Integration with DifferentialEquations.jl for particles with complex refractive index gradients, random inclusions, or extremely rough surfaces

#### L3: Aggregate Solvers (Multi-body Problems)

- [ ] **MSTM (Superposition T-matrix Method)** -- Multi-sphere aggregates (e.g., soot clusters) via Translation Addition Theorem for inter-body multiple scattering coupling; outputs global cluster T-matrix

#### Hybrid Strategies: Solver Fusion via Multiple Dispatch

Leverage multiple dispatch to pass boundary conditions or initial T-matrices between solvers, achieving 1+1>2.

- [ ] **SVM -> IITM (Regular Core + Rough Shell)** -- For particles with a perfect spheroidal core wrapped in a rough ice layer: first solve the core via SVM for the fast T-matrix, then pass it as ODE initial state `T(r0)` to IITM for the outer shell perturbation
- [ ] **MSTM -> EBCM (Aggregate Effective Medium)** -- Compute the exact T-matrix of a multi-sphere aggregate via MSTM, extract equivalent homogeneous parameters by inversion, then use EBCM for macroscopic averaging

---

### Phase 3: Post-Processing and Macroscopic Quantities

> **Status:** Not Started | **ETA:** 1 month

Transform raw T-matrices into experimentally measurable quantities.

**Deliverables:**

- [ ] **Fixed-Orientation Quantities**
  - Extinction, scattering, and absorption cross-sections (C_ext, C_sca, C_abs)
  - Amplitude scattering matrix S and Mueller matrix M at arbitrary scattering angles

- [ ] **Rotation and Orientation Averaging**
  - Fast Wigner D-matrix computation
  - Analytical random orientation averaging (Mishchenko's algorithm) for ensemble-averaged cross-sections and phase functions

---

### Phase 4: Ecosystem -- API, Optimization, and Open-Source Release

> **Status:** Not Started | **ETA:** Ongoing

**Deliverables:**

- [ ] **API Design** -- Idiomatic Julia interface:
  ```julia
  particle = Spheroid(a=1.0, c=2.0, refractive_index=1.5+0.01im)
  t_matrix = solve_tmatrix(particle, method=IITM(), N_max=15)
  cross_sections = calc_cross_sections(t_matrix, wavelength=0.5)
  ```

- [ ] **Performance Optimization**
  - `@code_warntype` audits for type stability
  - StaticArrays.jl for small matrix operations
  - `@threads` parallelism across wavelengths or matrix rows

- [ ] **Benchmarking** -- Rigorous accuracy and speed validation against NASA Fortran reference implementations (Mishchenko's T-matrix code)

---

## Architecture Guideline

Follow the **Traits / Abstract Type** design pattern:

- Define `AbstractTMatrixMethod` as the supertype
- Each solver (`EBCM`, `IITM`, `SVM`, `NFM_DS`, `MSTM`) implements this interface
- Post-processing code (cross-sections, Mueller matrices) operates generically on `AbstractTMatrixMethod` outputs, requiring zero duplication across solvers

This is the core architectural advantage of Julia for scientific computing: a shared post-processing pipeline that works with any solver through multiple dispatch.
