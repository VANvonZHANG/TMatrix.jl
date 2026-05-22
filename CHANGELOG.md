# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Comprehensive Documenter.jl documentation site with tutorials, API reference, theory guides, and developer docs
- Docstrings for all public functions across special_functions.jl, composed.jl, vswf.jl, and geometry.jl

## [0.1.0] - 2026-05-20

### Added

- Initial Julia package scaffolding with Project.toml, test suite, and CI workflows
- Phase 1 mathematical primitives: associated Legendre functions, spherical Bessel/Hankel functions, Wigner 3-j symbols, Gaunt coefficients
- Vector spherical wave functions (VSWF): M, N, L types with regular and outgoing variants
- Gauss-Legendre quadrature integration and azimuthal quadrature
- Spheroid surface parameterization and outward surface normals
- Defensive input validation for singular and edge-case inputs
- GitHub Actions CI (multi-version/platform testing), TagBot, CompatHelper, JuliaFormatter
- Aqua.jl code quality checks in test suite
