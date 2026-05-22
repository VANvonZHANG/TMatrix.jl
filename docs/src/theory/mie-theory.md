# Mie Theory

Lorenz-Mie scattering describes the exact solution for electromagnetic plane wave
scattering by a homogeneous sphere.

## Problem Setup

A plane wave propagating along the ``z``-axis illuminates a sphere of radius ``a``
and refractive index ``m = m' + i m''``, embedded in a medium with refractive index
``n_{med}`` (usually ``n_{med} = 1`` for vacuum/air).

The size parameter is:

```math
x = ka = \frac{2\pi a}{\lambda}
```

where ``\lambda`` is the wavelength in the surrounding medium.

## Mie Coefficients

The scattered field is expanded in outgoing VSWFs, and the internal field in regular VSWFs.
Matching boundary conditions at ``r = a`` yields the Mie coefficients:

```math
a_n = \frac{m \psi_n(mx) \psi_n'(x) - \psi_n(x) \psi_n'(mx)}{m \psi_n(mx) \xi_n'(x) - \xi_n(x) \psi_n'(mx)}
```

```math
b_n = \frac{\psi_n(mx) \psi_n'(x) - m \psi_n(x) \psi_n'(mx)}{\psi_n(mx) \xi_n'(x) - m \xi_n(x) \psi_n'(mx)}
```

where:
- ``\psi_n(\rho) = \rho j_n(\rho)`` (Riccati-Bessel function)
- ``\xi_n(\rho) = \rho h_n^{(1)}(\rho)`` (Riccati-Hankel function)

## T-matrix for a Sphere

For a sphere, the T-matrix is diagonal:

```math
T^{11}_{nmnm} = -b_n, \quad T^{22}_{nmnm} = -a_n
```

All off-diagonal elements are zero due to spherical symmetry.

## Cross-Sections

### Dimensional Cross-Sections

```math
C_{ext} = \frac{2\pi}{k^2} \sum_{n=1}^{\infty} (2n+1) \, \text{Re}(a_n + b_n)
```

```math
C_{sca} = \frac{2\pi}{k^2} \sum_{n=1}^{\infty} (2n+1) \left( |a_n|^2 + |b_n|^2 \right)
```

```math
C_{abs} = C_{ext} - C_{sca}
```

### Efficiencies

```math
Q_{ext} = \frac{C_{ext}}{\pi a^2}, \quad Q_{sca} = \frac{C_{sca}}{\pi a^2}, \quad Q_{abs} = \frac{C_{abs}}{\pi a^2}
```

### Asymmetry Parameter

```math
g = \frac{4\pi}{k^2 C_{sca}} \sum_{n=1}^{\infty} \frac{n(n+2)}{n+1} \, \text{Re}(a_n a_{n+1}^* + b_n b_{n+1}^*) + \frac{2n+1}{n(n+1)} \, \text{Re}(a_n b_n^*)
```

## Truncation

The infinite series is truncated at ``N_{\max}`` using the Wiscombe criterion:

```math
N_{\max} = \lceil x + 4 x^{1/3} + 2 \rceil
```

**Code mapping:**
- [`mie_ab`](@ref) → computes ``a_n``, ``b_n``
- [`mie_nmax`](@ref) → Wiscombe truncation
- [`solve_tmatrix`](@ref) → assembles the diagonal T-matrix
- [`calc_cross_sections`](@ref) → computes cross-sections from T-matrix
