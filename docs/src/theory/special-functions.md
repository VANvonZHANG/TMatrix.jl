# Special Functions

This page defines the special functions used as mathematical primitives in TMatrix.jl.

## Spherical Bessel Functions

### Definitions

The spherical Bessel functions of the first and second kind:

```math
j_n(x) = \sqrt{\frac{\pi}{2x}} J_{n+1/2}(x)
```

```math
y_n(x) = \sqrt{\frac{\pi}{2x}} Y_{n+1/2}(x)
```

### Spherical Hankel Functions

```math
h_n^{(1)}(x) = j_n(x) + i \, y_n(x) \quad \text{(outgoing)}
```

```math
h_n^{(2)}(x) = j_n(x) - i \, y_n(x) \quad \text{(incoming)}
```

### Derivatives

The derivative is computed via recurrence to avoid numerical differentiation:

```math
z_n'(x) = z_{n-1}(x) - \frac{n+1}{x} z_n(x)
```

For ``n = 0``: ``z_0'(x) = -z_1(x)``.

**Code mapping:**
- [`shankelh1`](@ref) → ``h_n^{(1)}(x)``
- [`sbesselj_deriv`](@ref) → ``j_n'(x)``
- [`shankelh1_deriv`](@ref) → ``{h_n^{(1)}}'(x)``

## Associated Legendre Polynomials

```math
P_n^m(\cos\theta) \quad \text{with Condon-Shortley phase } (-1)^m
```

### θ-Derivative

```math
\frac{dP_n^m}{d\theta} = \frac{n \cos\theta}{\sin\theta} P_n^m - \frac{n+m}{\sin\theta} P_{n-1}^m
```

**Code mapping:**
- [`dplmdtheta`](@ref) → ``dP_n^m/d\theta``
- [`plm_over_sintheta`](@ref) → ``P_n^m / \sin\theta``

## Gauss-Legendre Quadrature

For numerical integration over particle surfaces, Gauss-Legendre quadrature maps
nodes from ``[-1, 1]`` to the required angular intervals:

- ``\theta`` integration: ``[0, \pi]``
- ``\phi`` integration: ``[0, 2\pi]``

**Code mapping:**
- [`gl_theta`](@ref) → ``\theta`` quadrature
- [`gl_phi`](@ref) → ``\phi`` quadrature
