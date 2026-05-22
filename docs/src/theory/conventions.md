# Conventions

This page summarizes the conventions adopted throughout TMatrix.jl,
following **Sun et al. 2020, *T-matrix Concept***.

## Time-Harmonic Factor

The time dependence is ``e^{-i\omega t}``. This is the standard convention in optics
and determines that outgoing scattered waves use ``h_n^{(1)}(kr)``.

## Spherical Harmonics

### Condon-Shortley Phase

The associated Legendre polynomials ``P_n^m(x)`` include the Condon-Shortley phase
``(-1)^m``:

```math
P_n^m(x) = (-1)^m (1 - x^2)^{m/2} \frac{d^m}{dx^m} P_n(x)
```

This matches the default convention of `LegendrePolynomials.jl`.

### Unnormalized Scalar Harmonics

The scalar spherical harmonics used in T-matrix theory are unnormalized:

```math
Y_{mn}(\theta, \phi) = P_n^m(\cos\theta) \, e^{im\phi}
```

This differs from the standard physics normalization ``Y_{lm}`` by a factor of
``\sqrt{(2n+1)(n-m)! / (4\pi(n+m)!)}``.

## Vector Spherical Harmonic Normalization

Two normalization constants appear:

```math
\gamma'_{mn} = \sqrt{\frac{(2n+1)(n-m)!}{4\pi(n+m)!}}
```

```math
\gamma_{mn} = \sqrt{\frac{(2n+1)(n-m)!}{4\pi n(n+1)(n+m)!}}
```

## Mode Ordering

Mode indices:
- ``n = 1, 2, \ldots, N_{\max}``: Multipole order (degree)
- ``m = -n, \ldots, n``: Azimuthal mode number

## Spherical Coordinates

- ``r``: radial distance
- ``\theta``: polar angle from the ``+z`` axis, ``\theta \in [0, \pi]``
- ``\phi``: azimuthal angle in the ``xy``-plane from the ``+x`` axis, ``\phi \in [0, 2\pi]``

The unit vectors ``(\hat{r}, \hat{\theta}, \hat{\phi})`` form a right-handed system.
