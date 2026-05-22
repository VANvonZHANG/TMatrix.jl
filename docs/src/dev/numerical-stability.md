# Numerical Stability

This page documents numerical stability considerations and edge-case handling
in TMatrix.jl.

## Singularities at θ = 0 and θ = π

The vector spherical harmonics contain terms like ``P_n^m(\cos\theta) / \sin\theta``
and ``dP_n^m/d\theta``. Both become singular as ``\theta \to 0`` or ``\theta \to \pi``.

### P_n^m / sinθ

For ``m = 0``: genuinely singular, returns `Inf`.

For ``m = 1``: the limit is finite:

```math
\lim_{\theta \to 0} \frac{P_n^1(\cos\theta)}{\sin\theta} = -\frac{n(n+1)}{2}
```

For ``m \geq 2``: the limit is zero.

**Implementation:** [`plm_over_sintheta`](@ref) switches to `_plm_over_sintheta_limit`
when `abs(sin(theta)) < 1e-6`.

### dP_n^m/dθ

Handled via `_dplmdtheta_limit`, which computes a numerical derivative
with a small step size (``h = 10^{-7}``) when `abs(sin(theta)) < 1e-10`.

## Overflow in Associated Legendre Polynomials

For ``n > 100``, `Float64` overflow is possible in `Plm`. The `LegendrePolynomials.jl`
package supports `BigFloat` input for arbitrary precision.

**Current policy:** Document the limitation. Future work may add automatic `BigFloat`
fallback.

## Derivative Computation

All derivatives are computed analytically via recurrence relations, never by
numerical differentiation. This avoids catastrophic cancellation for large ``n``
or small ``\rho``.

Examples:
- ``z_n'(x) = z_{n-1}(x) - (n+1)/x \cdot z_n(x)``
- ``[\rho z_n(\rho)]'/\rho = z_{n-1}(\rho) - n/\rho \cdot z_n(\rho)``

## Gaunt Coefficient Selection Rules

The [`gaunt`](@ref) function returns `0.0` early when selection rules are violated:
- Triangular condition: `|l1-l2| > l3` or `l3 > l1+l2`
- `m1 + m2 + m3 != 0`
- `l1 + l2 + l3` is odd

This prevents unnecessary Wigner 3-j computation.
