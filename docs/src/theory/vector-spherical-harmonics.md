# Vector Spherical Harmonics

The vector spherical harmonics ``P_{mn}``, ``B_{mn}``, and ``C_{mn}`` form the
angular building blocks of the vector spherical wave functions.

## Definitions

### Scalar Harmonic

```math
Y_{mn}(\theta, \phi) = P_n^m(\cos\theta) \, e^{im\phi}
```

### Vector Harmonic P (Radial)

```math
P_{mn}(\theta, \phi) = \hat{r} \, Y_{mn}(\theta, \phi)
```

### Vector Harmonic B (Tangential Gradient)

```math
B_{mn}(\theta, \phi) = \left[ \hat{\theta} \frac{dP_n^m}{d\theta} + \hat{\phi} \frac{im}{\sin\theta} P_n^m \right] e^{im\phi}
```

### Vector Harmonic C (Curl)

```math
C_{mn}(\theta, \phi) = \left[ \hat{\theta} \frac{im}{\sin\theta} P_n^m - \hat{\phi} \frac{dP_n^m}{d\theta} \right] e^{im\phi}
```

## Relations

```math
B_{mn} = \hat{r} \times C_{mn}, \quad C_{mn} = B_{mn} \times \hat{r}
```

## Orthogonality

```math
\int B_{mn} \cdot C^*_{m'n'} \, d\Omega = 0
```

```math
\frac{1}{4\pi} \int B_{mn} \cdot B^*_{m'n'} \, d\Omega = \frac{\delta_{mm'} \delta_{nn'}}{\gamma_{mn}^2}
```

```math
\frac{1}{4\pi} \int P_{mn} \cdot P^*_{m'n'} \, d\Omega = \frac{\delta_{mm'} \delta_{nn'}}{\gamma_{mn}^{'2}}
```

**Code mapping:**
- [`P_vsh`](@ref) → ``P_{mn}``
- [`B_vsh`](@ref) → ``B_{mn}``
- [`C_vsh`](@ref) → ``C_{mn}``
