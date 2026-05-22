# Vector Spherical Wave Functions

The vector spherical wave functions (VSWFs) ``M_{mn}``, ``N_{mn}``, and ``L_{mn}``
are the complete solutions to the vector Helmholtz equation in spherical coordinates.

## Definitions

Let ``\rho = kr`` and ``z_n(\rho)`` be the radial function (``j_n`` for regular,
``h_n^{(1)}`` for outgoing).

### M-mode (Magnetic Multipole, TE)

```math
M_{mn}(\rho, \theta, \phi) = \gamma_{mn} \, z_n(\rho) \, C_{mn}(\theta, \phi)
```

No radial component. Represents transverse electric (TE) modes.

### N-mode (Electric Multipole, TM)

```math
N_{mn}(\rho, \theta, \phi) = \gamma_{mn} \left\{ \hat{r} \frac{n(n+1)}{\rho} z_n(\rho) P_{mn}(\theta, \phi) + \frac{[\rho z_n(\rho)]'}{\rho} B_{mn}(\theta, \phi) \right\}
```

Has all three components. Represents transverse magnetic (TM) modes.

The radial derivative term:

```math
\frac{[\rho z_n(\rho)]'}{\rho} = z_{n-1}(\rho) - \frac{n}{\rho} z_n(\rho)
```

### L-mode (Longitudinal)

```math
L_{mn}(\rho, \theta, \phi) = \frac{1}{k} \nabla \left[ z_n(kr) Y_{mn}(\theta, \phi) \right]
```

Used for plane wave expansion. Uses ``\gamma'_{mn}`` normalization.

## Properties

- Divergence-free: ``\nabla \cdot M_{mn} = 0``, ``\nabla \cdot N_{mn} = 0``
- Curl relation: ``\nabla \times M_{mn} = k \, N_{mn}``
- Far-field asymptotics:
  - ``M_{mn} \to (-i)^{n+1} \frac{e^{i\rho}}{\rho} \gamma_{mn} C_{mn}``
  - ``N_{mn} \to (-i)^n \frac{e^{i\rho}}{\rho} \gamma_{mn} B_{mn}``

## Regular vs. Outgoing

| Prefix | Radial function | Use case |
|--------|-----------------|----------|
| Regular (`Rg`) | ``j_n(kr)`` | Incident fields, internal fields |
| Outgoing | ``h_n^{(1)}(kr)`` | Scattered fields |

**Code mapping:**
- [`vswf_M`](@ref) → ``M_{mn}``
- [`vswf_N`](@ref) → ``N_{mn}``
- [`vswf_L`](@ref) → ``L_{mn}``
