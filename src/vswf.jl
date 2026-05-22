# Layer 2.1: Vector Spherical Wave Functions M, N, L
# (Sun §3.2, Eqs. 3.2.46–3.2.51)

# Select radial function: j_n for regular, h_n^(1) for outgoing
function _radial_func(n::Int, rho::Real; regular::Bool = true)
    regular && return sphericalbesselj(n, rho)
    return shankelh1(n, rho)
end

"""
    vswf_M(n::Int, m::Int, kr::Real, theta::Real, phi::Real; regular::Bool=true)

Vector spherical wave function ``M_{mn}(\\rho, \\theta, \\phi)`` — magnetic multipole (TE).

```
M_{mn}(\\rho, \\theta, \\phi) = \\gamma_{mn} \\cdot z_n(\\rho) \\cdot C_{mn}(\\theta, \\phi)
```

Has no radial component. `regular=true` uses ``j_n(kr)`` (regular at origin);
`regular=false` uses ``h_n^{(1)}(kr)`` (outgoing wave).

# Arguments
- `n::Int`: Multipole order (must be ≥ 1)
- `m::Int`: Azimuthal mode number
- `kr::Real`: Dimensionless radial coordinate ``\\rho = kr``
- `theta::Real`: Polar angle
- `phi::Real`: Azimuthal angle
- `regular::Bool=true`: Use regular (`j_n`) or outgoing (`h_n^{(1)}`) radial function

# Returns
`SVector{3, Complex}`: M-mode VSWF in spherical components `(r̂, θ̂, φ̂)`

# References
- Sun et al. 2020, Eq. 3.2.46
"""
function vswf_M(n::Int, m::Int, kr::Real, theta::Real, phi::Real; regular::Bool = true)
    n >= 1 || throw(ArgumentError("vswf_M requires n >= 1, got n=$n"))
    g = gamma_vsh(n, m)
    zn = _radial_func(n, kr; regular)
    return g * zn * C_vsh(n, m, theta, phi)
end

"""
    vswf_N(n::Int, m::Int, kr::Real, theta::Real, phi::Real; regular::Bool=true)

Vector spherical wave function ``N_{mn}(\\rho, \\theta, \\phi)`` — electric multipole (TM).

```
N_{mn} = \\gamma_{mn} \\{ \\hat{r} \\frac{n(n+1)}{\\rho} z_n(\\rho) P_{mn}(\\theta, \\phi)
         + B_{mn}(\\theta, \\phi) \\frac{[\\rho z_n(\\rho)]'}{\\rho} \\}
```

Has all three spherical components. `regular=true` uses ``j_n(kr)``;
`regular=false` uses ``h_n^{(1)}(kr)``.

# Arguments
- `n::Int`: Multipole order (must be ≥ 1)
- `m::Int`: Azimuthal mode number
- `kr::Real`: Dimensionless radial coordinate ``\\rho = kr``
- `theta::Real`: Polar angle
- `phi::Real`: Azimuthal angle
- `regular::Bool=true`: Use regular or outgoing radial function

# Returns
`SVector{3, Complex}`: N-mode VSWF in spherical components `(r̂, θ̂, φ̂)`

# References
- Sun et al. 2020, Eq. 3.2.46
"""
function vswf_N(n::Int, m::Int, kr::Real, theta::Real, phi::Real; regular::Bool = true)
    n >= 1 || throw(ArgumentError("vswf_N requires n >= 1, got n=$n"))
    g = gamma_vsh(n, m)
    zn = _radial_func(n, kr; regular)
    zn_prev = _radial_func(n - 1, kr; regular)
    # Radial component: γ_{mn} n(n+1)/(kr) z_n P_{mn}
    r_comp = g * n * (n + 1) / kr * zn * Plm(cos(theta), n, m) * exp(im * m * phi)
    # Tangential component: γ_{mn} [ρz_n]'/ρ · B_{mn}
    riccati = riccati_zn_prime(zn_prev, zn, n, kr)
    tang = g * riccati * B_vsh(n, m, theta, phi)
    T = complex(typeof(float(kr)))
    return SVector{3, T}(r_comp, tang[2], tang[3])
end

"""
    vswf_L(n::Int, m::Int, kr::Real, theta::Real, phi::Real; regular::Bool=true)

Longitudinal vector spherical wave function ``L_{mn}(\\rho, \\theta, \\phi)``.

Used for plane wave expansion. Uses ``\\gamma'_{mn}`` normalization instead of ``\\gamma_{mn}``:
```
L_{mn} = \\frac{1}{k} \\nabla [z_n(kr) Y_{mn}(\\theta, \\phi)]
```

In spherical components (Sun Eq. 3.2.50):
```
L_{mn} = \\gamma'_{mn} \\{ \\hat{r} z_n'(\\rho) P_n^m e^{im\\phi} + B_{mn} \\frac{z_n(\\rho)}{\\rho} \\}
```

# Arguments
- `n::Int`: Multipole order (must be ≥ 1)
- `m::Int`: Azimuthal mode number
- `kr::Real`: Dimensionless radial coordinate ``\\rho = kr``
- `theta::Real`: Polar angle
- `phi::Real`: Azimuthal angle
- `regular::Bool=true`: Use regular or outgoing radial function

# Returns
`SVector{3, Complex}`: L-mode VSWF in spherical components `(r̂, θ̂, φ̂)`

# References
- Sun et al. 2020, Eq. 3.2.50
"""
function vswf_L(n::Int, m::Int, kr::Real, theta::Real, phi::Real; regular::Bool = true)
    n >= 1 || throw(ArgumentError("vswf_L requires n >= 1, got n=$n"))
    gp = gamma_prime(n, m)
    zn = _radial_func(n, kr; regular)
    zn_deriv = regular ? sbesselj_deriv(n, kr) : shankelh1_deriv(n, kr)
    # Radial: γ'_{mn} z_n'(kr) P_n^m exp(imφ)
    r_comp = gp * zn_deriv * Plm(cos(theta), n, m) * exp(im * m * phi)
    # Tangential: γ'_{mn} z_n(kr)/(kr) B_{mn}(θ,φ)
    tang = gp * zn / kr * B_vsh(n, m, theta, phi)
    T = complex(typeof(float(kr)))
    return SVector{3, T}(r_comp, tang[2], tang[3])
end
