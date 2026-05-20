# Layer 2.1: Vector Spherical Wave Functions M, N, L
# (Sun §3.2, Eqs. 3.2.46–3.2.51)

# Select radial function: j_n for regular, h_n^(1) for outgoing
function _radial_func(n::Int, rho::Real; regular::Bool=true)
    regular && return sphericalbesselj(n, rho)
    return shankelh1(n, rho)
end

# M_{mn}(ρ,θ,φ) = γ_{mn} z_n(ρ) C_{mn}(θ,φ)
# Magnetic multipole (TE) — no radial component (Sun Eq. 3.2.46)
function vswf_M(n::Int, m::Int, kr::Real, theta::Real, phi::Real; regular::Bool=true)
    n >= 1 || throw(ArgumentError("vswf_M requires n >= 1, got n=$n"))
    g = gamma_vsh(n, m)
    zn = _radial_func(n, kr; regular)
    return g * zn * C_vsh(n, m, theta, phi)
end

# N_{mn}(ρ,θ,φ) = γ_{mn} { r̂ n(n+1)/ρ · z_n(ρ) P_{mn} + B_{mn} [ρz_n]'/ρ }
# Electric multipole (TM) — all three components (Sun Eq. 3.2.46)
function vswf_N(n::Int, m::Int, kr::Real, theta::Real, phi::Real; regular::Bool=true)
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

# L_{mn}(ρ,θ,φ) — longitudinal VSWF for plane wave expansion
# L = (1/k)∇[z_n(kr) Y_{mn}] with γ'_{mn} normalization (Sun Eq. 3.2.50)
function vswf_L(n::Int, m::Int, kr::Real, theta::Real, phi::Real; regular::Bool=true)
    n >= 1 || throw(ArgumentError("vswf_L requires n >= 1, got n=$n"))
    gp = gamma_prime(n, m)
    zn = _radial_func(n, kr; regular)
    if n == 0
        zn_deriv = regular ? sbesselj_deriv(0, kr) : shankelh1_deriv(0, kr)
    else
        zn_deriv = regular ? sbesselj_deriv(n, kr) : shankelh1_deriv(n, kr)
    end
    # Radial: γ'_{mn} z_n'(kr) P_n^m exp(imφ)
    r_comp = gp * zn_deriv * Plm(cos(theta), n, m) * exp(im * m * phi)
    # Tangential: γ'_{mn} z_n(kr)/(kr) B_{mn}(θ,φ)
    tang = gp * zn / kr * B_vsh(n, m, theta, phi)
    T = complex(typeof(float(kr)))
    return SVector{3, T}(r_comp, tang[2], tang[3])
end
