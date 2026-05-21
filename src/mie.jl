# Layer 0.5: Lorenz-Mie Coefficients (Sun §3.3, Eq. 3.3.41)

# ---- Solver dispatch tag ----------------------------------------------------

"""
    MieMethod

Solver dispatch tag for Lorenz-Mie (analytical sphere) T-matrix computation.
Passed as the second argument to `solve_tmatrix`.
"""
struct MieMethod end

# ---- Compact diagonal T-matrix storage --------------------------------------

"""
    MieTMatrix

Compact diagonal storage of the Lorenz-Mie T-matrix.

Fields
- `a::Vector{ComplexF64}`: T^22 diagonal, `T^22_{nmnm} = -a_n`
- `b::Vector{ComplexF64}`: T^11 diagonal, `T^11_{nmnm} = -b_n`
- `k::Float64`: wave number in the surrounding medium
- `radius::Float64`: sphere radius
- `N_max::Int`: maximum multipole order
"""
struct MieTMatrix
    a::Vector{ComplexF64}   # T^22 diagonal: T^22_{nmnm} = -a_n
    b::Vector{ComplexF64}   # T^11 diagonal: T^11_{nmnm} = -b_n
    k::Float64              # wave number in surrounding medium
    radius::Float64
    N_max::Int
end

# ---- Cross-sections container -----------------------------------------------

"""
    CrossSections

Scattering cross-sections and efficiencies for a Mie scatterer.

Fields
- `extinction::Float64`:  C_ext
- `scattering::Float64`:  C_sca
- `absorption::Float64`:  C_abs
- `Q_ext::Float64`:       extinction efficiency
- `Q_sca::Float64`:       scattering efficiency
- `Q_abs::Float64`:       absorption efficiency
- `asymmetry::Float64`:   g (asymmetry parameter)
"""
struct CrossSections
    extinction::Float64      # C_ext
    scattering::Float64      # C_sca
    absorption::Float64      # C_abs
    Q_ext::Float64           # extinction efficiency
    Q_sca::Float64           # scattering efficiency
    Q_abs::Float64           # absorption efficiency
    asymmetry::Float64       # g
end

function Base.show(io::IO, cs::CrossSections)
    println(io, "CrossSections:")
    println(io, "  C_ext = $(cs.extinction)")
    println(io, "  C_sca = $(cs.scattering)")
    println(io, "  C_abs = $(cs.absorption)")
    println(io, "  Q_ext = $(cs.Q_ext)")
    println(io, "  Q_sca = $(cs.Q_sca)")
    println(io, "  Q_abs = $(cs.Q_abs)")
    println(io, "  g     = $(cs.asymmetry)")
end

# ---- Truncation criteria ----------------------------------------------------

"""
    mie_nmax(x::Real)

Wiscombe truncation criterion for Mie series.
Returns `ceil(Int, x + 4 * x^(1/3) + 2)`.
"""
mie_nmax(x::Real) = ceil(Int, x + 4 * cbrt(x) + 2)

"""
    mie_nmax_conservative(x::Real)

Conservative variant of the Wiscombe truncation criterion.
Returns `ceil(Int, x + 4.05 * x^(1/3) + 8)`.
"""
mie_nmax_conservative(x::Real) = ceil(Int, x + 4.05 * cbrt(x) + 8)

# ---- Lorenz-Mie coefficients ------------------------------------------------

"""
    mie_ab(N_max::Int, x::Real, m::Complex)

Compute Lorenz-Mie coefficients `a_n` and `b_n` for `n = 1:N_max`.

Returns `(a, b)` where both are `Vector{Complex{Float64}}` (or promoted type).

# Arguments
- `N_max::Int`: maximum order to compute
- `x::Real`: size parameter `k·r`
- `m::Complex`: relative refractive index

# Mathematical formulas (Sun et al. 2020, Eq. 3.3.41)

Riccati-Bessel functions:
- `ζ_n(x) = x · j_n(x)`  (first kind)
- `ξ_n(x) = x · h_n^(1)(x)`  (third kind / Hankel)

Derivatives:
- `ζ_n'(x) = j_n(x) + x · j_n'(x)`
- `ξ_n'(x) = h_n^(1)(x) + x · h_n^(1)'(x)`

Lorenz-Mie coefficients:
```
b_n = [ζ_n(mx)·ζ_n'(x) - m·ζ_n(x)·ζ_n'(mx)] / [ξ_n(mx)·ζ_n'(x) - m·ζ_n(x)·ξ_n'(mx)]
a_n = [m·ζ_n(mx)·ζ_n'(x) - ζ_n(x)·ζ_n'(mx)] / [m·ξ_n(mx)·ζ_n'(x) - ζ_n(x)·ξ_n'(mx)]
```
"""
function mie_ab(N_max::Int, x::Real, m::Complex)
    T = promote_type(typeof(float(x)), typeof(m))
    a = Vector{T}(undef, N_max)
    b = Vector{T}(undef, N_max)

    mx = m * x

    for n in 1:N_max
        # Riccati-Bessel functions of first kind (ζ)
        zeta_x  = x  * sphericalbesselj(n, x)
        zeta_mx = mx * sphericalbesselj(n, mx)

        # Riccati-Bessel functions of third kind (ξ) — inside particle only
        # Sun Eq. 3.3.40a: ξ_n(x) = x·h_n^(1)(x)
        xi_mx = mx * shankelh1(n, mx)

        # Derivatives of Riccati-Bessel functions
        zeta_prime_x  = sphericalbesselj(n, x)  + x  * sbesselj_deriv(n, x)
        zeta_prime_mx = sphericalbesselj(n, mx) + mx * sbesselj_deriv(n, mx)

        xi_prime_mx = shankelh1(n, mx) + mx * shankelh1_deriv(n, mx)

        # Sun Eq. 3.3.41 — b_n
        num_b = zeta_mx * zeta_prime_x - m * zeta_x * zeta_prime_mx
        den_b = xi_mx  * zeta_prime_x - m * zeta_x * xi_prime_mx
        b[n] = num_b / den_b

        # Sun Eq. 3.3.41 — a_n
        num_a = m * zeta_mx * zeta_prime_x - zeta_x * zeta_prime_mx
        den_a = m * xi_mx  * zeta_prime_x - zeta_x * xi_prime_mx
        a[n] = num_a / den_a

        if !isfinite(a[n]) || !isfinite(b[n])
            @warn "mie_ab: non-finite coefficient at n=$n (x=$x, m=$m); setting to NaN"
            a[n] = T(NaN)
            b[n] = T(NaN)
        end
    end

    return a, b
end

# ---- High-level solver ------------------------------------------------------

"""
    solve_tmatrix(sphere::Sphere, ::MieMethod, wavelength, refractive_index; N_max)

Compute the Lorenz-Mie T-matrix for a sphere.

# Arguments
- `sphere::Sphere`: the spherical particle
- `::MieMethod`: solver dispatch tag
- `wavelength::Real`: wavelength in the surrounding medium
- `refractive_index::Complex`: relative refractive index of the sphere

# Keyword arguments
- `N_max::Int = mie_nmax(2π / wavelength * sphere.radius)`: maximum multipole order

# Returns
[`MieTMatrix`](@ref) containing the diagonal T-matrix coefficients.

# Validation
- `sphere.radius > 0`
- `wavelength > 0`
- `imag(refractive_index) >= 0` (passive medium)
- `N_max >= 1`
"""
function solve_tmatrix(
    sphere::Sphere,
    ::MieMethod,
    wavelength::Real,
    refractive_index::Complex;
    N_max::Int = mie_nmax(2π / wavelength * sphere.radius)
)
    sphere.radius > 0 || throw(ArgumentError("sphere radius must be positive, got $(sphere.radius)"))
    wavelength > 0 || throw(ArgumentError("wavelength must be positive, got $wavelength"))
    imag(refractive_index) >= 0 || throw(ArgumentError("imag(refractive_index) must be >= 0 for a passive medium, got $(imag(refractive_index))"))
    N_max >= 1 || throw(ArgumentError("N_max must be >= 1, got $N_max"))

    k = 2π / wavelength
    x = k * sphere.radius
    a, b = mie_ab(N_max, x, refractive_index)
    return MieTMatrix(a, b, k, sphere.radius, N_max)
end

# ---- Post-processing: cross-sections ----------------------------------------

"""
    calc_cross_sections(T::MieTMatrix)

Compute scattering cross-sections and efficiencies from a [`MieTMatrix`](@ref).

Formulas follow Bohren & Huffman (1983), §4.5.
"""
function calc_cross_sections(T::MieTMatrix)
    a = T.a
    b = T.b
    k = T.k
    r = T.radius
    N_max = T.N_max

    C_ext = 0.0
    C_sca = 0.0

    for n in 1:N_max
        coeff = 2n + 1
        C_ext += coeff * real(a[n] + b[n])
        C_sca += coeff * (abs2(a[n]) + abs2(b[n]))
    end

    prefactor = 2π / k^2
    C_ext *= prefactor
    C_sca *= prefactor
    C_abs = C_ext - C_sca

    area = π * r^2
    Q_ext = C_ext / area
    Q_sca = C_sca / area
    Q_abs = C_abs / area

    # Asymmetry parameter g  (Bohren & Huffman, Eq. 4.67)
    g = 0.0
    for n in 1:(N_max - 1)
        g += (n * (n + 2) / (n + 1)) * real(a[n] * conj(a[n+1]) + b[n] * conj(b[n+1]))
    end
    for n in 1:N_max
        g += ((2n + 1) / (n * (n + 1))) * real(a[n] * conj(b[n]))
    end
    g *= (4π / k^2) / C_sca

    return CrossSections(C_ext, C_sca, C_abs, Q_ext, Q_sca, Q_abs, g)
end
