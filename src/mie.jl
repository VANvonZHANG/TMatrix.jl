# Layer 0.5: Lorenz-Mie Coefficients (Sun §3.3, Eq. 3.3.41)

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
