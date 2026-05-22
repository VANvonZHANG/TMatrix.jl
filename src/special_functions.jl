# Layer 0.1: Spherical Bessel & Hankel Functions (Sun §3.2, Eqs. 3.2.27–3.2.33)

# Sun Eq. 3.2.28c: h_n^(1)(x) = j_n(x) + i y_n(x) — outgoing wave
"""
    shankelh1(n, x)

Spherical Hankel function of the first kind, ``h_n^{(1)}(x) = j_n(x) + i y_n(x)``.

Represents an outgoing spherical wave. Used for scattered fields.

# Arguments
- `n::Int`: Order
- `x`: Argument (real or complex)

# Returns
`Complex`: Outgoing spherical wave function value

# References
- Sun et al. 2020, Eq. 3.2.28c
"""
shankelh1(n, x) = sphericalbesselj(n, x) + im * sphericalbessely(n, x)

# Sun Eq. 3.2.28d: h_n^(2)(x) = j_n(x) - i y_n(x) — incoming wave
"""
    shankelh2(n, x)

Spherical Hankel function of the second kind, ``h_n^{(2)}(x) = j_n(x) - i y_n(x)``.

Represents an incoming spherical wave.

# Arguments
- `n::Int`: Order
- `x`: Argument (real or complex)

# Returns
`Complex`: Incoming spherical wave function value

# References
- Sun et al. 2020, Eq. 3.2.28d
"""
shankelh2(n, x) = sphericalbesselj(n, x) - im * sphericalbessely(n, x)

# Derivative via recurrence: z_n'(x) = z_{n-1}(x) - (n+1)/x z_n(x)
# For n=0: z_0'(x) = -z_1(x) (avoids negative-order Bessel)
"""
    sbesselj_deriv(n::Int, x)

Derivative of the spherical Bessel function of the first kind, ``j_n'(x)``.

Computed via recurrence relation: ``z_n'(x) = z_{n-1}(x) - (n+1)/x \\cdot z_n(x)``.
For ``n=0``, uses ``j_0'(x) = -j_1(x)``.

# Arguments
- `n::Int`: Order
- `x`: Argument

# Returns
Derivative value

# References
- Sun et al. 2020, Eq. 3.2.29
"""
function sbesselj_deriv(n::Int, x)
    n == 0 && return -sphericalbesselj(1, x)
    return sphericalbesselj(n - 1, x) - (n + 1) / x * sphericalbesselj(n, x)
end

"""
    sbessely_deriv(n::Int, x)

Derivative of the spherical Bessel function of the second kind, ``y_n'(x)``.

Computed via recurrence relation: ``z_n'(x) = z_{n-1}(x) - (n+1)/x \\cdot z_n(x)``.
For ``n=0``, uses ``y_0'(x) = -y_1(x)``.

# Arguments
- `n::Int`: Order
- `x`: Argument

# Returns
Derivative value

# References
- Sun et al. 2020, Eq. 3.2.29
"""
function sbessely_deriv(n::Int, x)
    n == 0 && return -sphericalbessely(1, x)
    return sphericalbessely(n - 1, x) - (n + 1) / x * sphericalbessely(n, x)
end

"""
    shankelh1_deriv(n::Int, x)

Derivative of the spherical Hankel function of the first kind, ``{h_n^{(1)}}'(x)``.

Computed via recurrence relation: ``z_n'(x) = z_{n-1}(x) - (n+1)/x \\cdot z_n(x)``.
For ``n=0``, uses ``{h_0^{(1)}}'(x) = -h_1^{(1)}(x)``.

# Arguments
- `n::Int`: Order
- `x`: Argument

# Returns
Derivative value

# References
- Sun et al. 2020, Eq. 3.2.29
"""
function shankelh1_deriv(n::Int, x)
    n == 0 && return -shankelh1(1, x)
    return shankelh1(n - 1, x) - (n + 1) / x * shankelh1(n, x)
end

# Riccati-Bessel derivative: [ρ z_n(ρ)]'/ρ = z_{n-1}(ρ) - n/ρ z_n(ρ)
# Used in N_{mn} construction (Sun Eq. 3.2.46)
riccati_zn_prime(zn_prev, zn, n, rho) = zn_prev - n / rho * zn

# Layer 0.2: Associated Legendre & θ-derivative (Sun §3.2, Eqs. 3.2.18–3.2.26)

# dP_n^m(cosθ)/dθ = n cosθ/sinθ · P_n^m(cosθ) - (n+m)/sinθ · P_{n-1}^m(cosθ)
# (Sun Eq. 3.2.24)
"""
    dplmdtheta(n::Int, m::Int, theta::Real)

Derivative of the associated Legendre polynomial with respect to ``\\theta``:
``dP_n^m(\\cos\\theta)/d\\theta``.

Computed via recurrence relation (Sun Eq. 3.2.24):
```
dP/dθ = n cosθ/sinθ · P_n^m(cosθ) - (n+m)/sinθ · P_{n-1}^m(cosθ)
```

Handles the removable singularity at ``\\theta = 0`` and ``\\theta = \\pi`` via
a finite-difference approximation.

# Arguments
- `n::Int`: Degree
- `m::Int`: Order
- `theta::Real`: Polar angle in radians

# Returns
`Real`: Derivative value

# References
- Sun et al. 2020, Eq. 3.2.24
"""
function dplmdtheta(n::Int, m::Int, theta::Real)
    abs(m) > n && return zero(float(theta))
    x = cos(theta)
    sth = sin(theta)
    if abs(sth) < 1e-10
        return _dplmdtheta_limit(n, m, theta)
    end
    pn = Plm(x, n, m)
    pn_prev = (n > 0 && m <= n - 1) ? Plm(x, n - 1, m) : zero(x)
    return n * x / sth * pn - (n + m) / sth * pn_prev
end

# P_n^m(cosθ)/sinθ with removable singularity handling
"""
    plm_over_sintheta(n::Int, m::Int, theta::Real)

Associated Legendre polynomial divided by ``\\sin\\theta``:
``P_n^m(\\cos\\theta) / \\sin\\theta``.

This ratio appears in the vector spherical harmonics ``B_{mn}`` and ``C_{mn}``.
For ``m = 0``, the function is genuinely singular and returns `Inf`.
For ``m \\geq 2``, the limit is zero.
For ``m = 1``, the limit is ``-n(n+1)/2`` at ``\\theta = 0`` and ``(-1)^n n(n+1)/2`` at ``\\theta = \\pi``.

# Arguments
- `n::Int`: Degree
- `m::Int`: Order
- `theta::Real`: Polar angle in radians

# Returns
`Real`: The ratio value, or limiting value near singularities

# References
- Sun et al. 2020, Eq. 3.2.41
"""
function plm_over_sintheta(n::Int, m::Int, theta::Real)
    abs(m) > n && return zero(float(theta))
    sth = sin(theta)
    # Threshold must be large enough that cos(theta) hasn't rounded to ±1
    if abs(sth) < 1e-6
        return _plm_over_sintheta_limit(n, m, theta)
    end
    return Plm(cos(theta), n, m) / sth
end

# Limiting value of P_n^m(cosθ)/sinθ as θ → 0 or π
function _plm_over_sintheta_limit(n::Int, m::Int, theta::Real)
    m == 0 && return oftype(float(theta), Inf)  # genuinely singular for m=0
    m >= 2 && return zero(float(theta))
    # m == 1: P_n^1(x)/sinθ → -n(n+1)/2 as θ→0 (Condon-Shortley phase)
    x = cos(theta)
    sign = x > 0 ? -1.0 : (-1.0)^(n + 2)
    return sign * n * (n + 1) / 2
end

# Limiting value of dP_n^m(cosθ)/dθ at θ → 0 or π
function _dplmdtheta_limit(n::Int, m::Int, theta::Real)
    h = oftype(float(theta), 1e-7)
    p_plus = Plm(cos(theta + h), n, m)
    p_minus = Plm(cos(theta - h), n, m)
    return (p_plus - p_minus) / (2h)
end

# Layer 0.3: Wigner 3-j Symbols — no thin compositions needed, direct use of WignerSymbols.jl

# Layer 0.4: Gauss-Legendre Quadrature with interval mapping

# Map Gauss-Legendre nodes from [-1,1] to [0,π] for θ integration
"""
    gl_theta(N::Int)

Gauss-Legendre quadrature nodes and weights mapped to ``[0, \\pi]`` for ``\\theta`` integration.

Transforms the standard `[-1, 1]` interval via ``\\theta = \\pi/2 \\cdot (x + 1)``.

# Arguments
- `N::Int`: Number of quadrature points

# Returns
`(theta, w_theta)` where `theta` and `w_theta` are vectors of length `N`

# References
- Sun et al. 2020, used in surface integrals for EBCM
"""
function gl_theta(N::Int)
    x, w = gausslegendre(N)
    theta = (pi / 2) .* (x .+ 1)
    w_theta = (pi / 2) .* w
    return theta, w_theta
end

# Map Gauss-Legendre nodes from [-1,1] to [0,2π] for φ integration
"""
    gl_phi(N::Int)

Gauss-Legendre quadrature nodes and weights mapped to ``[0, 2\\pi]`` for ``\\phi`` integration.

Transforms the standard `[-1, 1]` interval via ``\\phi = \\pi \\cdot (x + 1)``.

# Arguments
- `N::Int`: Number of quadrature points

# Returns
`(phi, w_phi)` where `phi` and `w_phi` are vectors of length `N`

# References
- Sun et al. 2020, used in surface integrals for EBCM
"""
function gl_phi(N::Int)
    x, w = gausslegendre(N)
    phi = pi .* (x .+ 1)
    w_phi = pi .* w
    return phi, w_phi
end
