# Layer 0.1: Spherical Bessel & Hankel Functions (Sun §3.2, Eqs. 3.2.27–3.2.33)

# Sun Eq. 3.2.28c: h_n^(1)(x) = j_n(x) + i y_n(x) — outgoing wave
shankelh1(n, x) = sphericalbesselj(n, x) + im * sphericalbessely(n, x)

# Sun Eq. 3.2.28d: h_n^(2)(x) = j_n(x) - i y_n(x) — incoming wave
shankelh2(n, x) = sphericalbesselj(n, x) - im * sphericalbessely(n, x)

# Derivative via recurrence: z_n'(x) = z_{n-1}(x) - (n+1)/x z_n(x)
# For n=0: z_0'(x) = -z_1(x) (avoids negative-order Bessel)
function sbesselj_deriv(n::Int, x)
    n == 0 && return -sphericalbesselj(1, x)
    return sphericalbesselj(n - 1, x) - (n + 1) / x * sphericalbesselj(n, x)
end

function sbessely_deriv(n::Int, x)
    n == 0 && return -sphericalbessely(1, x)
    return sphericalbessely(n - 1, x) - (n + 1) / x * sphericalbessely(n, x)
end

function shankelh1_deriv(n::Int, x)
    n == 0 && return -shankelh1(1, x)
    return shankelh1(n - 1, x) - (n + 1) / x * shankelh1(n, x)
end

function shankelh2_deriv(n::Int, x)
    n == 0 && return -shankelh2(1, x)
    return shankelh2(n - 1, x) - (n + 1) / x * shankelh2(n, x)
end

# Riccati-Bessel derivative: [ρ z_n(ρ)]'/ρ = z_{n-1}(ρ) - n/ρ z_n(ρ)
# Used in N_{mn} construction (Sun Eq. 3.2.46)
riccati_zn_prime(zn_prev, zn, n, rho) = zn_prev - n / rho * zn

# Layer 0.2: Associated Legendre & θ-derivative (Sun §3.2, Eqs. 3.2.18–3.2.26)

# dP_n^m(cosθ)/dθ = n cosθ/sinθ · P_n^m(cosθ) - (n+m)/sinθ · P_{n-1}^m(cosθ)
# (Sun Eq. 3.2.24)
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
function gl_theta(N::Int)
    x, w = gausslegendre(N)
    theta = (pi / 2) .* (x .+ 1)
    w_theta = (pi / 2) .* w
    return theta, w_theta
end

# Map Gauss-Legendre nodes from [-1,1] to [0,2π] for φ integration
function gl_phi(N::Int)
    x, w = gausslegendre(N)
    phi = pi .* (x .+ 1)
    w_phi = pi .* w
    return phi, w_phi
end
