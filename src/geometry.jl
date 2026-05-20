# Layer 2.2: Geometry & Surface Parameterization

using StaticArrays: SVector

abstract type AbstractParticle end

# ===== Sphere =====

struct Sphere <: AbstractParticle
    radius::Float64
end

surface_radius(s::Sphere, theta::Real) = s.radius
surface_dradius_dtheta(s::Sphere, theta::Real) = 0.0

function surface_element(s::Sphere, theta::Real)
    return s.radius^2 * sin(theta)
end

function surface_normal(s::Sphere, theta::Real)
    return SVector{3, Float64}(1.0, 0.0, 0.0)
end

# ===== Spheroid =====
# Ellipsoid: (x²+y²)/a² + z²/b² = 1  →  r(θ) = ab / √(b²sin²θ + a²cos²θ)
# a=equatorial (xy-plane), b=polar (z-axis)

struct Spheroid <: AbstractParticle
    a::Float64  # equatorial radius
    b::Float64  # polar radius
end

function surface_radius(sp::Spheroid, theta::Real)
    a, b = sp.a, sp.b
    sth, cth = sin(theta), cos(theta)
    return a * b / sqrt(b^2 * sth^2 + a^2 * cth^2)
end

function surface_dradius_dtheta(sp::Spheroid, theta::Real)
    a, b = sp.a, sp.b
    sth, cth = sin(theta), cos(theta)
    denom = b^2 * sth^2 + a^2 * cth^2
    return a * b * (a^2 - b^2) * sth * cth / denom^1.5
end

function surface_element(sp::Spheroid, theta::Real)
    r = surface_radius(sp, theta)
    dr = surface_dradius_dtheta(sp, theta)
    return r * sin(theta) * sqrt(r^2 + dr^2)
end

function surface_normal(sp::Spheroid, theta::Real)
    r = surface_radius(sp, theta)
    dr = surface_dradius_dtheta(sp, theta)
    n_r = r
    n_theta = -dr
    norm_val = sqrt(n_r^2 + n_theta^2)
    return SVector{3, Float64}(n_r / norm_val, n_theta / norm_val, 0.0)
end

# ===== Chebyshev Particle =====
# r(θ) = a[1 + ε T_n(cosθ)]

struct ChebyshevParticle <: AbstractParticle
    a::Float64
    epsilon::Float64
    n_cheb::Int
end

function _chebyshev_T(n::Int, x::Real)
    n == 0 && return one(x)
    n == 1 && return x
    T_prev2, T_prev1 = one(x), x
    for _ in 2:n
        T_prev2, T_prev1 = T_prev1, 2x * T_prev1 - T_prev2
    end
    return T_prev1
end

function _chebyshev_T_deriv(n::Int, x::Real)
    n == 0 && return zero(x)
    n == 1 && return one(x)
    # dT_n/dx = n U_{n-1}(x), U = Chebyshev of second kind
    U_prev2, U_prev1 = one(x), 2x
    for _ in 2:(n-1)
        U_prev2, U_prev1 = U_prev1, 2x * U_prev1 - U_prev2
    end
    return n * U_prev1
end

function surface_radius(cp::ChebyshevParticle, theta::Real)
    return cp.a * (1 + cp.epsilon * _chebyshev_T(cp.n_cheb, cos(theta)))
end

function surface_dradius_dtheta(cp::ChebyshevParticle, theta::Real)
    return -cp.a * cp.epsilon * _chebyshev_T_deriv(cp.n_cheb, cos(theta)) * sin(theta)
end

function surface_element(cp::ChebyshevParticle, theta::Real)
    r = surface_radius(cp, theta)
    dr = surface_dradius_dtheta(cp, theta)
    return r * sin(theta) * sqrt(r^2 + dr^2)
end

function surface_normal(cp::ChebyshevParticle, theta::Real)
    r = surface_radius(cp, theta)
    dr = surface_dradius_dtheta(cp, theta)
    n_r, n_theta = r, -dr
    norm_val = sqrt(n_r^2 + n_theta^2)
    return SVector{3, Float64}(n_r / norm_val, n_theta / norm_val, 0.0)
end
