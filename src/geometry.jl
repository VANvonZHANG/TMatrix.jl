# Layer 2.2: Geometry & Surface Parameterization

using StaticArrays: SVector

"""
    AbstractParticle

Abstract supertype for all particle geometries used in T-matrix computations.

Concrete subtypes: [`Sphere`](@ref), [`Spheroid`](@ref), [`ChebyshevParticle`](@ref).

Each subtype must implement:
- `surface_radius(p, θ)` — surface-to-origin distance
- `surface_dradius_dtheta(p, θ)` — derivative for normal computation
- `surface_element(p, θ)` — Jacobian × sinθ for integration
- `surface_normal(p, θ)` — outward unit normal vector
"""
abstract type AbstractParticle end

# ===== Sphere =====

"""
    Sphere(radius::Float64)

Spherical particle geometry.

# Fields
- `radius::Float64`: Sphere radius

# Example
```julia
s = Sphere(1.0)
```
"""
struct Sphere <: AbstractParticle
    radius::Float64
end

"""
    surface_radius(p::AbstractParticle, theta::Real)

Surface-to-origin distance ``r(\\theta)`` for a particle.

# Arguments
- `p::AbstractParticle`: Particle geometry
- `theta::Real`: Polar angle

# Returns
`Float64`: Radial distance from origin to surface
"""
surface_radius(s::Sphere, theta::Real) = s.radius
"""
    surface_dradius_dtheta(p::AbstractParticle, theta::Real)

Derivative of the surface radius with respect to ``\\theta``:
``dr(\\theta)/d\\theta``.

Needed for surface normal computation.

# Arguments
- `p::AbstractParticle`: Particle geometry
- `theta::Real`: Polar angle

# Returns
`Float64`: Derivative value
"""
surface_dradius_dtheta(s::Sphere, theta::Real) = 0.0

"""
    surface_element(p::AbstractParticle, theta::Real)

Surface element ``dS = r(\\theta) \\sin\\theta \\sqrt{r^2 + (dr/d\\theta)^2} \\, d\\theta \\, d\\phi``
for axisymmetric particles.

Returns the ``\\theta``-dependent part (the integrand before ``d\\theta \\, d\\phi``).

# Arguments
- `p::AbstractParticle`: Particle geometry
- `theta::Real`: Polar angle

# Returns
`Float64`: Surface element value
"""
function surface_element(s::Sphere, theta::Real)
    return s.radius^2 * sin(theta)
end

"""
    surface_normal(p::AbstractParticle, theta::Real)

Outward-pointing unit surface normal vector in spherical coordinates.

For axisymmetric particles, the normal lies in the ``(\\hat{r}, \\hat{\\theta})`` plane
(no ``\\hat{\\phi}`` component).

# Arguments
- `p::AbstractParticle`: Particle geometry
- `theta::Real`: Polar angle

# Returns
`SVector{3, Float64}`: Unit normal vector `(n_r, n_θ, 0.0)`
"""
function surface_normal(s::Sphere, theta::Real)
    return SVector{3, Float64}(1.0, 0.0, 0.0)
end

# ===== Spheroid =====
# Ellipsoid: (x²+y²)/a² + z²/b² = 1  →  r(θ) = ab / √(b²sin²θ + a²cos²θ)
# a=equatorial (xy-plane), b=polar (z-axis)

"""
    Spheroid(a::Float64, b::Float64)

Spheroidal (ellipsoid of revolution) particle geometry.

Surface equation: ``(x^2+y^2)/a^2 + z^2/b^2 = 1``

# Fields
- `a::Float64`: Equatorial radius (xy-plane)
- `b::Float64`: Polar radius (z-axis)

`a > b` gives an oblate spheroid; `a < b` gives a prolate spheroid.

# Example
```julia
sp = Spheroid(1.0, 2.0)  # prolate spheroid
```
"""
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

"""
    ChebyshevParticle(a::Float64, epsilon::Float64, n_cheb::Int)

Chebyshev-deformed sphere particle geometry.

Surface equation: ``r(\\theta) = a [1 + \\varepsilon T_n(\\cos\\theta)]``

# Fields
- `a::Float64`: Mean radius
- `epsilon::Float64`: Deformation amplitude
- `n_cheb::Int`: Chebyshev polynomial degree

# Example
```julia
cp = ChebyshevParticle(1.0, 0.1, 4)
```
"""
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
    for _ in 2:(n - 1)
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
