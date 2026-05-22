# Layer 1.1: Gaunt Coefficients
# G(l1,l2,l3; m1,m2,m3) = √[(2l1+1)(2l2+1)(2l3+1)/(4π)] × 3j(l1,l2,l3; 0,0,0) × 3j(l1,l2,l3; m1,m2,m3)

function gaunt(l1::Int, l2::Int, l3::Int, m1::Int, m2::Int, m3::Int)
    # Selection rules for early return
    abs(l1 - l2) > l3 && return 0.0
    l3 > l1 + l2 && return 0.0
    m1 + m2 + m3 != 0 && return 0.0
    isodd(l1 + l2 + l3) && return 0.0

    prefactor = sqrt((2l1 + 1) * (2l2 + 1) * (2l3 + 1) / (4pi))
    w3j_zero = wigner3j(Float64, l1, l2, l3, 0, 0, 0)
    w3j_m = wigner3j(Float64, l1, l2, l3, m1, m2, m3)
    return prefactor * w3j_zero * w3j_m
end

# Layer 1.2: Vector Spherical Harmonics P, B, C (Sun §3.2, Eqs. 3.2.39–3.2.45)

# Normalization constant γ'_{mn} = √[(2n+1)(n-m)! / (4π(n+m)!)]
function gamma_prime(n::Int, m::Int)
    return sqrt((2n + 1) / (4pi) * exp(loggamma(n - m + 1) - loggamma(n + m + 1)))
end

# Normalization constant γ_{mn} = √[(2n+1)(n-m)! / (4πn(n+1)(n+m)!)]
function gamma_vsh(n::Int, m::Int)
    n == 0 && return oftype(1.0, Inf)
    return sqrt((2n + 1) / (4pi * n * (n + 1)) *
                exp(loggamma(n - m + 1) - loggamma(n + m + 1)))
end

# Scalar harmonic Y_{mn}(θ,φ) = P_n^m(cosθ) exp(imφ) (Sun Eq. 3.2.40)
Y_mn(n::Int, m::Int, theta::Real, phi::Real) = Plm(cos(theta), n, m) * exp(im * m * phi)

# Vector harmonic P_{mn}: radial component only (Sun Eq. 3.2.41a)
function P_vsh(n::Int, m::Int, theta::Real, phi::Real)
    T = complex(typeof(float(theta)))
    val = Plm(cos(theta), n, m) * exp(im * m * phi)
    return SVector{3, T}(val, zero(T), zero(T))
end

# Vector harmonic B_{mn}: tangential gradient type (Sun Eq. 3.2.41b)
# B_{mn} = [θ̂ dP/dθ + φ̂ im/sinθ P_n^m] exp(imφ)
function B_vsh(n::Int, m::Int, theta::Real, phi::Real)
    T = complex(typeof(float(theta)))
    eiphi = exp(im * m * phi)
    dpdt = dplmdtheta(n, m, theta)
    p_over_s = plm_over_sintheta(n, m, theta)
    return SVector{3, T}(zero(T), dpdt * eiphi, im * m * p_over_s * eiphi)
end

# Vector harmonic C_{mn}: curl type (Sun Eq. 3.2.41c)
# C_{mn} = [θ̂ im/sinθ P_n^m - φ̂ dP/dθ] exp(imφ)
function C_vsh(n::Int, m::Int, theta::Real, phi::Real)
    T = complex(typeof(float(theta)))
    eiphi = exp(im * m * phi)
    dpdt = dplmdtheta(n, m, theta)
    p_over_s = plm_over_sintheta(n, m, theta)
    return SVector{3, T}(zero(T), im * m * p_over_s * eiphi, -dpdt * eiphi)
end
