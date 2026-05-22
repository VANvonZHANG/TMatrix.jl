# Layer 1.1: Gaunt Coefficients
# G(l1,l2,l3; m1,m2,m3) = √[(2l1+1)(2l2+1)(2l3+1)/(4π)] × 3j(l1,l2,l3; 0,0,0) × 3j(l1,l2,l3; m1,m2,m3)

"""
    gaunt(l1, l2, l3, m1, m2, m3)

Gaunt coefficient ``G(l_1, l_2, l_3; m_1, m_2, m_3)``.

Integral of the triple product of spherical harmonics over the full sphere.
Decomposes into Wigner 3-j symbols:
```
G = sqrt[(2l1+1)(2l2+1)(2l3+1)/(4pi)] * 3j(l1,l2,l3; 0,0,0) * 3j(l1,l2,l3; m1,m2,m3)
```

Selection rules return zero early when violated:
- Triangular condition: `|l1-l2| <= l3 <= l1+l2`
- `m1 + m2 + m3 = 0`
- `l1 + l2 + l3` must be even

# Arguments
- `l1, l2, l3::Int`: Angular momenta
- `m1, m2, m3::Int`: Magnetic quantum numbers

# Returns
`Float64`: Gaunt coefficient value

# References
- Sun et al. 2020, mode coupling in surface integrals
"""
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

"""
    gamma_prime(n::Int, m::Int)

Normalization constant ``\\gamma'_{mn}`` for scalar spherical harmonics.

```
\\gamma'_{mn} = \\sqrt{\\frac{(2n+1)(n-m)!}{4\\pi(n+m)!}}
```

# Arguments
- `n::Int`: Degree
- `m::Int`: Order

# Returns
`Float64`: Normalization constant

# References
- Sun et al. 2020, Eq. 3.2.39a
"""
function gamma_prime(n::Int, m::Int)
    return sqrt((2n + 1) / (4pi) * exp(loggamma(n - m + 1) - loggamma(n + m + 1)))
end

"""
    gamma_vsh(n::Int, m::Int)

Normalization constant ``\\gamma_{mn}`` for vector spherical harmonics.

```
\\gamma_{mn} = \\sqrt{\\frac{(2n+1)(n-m)!}{4\\pi n(n+1)(n+m)!}}
```

Returns `Inf` for `n = 0` since the denominator vanishes.

# Arguments
- `n::Int`: Degree
- `m::Int`: Order

# Returns
`Float64`: Normalization constant

# References
- Sun et al. 2020, Eq. 3.2.39b
"""
function gamma_vsh(n::Int, m::Int)
    n == 0 && return oftype(1.0, Inf)
    return sqrt((2n + 1) / (4pi * n * (n + 1)) *
                exp(loggamma(n - m + 1) - loggamma(n + m + 1)))
end

# Scalar harmonic Y_{mn}(θ,φ) = P_n^m(cosθ) exp(imφ) (Sun Eq. 3.2.40)
Y_mn(n::Int, m::Int, theta::Real, phi::Real) = Plm(cos(theta), n, m) * exp(im * m * phi)

"""
    P_vsh(n::Int, m::Int, theta::Real, phi::Real)

Vector spherical harmonic ``P_{mn}(\\theta, \\phi)`` — radial component only.

```
P_{mn} = \\hat{r} \\cdot Y_{mn}(\\theta, \\phi)
```

Returns an `SVector{3, Complex}` with only the radial (first) component non-zero.

# Arguments
- `n::Int`: Degree
- `m::Int`: Order
- `theta::Real`: Polar angle
- `phi::Real`: Azimuthal angle

# Returns
`SVector{3, Complex}`: Radial vector harmonic

# References
- Sun et al. 2020, Eq. 3.2.41a
"""
function P_vsh(n::Int, m::Int, theta::Real, phi::Real)
    T = complex(typeof(float(theta)))
    val = Plm(cos(theta), n, m) * exp(im * m * phi)
    return SVector{3, T}(val, zero(T), zero(T))
end

"""
    B_vsh(n::Int, m::Int, theta::Real, phi::Real)

Vector spherical harmonic ``B_{mn}(\\theta, \\phi)`` — tangential gradient type.

```
B_{mn} = [\\hat{\\theta} \\frac{dP}{d\\theta} + \\hat{\\phi} \\frac{im}{\\sin\\theta} P_n^m] e^{im\\phi}
```

# Arguments
- `n::Int`: Degree
- `m::Int`: Order
- `theta::Real`: Polar angle
- `phi::Real`: Azimuthal angle

# Returns
`SVector{3, Complex}`: Tangential gradient vector harmonic

# References
- Sun et al. 2020, Eq. 3.2.41b
"""
function B_vsh(n::Int, m::Int, theta::Real, phi::Real)
    T = complex(typeof(float(theta)))
    eiphi = exp(im * m * phi)
    dpdt = dplmdtheta(n, m, theta)
    p_over_s = plm_over_sintheta(n, m, theta)
    return SVector{3, T}(zero(T), dpdt * eiphi, im * m * p_over_s * eiphi)
end

"""
    C_vsh(n::Int, m::Int, theta::Real, phi::Real)

Vector spherical harmonic ``C_{mn}(\\theta, \\phi)`` — curl type.

```
C_{mn} = [\\hat{\\theta} \\frac{im}{\\sin\\theta} P_n^m - \\hat{\\phi} \\frac{dP}{d\\theta}] e^{im\\phi}
```

Related to ``B_{mn}`` via ``C_{mn} = B_{mn} \\times \\hat{r}``.
(Some texts define ``C_{mn} = \\hat{r} \\times B_{mn}``, which differs by a sign.)

# Arguments
- `n::Int`: Degree
- `m::Int`: Order
- `theta::Real`: Polar angle
- `phi::Real`: Azimuthal angle

# Returns
`SVector{3, Complex}`: Curl-type vector harmonic

# References
- Sun et al. 2020, Eq. 3.2.41c
"""
function C_vsh(n::Int, m::Int, theta::Real, phi::Real)
    T = complex(typeof(float(theta)))
    eiphi = exp(im * m * phi)
    dpdt = dplmdtheta(n, m, theta)
    p_over_s = plm_over_sintheta(n, m, theta)
    return SVector{3, T}(zero(T), im * m * p_over_s * eiphi, -dpdt * eiphi)
end
