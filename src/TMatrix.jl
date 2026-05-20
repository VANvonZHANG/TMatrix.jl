module TMatrix

using SpecialFunctions: sphericalbesselj, sphericalbessely, loggamma
using LegendrePolynomials: Plm, Pl
using WignerSymbols: wigner3j
using FastGaussQuadrature: gausslegendre
using StaticArrays: SVector
using LinearAlgebra: dot, cross

# Layer 0: Direct dependencies & thin compositions
include("special_functions.jl")

# Layer 1: Composed functions
include("composed.jl")

# Layer 2: VSWFs & Geometry
include("vswf.jl")
include("geometry.jl")

# Layer 0 exports
export shankelh1, shankelh2, sbesselj_deriv, sbessely_deriv, shankelh1_deriv
export dplmdtheta, plm_over_sintheta, gl_theta, gl_phi

# Layer 1 exports
export gaunt, gamma_prime, gamma_vsh
export P_vsh, B_vsh, C_vsh

# Layer 2 exports
export vswf_M, vswf_N, vswf_L
export Sphere, Spheroid, ChebyshevParticle
export surface_radius, surface_dradius_dtheta, surface_element, surface_normal

end
