using TMatrix
using Test
using SpecialFunctions: sphericalbesselj, sphericalbessely
using LegendrePolynomials: Plm, Pl
using StaticArrays
using WignerSymbols: wigner3j
using FastGaussQuadrature: gausslegendre

@testset "TMatrix.jl Phase 1" begin
    include("test_special_functions.jl")
    include("test_composed.jl")
    include("test_vswf.jl")
    include("test_geometry.jl")
end
