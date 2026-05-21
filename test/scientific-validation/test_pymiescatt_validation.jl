using TMatrix
using Test
using JSON

# Helper for comparing near-zero values
function approx_zero_or_rel(actual, expected; rtol=1e-6, atol=1e-10)
    if abs(expected) < atol
        return abs(actual) < atol
    else
        return isapprox(actual, expected; rtol=rtol)
    end
end

@testset "Mie vs PyMieScatt" begin
    data = JSON.parsefile(joinpath(@__DIR__, "reference_data.json"))

    for case in data
        m = case["m"]["real"] + case["m"]["imag"] * im
        λ = case["wavelength_nm"] * 1e-9  # nm → m
        r = case["diameter_nm"] * 1e-9 / 2  # nm → m, diameter → radius

        T = solve_tmatrix(Sphere(r), MieMethod(), λ, m)
        C = calc_cross_sections(T)

        # Compare efficiencies
        @test C.Q_ext ≈ case["Qext"] rtol=1e-6
        @test C.Q_sca ≈ case["Qsca"] rtol=1e-6
        @test approx_zero_or_rel(C.Q_abs, case["Qabs"]; rtol=1e-6, atol=1e-10)
        @test C.asymmetry ≈ case["g"] rtol=1e-6

        # Compare cross-sections (in nm²)
        @test C.extinction * 1e18 ≈ case["Cext"] rtol=1e-6
        @test C.scattering * 1e18 ≈ case["Csca"] rtol=1e-6
        @test approx_zero_or_rel(C.absorption * 1e18, case["Cabs"]; rtol=1e-6, atol=1e-6)
    end
end
