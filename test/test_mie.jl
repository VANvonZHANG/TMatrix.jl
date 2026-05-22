using LinearAlgebra

@testset "Mie Coefficients - Analytical Limits" begin
    # No scattering: m = 1.0 → all a_n = b_n = 0
    a, b = TMatrix.mie_ab(10, 5.0, 1.0 + 0.0im)
    @test all(abs.(a) .< 1e-14)
    @test all(abs.(b) .< 1e-14)

    # Rayleigh limit: x ≪ 1, a_1 ≈ -i·(2/3)·(m²-1)/(m²+2)·x³
    m = 1.5 + 0.0im
    x = 0.01
    a, b = TMatrix.mie_ab(5, x, m)
    a1_expected = -im * (2/3) * (m^2 - 1) / (m^2 + 2) * x^3
    @test a[1] ≈ a1_expected rtol=1e-3
    @test abs(b[1]) < abs(a1_expected) * 1e-2  # b_1 is O(x⁵)

    # Energy conservation for real m: Q_abs ≈ 0, Q_ext ≈ Q_sca
    T = solve_tmatrix(Sphere(1.0), MieMethod(), 0.5, 1.5 + 0.0im; N_max = 20)
    C = calc_cross_sections(T)
    @test C.Q_abs ≈ 0.0 atol=1e-12
    @test C.Q_ext ≈ C.Q_sca rtol=1e-12
end

@testset "Mie Coefficients - Decay Verification" begin
    # For n > x, coefficients should decay rapidly
    x = 5.0
    m = 1.5 + 0.01im
    N_max = TMatrix.mie_nmax(x)
    a, b = TMatrix.mie_ab(N_max, x, m)
    @test abs(a[N_max]) < 1e-10
    @test abs(b[N_max]) < 1e-10
end

@testset "solve_tmatrix Interface" begin
    # Valid inputs
    T = solve_tmatrix(Sphere(0.5), MieMethod(), 0.532, 1.5 + 0.01im)
    @test T isa TMatrix.MieTMatrix
    @test T.N_max == TMatrix.mie_nmax(2π / 0.532 * 0.5)
    @test length(T.a) == T.N_max
    @test length(T.b) == T.N_max

    # Invalid inputs
    @test_throws ArgumentError solve_tmatrix(Sphere(-1.0), MieMethod(), 0.532, 1.5 + 0.01im)
    @test_throws ArgumentError solve_tmatrix(Sphere(0.5), MieMethod(), -0.532, 1.5 + 0.01im)
    @test_throws ArgumentError solve_tmatrix(Sphere(0.5), MieMethod(), 0.532, 1.5 - 0.01im)
end

@testset "Cross Sections - Properties" begin
    m = 1.5 + 0.1im
    x = 2.0
    r = 0.5
    λ = 2π * r / x
    T = solve_tmatrix(Sphere(r), MieMethod(), λ, m)
    C = calc_cross_sections(T)

    # Non-negative cross-sections
    @test C.extinction >= 0
    @test C.scattering >= 0
    @test C.absorption >= 0

    # Energy conservation: Q_ext = Q_sca + Q_abs
    @test C.Q_ext ≈ C.Q_sca + C.Q_abs rtol=1e-12

    # Efficiencies should be non-negative
    @test C.Q_ext >= 0
    @test C.Q_sca >= 0
    @test C.Q_abs >= 0

    # g should be in [-1, 1]
    @test -1.0 <= C.asymmetry <= 1.0
end

@testset "Diagonal Conversion" begin
    T = solve_tmatrix(Sphere(0.5), MieMethod(), 0.532, 1.5 + 0.01im; N_max = 5)
    D = Diagonal(T)
    @test D isa LinearAlgebra.Diagonal
    N_total = 2 * 5 * (5 + 2)  # 2 * N_max * (N_max + 2)
    @test size(D) == (N_total, N_total)
end
