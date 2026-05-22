using FastGaussQuadrature: gausslegendre

@testset "Sphere Geometry" begin
    s = TMatrix.Sphere(1.0)
    @test TMatrix.surface_radius(s, 1.0) ≈ 1.0
    @test TMatrix.surface_dradius_dtheta(s, 1.0) ≈ 0.0 atol=1e-15

    # Surface area = 4πa² via quadrature
    xg, wg = gausslegendre(60)
    area = sum(wg .* [(pi / 2) * TMatrix.surface_element(s, (pi / 2) * (xi + 1)) * 2pi
                for xi in xg])
    @test area ≈ 4pi rtol=1e-10

    @test TMatrix.surface_normal(s, 1.0)[1] > 0
end

@testset "Spheroid Geometry" begin
    sp = TMatrix.Spheroid(1.0, 2.0)
    @test TMatrix.surface_radius(sp, 1e-10) ≈ 2.0 rtol=1e-6   # pole: r = b
    @test TMatrix.surface_radius(sp, pi / 2) ≈ 1.0 rtol=1e-14  # equator: r = a

    # a=b sphere case: area = 4πa²
    sp_s = TMatrix.Spheroid(1.0, 1.0)
    xg, wg = gausslegendre(60)
    area = sum(wg .* [(pi / 2) * TMatrix.surface_element(sp_s, (pi / 2) * (xi + 1)) * 2pi
                for xi in xg])
    @test area ≈ 4pi rtol=1e-8
end

@testset "Chebyshev Particle" begin
    cp = TMatrix.ChebyshevParticle(1.0, 0.0, 2)
    @test TMatrix.surface_radius(cp, 1.0) ≈ 1.0 rtol=1e-14  # ε=0: sphere

    cp2 = TMatrix.ChebyshevParticle(1.0, 0.1, 2)
    @test TMatrix.surface_radius(cp2, 0.0) ≈ 1.1 rtol=1e-14  # r(0) = a(1+ε)
end
