@testset "Gaunt Coefficients" begin
    # G(0,0,0; 0,0,0) = 1/√(4π)
    @test TMatrix.gaunt(0, 0, 0, 0, 0, 0) ≈ 1 / sqrt(4pi) rtol=1e-14

    # Selection rules
    @test TMatrix.gaunt(1, 1, 3, 0, 0, 0) ≈ 0.0 atol=1e-15  # triangular
    @test TMatrix.gaunt(1, 1, 2, 0, 0, 1) ≈ 0.0 atol=1e-15  # m-sum
    @test TMatrix.gaunt(1, 2, 2, 0, 0, 0) ≈ 0.0 atol=1e-15  # parity

    # Closed form: G(l,l,0; m,-m,0) = (-1)^m / √(4π)
    @test TMatrix.gaunt(2, 2, 0, 1, -1, 0) ≈ -1 / sqrt(4pi) rtol=1e-12
    @test TMatrix.gaunt(3, 3, 0, 2, -2, 0) ≈ 1 / sqrt(4pi) rtol=1e-12
end

@testset "Normalization Constants" begin
    @test TMatrix.gamma_prime(0, 0) ≈ 1 / sqrt(4pi) rtol=1e-14
    @test TMatrix.gamma_prime(1, 0) ≈ sqrt(3 / (4pi)) rtol=1e-14
    @test TMatrix.gamma_vsh(1, 0) ≈ sqrt(3 / (8pi)) rtol=1e-14
    @test TMatrix.gamma_vsh(1, 1) ≈ sqrt(3 / (16pi)) rtol=1e-14
end

@testset "Vector Spherical Harmonics" begin
    # P_vsh: radial only
    let v = TMatrix.P_vsh(2, 1, 0.5, 0.3)
        @test v[1] ≈ Plm(cos(0.5), 2, 1) * exp(im * 0.3) rtol=1e-14
        @test abs(v[2]) < 1e-14 && abs(v[3]) < 1e-14
    end

    # B, C: no radial component; B = r̂ × C (Sun Eq. 3.2.41b-c)
    let vB = TMatrix.B_vsh(3, 1, 0.7, 0.4), vC = TMatrix.C_vsh(3, 1, 0.7, 0.4)
        @test abs(vB[1]) < 1e-14 && abs(vC[1]) < 1e-14
        @test vB[2] ≈ -vC[3] rtol=1e-12
        @test vB[3] ≈ vC[2] rtol=1e-12
    end

    # Parity: C_θ(π-θ, π+φ) = (-1)^n C_θ(θ,φ) (Sun Eq. 3.2.43c)
    let n = 3, m = 1, theta = 0.7, phi = 0.4
        v1 = TMatrix.C_vsh(n, m, theta, phi)
        v2 = TMatrix.C_vsh(n, m, pi - theta, pi + phi)
        abs(v1[2]) > 1e-10 && @test v2[2] ≈ (-1.0)^n * v1[2] rtol=1e-10
    end
end
