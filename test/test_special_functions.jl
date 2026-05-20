@testset "Spherical Bessel & Hankel" begin
    # Closed-form: j_0(x) = sin(x)/x, y_0(x) = -cos(x)/x (Sun Eq. 3.2.33)
    @test sphericalbesselj(0, 1.0) ≈ sin(1.0) / 1.0 rtol=1e-14
    @test sphericalbessely(0, 1.0) ≈ -cos(1.0) / 1.0 rtol=1e-14

    # h_0^(1)(x) = exp(ix)/(ix), h_0^(2)(x) = exp(-ix)/(-ix) (Sun Eq. 3.2.28c-d)
    @test TMatrix.shankelh1(0, 1.0) ≈ exp(im) / im rtol=1e-14
    @test TMatrix.shankelh2(0, 1.0) ≈ exp(-im) / (-im) rtol=1e-14

    # Wronskian j_n y_n' - j_n' y_n = 1/x² (Sun Eq. 3.2.32a)
    let n = 5, x = 3.0
        @test sphericalbesselj(n, x) * TMatrix.sbessely_deriv(n, x) -
              TMatrix.sbesselj_deriv(n, x) * sphericalbessely(n, x) ≈ 1 / x^2 rtol=1e-12
    end

    # Wronskian j_n h_n' - j_n' h_n = i/x² (Sun Eq. 3.2.32b)
    let n = 5, x = 3.0
        @test sphericalbesselj(n, x) * TMatrix.shankelh1_deriv(n, x) -
              TMatrix.sbesselj_deriv(n, x) * TMatrix.shankelh1(n, x) ≈ im / x^2 rtol=1e-12
    end

    # Recurrence: (2n+1)/x z_n = z_{n-1} + z_{n+1} (Sun Eq. 3.2.29)
    let n = 10, x = 2.0
        @test (2n + 1) / x * sphericalbesselj(n, x) ≈
              sphericalbesselj(n - 1, x) + sphericalbesselj(n + 1, x) rtol=1e-12
    end
end

@testset "Associated Legendre" begin
    # P_n(1) = 1, P_n(-1) = (-1)^n (Sun Eq. 3.2.15)
    @test Plm(1.0, 5, 0) ≈ 1.0 rtol=1e-14
    @test Plm(-1.0, 5, 0) ≈ -1.0 rtol=1e-14

    # P_2^2(cosθ) = 3 sin²θ (Sun Eq. 3.2.22c)
    @test Plm(cos(1.0), 2, 2) ≈ 3 * sin(1.0)^2 rtol=1e-14

    # Orthogonality: ∫ P_3^1 P_3^1 sinθ dθ = (3+1)!/(3-1)! · 2/7 (Sun Eq. 3.2.23)
    let x_gl = [-0.5, 0.0, 0.5], w_gl = [0.5, 0.6666667, 0.5] # rough check
        using FastGaussQuadrature: gausslegendre
        xg, wg = gausslegendre(60)
        @test sum(wg .* [Plm(xi, 3, 1) * Plm(xi, 3, 1) for xi in xg]) ≈
              float(factorial(4)) / factorial(2) * 2 / 7 rtol=1e-8
    end
end

@testset "Legendre θ-derivative" begin
    # dP/dθ vs finite difference
    let n = 3, m = 1, theta = 0.8, h = 1e-7
        fd = (Plm(cos(theta + h), n, m) - Plm(cos(theta - h), n, m)) / 2h
        @test TMatrix.dplmdtheta(n, m, theta) ≈ fd rtol=1e-6
    end

    # plm_over_sintheta matches direct division
    @test TMatrix.plm_over_sintheta(3, 1, 1.0) ≈ Plm(cos(1.0), 3, 1) / sin(1.0) rtol=1e-13

    # Limit at θ→0 for m=1: P_n^1/sinθ → -n(n+1)/2 (Condon-Shortley)
    @test TMatrix.plm_over_sintheta(3, 1, 1e-8) ≈ -3 * 4 / 2 rtol=1e-6
end

@testset "Wigner 3-j Symbols" begin
    @test wigner3j(Float64, 1, 1, 3, 0, 0, 0) ≈ 0.0 atol=1e-15  # triangular rule
    @test wigner3j(Float64, 1, 1, 1, 0, 0, 1) ≈ 0.0 atol=1e-15  # m-sum rule
    @test wigner3j(Float64, 1//2, 1//2, 1, -1//2, 1//2, 0) ≈ 1 / sqrt(6) rtol=1e-14
    @test wigner3j(Float64, 1, 1, 0, -1, 1, 0) ≈ 1 / sqrt(3) rtol=1e-14
end

@testset "Gauss-Legendre Quadrature" begin
    x, w = gausslegendre(50)
    @test sum(w) ≈ 2.0 rtol=1e-14
    @test sum(w .* x .^ 2) ≈ 2 / 3 rtol=1e-12  # exact for x^2

    # Interval mappings
    tn, tw = TMatrix.gl_theta(50)
    @test sum(tw) ≈ Float64(pi) rtol=1e-14
    pn, pw = TMatrix.gl_phi(50)
    @test sum(pw) ≈ Float64(2pi) rtol=1e-14
end
