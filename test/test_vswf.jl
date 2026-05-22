@testset "VSWFs M, N" begin
    # M has no radial component; N has radial component (Sun Eq. 3.2.46)
    @test abs(TMatrix.vswf_M(2, 1, 3.0, 0.7, 0.4; regular = true)[1]) < 1e-14
    @test abs(TMatrix.vswf_N(2, 1, 3.0, 0.7, 0.4; regular = true)[1]) > 1e-15

    # Far-field M asymptotic: M → (-i)^{n+1} e^{iρ}/ρ · γ · C (Sun Eq. 3.2.46b)
    let n = 3, m = 1, theta = 0.5, phi = 0.3, rho = 200000.0
        vM = TMatrix.vswf_M(n, m, rho, theta, phi; regular = false)
        g = TMatrix.gamma_vsh(n, m)
        vC = TMatrix.C_vsh(n, m, theta, phi)
        asymp = (-im)^(n + 1) * exp(im * rho) / rho * g * vC
        @test vM ≈ asymp rtol=1e-4
    end

    # Far-field N tangential: N_tang → (-i)^n e^{iρ}/ρ · γ · B (Sun Eq. 3.2.46c)
    let n = 3, m = 1, theta = 0.5, phi = 0.3, rho = 200000.0
        vN = TMatrix.vswf_N(n, m, rho, theta, phi; regular = false)
        g = TMatrix.gamma_vsh(n, m)
        vB = TMatrix.B_vsh(n, m, theta, phi)
        asymp = (-im)^n * exp(im * rho) / rho * g * vB
        @test vN[2:end] ≈ asymp[2:end] rtol=1e-3
    end

    # Divergence-free M: numerical ∇·M ≈ 0
    let n = 2, m = 1, kr = 3.0, theta = 0.7, phi = 0.4, h = 1e-6
        vM_p = TMatrix.vswf_M(n, m, kr, theta + h, phi)
        vM_m = TMatrix.vswf_M(n, m, kr, theta - h, phi)
        vM_pp = TMatrix.vswf_M(n, m, kr, theta, phi + h)
        vM_pm = TMatrix.vswf_M(n, m, kr, theta, phi - h)
        div_theta = (sin(theta + h) * vM_p[2] - sin(theta - h) * vM_m[2]) /
                    (2h * sin(theta))
        div_phi = (vM_pp[3] - vM_pm[3]) / (2h * sin(theta))
        @test abs(div_theta + div_phi) < 1e-4
    end

    # L uses gamma_prime normalization (not gamma_vsh), has all three components
    let vL = TMatrix.vswf_L(2, 1, 3.0, 0.7, 0.4; regular = true)
        @test abs(vL[1]) > 1e-15  # radial component non-zero
        @test abs(vL[2]) > 1e-15  # tangential component non-zero
    end
end
