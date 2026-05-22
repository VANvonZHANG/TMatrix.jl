# Mie Scattering Tutorial

This tutorial demonstrates a parametric scan over the size parameter and
visualizes how Mie cross-sections vary with particle size.

## Size Parameter Scan

The size parameter ``x = ka = 2\pi a / \lambda`` controls the scattering regime.
Small ``x`` (``x \ll 1``) is the Rayleigh regime; large ``x`` is the geometric-optics regime.

### Setup

```julia
using TMatrix

m = 1.5 + 0.01im       # refractive index
k = 1.0                # set k = 1, scan radius instead
radii = range(0.01, 10.0, length=200)

C_ext = Float64[]
C_sca = Float64[]
C_abs = Float64[]

for a in radii
    particle = Sphere(a)
    x = k * a
    N_max = mie_nmax(x)
    t_matrix = solve_tmatrix(particle, MieMethod(), k, m, N_max)
    cs = calc_cross_sections(t_matrix)
    push!(C_ext, cs.extinction)
    push!(C_sca, cs.scattering)
    push!(C_abs, cs.absorption)
end
```

### Plotting

With your favorite plotting library (e.g., `Plots.jl`):

```julia
using Plots

x_vals = k .* radii
plot(x_vals, [C_ext C_sca C_abs],
     label=["C_ext" "C_sca" "C_abs"],
     xlabel="Size parameter x = ka",
     ylabel="Cross-section",
     title="Mie Cross-sections vs. Size Parameter")
```

## Validation

TMatrix.jl has been validated against [PyMieScatt](https://github.com/bsunwar/pyMieScatt)
to machine precision for efficiency quantities (``Q_ext``, ``Q_sca``, ``Q_abs``).

See the `validation/mie/` directory in the repository for the validation scripts.

## Interpreting the Results

- **Small ``x``**: In the Rayleigh regime, ``C_sca \propto x^4`` and ``C_{abs} \propto x``.
- **Resonances**: For moderate ``x``, cross-sections show oscillatory behavior due to
  interference between different multipole orders (Mie resonances).
- **Large ``x``**: Cross-sections approach the geometric-optics limit
  ``C_{ext} \to 2 \pi a^2`` (the extinction paradox).

## Exercises

1. Vary the refractive index `m` and observe how resonance positions shift.
2. Compute the asymmetry parameter `g` as a function of `x`.
3. Compare `mie_nmax` with a manually chosen truncation order.
