"""Generate TMatrix.jl output from input CSV.

Usage:
    julia --project=../.. generate_julia_output.jl

Reads:  ../inputs/test_cases.csv
Writes: ../outputs/julia_results.csv
"""

using CSV
using DataFrames
using TMatrix

# Paths
base_dir = dirname(dirname(@__FILE__))
input_csv = joinpath(base_dir, "inputs", "test_cases.csv")
output_csv = joinpath(base_dir, "outputs", "julia_results.csv")

function main()
    if !isfile(input_csv)
        error("Input file not found: $input_csv")
    end

    mkpath(dirname(output_csv))

    df = CSV.read(input_csv, DataFrame)
    results = DataFrame(
        name = String[],
        m_real = Float64[],
        m_imag = Float64[],
        wavelength_nm = Float64[],
        diameter_nm = Float64[],
        Qext = Float64[],
        Qsca = Float64[],
        Qabs = Float64[],
        g = Float64[],
        Cext = Float64[],
        Csca = Float64[],
        Cabs = Float64[],
    )

    for row in eachrow(df)
        name = row.name
        m = row.m_real + row.m_imag * im
        λ = row.wavelength_nm * 1e-9  # nm → m
        r = row.diameter_nm * 1e-9 / 2  # nm → m, diameter → radius

        T = solve_tmatrix(Sphere(r), MieMethod(), λ, m)
        C = calc_cross_sections(T)

        # Convert cross-sections to nm² for comparison
        push!(results, (
            name,
            row.m_real,
            row.m_imag,
            row.wavelength_nm,
            row.diameter_nm,
            C.Q_ext,
            C.Q_sca,
            C.Q_abs,
            C.asymmetry,
            C.extinction * 1e18,
            C.scattering * 1e18,
            C.absorption * 1e18,
        ))
    end

    CSV.write(output_csv, results)
    println("TMatrix.jl results written to $output_csv")
    for row in eachrow(results)
        println("  $(row.name): Qext=$(round(row.Qext, digits=6)), Qsca=$(round(row.Qsca, digits=6)), g=$(round(row.g, digits=6))")
    end
end

main()
