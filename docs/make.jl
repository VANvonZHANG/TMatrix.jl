using Documenter
using TMatrix

makedocs(
    sitename = "TMatrix.jl",
    format = Documenter.HTML(
        prettyurls = false,
        canonical = "https://vanvonzhang.github.io/T-Matrix.jl",
    ),
    modules = [TMatrix],
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Getting Started" => "tutorials/getting-started.md",
            "Mie Scattering" => "tutorials/mie-scattering.md",
        ],
        "API Reference" => [
            "Special Functions" => "api/special-functions.md",
            "Composed Functions" => "api/composed-functions.md",
            "Vector Spherical Wave Functions" => "api/vswf.md",
            "Geometry" => "api/geometry.md",
            "Mie Solver" => "api/mie-solver.md",
        ],
        "Theory" => [
            "Conventions" => "theory/conventions.md",
            "Special Functions" => "theory/special-functions.md",
            "Vector Spherical Harmonics" => "theory/vector-spherical-harmonics.md",
            "VSWFs" => "theory/vswf.md",
            "Mie Theory" => "theory/mie-theory.md",
        ],
        "Developer Docs" => [
            "Architecture" => "dev/architecture.md",
            "Numerical Stability" => "dev/numerical-stability.md",
            "Solver Interface" => "dev/solver-interface.md",
        ],
        "References" => "references.md",
    ],
    warnonly = [:missing_docs],
)
