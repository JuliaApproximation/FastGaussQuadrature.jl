using Documenter
using FastGaussQuadrature

# Setup for doctests in docstrings
DocMeta.setdocmeta!(FastGaussQuadrature, :DocTestSetup, :(using LinearAlgebra, SpecialFunctions, FastGaussQuadrature))

makedocs(;
    modules = [FastGaussQuadrature],
    format = Documenter.HTML(
        canonical = "https://JuliaApproximation.github.io/FastGaussQuadrature.jl/stable/",
        assets = ["assets/favicon.ico"],
        repolink = "https://github.com/JuliaApproximation/FastGaussQuadrature.jl"
    ),
    pages = [
        "Home" => "index.md",
        "Gaussian Quadrature" => "gaussquadrature.md",
        "Benchmark" => "benchmark.md",
        "Roots of Bessel function" => "besselroots.md",
        "Misc." => "misc.md",
        "References" => "reference.md",
    ],
    repo = "https://github.com/JuliaApproximation/FastGaussQuadrature.jl/blob/{commit}{path}#L{line}",
    sitename = "FastGaussQuadrature.jl",
    authors = "Alex Townsend, Sheehan Olver, Peter Opsomer, and contributors.",
)

deploydocs(; repo = "github.com/JuliaApproximation/FastGaussQuadrature.jl")
