using Documenter
using GasChromatographySystems

makedocs(
    sitename = "GasChromatographySystems",
    format = Documenter.HTML(),
    modules = [GasChromatographySystems]
    pages = Any[
                "Home" => "index.md",
                "Instalation" => "installation.md",
                "Definition of a GC system" => "definition_sys.md",
                "Flow calculation" => "flowcalc.md",
                "Simulation" => "simulation.md",
                "Examples" => "examples.md",
                "Docstrings" => "docstrings.md"
            ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/JanLeppert/GasChromatographySystems.jl"#,
    #devbranch = "main"
)
