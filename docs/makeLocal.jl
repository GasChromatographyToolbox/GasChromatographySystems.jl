# # make file only for local make of the document.
# this result in errors in Travis
#
# ## To build the documentation locally
#
# ### Running inside Julia REPL:
#  cd to docs folder `cd("path to docs")` and run the following command:
# ```
# include("makeLocal.jl")
# ```
#
# ### Running inside OS Terminal:
# cd to docs folder using OS terminal and run the following command (julia path should be added to OS path):
# ```
# julia --color=yes makeLocal.jl
# ```


using Pkg
pkg"activate .."
push!(LOAD_PATH,"../src/")

#

using Documenter, GasChromatographySystems
using Graphs, CairoMakie

makedocs(;
    modules=[GasChromatographySystems],
	format = Documenter.HTML(
        prettyurls = prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    pages=[
        "Home" => "index.md",
        "Instalation" => "installation.md",
        "Definition of a GC system" => "definition_sys.md",
        "Flow calculation" => "flowcalc.md",
        "Simulation" => "simulation.md",
        "Examples" => "examples.md"#,
        #"Docstrings" => "docstrings.md"
    ],
    repo="github.com/JanLeppert/GasChromatographySystems.jl",
    sitename="GasChromatographySystems",
    authors="Jan Leppert",
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
# deploydocs(;
#     repo="github.com/juliamatlab/MatLang.git",
# )