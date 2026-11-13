module GasChromatographySystems

using Reexport
@reexport using Graphs
using GraphMakie
@reexport using Symbolics
using NetworkLayout
using CairoMakie
using GasChromatographySimulator
using Interpolations
@reexport using CSV
@reexport using DataFrames
using Plots
using SpecialFunctions
using OrdinaryDiffEq
using Integrals
using LsqFit

# some constants
const Tst = 273.15 # K
const R = 8.31446261815324 # J mol⁻¹ K⁻¹
const Tn = 25.0 + Tst # K
const pn = 101300 # Pa

include("./Structures.jl")
include("./Flowcalc.jl")
include("./Systems.jl")
include("./SystemToParameters.jl")
include("./SolvingSystems.jl")
include("./ThermalModulator.jl")
include("./Misc.jl")
include("./Plot.jl")
include("./GCxGC.jl")

end # module
