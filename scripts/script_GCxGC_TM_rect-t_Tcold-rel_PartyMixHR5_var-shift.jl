### A Pluto.jl notebook ###
# v0.19.26

#using Markdown
#using InteractiveUtils

# ╔═╡ 113c6e70-2168-11ee-3f7f-2775a67bab90
#begin
#	import Pkg
    # activate the shared project environment
#	Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
	#Pkg.upgrade_manifest()
	#Pkg.precompile()
 #   Pkg.instantiate()

	using CSV
    using DataFrames
	#using Plots
    #using CairoMakie
    #using GraphMakie
	#using Graphs
    #using NetworkLayout
    #using Symbolics
	using GasChromatographySimulator
	#using PlutoUI
	using OrdinaryDiffEq
	#using LsqFit
	#using Interpolations
	using UrlDownload
	include(joinpath(pwd(), "src", "GasChromatographySystems.jl"))

#=
- smooth rectangle temperature function over time at modulation point (`spatial=false`, `tflank=20`, `algTM=Vern9` and `ng=true`)
- Tcold is relative (`Tcold_abs=false`)
- variation of the modulator shift and its influence on the retention times
=#

optTM = GasChromatographySystems.ModuleTMOptions(Tcold_abs=false, ng=true, alg="simplifiedTM", tflank=Inf, sflank=Inf, dtinit=0.5e-7)
optCol = GasChromatographySystems.ModuleColumnOptions(ng=true)

## Definition System
### Temperature program
tsteps_, Tsteps_ = GasChromatographySimulator.conventional_program([50.0, 2.0, 5.0, 200.0, 0.0, 5.0, 280.0, 10.0]) 
GCxGC_TP = GasChromatographySystems.TemperatureProgram(tsteps_, Tsteps_)
TP1 = GCxGC_TP
TPM = GCxGC_TP
TP2 = GCxGC_TP
TTL = 280.0

### Columns
# 1st D 
L1 = 29.5 # m
d1 = 0.25 # mm
df1 = 0.25 # µm
sp1 = "Rxi5SilMS"
# modulator
LM = [0.30, 0.005, 0.9, 0.005, 0.3]#[0.3, 0.005, 0.9, 0.005, 0.3]
dM = 0.25
dfM = 0.25
spM = "Rxi17SilMS"
# 2nd D
L2 = 0.245#0.56 # m
d2 = 0.25 # mm
df2 = 0.25 # µm
sp2 = "Rxi17SilMS"
# TL
LTL = 0.235
dTL = 0.25
dfTL = 0.25
spTL = "Rxi17SilMS"

### Modulator
PM = 4.0
shift = 0.0:0.1:4.0#3.35#(PM-0.35)#0.0#-0.35*2
thot = 0.35
ratio = (PM-thot)/PM
Thot = 25.0
Tcold = -130.0

### Flows and pressures
F = 0.8
pin = NaN
pout = 0.0

### System
sys = Array{GasChromatographySystems.System}(undef, length(shift))
Threads.@threads for i=1:length(shift)
	sys[i] = GasChromatographySystems.GCxGC_TM(L1, d1, df1, sp1, TP1, L2, d2, df2, sp2, TP2, LTL, dTL, dfTL, spTL, TTL, LM, dM, dfM, spM, shift[i], PM, ratio, Thot, Tcold, TPM, F, pin, pout; optTM=optTM, optCol=optCol)
end

### Measurement
meas = DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/GCxGC/Simulation/PartyMixHR5_RT.csv", header=1, silencewarnings=true, stringtype=String))
filter!([:CAS] => x -> !ismissing(x), meas)
### Substances
db = DataFrame(urldownload("https://raw.githubusercontent.com/JanLeppert/RetentionData/main/Databases/GCSim_database_nonflag.csv"))
insertcols!(db, 1, :No => collect(1:length(db.Name)))
filter!([:phi0, :CAS, :Phase] => (x, y, z) -> x == 0.001 && y in meas.CAS && z in [sp1, sp2], db)

selected_solutes_ = GasChromatographySystems.common_solutes(db, sys[1]).Name
selected_solutes = selected_solutes_[17:end]

## Simulation
par = Array{Array{GasChromatographySimulator.Parameters,1}}(undef, length(shift))
Threads.@threads for i=1:length(shift)
    par[i] = GasChromatographySystems.graph_to_parameters(sys[i], db, selected_solutes)
end
paths = GasChromatographySystems.all_paths(sys[1].g, 1)[2]
sim = Array{Tuple{Vector{String}, Vector{Vector{DataFrame}}, Vector{Vector{Any}}, Vector{GasChromatographySimulator.Parameters}}}(undef, length(shift))
Threads.@threads for i=1:length(shift)
#for i=1:length(shift)
    #println("i: $(i)")
    sim[i] = GasChromatographySystems.simulate_along_paths(sys[i], paths, par[i])
end

## Peaklist
pl_GCxGC = Array{DataFrame}(undef, length(shift))
Threads.@threads for i=1:length(shift)
	pl_GCxGC[i] = GasChromatographySystems.peaklist_GCxGC(sim[i][2][1][8], PM)
end

## Result
Names = String[]
shifts = Float64[]
tR1s = Float64[]
tR2s = Float64[]
for j=1:length(shift)
    for i=1:length(selected_solutes)
        push!(Names, pl_GCxGC[j].Name[i])
        push!(shifts, shift[j])
        push!(tR1s, pl_GCxGC[j].tR1[i])
        push!(tR2s, pl_GCxGC[j].tR2[i])
    end
end
res = DataFrame(Name=Names, shift=shifts, tR1=tR1s, tR2=tR2s)

## Save Result
CSV.write(joinpath(pwd(), "scripts", "var-shift_PartyMixHR5_tflankInf.csv"), res)

settings = Dict(fieldnames(GasChromatographySystems.System) .=> getfield.(Ref(sys[1]), fieldnames(GasChromatographySystems.System)))
CSV.write(joinpath(pwd(), "scripts", "var-shift_PartyMixHR5_settings_tflankInf.csv"), settings)