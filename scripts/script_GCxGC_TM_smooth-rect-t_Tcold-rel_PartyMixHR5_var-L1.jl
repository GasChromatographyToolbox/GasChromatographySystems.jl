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
- variation of the length of 1st D column and its influence on the retention times
- modulation shift set to 3.50 s
=#

optTM = GasChromatographySystems.ModuleTMopt(Tcold_abs=false, ng=true, alg=Vern9(), tflank=20.0, sflank=Inf, dtinit=0.5e-7)
optCol = GasChromatographySystems.ModuleColumnOpt(ng=true)

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
L1 = 28.0:0.5:32.0 # m
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
shift = 0.0
thot = 0.35
ratio = (PM-thot)/PM
Thot = 25.0
Tcold = -130.0

### Flows and pressures
F = 0.8
pin = NaN
pout = 0.0

### System
sys = Array{GasChromatographySystems.System}(undef, length(L1))
Threads.@threads for i=1:length(L1)
	sys[i] = GasChromatographySystems.GCxGC_TM(L1[i], d1, df1, sp1, TP1, L2, d2, df2, sp2, TP2, LTL, dTL, dfTL, spTL, TTL, LM, dM, dfM, spM, shift, PM, ratio, Thot, Tcold, TPM, F, pin, pout; optTM=optTM, optCol=optCol)
end

### Measurement
meas = DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/GCxGC/Simulation/PartyMixHR5_RT.csv", header=1, silencewarnings=true, stringtype=String))
filter!([:CAS] => x -> !ismissing(x), meas)
### Substances
db = DataFrame(urldownload("https://raw.githubusercontent.com/JanLeppert/RetentionData/main/Databases/GCSim_database_nonflag.csv"))
insertcols!(db, 1, :No => collect(1:length(db.Name)))
filter!([:phi0, :CAS, :Phase] => (x, y, z) -> x == 0.001 && y in meas.CAS && z in [sp1, sp2], db)

selected_solutes_ = GasChromatographySystems.common_solutes(db, sys[1]).Name
selected_solutes = selected_solutes_[17:end]#[Not([13, 15, 16])]
# 85 -> multiple entries in database for one of the stat. phases

## Simulation
par = Array{Array{GasChromatographySimulator.Parameters,1}}(undef, length(L1))
Threads.@threads for i=1:length(L1)
    par[i] = GasChromatographySystems.graph_to_parameters(sys[i], db, selected_solutes)
end
paths = GasChromatographySystems.all_paths(sys[1].g, 1)[2]
sim = Array{Tuple{Vector{String}, Vector{Vector{DataFrame}}, Vector{Vector{Any}}, Vector{GasChromatographySimulator.Parameters}}}(undef, length(L1))
Threads.@threads for i=1:length(L1)
#for i=1:length(shift)
    #println("i: $(i)")
    try
        sim[i] = GasChromatographySystems.simulate_along_paths(sys[i], paths, par[i])
    catch
        sim[i] = undef
    end
end

## Peaklist
pl_GCxGC = Array{DataFrame}(undef, length(L1))
Threads.@threads for i=1:length(L1)
	pl_GCxGC[i] = GasChromatographySystems.peaklist_GCxGC(sim[i][2][1][8], PM)
end

## Result
Names = String[]
L1s = Float64[]
tR1s = Float64[]
tR2s = Float64[]
for j=1:length(L1)
    for i=1:length(selected_solutes)
        push!(Names, pl_GCxGC[j].Name[i])
        push!(L1s, L1[j])
        push!(tR1s, pl_GCxGC[j].tR1[i])
        push!(tR2s, pl_GCxGC[j].tR2[i])
    end
end
res = DataFrame(Name=Names, L1=L1s, tR1=tR1s, tR2=tR2s)

## Save Result
CSV.write(joinpath(pwd(), "scripts", "var-L1_PartyMixHR5_shift0.csv"), res)

settings = Dict(fieldnames(GasChromatographySystems.System) .=> getfield.(Ref(sys[1]), fieldnames(GasChromatographySystems.System)))
CSV.write(joinpath(pwd(), "scripts", "var-L1_PartyMixHR5_shift0_settings.csv"), settings)