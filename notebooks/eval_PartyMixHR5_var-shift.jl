### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 7656e75e-2e47-11ee-31b0-f9a794c245bb
begin
	import Pkg
    # activate the shared project environment
	Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
	#Pkg.upgrade_manifest()
	#Pkg.precompile()
    Pkg.instantiate()

	using CSV, DataFrames
	using Plots, CairoMakie, GraphMakie
	using Graphs, NetworkLayout, Symbolics
	using GasChromatographySimulator
	using PlutoUI
	using OrdinaryDiffEq
	using LsqFit
	using Interpolations
	using UrlDownload
	include(joinpath(dirname(pwd()), "src", "GasChromatographySystems.jl"))
	TableOfContents(depth=4)
end

# ╔═╡ dc716bc4-153c-4aef-9d00-55f99ab3f589
md"""
# Evaluation of variation of modulation shift

Using the simulation of GCxGC with looped thermal modulation with a 
- smoothed rectangular temperature function over time at the modulation points
- uniform temperature over the spot
- the cooling temperature is relative (reduction of the programmed oven temperature by the value of the cooling temperature)
"""

# ╔═╡ c7b83c87-2a13-487c-bc45-6915317201d0
md"""
## Measurement
"""

# ╔═╡ 10353347-1f2c-4f00-9a7f-cf5518eec85d
begin
	meas = DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/GCxGC/Simulation/PartyMixHR5_RT.csv", header=1, silencewarnings=true, stringtype=String))
	filter!([:CAS] => x -> !ismissing(x), meas)
end

# ╔═╡ 4c15c7db-01dc-4de5-a60b-5fd56c503926
md"""
## Simulation results
"""

# ╔═╡ 41d3693b-bfb4-4732-9a59-d92640f1ab79
sim = DataFrame(CSV.File("/Users/janleppert/Documents/GitHub/GasChromatographySystems/scripts/var-shift_PartyMixHR5.csv", header=1, silencewarnings=true, stringtype=String))

# ╔═╡ ef286a7a-246c-4e2c-a257-6ab301ad8590
function eval_varshift(sim)
	name = unique(sim.Name)
	shift = unique(sim.shift)
	results = Array{DataFrame}(undef, length(name))
	p_tR1 = Array{Plots.Plot}(undef, length(name))
	p_tR2 = Array{Plots.Plot}(undef, length(name))
	p_tR1s = Plots.plot(xlabel="modulation shift in s", ylabel="1st dim retention time tR1 in s")
	p_tR2s = Plots.plot(xlabel="modulation shift in s", ylabel="2nd dim retention time tR2 in s")
	p_ΔtR1s = Plots.plot(xlabel="modulation shift in s", ylabel="1st dim difference meas-sim ΔtR1 in s")
	p_ΔtR2s = Plots.plot(xlabel="modulation shift in s", ylabel="2nd dim difference meas-sim ΔtR2 in s")
	shift_opt_tR2 = fill(NaN, length(name))
	for i=1:length(name)
		results[i] = sort(filter([:Name] => x -> x == name[i], sim), :shift)
		p_tR1[i] = Plots.plot(results[i].shift, results[i].tR1, xlabel="modulation shift in s", ylabel="1st dim retention time tR1 in s", title=results[i].Name[1], label="")
		p_tR2[i] = Plots.plot(results[i].shift, results[i].tR2, xlabel="modulation shift in s", ylabel="2nd dim retention time tR2 in s", title=results[i].Name[1], label="")
		Plots.plot!(p_tR1s, results[i].shift, results[i].tR1, label=results[i].Name[1])
		Plots.plot!(p_tR2s, results[i].shift, results[i].tR2, label=results[i].Name[1])
		if !ismissing(meas.tR1[findfirst(name[i].==meas.Name)]) && !isnothing(meas.tR1[findfirst(name[i].==meas.Name)])
			Plots.plot!(p_ΔtR1s, results[i].shift, meas.tR1[findfirst(name[i].==meas.Name)] .- results[i].tR1, label=results[i].Name[1])
			Plots.plot!(p_ΔtR2s, results[i].shift, meas.tR2[findfirst(name[i].==meas.Name)] .- results[i].tR2, label=results[i].Name[1])

			itp = LinearInterpolation((meas.tR2[findfirst(name[i].==meas.Name)] .- results[i].tR2[21:41], ), shift[21:41], extrapolation_bc=Flat())
			shift_opt_tR2[i] = itp(0.0)
		end
	end
	return p_tR1, p_tR2, p_tR1s, p_tR2s, p_ΔtR1s, p_ΔtR2s, shift_opt_tR2
end

# ╔═╡ 74e6708e-985b-432d-87d6-9354ce903a31
function optimal_shift(res)
	opt_shift = sum(res[7][isnan.(res[7]).==false])/length(res[7][isnan.(res[7]).==false]), maximum(res[7][isnan.(res[7]).==false]), minimum(res[7][isnan.(res[7]).==false])
	return opt_shift
end

# ╔═╡ 0d8cddb6-1b32-4abb-a71d-4063bd0a7c0d
res = eval_varshift(sim)

# ╔═╡ 7feb935b-aba6-4c7b-8138-a45777978c55
Plots.plot!(res[6], legend=:bottomright)

# ╔═╡ e281594a-0faf-4f50-915c-5bca8abb6969
opt_shift_tR2 = optimal_shift(res)

# ╔═╡ 9bd6dd6a-f8ed-476d-98e6-0ff9da5288f1
md"""
### Optima regarding tR2

The the optima is where the difference between measured and simulated 2nd dimension retention time is 0.

Mean value of the optima: $(opt_shift_tR2[1]) s

Minima: $(opt_shift_tR2[3]) s

Maxima: $(opt_shift_tR2[2]) s

The optima is different for different solutes. 

The edge, where ΔtR2 changes abruptly is different for different solutes. With this the maxima ΔtR2 values are different for different solutes.
"""

# ╔═╡ ef1942c3-8076-4830-8c5c-aee353ca1add
md"""
### Optima regarding tR1

The the optima is where the difference between measured and simulated 2nd dimension retention time is 0.

No general optima regarding tR1 can be found. Only few solutes have a ΔtR1 value of 0 inside the tested range. (this is periodically, therefore the value 0 is never reached for the other solutes)

An edge with a larger change of ΔtR1 is observed in the same range of the shift as for ΔtR1. 

Occasionally jumps in ΔtR1 for single values of shift are observed. This could be an error/artefact from the simulation?
"""

# ╔═╡ d59520b4-5707-4236-b91a-56350025b0e5
res[5]

# ╔═╡ 97dc29d4-c0dd-4dac-a0d2-805e0eefc13c
md"""
## tflank = 12
"""

# ╔═╡ adeaf86a-4813-4a87-b8c3-35607d17158d
sim_tflank12 = DataFrame(CSV.File("/Users/janleppert/Documents/GitHub/GasChromatographySystems/scripts/var-shift_PartyMixHR5_tflank12.csv", header=1, silencewarnings=true, stringtype=String))

# ╔═╡ c7138529-ddd5-43a0-815b-64a02c3c1140
res_tflank12 = eval_varshift(sim_tflank12)

# ╔═╡ 11ebcf34-9d24-41ac-bc34-c9b114592b1f
opt_shift_tR2_tflank12 = optimal_shift(res_tflank12)

# ╔═╡ 28d4fdc6-b68d-442d-88cc-2270502ae1ba
Plots.plot!(res_tflank12[6], legend=:bottomright)

# ╔═╡ 83c54269-7275-4090-86b7-fc02a6cb9752
res_tflank12[5]

# ╔═╡ 127dfbeb-3857-4843-8e18-48e3a94eb82e
md"""
## tflank = 40
"""

# ╔═╡ 0c192eb2-0b58-47cd-a649-99a54d3d51ea
sim_tflank40 = DataFrame(CSV.File("/Users/janleppert/Documents/GitHub/GasChromatographySystems/scripts/var-shift_PartyMixHR5_tflank40.csv", header=1, silencewarnings=true, stringtype=String))

# ╔═╡ efd7f743-5233-498f-a558-c3a80a3a9e1a
res_tflank40 = eval_varshift(sim_tflank40)

# ╔═╡ a44a3db1-a0e2-4f83-bcc1-266098401db3
opt_shift_tR2_tflank40 = optimal_shift(res_tflank40)

# ╔═╡ a8381176-407b-419c-9d18-e01ae6fe9f69
Plots.plot!(res_tflank40[6], legend=:bottomright)

# ╔═╡ 62f84322-acde-49de-8a19-63deaa75a003
res_tflank40[5]

# ╔═╡ a61089be-7986-4fd7-a8f9-17474f2d5a8a
md"""
## tflank = 100
"""

# ╔═╡ 371e047b-ef3e-407d-bf52-c8f8abacc731
sim_tflank100 = DataFrame(CSV.File("/Users/janleppert/Documents/GitHub/GasChromatographySystems/scripts/var-shift_PartyMixHR5_tflank100.csv", header=1, silencewarnings=true, stringtype=String))

# ╔═╡ 86a1b245-eb20-478f-90f0-619eed4d651c
res_tflank100 = eval_varshift(sim_tflank100)

# ╔═╡ b6e5e8b4-85ed-4ceb-bc3f-3369f649c448
opt_shift_tR2_tflank100 = optimal_shift(res_tflank100)

# ╔═╡ 6fb4587d-8ad1-43e9-be90-1e66eeddf9ed
Plots.plot!(res_tflank100[6], legend=:bottomright)

# ╔═╡ 574a8271-3e2b-4f47-b2e6-2c9af66024c8
res_tflank100[5]

# ╔═╡ 69cbba05-3f41-4020-a0af-7ae643580fbb
md"""
## tflank = 1000
"""

# ╔═╡ 565b6767-1e6b-4c43-be9c-d51efe19a2ca
sim_tflank1000 = DataFrame(CSV.File("/Users/janleppert/Documents/GitHub/GasChromatographySystems/scripts/var-shift_PartyMixHR5_tflank1000.csv", header=1, silencewarnings=true, stringtype=String))

# ╔═╡ 3a3a313c-27ba-47ca-bd6e-be1330d646be
res_tflank1000 = eval_varshift(sim_tflank1000)

# ╔═╡ 4f93d29f-a6c1-466f-8758-eb25564b1227
opt_shift_tR2_tflank1000 = optimal_shift(res_tflank1000)

# ╔═╡ 6bed4880-b653-42ff-95d2-c6f0a771df22
Plots.plot!(res_tflank1000[6], legend=:bottomright)

# ╔═╡ 18fcb267-69c3-4e13-9e8f-fe754e8f003a
res_tflank1000[5]

# ╔═╡ 8416b15b-9f8e-41e1-ae2b-476096cc26fc
md"""
## tflank = Inf

`alg="simplifiedTM"`
"""

# ╔═╡ 5193b7d5-05bc-40b8-9f9f-cf549f4e57a6
sim_tflankInf = DataFrame(CSV.File("/Users/janleppert/Documents/GitHub/GasChromatographySystems/scripts/var-shift_PartyMixHR5_tflankInf.csv", header=1, silencewarnings=true, stringtype=String))

# ╔═╡ 0655003d-df28-4115-b24f-cd8124b1fff9
res_tflankInf = eval_varshift(sim_tflankInf)

# ╔═╡ 06f27a3e-2a59-4815-b97d-d974d288a463
opt_shift_tR2_tflankInf = optimal_shift(res_tflankInf)

# ╔═╡ fd2c66cf-6acd-4f66-aee2-7231637c4c27
Plots.plot!(res_tflankInf[6], legend=:bottomright)

# ╔═╡ 1ff8bab6-24c5-4943-89e2-35ce3adbb01f
res_tflankInf[5]

# ╔═╡ 458dd99e-ae67-4fb1-a62e-efedc7f01772
md"""
## sflank = 40, tflank = 20
"""

# ╔═╡ 971faee5-4ace-4918-9602-2bf922ff2bba
sim_sflank40 = DataFrame(CSV.File("/Users/janleppert/Documents/GitHub/GasChromatographySystems/scripts/var-shift_rect-xt_PartyMixHR5.csv", header=1, silencewarnings=true, stringtype=String))

# ╔═╡ 49964c8f-8022-46a8-aeb2-0a20e0582f61
res_sflank40 = eval_varshift(sim_sflank40)

# ╔═╡ 9362201e-72a5-4f23-9895-4f0fa6764473
opt_shift_tR2_sflank40 = optimal_shift(res_sflank40)

# ╔═╡ a35960ad-451b-468f-82c2-b8c6737cbc0c
Plots.plot!(res_sflank40[6], legend=:bottomright)

# ╔═╡ 76b220e5-4664-4265-8734-9f777f87b449
res_sflank40[5]

# ╔═╡ 97b4b763-44fd-4f74-89de-78547e6f6858
md"""
## Plot Modulation
"""

# ╔═╡ ee59212e-9c44-4fbc-913a-bff91845286c
begin
	tsteps_, Tsteps_ = GasChromatographySimulator.conventional_program([50.0, 2.0, 5.0, 200.0, 0.0, 5.0, 280.0, 10.0]) 
	GCxGC_TP = GasChromatographySystems.TemperatureProgram(tsteps_, Tsteps_)
end

# ╔═╡ 11eb9a57-1ee5-40ae-9c7d-50ebd3d5d1ed
sys_tflank20 = GasChromatographySystems.GCxGC_TM(29.5, 0.25, 0.25, "Rxi5SilMS", GCxGC_TP, 0.245, 0.25, 0.25, "Rxi17SilMS", GCxGC_TP, 0.235, 0.25, 0.25, "Rxi17SilMS", 280.0, [0.3, 0.005, 0.9, 0.005, 0.3], 0.25, 0.25, "Rxi17SilMS", opt_shift_tR2[1], 4.0, (4.0-0.35)/4.0, 25.0, -130.0, GCxGC_TP, 0.8, NaN, 0.0; optTM=GasChromatographySystems.ModuleTMOptions(Tcold_abs=false, ng=true, alg=Vern9(), tflank=20.0, sflank=Inf, dtinit=0.5e-7), optCol=GasChromatographySystems.ModuleColumnOptions(ng=true))

# ╔═╡ 97fa54b4-a3d7-402c-9b44-8ac75211562d
sys_tflank12 = GasChromatographySystems.GCxGC_TM(29.5, 0.25, 0.25, "Rxi5SilMS", GCxGC_TP, 0.245, 0.25, 0.25, "Rxi17SilMS", GCxGC_TP, 0.235, 0.25, 0.25, "Rxi17SilMS", 280.0, [0.3, 0.005, 0.9, 0.005, 0.3], 0.25, 0.25, "Rxi17SilMS", opt_shift_tR2_tflank12[1], 4.0, (4.0-0.35)/4.0, 25.0, -130.0, GCxGC_TP, 0.8, NaN, 0.0; optTM=GasChromatographySystems.ModuleTMOptions(Tcold_abs=false, ng=true, alg=Vern9(), tflank=12.0, sflank=Inf, dtinit=0.5e-7), optCol=GasChromatographySystems.ModuleColumnOptions(ng=true))

# ╔═╡ feb9887a-2c73-4dfe-a823-6c12daa6cd98
sys_tflank40 = GasChromatographySystems.GCxGC_TM(29.5, 0.25, 0.25, "Rxi5SilMS", GCxGC_TP, 0.245, 0.25, 0.25, "Rxi17SilMS", GCxGC_TP, 0.235, 0.25, 0.25, "Rxi17SilMS", 280.0, [0.3, 0.005, 0.9, 0.005, 0.3], 0.25, 0.25, "Rxi17SilMS", opt_shift_tR2_tflank40[1], 4.0, (4.0-0.35)/4.0, 25.0, -130.0, GCxGC_TP, 0.8, NaN, 0.0; optTM=GasChromatographySystems.ModuleTMOptions(Tcold_abs=false, ng=true, alg=Vern9(), tflank=40.0, sflank=Inf, dtinit=0.5e-7), optCol=GasChromatographySystems.ModuleColumnOptions(ng=true))

# ╔═╡ 5c2b81c6-6e0f-4d5c-a7e4-7da647c3e4c5
sys_tflank100 = GasChromatographySystems.GCxGC_TM(29.5, 0.25, 0.25, "Rxi5SilMS", GCxGC_TP, 0.245, 0.25, 0.25, "Rxi17SilMS", GCxGC_TP, 0.235, 0.25, 0.25, "Rxi17SilMS", 280.0, [0.3, 0.005, 0.9, 0.005, 0.3], 0.25, 0.25, "Rxi17SilMS", opt_shift_tR2_tflank100[1], 4.0, (4.0-0.35)/4.0, 25.0, -130.0, GCxGC_TP, 0.8, NaN, 0.0; optTM=GasChromatographySystems.ModuleTMOptions(Tcold_abs=false, ng=true, alg=Vern9(), tflank=100.0, sflank=Inf, dtinit=0.5e-7), optCol=GasChromatographySystems.ModuleColumnOptions(ng=true))

# ╔═╡ 9bb3d0d9-71ac-4f56-85b1-ab6b50c959c8
sys_tflank1000 = GasChromatographySystems.GCxGC_TM(29.5, 0.25, 0.25, "Rxi5SilMS", GCxGC_TP, 0.245, 0.25, 0.25, "Rxi17SilMS", GCxGC_TP, 0.235, 0.25, 0.25, "Rxi17SilMS", 280.0, [0.3, 0.005, 0.9, 0.005, 0.3], 0.25, 0.25, "Rxi17SilMS", opt_shift_tR2_tflank1000[1], 4.0, (4.0-0.35)/4.0, 25.0, -130.0, GCxGC_TP, 0.8, NaN, 0.0; optTM=GasChromatographySystems.ModuleTMOptions(Tcold_abs=false, ng=true, alg=Vern9(), tflank=1000.0, sflank=Inf, dtinit=0.5e-7), optCol=GasChromatographySystems.ModuleColumnOptions(ng=true))

# ╔═╡ 54d7c1eb-f9c4-49d9-a4bd-f06fd41d0ab1
sys_tflankInf = GasChromatographySystems.GCxGC_TM(29.5, 0.25, 0.25, "Rxi5SilMS", GCxGC_TP, 0.245, 0.25, 0.25, "Rxi17SilMS", GCxGC_TP, 0.235, 0.25, 0.25, "Rxi17SilMS", 280.0, [0.3, 0.005, 0.9, 0.005, 0.3], 0.25, 0.25, "Rxi17SilMS", opt_shift_tR2_tflankInf[1], 4.0, (4.0-0.35)/4.0, 25.0, -130.0, GCxGC_TP, 0.8, NaN, 0.0; optTM=GasChromatographySystems.ModuleTMOptions(Tcold_abs=false, ng=true, alg="simplifiedTM", tflank=Inf, sflank=Inf, dtinit=0.5e-7), optCol=GasChromatographySystems.ModuleColumnOptions(ng=true))

# ╔═╡ 938d7398-484e-461e-97d4-fd15847613f6
sys_sflank40 = GasChromatographySystems.GCxGC_TM(29.5, 0.25, 0.25, "Rxi5SilMS", GCxGC_TP, 0.245, 0.25, 0.25, "Rxi17SilMS", GCxGC_TP, 0.235, 0.25, 0.25, "Rxi17SilMS", 280.0, [0.3, 0.005, 0.9, 0.005, 0.3], 0.25, 0.25, "Rxi17SilMS", opt_shift_tR2_sflank40[1], 4.0, (4.0-0.35)/4.0, 25.0, -130.0, GCxGC_TP, 0.8, NaN, 0.0; optTM=GasChromatographySystems.ModuleTMOptions(Tcold_abs=false, ng=false, alg=Vern9(), tflank=20, sflank=40, dtinit=0.5e-7), optCol=GasChromatographySystems.ModuleColumnOptions(ng=true))

# ╔═╡ 0cbf1231-69e1-4ead-b231-9ef51eba222e
sys = [sys_tflank12, sys_tflank20, sys_tflank40, sys_tflank100, sys_tflank1000, sys_tflankInf, sys_sflank40]

# ╔═╡ 9723cf14-2e95-427e-b346-597cec0f516b
begin
	tflank = [12, 20, 40, 100, 1000, Inf, 20]
	opt_shift = [opt_shift_tR2_tflank12[1], opt_shift_tR2[1], opt_shift_tR2_tflank40[1], opt_shift_tR2_tflank100[1], opt_shift_tR2_tflank1000[1], opt_shift_tR2_tflankInf[1], opt_shift_tR2_sflank40[1]]
	Plots.scatter(tflank, opt_shift, xlabel="tflank factor", ylabel="optimal shift in s", legend=false)
end

# ╔═╡ de46a610-2999-4a62-b8cd-959f4b785d53
begin
	plotly()
	p_TM = Plots.plot(title="temperature at modulator point", xlabel="time in s")

	for i=1:length(sys)
		TM_itp = GasChromatographySystems.module_temperature(sys[i].modules[3], sys[i])[5]
		Plots.plot!(p_TM, 0.0:0.001:sys[i].modules[3].PM, TM_itp.(sys[i].modules[3].L/2, 0.0:0.001:sys[i].modules[3].PM).-273.115, label="tflank = $(tflank[i])")
	end
	p_TM
end

# ╔═╡ c818c94c-9278-4f66-976a-b7e817847041
opt_shift

# ╔═╡ f13f2720-45cb-443c-8c5c-75691ae41001
sims = [sim_tflank12, sim, sim_tflank40, sim_tflank100, sim_tflank1000, sim_tflankInf]

# ╔═╡ c5aeca56-3377-4e9e-ab65-86e860fcf5f3
function plot_tR1_tR2(sims, select_solute)
	#select_solute = "Valerophenone"
	shifts = unique(sims[1].shift)
	tR1 = Array{Float64}(undef, length(sims), length(shifts))
	tR2 = Array{Float64}(undef, length(sims), length(shifts))
	for j=1:length(shifts)
		for i=1:length(sims)
			df_filter = filter([:Name, :shift] => (x, y) -> x == select_solute && y == shifts[j], sims[i])
			tR1[i,j] = df_filter.tR1[1]
			tR2[i,j] = df_filter.tR2[1]
		end
	end
	p_tR1 = Plots.plot(xlabel="shift in s", ylabel="tR1 in s", title=select_solute)
	for i=1:length(sims)
		Plots.plot!(shifts, tR1[i,:], label="tflank = $(tflank[i])")
	end
	p_tR2 = Plots.plot(xlabel="shift in s", ylabel="tR2 in s", title=select_solute)
	for i=1:length(sims)
		Plots.plot!(shifts, tR2[i,:], label="tflank = $(tflank[i])")
	end
	return p_tR1, p_tR2
end

# ╔═╡ 074d1d38-1aff-497b-b8ad-0ab9bdc55214
begin
	names = unique(sims[1].Name)
	p_tR1s = Array{Plots.Plot}(undef, length(names))
	p_tR2s = Array{Plots.Plot}(undef, length(names))
	for i=1:length(names)
		p_tR1s[i], p_tR2s[i] = plot_tR1_tR2(sims, names[i])
	end
end

# ╔═╡ 9be1dfaf-d06d-4275-8884-bdfa6c9c4b90
p_tR1s[5]

# ╔═╡ 86b7acf1-42db-4481-8afe-0b30e202ea15
p_tR2s[5]

# ╔═╡ d995457c-6dae-4601-a68e-5e3ce1a98f8b
md"""
## Simulations at optimal shift
"""

# ╔═╡ e2c03d8b-20e3-4780-b56f-0717aad1ebf8
opt_shift

# ╔═╡ 53a76ef7-65a4-4498-b464-48c236a057e4
begin
	db = DataFrame(urldownload("https://raw.githubusercontent.com/JanLeppert/RetentionData/main/Databases/GCSim_database_nonflag.csv"))
	insertcols!(db, 1, :No => collect(1:length(db.Name)))
	filter!([:phi0, :CAS, :Phase] => (x, y, z) -> x == 0.001 && y in meas.CAS && z in [sys[1].modules[1].sp, sys[1].modules[end].sp], db)
	
	selected_solutes_ = GasChromatographySystems.common_solutes(db, sys[1]).Name
	selected_solutes = selected_solutes_[17:end]
end

# ╔═╡ 0528f27c-f191-4347-b68b-8812b3f4f434
begin
	par_opt = Array{Array{GasChromatographySimulator.Parameters,1}}(undef, length(sys))
	Threads.@threads for i=1:length(sys)
	    par_opt[i] = GasChromatographySystems.graph_to_parameters(sys[i], db, selected_solutes)
	end
	paths = GasChromatographySystems.all_paths(sys[1].g, 1)[2]
	sim_opt = Array{Tuple{Vector{String}, Vector{Vector{DataFrame}}, Vector{Vector{Any}}, Vector{GasChromatographySimulator.Parameters}}}(undef, length(sys))
	Threads.@threads for i=1:length(sys)
	#for i=1:length(shift)
	    #println("i: $(i)")
	    sim_opt[i] = GasChromatographySystems.simulate_along_paths(sys[i], paths, par_opt[i])
	end
end

# ╔═╡ b00f669a-1ef2-421a-8169-6d3167c84957
sim_opt[1][2][1][end]

# ╔═╡ 72c2eeca-a1b6-46ab-85d4-622ea1193d54
begin
	pls_GCxGC = Array{DataFrame}(undef, length(sys))
	for i=1:length(sys)
		pls_GCxGC[i] = GasChromatographySystems.peaklist_GCxGC(sim_opt[i][2][1][end], 4.0)
	end
end

# ╔═╡ cb64fb1d-659e-4285-83ce-f60d8e5cd434
begin
	compare = Array{DataFrame}(undef, length(sys))
	avreltR1 = Array{Float64}(undef, length(sys))
	avreltR2 = Array{Float64}(undef, length(sys))
	for i=1:length(sys)
		compare[i] = GasChromatographySystems.comparison_meas_sim(meas, pls_GCxGC[i])
		avreltR1[i] = sum(abs.(collect(skipmissing(compare[i].relΔtR1_percent))))/length(collect(skipmissing(compare[i].relΔtR1_percent)))
		avreltR2[i] = sum(abs.(collect(skipmissing(compare[i].relΔtR2_percent))))/length(collect(skipmissing(compare[i].relΔtR2_percent)))
	end	
end

# ╔═╡ 372d9292-21ee-495b-bcd1-aa9e008eefc9
compare

# ╔═╡ 273722f3-44e8-4343-8274-8fe4a3fcc6fd
DataFrame(tflank=tflank, avreltR1=avreltR1, avreltR2=avreltR2)

# ╔═╡ 0b70e24c-5187-4350-a340-c74566f4d1a3
Plots.scatter(tflank, avreltR1)

# ╔═╡ 931450ba-18a3-41e3-840b-d0dd3ada6d29
Plots.scatter(tflank, avreltR2)

# ╔═╡ 3dc26e76-23fe-4a64-bc4a-21f4528cc5d4
meas_chrom = sort!(DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/GCxGC/Simulation/PartyMixHR5.csv", header=1, silencewarnings=true)), :RT)

# ╔═╡ dd03770d-c853-4617-91b3-52060d208c18
begin
	c_slices_m, t_D1_m, t_D2_m = GasChromatographySystems.chrom_slicing(meas_chrom.RT.*60, meas_chrom.TIC, 4.0)

	# interpolation in tR2 for the same time basis
	itps = Array{Interpolations.Extrapolation}(undef, length(c_slices_m))
	for i=1:length(c_slices_m)
		itps[i] = LinearInterpolation((t_D2_m[i],), c_slices_m[i], extrapolation_bc=Flat())
	end
	chrom_mat = Array{Float64}(undef, length(c_slices_m), length(0.0:0.01:4.0))
	for i=1:length(c_slices_m)
		chrom_mat[i,:] = itps[i].(0.0:0.01:4.0)
	end
	chrom_mat
end;

# ╔═╡ 761c0b2d-a692-4d8c-b34d-848df79c378f
begin
	plotly()
	chrom_mat_mod = Array{Float64}(undef, size(chrom_mat)[1], size(chrom_mat)[2])
	for j=1:size(chrom_mat)[2]
		for i=1:size(chrom_mat)[1]
			if chrom_mat[i,j] > 25000.0
				chrom_mat_mod[i,j] = 25000.0
			else
				chrom_mat_mod[i,j] = chrom_mat[i,j]
			end
		end
	end
	#gr()
	p_meas = Plots.heatmap(0.0:4.0:(size(chrom_mat)[1]-1)*4.0, 0.0:0.01:4.0, chrom_mat_mod', c=:jet1, colorbar=false, legend=:bottomright);


	# simulation
	# spot marker
	for j=1:length(sys)
		Plots.scatter!(pls_GCxGC[j].tR1, pls_GCxGC[j].tR2, c=j, m=:diamond, markeralpha=1, msize=3, label="tflank=$(tflank[j])", xlims=(0.0, 3500.0), ylims=(0.0, 4.0), xlabel="tR1 in s", ylabel="tR2 in s");
		# line marker
		RT1 = Tuple{Float64, Float64}[]
		RT2 = Tuple{Float64, Float64}[]
		for i=1:length(compare[j].tR1_meas)
			if ismissing(compare[j].tR1_meas[i])==false && 	ismissing(compare[j].tR1_sim[i])==false
				push!(RT1, (compare[j].tR1_meas[i], compare[j].tR1_sim[i]))
				push!(RT2, (compare[j].tR2_meas[i], compare[j].tR2_sim[i]))
			end
		end
		Plots.plot!(RT1, RT2, c=j, linewidth=2, label="")
	end
	p_meas
end

# ╔═╡ Cell order:
# ╠═7656e75e-2e47-11ee-31b0-f9a794c245bb
# ╠═dc716bc4-153c-4aef-9d00-55f99ab3f589
# ╠═c7b83c87-2a13-487c-bc45-6915317201d0
# ╠═10353347-1f2c-4f00-9a7f-cf5518eec85d
# ╠═4c15c7db-01dc-4de5-a60b-5fd56c503926
# ╠═41d3693b-bfb4-4732-9a59-d92640f1ab79
# ╠═ef286a7a-246c-4e2c-a257-6ab301ad8590
# ╠═74e6708e-985b-432d-87d6-9354ce903a31
# ╠═0d8cddb6-1b32-4abb-a71d-4063bd0a7c0d
# ╠═9bd6dd6a-f8ed-476d-98e6-0ff9da5288f1
# ╠═7feb935b-aba6-4c7b-8138-a45777978c55
# ╠═e281594a-0faf-4f50-915c-5bca8abb6969
# ╟─ef1942c3-8076-4830-8c5c-aee353ca1add
# ╠═d59520b4-5707-4236-b91a-56350025b0e5
# ╠═97dc29d4-c0dd-4dac-a0d2-805e0eefc13c
# ╠═adeaf86a-4813-4a87-b8c3-35607d17158d
# ╠═c7138529-ddd5-43a0-815b-64a02c3c1140
# ╠═11ebcf34-9d24-41ac-bc34-c9b114592b1f
# ╠═28d4fdc6-b68d-442d-88cc-2270502ae1ba
# ╠═83c54269-7275-4090-86b7-fc02a6cb9752
# ╠═127dfbeb-3857-4843-8e18-48e3a94eb82e
# ╠═0c192eb2-0b58-47cd-a649-99a54d3d51ea
# ╠═efd7f743-5233-498f-a558-c3a80a3a9e1a
# ╠═a44a3db1-a0e2-4f83-bcc1-266098401db3
# ╠═a8381176-407b-419c-9d18-e01ae6fe9f69
# ╠═62f84322-acde-49de-8a19-63deaa75a003
# ╠═a61089be-7986-4fd7-a8f9-17474f2d5a8a
# ╠═371e047b-ef3e-407d-bf52-c8f8abacc731
# ╠═86a1b245-eb20-478f-90f0-619eed4d651c
# ╠═b6e5e8b4-85ed-4ceb-bc3f-3369f649c448
# ╠═6fb4587d-8ad1-43e9-be90-1e66eeddf9ed
# ╠═574a8271-3e2b-4f47-b2e6-2c9af66024c8
# ╠═69cbba05-3f41-4020-a0af-7ae643580fbb
# ╠═565b6767-1e6b-4c43-be9c-d51efe19a2ca
# ╠═3a3a313c-27ba-47ca-bd6e-be1330d646be
# ╠═4f93d29f-a6c1-466f-8758-eb25564b1227
# ╠═6bed4880-b653-42ff-95d2-c6f0a771df22
# ╠═18fcb267-69c3-4e13-9e8f-fe754e8f003a
# ╠═8416b15b-9f8e-41e1-ae2b-476096cc26fc
# ╠═5193b7d5-05bc-40b8-9f9f-cf549f4e57a6
# ╠═0655003d-df28-4115-b24f-cd8124b1fff9
# ╠═06f27a3e-2a59-4815-b97d-d974d288a463
# ╠═fd2c66cf-6acd-4f66-aee2-7231637c4c27
# ╠═1ff8bab6-24c5-4943-89e2-35ce3adbb01f
# ╠═458dd99e-ae67-4fb1-a62e-efedc7f01772
# ╠═971faee5-4ace-4918-9602-2bf922ff2bba
# ╠═49964c8f-8022-46a8-aeb2-0a20e0582f61
# ╠═9362201e-72a5-4f23-9895-4f0fa6764473
# ╠═a35960ad-451b-468f-82c2-b8c6737cbc0c
# ╠═76b220e5-4664-4265-8734-9f777f87b449
# ╠═97b4b763-44fd-4f74-89de-78547e6f6858
# ╠═ee59212e-9c44-4fbc-913a-bff91845286c
# ╠═11eb9a57-1ee5-40ae-9c7d-50ebd3d5d1ed
# ╠═97fa54b4-a3d7-402c-9b44-8ac75211562d
# ╠═feb9887a-2c73-4dfe-a823-6c12daa6cd98
# ╠═5c2b81c6-6e0f-4d5c-a7e4-7da647c3e4c5
# ╠═9bb3d0d9-71ac-4f56-85b1-ab6b50c959c8
# ╠═54d7c1eb-f9c4-49d9-a4bd-f06fd41d0ab1
# ╠═938d7398-484e-461e-97d4-fd15847613f6
# ╠═0cbf1231-69e1-4ead-b231-9ef51eba222e
# ╠═de46a610-2999-4a62-b8cd-959f4b785d53
# ╠═9723cf14-2e95-427e-b346-597cec0f516b
# ╠═c818c94c-9278-4f66-976a-b7e817847041
# ╠═f13f2720-45cb-443c-8c5c-75691ae41001
# ╠═c5aeca56-3377-4e9e-ab65-86e860fcf5f3
# ╠═074d1d38-1aff-497b-b8ad-0ab9bdc55214
# ╠═9be1dfaf-d06d-4275-8884-bdfa6c9c4b90
# ╠═86b7acf1-42db-4481-8afe-0b30e202ea15
# ╠═d995457c-6dae-4601-a68e-5e3ce1a98f8b
# ╠═e2c03d8b-20e3-4780-b56f-0717aad1ebf8
# ╠═53a76ef7-65a4-4498-b464-48c236a057e4
# ╠═0528f27c-f191-4347-b68b-8812b3f4f434
# ╠═b00f669a-1ef2-421a-8169-6d3167c84957
# ╠═72c2eeca-a1b6-46ab-85d4-622ea1193d54
# ╠═cb64fb1d-659e-4285-83ce-f60d8e5cd434
# ╠═372d9292-21ee-495b-bcd1-aa9e008eefc9
# ╠═273722f3-44e8-4343-8274-8fe4a3fcc6fd
# ╠═0b70e24c-5187-4350-a340-c74566f4d1a3
# ╠═931450ba-18a3-41e3-840b-d0dd3ada6d29
# ╠═3dc26e76-23fe-4a64-bc4a-21f4528cc5d4
# ╠═dd03770d-c853-4617-91b3-52060d208c18
# ╠═761c0b2d-a692-4d8c-b34d-848df79c378f
