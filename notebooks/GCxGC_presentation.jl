### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ e9d64f0e-95b7-11ed-2fee-ab579cd02ebc
begin
	import Pkg
    # activate the shared project environment
	Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()

	using CSV, DataFrames
	using Plots, CairoMakie, GraphMakie
	using Graphs, NetworkLayout, Symbolics
	using GasChromatographySimulator
	using PlutoUI
	include(joinpath(dirname(pwd()), "src", "GasChromatographySystems.jl"))
	TableOfContents()
end

# ╔═╡ 139bfe64-17ac-455d-84a2-09c9564156cd
md"""
# Simple GCxGC system
"""

# ╔═╡ f6c6443d-0e39-4e82-b2c4-49f2954e5b49
md"""
## Definition of system
"""

# ╔═╡ 719dc69d-2692-4d41-868e-14b664019d57
md"""
### Temperature programs
"""

# ╔═╡ b6a40aab-fd24-471f-a0d5-72c06c7519ad
(300-40)/(3120/60)

# ╔═╡ 0235c359-e59e-4c20-905c-67e9973c1530
# rate 1°C/min [0.0, 120.0, 15600.0, 300.0], [40.0, 40.0, 300.0, 300.0]
# rate 3°C/min [0.0, 120.0, 5200.0, 300.0], [40.0, 40.0, 300.0, 300.0]
# rate 5°C/min [0.0, 120.0, 3120.0, 300.0], [40.0, 40.0, 300.0, 300.0]
# rate 8°C/min [0.0, 120.0, 1950.0, 300.0], [40.0, 40.0, 300.0, 300.0]
# rate 10°C/min [0.0, 120.0, 1560.0, 300.0], [40.0, 40.0, 300.0, 300.0]

# ╔═╡ 186e6bb5-14b4-4b5d-a057-1e3c649262cb
(300-40)/8*60

# ╔═╡ cec32771-d528-45ce-9bb4-59c6b3e95e15
tsteps = [0.0, 120.0, 3120.0, 300.0]

# ╔═╡ c9941593-0311-473a-9c0c-f764d41ac7e7
GCxGC_TP = GasChromatographySystems.TemperatureProgram(tsteps, [40.0, 40.0, 300.0, 300.0])

# ╔═╡ f9041774-42ff-4f86-993a-e3df618cc36e
begin
	T_off = 20.0
	offset_TP = GasChromatographySystems.TemperatureProgram(tsteps, [40.0, 40.0, 300.0, 300.0].+T_off)
end

# ╔═╡ 3444976f-7687-4a63-9810-45f74406f6f6
begin # wobble (prog higher heating rate in 2nd dimension, short cooldown but not down to start temperature)
end

# ╔═╡ 7a47ad5b-23ad-4385-abc6-09b92cecf09a
begin # thermal gradient in 2nd dimension
end

# ╔═╡ 03190a64-09cf-4a9b-9b86-e2b50392481c
md"""
### 1. GCxGC common temperature program
"""

# ╔═╡ a5dc0d72-03cd-4a54-9420-1b237171a279
# F ≈ 1mL/min [159550.0, 159550.0, 265600.0, 265600.0]
# F ≈ 1.5mL/min [190720.0, 190720.0, 317515.0, 317515.0]
# F ≈ 2mL/min [220235.0, 220235.0, 366655.0, 366655.0]

# ╔═╡ d571712f-981e-4d1e-95ce-b4dbf7de50c6
begin
	g = SimpleDiGraph(3)
	add_edge!(g, 1, 2) # 1st-D GC
	add_edge!(g, 2, 3) # 2nd-D GC
	#add_edge!(g2, 3, 4) # 2nd TL column inlet -> TL column -> Det
	# pressure points
	pp = Array{GasChromatographySystems.PressurePoint}(undef, nv(g))
	pp[1] = GasChromatographySystems.PressurePoint("p₁", tsteps, [159550.0, 159550.0, 265600.0, 265600.0]) # inlet 
	pp[2] = GasChromatographySystems.PressurePoint("p₂", tsteps, [NaN, NaN, NaN, NaN]) 
	pp[3] = GasChromatographySystems.PressurePoint("p₃", tsteps, [eps(Float64), eps(Float64), eps(Float64), eps(Float64)]) 
	#pp2[4] = GasChromatographySystems.PressurePoint("p₄", [0.0, 1800.0], [101300.0, 101300.0]) # outlet
	# modules
	modules = Array{GasChromatographySystems.AbstractModule}(undef, ne(g))
	modules[1] = GasChromatographySystems.ModuleColumn("GC column 1", 30.0, 0.25e-3, 0.25e-6, "Rxi5ms", GCxGC_TP)
	modules[2] = GasChromatographySystems.ModuleColumn("GC column 2", 2.0, 0.25e-3, 0.25e-6, "Rxi17SilMS", GCxGC_TP)
	#modules2[3] = GasChromatographySystems.ModuleColumn("GC column 3", 10.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	# system
	sys = GasChromatographySystems.update_system(GasChromatographySystems.System(g, pp, modules, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ af30bbcd-d6e4-4bd7-97f0-3243c514194a
GasChromatographySystems.plot_flow_over_time(sys)

# ╔═╡ f04130dc-e108-4e82-a2a2-a2c123b99df1
md"""
### 1a. Constant flow
"""

# ╔═╡ 5a250aed-458e-443e-849c-5e7580d260f1
begin
	g_ = SimpleDiGraph(3)
	add_edge!(g_, 1, 2) # 1st-D GC
	add_edge!(g_, 2, 3) # 2nd-D GC
	# pressure points
	pp_ = Array{GasChromatographySystems.PressurePoint}(undef, nv(g_))
	pp_[1] = GasChromatographySystems.PressurePoint("p₁", tsteps, [NaN, NaN, NaN, NaN]) # inlet 
	pp_[2] = GasChromatographySystems.PressurePoint("p₂", tsteps, [NaN, NaN, NaN, NaN]) # mid-point 
	pp_[3] = GasChromatographySystems.PressurePoint("p₃", tsteps, [eps(Float64), eps(Float64), eps(Float64), eps(Float64)]) # outlet
	# modules
	modules_ = Array{GasChromatographySystems.AbstractModule}(undef, ne(g_))
	modules_[1] = GasChromatographySystems.ModuleColumn("GC column 1", 30.0, 0.25e-3, 0.25e-6, "Rxi5ms", GCxGC_TP, 1.0/60e6)
	modules_[2] = GasChromatographySystems.ModuleColumn("GC column 2", 2.0, 0.25e-3, 0.25e-6, "Rxi17SilMS", GCxGC_TP, NaN) # defining here also the flow results in an error
	# system
	sys_1a = GasChromatographySystems.update_system(GasChromatographySystems.System(g_, pp_, modules_, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ cfb715f1-8ccd-4035-ab6a-cb2746e0e6bf
GasChromatographySystems.plot_graph_with_flow(sys_1a, 0; lay=SquareGrid(cols=3), node_size=90, nlabels_fontsize=20, elabels_fontsize=20, elabels_distance=28)

# ╔═╡ 7ee6e1c0-8102-49dd-a4d9-0a84d3461e06
GasChromatographySystems.plot_flow_over_time(sys_1a)

# ╔═╡ 80df8382-bb80-4937-804b-58c861b76c10
GasChromatographySystems.plot_pressure_over_time(sys_1a)

# ╔═╡ 6cfda8f9-f918-4451-97c3-084eb1d9d339
md"""
### 1b. reverse phases
"""

# ╔═╡ a8971018-3f53-4276-b004-9245985ca490
begin
	g_rev = SimpleDiGraph(3)
	add_edge!(g_rev, 1, 2) # 1st-D GC
	add_edge!(g_rev, 2, 3) # 2nd-D GC
	#add_edge!(g2, 3, 4) # 2nd TL column inlet -> TL column -> Det
	# pressure points
	pp_rev = Array{GasChromatographySystems.PressurePoint}(undef, nv(g_rev))
	pp_rev[1] = GasChromatographySystems.PressurePoint("p₁", tsteps, [190720.0, 190720.0, 317515.0, 317515.0]) # inlet 
	pp_rev[2] = GasChromatographySystems.PressurePoint("p₂", tsteps, [NaN, NaN, NaN, NaN]) 
	pp_rev[3] = GasChromatographySystems.PressurePoint("p₃", tsteps, [eps(Float64), eps(Float64), eps(Float64), eps(Float64)]) 
	#pp2[4] = GasChromatographySystems.PressurePoint("p₄", [0.0, 1800.0], [101300.0, 101300.0]) # outlet
	# modules
	modules_rev = Array{GasChromatographySystems.AbstractModule}(undef, ne(g_rev))
	modules_rev[1] = GasChromatographySystems.ModuleColumn("GC column 1", 30.0, 0.25e-3, 0.25e-6, "Rxi17SilMS", GCxGC_TP)
	modules_rev[2] = GasChromatographySystems.ModuleColumn("GC column 2", 2.0, 0.25e-3, 0.25e-6, "Rxi5ms", GCxGC_TP)
	#modules2[3] = GasChromatographySystems.ModuleColumn("GC column 3", 10.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	# system
	sys_rev = GasChromatographySystems.update_system(GasChromatographySystems.System(g_rev, pp_rev, modules_rev, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ aa2e85ad-4b19-4f90-9650-1dadb47306ae
md"""
### 1c. Flow modulator
"""

# ╔═╡ 9509bc4d-2906-4a4c-919b-99412f1e0985
begin
	g_1c = SimpleDiGraph(4)
	add_edge!(g_1c, 1, 3) # 1st-D GC
	add_edge!(g_1c, 2, 3) # Flow modulator
	add_edge!(g_1c, 3, 4) # 2nd-D GC
	# pressure points
	pp_1c = Array{GasChromatographySystems.PressurePoint}(undef, nv(g_1c))
	pp_1c[1] = GasChromatographySystems.PressurePoint("p₁", tsteps, [NaN, NaN, NaN, NaN]) # inlet 
	pp_1c[2] = GasChromatographySystems.PressurePoint("p₂", tsteps, [NaN, NaN, NaN, NaN]) # flow modulator inlet
	pp_1c[3] = GasChromatographySystems.PressurePoint("p₃", tsteps, [NaN, NaN, NaN, NaN]) # mid-point 
	pp_1c[4] = GasChromatographySystems.PressurePoint("p₄", tsteps, [eps(Float64), eps(Float64), eps(Float64), eps(Float64)])# outlet
	# modules
	modules_1c = Array{GasChromatographySystems.AbstractModule}(undef, ne(g_1c))
	modules_1c[1] = GasChromatographySystems.ModuleColumn("GC column 1", 30.0, 0.25e-3, 0.25e-6, "Rxi5ms", GCxGC_TP, 1.0/60e6)
	modules_1c[2] = GasChromatographySystems.ModuleColumn("FM", 0.1, 0.5e-3, 0.0, "", GCxGC_TP, NaN) 
	modules_1c[3] = GasChromatographySystems.ModuleColumn("GC column 2", 2.0, 0.25e-3, 0.25e-6, "Rxi17SilMS", GCxGC_TP, 2.0/60e6) 
	# system
	sys_1c = GasChromatographySystems.update_system(GasChromatographySystems.System(g_1c, pp_1c, modules_1c, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ 1b60cf45-d8db-4e28-b287-2ac6ff2cd1ea
GasChromatographySystems.substitute_unknown_flows(sys_1c)

# ╔═╡ a72f4c81-0a48-4343-8840-acb0514fc2ee
GasChromatographySystems.plot_graph_with_flow(sys_1c, 0; node_size=90, nlabels_fontsize=20, elabels_fontsize=20, elabels_distance=28)

# ╔═╡ 22674e61-db92-4181-ad84-c36955ac77a7
GasChromatographySystems.plot_flow_over_time(sys_1c)

# ╔═╡ 56ac1bb7-7e33-45ad-965c-d9b6e8d29bf3
GasChromatographySystems.plot_pressure_over_time(sys_1c)

# ╔═╡ e4fa29f4-fc83-4aa4-a548-f2683ad5d335
md"""
### 2. GCxGC offset temperature program in 2nd dimension
"""

# ╔═╡ 38a4081d-49ed-4688-84ba-88218f535911
begin
	g_offset = SimpleDiGraph(3)
	add_edge!(g_offset, 1, 2) # 1st-D GC
	add_edge!(g_offset, 2, 3) # 2nd-D GC
	# pressure points
	pp_offset = Array{GasChromatographySystems.PressurePoint}(undef, nv(g_offset))
	pp_offset[1] = GasChromatographySystems.PressurePoint("p₁", tsteps, [159550.0, 159550.0, 265600.0, 265600.0]) # inlet 
	pp_offset[2] = GasChromatographySystems.PressurePoint("p₂", tsteps, [NaN, NaN, NaN, NaN]) 
	pp_offset[3] = GasChromatographySystems.PressurePoint("p₃", tsteps, [eps(Float64), eps(Float64), eps(Float64), eps(Float64)])
	# modules
	modules_offset = Array{GasChromatographySystems.AbstractModule}(undef, ne(g_offset))
	modules_offset[1] = GasChromatographySystems.ModuleColumn("GC column 1", 30.0, 0.25e-3, 0.25e-6, "Rxi5ms", GCxGC_TP)
	modules_offset[2] = GasChromatographySystems.ModuleColumn("GC column 2", 2.0, 0.25e-3, 0.25e-6, "Rxi17SilMS", offset_TP)
	# system
	sys_offset = GasChromatographySystems.update_system(GasChromatographySystems.System(g_offset, pp_offset, modules_offset, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ 4b379226-1646-4b13-b64a-7b55b66442f5
md"""
## Pressure and flow calculation
"""

# ╔═╡ 06b9ba3f-4dff-4fba-afce-923b95fae9c3
GasChromatographySystems.plot_graph_with_flow(sys, 0; lay=SquareGrid(cols=3), node_size=90, nlabels_fontsize=20, elabels_fontsize=20, elabels_distance=28)

# ╔═╡ 9c9aed82-3015-419b-ae6a-aa6f97a9d753
GasChromatographySystems.plot_pressure_over_time(sys)

# ╔═╡ a65c868c-fa26-4a7f-be94-2d86ee9518fc
md"""
## Selection of solutes
"""

# ╔═╡ b424f5f1-aca8-40e5-aa0d-32fc7c4d4ee6
#db = filter([:CAS, :phi0] => (x, y) -> !ismissing(x) && y == 0.001, DataFrame(CSV.File(joinpath(dirname(pwd()), "data", "database_Rxi17SilMS_Rxi5SilMS_.csv"), header=1, silencewarnings=true)))
db = filter([:CAS] => x -> !ismissing(x), DataFrame(CSV.File(joinpath(dirname(pwd()), "data", "database_Rxi17SilMS_Rxi5SilMS_.csv"), header=1, silencewarnings=true)))

# ╔═╡ 0b0c0577-ec2e-450f-9f19-4a29da213efe
unique(GasChromatographySystems.all_stationary_phases(sys))

# ╔═╡ 3cd3ea64-a45c-49fb-90bb-c807c47eefc0
selected_solutes = GasChromatographySystems.common_solutes(db, sys).Name

# ╔═╡ e935ae05-c654-4365-a569-fed418f35fe1
md"""
## Graph to Parameters
"""

# ╔═╡ d060c9a3-d760-472f-b52b-3c481b39a85c
par = GasChromatographySystems.graph_to_parameters(sys, db, selected_solutes)

# ╔═╡ 79bc683a-8d51-4d81-a27c-e42663720917
md"""
## Simulation
"""

# ╔═╡ 1436c866-9124-4ab4-8ab7-0b04f4de896c
md"""
### 1. GCxGC common temperature program
"""

# ╔═╡ 457df2b2-8bfa-4e04-b4f2-ae1ca8902a4e
paths = GasChromatographySystems.all_paths(sys.g, 1)[2]

# ╔═╡ f106b5fd-e17a-4876-8562-e82670a1ab25
sim_res = GasChromatographySystems.simulate_along_paths(sys, paths, par)

# ╔═╡ d8c30c13-7209-43ea-93c8-ac0310b0a8e5
begin
	par_ = GasChromatographySystems.graph_to_parameters(sys, db, selected_solutes)
	paths_ = GasChromatographySystems.all_paths(sys.g, 1)[2]
	sim_res_ = GasChromatographySystems.simulate_along_paths(sys, paths_, par_)
end

# ╔═╡ 7dac264f-2e91-44bd-b96d-207f072cba9c
sim_res_foc = GasChromatographySystems.simulate_along_paths(sys, paths, par; refocus=[true, true])

# ╔═╡ 5373eaae-4372-4561-83cf-b62664bf28af
sim_res[2][1][1]

# ╔═╡ c5888be4-c641-4871-81a5-11d01ace12ed
#savefig(p_GCxGC_cats, "chrom_GCxGC_10.svg")

# ╔═╡ 987efbb5-8a48-4c83-9a20-50e5fbd0a16d
function plot_GCxGC_contour_(pl_GCxGC; split_ratio=ones(length(pl_GCxGC.Name)), mod_τ1=1.0)
	t¹ = 0.9*minimum(pl_GCxGC.tR1):(1.1*maximum(pl_GCxGC.tR1)-0.9*minimum(pl_GCxGC.tR1))/1000:1.1*maximum(pl_GCxGC.tR1)
	t² = 0.9*minimum(pl_GCxGC.tR2.-pl_GCxGC.tR1):(1.1*maximum(pl_GCxGC.tR2.-pl_GCxGC.tR1)-0.9*minimum(pl_GCxGC.tR2.-pl_GCxGC.tR1))/1000:1.1*maximum(pl_GCxGC.tR2.-pl_GCxGC.tR1)
	chrom2D = Array{Float64}(undef, length(t¹), length(t²), length(pl_GCxGC.Name))
	for j=1:length(t²)
		for i=1:length(t¹)
			for k=1:length(pl_GCxGC.Name)
				chrom2D[i,j,k] = split_ratio[k]*exp(-((t¹[i]-pl_GCxGC.tR1[k])^2/(2*(mod_τ1*pl_GCxGC.τR1[k])^2)+(t²[j]-(pl_GCxGC.tR2[k]-pl_GCxGC.tR1[k]))^2/(2*pl_GCxGC.τR2[k]^2)))
			end
		end
	end
	#p_2D = Plots.contourf(t¹, t², sum(chrom2D, dims=3)[:,:,1]', fill=true, legend=false, xlabel="tR1 in s", ylabel="tR2 in s", c=:turbo, levels=40)
	return t¹, t², sum(chrom2D, dims=3)[:,:,1]'
end

# ╔═╡ fba31c65-0465-458a-8fa5-7df1abdf179c
# plot 1d chromatogram measured at detector

# ╔═╡ 88bf3e84-6772-453d-8bee-2ccbbcf7cb34
begin
	p_chrom = GasChromatographySimulator.plot_chromatogram(sim_res[2][1][1], (0, sum(GasChromatographySystems.common_timesteps(sys))))[1]
	GasChromatographySimulator.plot_chromatogram!(p_chrom, sim_res[2][1][2], (0, sum(GasChromatographySystems.common_timesteps(sys))); mirror=true)
	p_chrom
end

# ╔═╡ 0cc7f988-bf98-4e45-93d5-7e3dea48c4ed
begin
	p_chrom_foc = GasChromatographySimulator.plot_chromatogram(sim_res_foc[2][1][1], (0, sum(GasChromatographySystems.common_timesteps(sys))))[1]
	GasChromatographySimulator.plot_chromatogram!(p_chrom_foc, sim_res_foc[2][1][2], (0, sum(GasChromatographySystems.common_timesteps(sys))); mirror=true)
	p_chrom_foc
end

# ╔═╡ a3492676-6d52-49ef-9623-ffb3e5bda682
md"""
### 1a. Constant flow
"""

# ╔═╡ 81f7d27c-4048-4c06-acc7-8c01ab580c58
par_1a = GasChromatographySystems.graph_to_parameters(sys_1a, db, selected_solutes)

# ╔═╡ cce329b8-234a-45a9-aa6d-b4a57a0f4627
paths_1a = GasChromatographySystems.all_paths(sys_1a.g, 1)[2]

# ╔═╡ dfc605b8-bc23-4456-8dbe-6cb83f340785
sim_res_1a = GasChromatographySystems.simulate_along_paths(sys_1a, paths_1a, par_1a; refocus=[true, true])

# ╔═╡ e7956916-dccb-46c1-a30b-ef77650e54ce
md"""
### 1b. reverse phases
"""

# ╔═╡ 2c45c22a-842f-4d49-8b0b-d6f9546932c9
par_rev = GasChromatographySystems.graph_to_parameters(sys_rev, db, selected_solutes)

# ╔═╡ eccd7341-eae1-4a37-9cfd-f1817857e572
paths_rev = GasChromatographySystems.all_paths(sys_rev.g, 1)[2]

# ╔═╡ aa0d4fda-3f1b-468c-9568-acadf61b8b95
sim_res_rev = GasChromatographySystems.simulate_along_paths(sys_rev, paths_rev, par_rev)

# ╔═╡ 834739e9-fa7f-43f0-beac-ae5790575255
#savefig(p_GCxGC_cats_rev, "chrom_GCxGC_8_rev.svg")

# ╔═╡ 28316a84-710b-443e-a136-43fcfd821afa
md"""
### 1c. Flow modulator
"""

# ╔═╡ fa011932-61e2-4377-bbf9-eb0bcd8c1b22
par_1c = GasChromatographySystems.graph_to_parameters(sys_1c, db, selected_solutes; interp=true, dt=10) # using the calculated functions for the pressures
# -> in this case the simulation takes to long
# interpolation for an increased stepsize seems workable

# ╔═╡ 1783d8de-90df-434d-b71d-43c39a0d3579
paths_1c = GasChromatographySystems.all_paths(sys_1c.g, 1)[2]

# ╔═╡ fcba3e01-9464-4f95-8144-d8d1671b265e
sim_res_1c = GasChromatographySystems.simulate_along_paths(sys_1c, paths_1c, par_1c)

# ╔═╡ 679140af-7598-41fe-aad1-90c50e503218
md"""
### 2. GCxGC offset temperature program in 2nd dimension
"""

# ╔═╡ 1f92dabf-8597-4ef8-9e05-628e922defbf
par_offset = GasChromatographySystems.graph_to_parameters(sys_offset, db, selected_solutes)

# ╔═╡ 6efdb3ef-4f6e-4073-bfa1-501732088e5c
paths_offset = GasChromatographySystems.all_paths(sys_offset.g, 1)[2]

# ╔═╡ 097f05ad-8723-4693-bfba-c350cd478f54
sim_res_offset = GasChromatographySystems.simulate_along_paths(sys_offset, paths_offset, par_offset)

# ╔═╡ 8d14ae41-a062-4d89-8e5a-e78b62cc0564
#savefig(p_GCxGC_cats_offset, "chrom_GCxGC_8_offset.svg")

# ╔═╡ 67c99b55-6874-4b3f-adea-caf12504bf46
md"""
# End
"""

# ╔═╡ 127069e6-6da1-4742-88b4-a0d398a88b24
function plot_GCxGC(pl_GCxGC, sys; categories = String[])
	t1end = sum(GasChromatographySystems.common_timesteps(sys))
	t2end = ceil(maximum(pl_GCxGC.tR2.-pl_GCxGC.tR1))
	p_gcxgc = Plots.scatter(pl_GCxGC.tR1, pl_GCxGC.tR2.-pl_GCxGC.tR1, xlims=(0.0,t1end), ylims=(0.0, t2end), label="", xlabel=string(sys.modules[1].stationary_phase, ", t¹ in s"), ylabel=string(sys.modules[2].stationary_phase, ", t² in s"))
	if !isempty(categories)
		for i=1:length(categories) 
			pl_f = filter([:Cat] => x -> occursin(categories[i], x), pl_GCxGC)
			Plots.scatter!(p_gcxgc, pl_f.tR1, pl_f.tR2.-pl_f.tR1, label=categories[i], legend=:bottomright)
		end
	end
	return p_gcxgc
end

# ╔═╡ 417fa98b-38d1-434a-a44c-2948d71bec42
function peaklist_GCxGC_with_categories(pl_D1, pl_D2, db)
	pl = GasChromatographySystems.peaklist_GCxGC(pl_D1, pl_D2)
	cat = Array{String}(undef, length(pl.Name))
	for i=1:length(pl.Name)
		ii = findfirst(pl.Name[i].==db.Name)
		cat[i] = ""
		for j=1:length(findall(occursin.("Cat", names(db))))
			if !ismissing(db[ii, Symbol("Cat_$(j)")])
				if j==1
					cat[i] = cat[i]*db[ii, Symbol("Cat_$(j)")]
				else
					cat[i] = cat[i]*", "*db[ii, Symbol("Cat_$(j)")]
				end
			end
		end
	end
	pl[!, :Cat] = cat
	return pl
end

# ╔═╡ 0bc7a04d-7534-4e0b-993c-a6ca65f8d363
pl = peaklist_GCxGC_with_categories(sim_res[2][1][1], sim_res[2][1][2], db)

# ╔═╡ 3986a24a-689c-48b8-8e78-ec5d3842d0d7
plot_GCxGC(pl, sys)

# ╔═╡ 938d8677-7898-43c8-b07c-6e9d7f3c7e2e
begin
	gr()
	p_GCxGC_cats = plot_GCxGC(pl, sys; categories = ["saturated FAME", "unsaturated FAME", "aromatic", "PCB", "alkanes", "alcohols", "alkanones", "phenones"])
	p_GCxGC_cats
end

# ╔═╡ 5429291b-8efd-4ebb-9d4d-1be143073482
pl_foc =peaklist_GCxGC_with_categories(sim_res_foc[2][1][1], sim_res_foc[2][1][2], db)

# ╔═╡ 939e3a45-03cd-4e07-a672-10c32cf2637e
t¹, t², chrom = plot_GCxGC_contour_(pl_foc)

# ╔═╡ 7100eeed-5604-4d81-855e-b08cd40135a6
Plots.contour(t¹, t², chrom, legend=false, xlabel="tR1 in s", ylabel="tR2 in s", c=:turbo, levels=100, xlims=(0.0, 1700.0), ylims=(0.0, 5.0), grid=false)

# ╔═╡ 00fbcd7c-71a1-43f7-b358-abfde851a237
begin
	plotly()
	p_GCxGC = Plots.heatmap(t¹, t², chrom.^(1//3), legend=false, xlabel="tR1 in s", ylabel="tR2 in s", c=:jet1, dpi=500)#, xlims=(0.0, 1700.0), ylims=(0.0, 5.0), grid=false)
	#gui()
end

# ╔═╡ 6a1e1c53-93da-4b2d-a65d-dd13fd86c174
pl_1a = peaklist_GCxGC_with_categories(sim_res_1a[2][1][1], sim_res_1a[2][1][2], db)

# ╔═╡ 9976debc-cbbf-462a-a1f6-b4f9aa468cd2
plot_GCxGC(pl_1a, sys_1a; categories = ["saturated FAME", "unsaturated FAME", "aromatic", "PCB", "alkanes", "alcohols", "alkanones", "phenones"])

# ╔═╡ 1143dfde-bc45-4c2f-b90c-be6afe7bb35a
t¹_1a, t²_1a, chrom_1a = plot_GCxGC_contour_(pl_1a; mod_τ1=1.5)

# ╔═╡ 5771318a-215a-440d-bc13-d71e042a5c03
begin
	plotly()
	#gr()
	p_GCxGC_1a = Plots.heatmap(t¹_1a, t²_1a, chrom_1a.^(1//3), legend=false, xlabel="tR1 in s", ylabel="tR2 in s", c=:jet1, dpi=500)#, xlims=(0.0, 1700.0), ylims=(0.0, 5.0), grid=false)
	#gui()
end

# ╔═╡ 397a913b-fa6c-4eba-b8b7-7571256a001b
pl_rev = peaklist_GCxGC_with_categories(sim_res_rev[2][1][1], sim_res_rev[2][1][2], db)

# ╔═╡ c3f16558-43b8-4822-98af-daf8e06cd1df
p_GCxGC_cats_rev = plot_GCxGC(pl_rev, sys_rev; categories = ["saturated FAME", "unsaturated FAME", "aromatic", "PCB", "alkanes", "alcohols", "alkanones", "phenones"])

# ╔═╡ 7bcca943-1b5a-4eb5-a419-30a1959c53aa
pl_1c = peaklist_GCxGC_with_categories(sim_res_1c[2][1][1], sim_res_1c[2][1][2], db)

# ╔═╡ ee516d3e-d60b-49c6-aa2f-eb7cb3a9acd1
plot_GCxGC(pl_1c, sys_1c; categories = ["saturated FAME", "unsaturated FAME", "aromatic", "PCB", "alkanes", "alcohols", "alkanones", "phenones"])

# ╔═╡ ec9dcf65-93ca-4299-a5a7-578aca60e906
t¹_1c, t²_1c, chrom_1c = plot_GCxGC_contour_(pl_1c)

# ╔═╡ a6b8b92d-e7a7-4b1a-9a12-1db83b530e13
begin
	plotly()
	p_GCxGC_1c = Plots.heatmap(t¹_1c, t²_1c, chrom_1c.^(1//3), legend=false, xlabel="tR1 in s", ylabel="tR2 in s", c=:jet1, dpi=500)#, xlims=(0.0, 1700.0), ylims=(0.0, 5.0), grid=false)
	#gui()
end

# ╔═╡ a0d12beb-a1e8-4f32-902d-233308f0a783
pl_offset = peaklist_GCxGC_with_categories(sim_res_offset[2][1][1], sim_res_offset[2][1][2], db)

# ╔═╡ 28fbe11a-4132-41e1-8ea8-3b65749f41b1
p_GCxGC_cats_offset = plot_GCxGC(pl_offset, sys_offset; categories = ["saturated FAME", "unsaturated FAME", "aromatic", "PCB", "alkanes", "alcohols", "alkanones", "phenones"])

# ╔═╡ Cell order:
# ╠═e9d64f0e-95b7-11ed-2fee-ab579cd02ebc
# ╠═139bfe64-17ac-455d-84a2-09c9564156cd
# ╠═f6c6443d-0e39-4e82-b2c4-49f2954e5b49
# ╠═719dc69d-2692-4d41-868e-14b664019d57
# ╠═b6a40aab-fd24-471f-a0d5-72c06c7519ad
# ╠═0235c359-e59e-4c20-905c-67e9973c1530
# ╠═186e6bb5-14b4-4b5d-a057-1e3c649262cb
# ╠═cec32771-d528-45ce-9bb4-59c6b3e95e15
# ╠═c9941593-0311-473a-9c0c-f764d41ac7e7
# ╠═f9041774-42ff-4f86-993a-e3df618cc36e
# ╠═3444976f-7687-4a63-9810-45f74406f6f6
# ╠═7a47ad5b-23ad-4385-abc6-09b92cecf09a
# ╠═03190a64-09cf-4a9b-9b86-e2b50392481c
# ╠═a5dc0d72-03cd-4a54-9420-1b237171a279
# ╠═af30bbcd-d6e4-4bd7-97f0-3243c514194a
# ╠═d571712f-981e-4d1e-95ce-b4dbf7de50c6
# ╠═f04130dc-e108-4e82-a2a2-a2c123b99df1
# ╠═5a250aed-458e-443e-849c-5e7580d260f1
# ╠═cfb715f1-8ccd-4035-ab6a-cb2746e0e6bf
# ╠═7ee6e1c0-8102-49dd-a4d9-0a84d3461e06
# ╠═80df8382-bb80-4937-804b-58c861b76c10
# ╠═6cfda8f9-f918-4451-97c3-084eb1d9d339
# ╠═a8971018-3f53-4276-b004-9245985ca490
# ╠═aa2e85ad-4b19-4f90-9650-1dadb47306ae
# ╠═9509bc4d-2906-4a4c-919b-99412f1e0985
# ╠═1b60cf45-d8db-4e28-b287-2ac6ff2cd1ea
# ╠═a72f4c81-0a48-4343-8840-acb0514fc2ee
# ╠═22674e61-db92-4181-ad84-c36955ac77a7
# ╠═56ac1bb7-7e33-45ad-965c-d9b6e8d29bf3
# ╠═e4fa29f4-fc83-4aa4-a548-f2683ad5d335
# ╠═38a4081d-49ed-4688-84ba-88218f535911
# ╠═4b379226-1646-4b13-b64a-7b55b66442f5
# ╠═06b9ba3f-4dff-4fba-afce-923b95fae9c3
# ╠═9c9aed82-3015-419b-ae6a-aa6f97a9d753
# ╠═a65c868c-fa26-4a7f-be94-2d86ee9518fc
# ╠═b424f5f1-aca8-40e5-aa0d-32fc7c4d4ee6
# ╠═0b0c0577-ec2e-450f-9f19-4a29da213efe
# ╠═3cd3ea64-a45c-49fb-90bb-c807c47eefc0
# ╠═e935ae05-c654-4365-a569-fed418f35fe1
# ╠═d060c9a3-d760-472f-b52b-3c481b39a85c
# ╠═79bc683a-8d51-4d81-a27c-e42663720917
# ╟─1436c866-9124-4ab4-8ab7-0b04f4de896c
# ╠═457df2b2-8bfa-4e04-b4f2-ae1ca8902a4e
# ╠═f106b5fd-e17a-4876-8562-e82670a1ab25
# ╠═d8c30c13-7209-43ea-93c8-ac0310b0a8e5
# ╠═7dac264f-2e91-44bd-b96d-207f072cba9c
# ╠═5373eaae-4372-4561-83cf-b62664bf28af
# ╠═0bc7a04d-7534-4e0b-993c-a6ca65f8d363
# ╠═5429291b-8efd-4ebb-9d4d-1be143073482
# ╠═3986a24a-689c-48b8-8e78-ec5d3842d0d7
# ╠═938d8677-7898-43c8-b07c-6e9d7f3c7e2e
# ╠═c5888be4-c641-4871-81a5-11d01ace12ed
# ╠═987efbb5-8a48-4c83-9a20-50e5fbd0a16d
# ╠═939e3a45-03cd-4e07-a672-10c32cf2637e
# ╠═7100eeed-5604-4d81-855e-b08cd40135a6
# ╠═00fbcd7c-71a1-43f7-b358-abfde851a237
# ╠═fba31c65-0465-458a-8fa5-7df1abdf179c
# ╠═88bf3e84-6772-453d-8bee-2ccbbcf7cb34
# ╠═0cc7f988-bf98-4e45-93d5-7e3dea48c4ed
# ╠═a3492676-6d52-49ef-9623-ffb3e5bda682
# ╠═81f7d27c-4048-4c06-acc7-8c01ab580c58
# ╠═cce329b8-234a-45a9-aa6d-b4a57a0f4627
# ╠═dfc605b8-bc23-4456-8dbe-6cb83f340785
# ╠═6a1e1c53-93da-4b2d-a65d-dd13fd86c174
# ╠═9976debc-cbbf-462a-a1f6-b4f9aa468cd2
# ╠═1143dfde-bc45-4c2f-b90c-be6afe7bb35a
# ╠═5771318a-215a-440d-bc13-d71e042a5c03
# ╠═e7956916-dccb-46c1-a30b-ef77650e54ce
# ╠═2c45c22a-842f-4d49-8b0b-d6f9546932c9
# ╠═eccd7341-eae1-4a37-9cfd-f1817857e572
# ╠═aa0d4fda-3f1b-468c-9568-acadf61b8b95
# ╠═397a913b-fa6c-4eba-b8b7-7571256a001b
# ╠═c3f16558-43b8-4822-98af-daf8e06cd1df
# ╠═834739e9-fa7f-43f0-beac-ae5790575255
# ╠═28316a84-710b-443e-a136-43fcfd821afa
# ╠═fa011932-61e2-4377-bbf9-eb0bcd8c1b22
# ╠═1783d8de-90df-434d-b71d-43c39a0d3579
# ╠═fcba3e01-9464-4f95-8144-d8d1671b265e
# ╠═7bcca943-1b5a-4eb5-a419-30a1959c53aa
# ╠═ee516d3e-d60b-49c6-aa2f-eb7cb3a9acd1
# ╠═ec9dcf65-93ca-4299-a5a7-578aca60e906
# ╠═a6b8b92d-e7a7-4b1a-9a12-1db83b530e13
# ╠═679140af-7598-41fe-aad1-90c50e503218
# ╠═1f92dabf-8597-4ef8-9e05-628e922defbf
# ╠═6efdb3ef-4f6e-4073-bfa1-501732088e5c
# ╠═097f05ad-8723-4693-bfba-c350cd478f54
# ╠═a0d12beb-a1e8-4f32-902d-233308f0a783
# ╠═28fbe11a-4132-41e1-8ea8-3b65749f41b1
# ╠═8d14ae41-a062-4d89-8e5a-e78b62cc0564
# ╠═67c99b55-6874-4b3f-adea-caf12504bf46
# ╠═127069e6-6da1-4742-88b4-a0d398a88b24
# ╠═417fa98b-38d1-434a-a44c-2948d71bec42
