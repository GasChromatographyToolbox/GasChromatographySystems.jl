### A Pluto.jl notebook ###
# v0.19.19

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
# Simple GCxGCxGC system
"""

# ╔═╡ f6c6443d-0e39-4e82-b2c4-49f2954e5b49
md"""
## Definition of system
"""

# ╔═╡ 719dc69d-2692-4d41-868e-14b664019d57
md"""
### Temperature programs
"""

# ╔═╡ c9941593-0311-473a-9c0c-f764d41ac7e7
GCxGC_TP = GasChromatographySystems.TemperatureProgram([0.0, 120.0, 3120.0, 300.0], [40.0, 40.0, 300.0, 300.0])

# ╔═╡ f9041774-42ff-4f86-993a-e3df618cc36e
begin
	T_off = 20.0
	offset_TP = GasChromatographySystems.TemperatureProgram([0.0, 120.0, 3120.0, 300.0], [40.0, 40.0, 300.0, 300.0].+T_off)
end

# ╔═╡ 3444976f-7687-4a63-9810-45f74406f6f6
begin # wobble (prog higher heating rate in 2nd dimension, short cooldown but not down to start temperature)
end

# ╔═╡ 7a47ad5b-23ad-4385-abc6-09b92cecf09a
begin # thermal gradient in 2nd dimension
end

# ╔═╡ 03190a64-09cf-4a9b-9b86-e2b50392481c
md"""
### 1. GCxGCxGC common temperature program
"""

# ╔═╡ d571712f-981e-4d1e-95ce-b4dbf7de50c6
begin
	g = SimpleDiGraph(4)
	add_edge!(g, 1, 2) # 1st-D GC
	add_edge!(g, 2, 3) # 2nd-D GC
	add_edge!(g, 3, 4) # 3rd-D GC
	# pressure points
	pp = Array{GasChromatographySystems.PressurePoint}(undef, nv(g))
	pp[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 120.0, 3120.0, 300.0], [159550.0, 159550.0, 265600.0, 265600.0]) # inlet 
	pp[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 120.0, 3120.0, 300.0], [NaN, NaN, NaN, NaN]) 
	pp[3] = GasChromatographySystems.PressurePoint("p₂", [0.0, 120.0, 3120.0, 300.0], [NaN, NaN, NaN, NaN]) 
	pp[4] = GasChromatographySystems.PressurePoint("p₃", [0.0, 120.0, 3120.0, 300.0], [eps(Float64), eps(Float64), eps(Float64), eps(Float64)]) # outlet
	# modules
	modules = Array{GasChromatographySystems.AbstractModule}(undef, ne(g))
	modules[1] = GasChromatographySystems.ModuleColumn("GC column 1", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", GCxGC_TP)
	modules[2] = GasChromatographySystems.ModuleColumn("GC column 2", 2.0, 0.25e-3, 0.25e-6, "SPB50", GCxGC_TP)
	modules[2] = GasChromatographySystems.ModuleColumn("GC column 3", 0.5, 0.1e-3, 0.1e-6, "Wax", GCxGC_TP)
	# system
	sys = GasChromatographySystems.update_system(GasChromatographySystems.System(g, pp, modules, GasChromatographySystems.Options(ng=true)))
end

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
	pp_offset[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 120.0, 3120.0, 300.0], [159550.0, 159550.0, 265600.0, 265600.0]) # inlet 
	pp_offset[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 120.0, 3120.0, 300.0], [NaN, NaN, NaN, NaN]) 
	pp_offset[3] = GasChromatographySystems.PressurePoint("p₃", [0.0, 120.0, 3120.0, 300.0], [eps(Float64), eps(Float64), eps(Float64), eps(Float64)])
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

# ╔═╡ 7783216f-6a73-44aa-98f7-e0ca628ba145
function plot_flow_over_time(sys)
	plotly()
	F_func = GasChromatographySystems.flow_functions(sys)
	com_timesteps = GasChromatographySystems.common_timesteps(sys)
	trange = 0:sum(com_timesteps)/100:sum(com_timesteps)
	p_flow = Plots.plot(xlabel="time in s", ylabel="flow in mL/min")
	for i=1:ne(sys.g)
		Plots.plot!(p_flow, trange, F_func[i].(trange)*60*1e6, label="F_$(sys.modules[i].name)")
	end
	return p_flow
end

# ╔═╡ b51877ef-a8a6-48d4-a521-5b1e01d2be40
function plot_pressure_over_time(sys)
	plotly()
	p_func = GasChromatographySystems.pressure_functions(sys)
	com_timesteps = GasChromatographySystems.common_timesteps(sys)
	trange = 0:sum(com_timesteps)/100:sum(com_timesteps)
	p_pres = Plots.plot(xlabel="time in s", ylabel="pressure in Pa(a)", legend=:topleft)
	for i=1:nv(sys.g)
		Plots.plot!(p_pres, trange, p_func[i].(trange), label="$(sys.pressurepoints[i].name)")
	end
	return p_pres
end

# ╔═╡ 06b9ba3f-4dff-4fba-afce-923b95fae9c3
GasChromatographySystems.plot_graph_with_flow(sys, 0)

# ╔═╡ af30bbcd-d6e4-4bd7-97f0-3243c514194a
plot_flow_over_time(sys)

# ╔═╡ 9c9aed82-3015-419b-ae6a-aa6f97a9d753
plot_pressure_over_time(sys)

# ╔═╡ a65c868c-fa26-4a7f-be94-2d86ee9518fc
md"""
## Selection of solutes
"""

# ╔═╡ b424f5f1-aca8-40e5-aa0d-32fc7c4d4ee6
db = filter([:CAS, :phi0] => (x, y) -> !ismissing(x) && y == 0.001, DataFrame(CSV.File(joinpath(dirname(pwd()), "data", "database_Blumberg2017.csv"), header=1, silencewarnings=true)))

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
sim_res = GasChromatographySystems.simulate_along_paths(sys, paths, db, selected_solutes, par)

# ╔═╡ 5373eaae-4372-4561-83cf-b62664bf28af
sim_res[2][1][1]

# ╔═╡ e57ab914-8db4-41ca-9950-1134e35e8c40
function peaklist_GCxGC(pl_1, pl_2)
	CAS1 = GasChromatographySimulator.CAS_identification(pl_1.Name).CAS
	pl_GCxGC = DataFrame(Name = pl_1.Name, CAS = CAS1, tR1 = pl_1.tR, τR1 = pl_1.τR)
	tR2 = Array{Float64}(undef, length(pl_GCxGC.Name))
	τR2 = Array{Float64}(undef, length(pl_GCxGC.Name))
	Cat = Array{Array{String}}(undef, length(pl_GCxGC.Name))
	#Cat_2 = Array{Union{Missing, String}}(undef, length(pl_GCxGC.Name))
	#Cat_3 = Array{Union{Missing, String}}(undef, length(pl_GCxGC.Name))
	#Cat_4 = Array{Union{Missing, String}}(undef, length(pl_GCxGC.Name))
	CAS2 = GasChromatographySimulator.CAS_identification(pl_2.Name).CAS
	for i=1:length(pl_GCxGC.Name)
		ii = findfirst(CAS1[i].==CAS2)
		tR2[i] = pl_2.tR[ii]
		τR2[i] = pl_2.τR[ii]
		i_db = findfirst(CAS1[i].==db.CAS)
		Cat[i] = unique([string(db.Cat_1[i_db]), string(db.Cat_2[i_db]), string(db.Cat_3[i_db]), string(db.Cat_4[i_db])])
	end
	pl_GCxGC[!, :tR2] = tR2
	pl_GCxGC[!, :τR2] = τR2
	pl_GCxGC[!, :Cat] = Cat
	return pl_GCxGC
end

# ╔═╡ 0bc7a04d-7534-4e0b-993c-a6ca65f8d363
pl = peaklist_GCxGC(sim_res[2][1][1], sim_res[2][1][2])

# ╔═╡ ab55966a-ce1c-440d-b7d3-d6f902332e10
["test", "text", "test", missing]

# ╔═╡ 0935e104-bc25-4a7d-a430-18474d6dc332
unique(skipmissing(["test", "text", "test", missing]))

# ╔═╡ 5b0e7e1f-bf59-4500-ada6-fef3edbda03f
unique([skipmissing("test"), skipmissing("text"), skipmissing("test"), skipmissing(missing)])

# ╔═╡ 024984d7-1771-4bbb-8d58-be43c7793710
Plots.scatter(pl.tR1, pl.tR2.-pl.tR1)

# ╔═╡ 1a92f83f-d67d-461a-92d7-8ffe2cae850d
in("FAME", pl.Cat[1])

# ╔═╡ 9bf9cf0e-0be9-4d3f-93d7-22da4f1a9d1a
pl.Cat[1]

# ╔═╡ 6e48e199-699c-464f-940e-6b78b67d80a4
pl_ = DataFrame()

# ╔═╡ 8c676a07-f179-4a71-ae16-f7791fff054f
pl[1,:]

# ╔═╡ 9d5ea53f-f1a6-4746-9000-caf6c3d8d9d0
for i=1:length(pl.Name)
	if in("FAME", pl.Cat[i])
		append!(pl_, pl[i,:])
	end
end

# ╔═╡ b19c9703-d57c-4ac8-8b50-7314e0218f2f
pl_

# ╔═╡ 4e523bf4-34cf-4a1c-b08b-8bcf8f4119d8


# ╔═╡ 92cdb34f-a791-4352-9855-29e7dafa6741
begin
	cats = unique([db.Cat_1; db.Cat_2; db.Cat_3; db.Cat_4])
	p_gcxgc = Plots.scatter(pl.tR1, pl.tR2.-pl.tR1, xlims=(0.0,3000.0), ylims=(0.0, 9.0))
	#cats_ = ["alkanes", "saturated FAME", "unsaturated FAME", ]
	#for i=1:length(cats) 
end

# ╔═╡ 1e102fe9-8c3f-4392-9eb7-c9d265876574
cats[1]

# ╔═╡ 9bed6c0b-7a48-4a30-ae15-4ce4784570ef
filter([:Cat] => x -> in(cats, x), pl)

# ╔═╡ df791608-5c03-42d6-b8ef-205882a18553
(pl.Cat_1 .== cats[1])

# ╔═╡ 7c5b1955-bb74-430c-a80b-284aecabd3a6
(skipmissing(pl.Cat_2) .== cats[1])

# ╔═╡ fba31c65-0465-458a-8fa5-7df1abdf179c
# plot 1d chromatogram measured at detector

# ╔═╡ 679140af-7598-41fe-aad1-90c50e503218
md"""
### 2. GCxGC offset temperature program in 2nd dimension
"""

# ╔═╡ 1f92dabf-8597-4ef8-9e05-628e922defbf
par_offset = GasChromatographySystems.graph_to_parameters(sys_offset, db, selected_solutes)

# ╔═╡ 6efdb3ef-4f6e-4073-bfa1-501732088e5c
paths_offset = GasChromatographySystems.all_paths(sys_offset.g, 1)[2]

# ╔═╡ 097f05ad-8723-4693-bfba-c350cd478f54
sim_res_offset = GasChromatographySystems.simulate_along_paths(sys_offset, paths_offset, db, selected_solutes, par_offset)

# ╔═╡ 417fa98b-38d1-434a-a44c-2948d71bec42
pl_offset = peaklist_GCxGC(sim_res_offset[2][1][1], sim_res_offset[2][1][2])

# ╔═╡ 28fbe11a-4132-41e1-8ea8-3b65749f41b1
p_offset = Plots.scatter(pl_offset.tR1, pl_offset.tR2.-pl_offset.tR1, xlims=(0.0,3000.0), ylims=(0.0, 9.0))

# ╔═╡ 67c99b55-6874-4b3f-adea-caf12504bf46
md"""
# End
"""

# ╔═╡ Cell order:
# ╠═e9d64f0e-95b7-11ed-2fee-ab579cd02ebc
# ╠═139bfe64-17ac-455d-84a2-09c9564156cd
# ╠═f6c6443d-0e39-4e82-b2c4-49f2954e5b49
# ╠═719dc69d-2692-4d41-868e-14b664019d57
# ╠═c9941593-0311-473a-9c0c-f764d41ac7e7
# ╠═f9041774-42ff-4f86-993a-e3df618cc36e
# ╠═3444976f-7687-4a63-9810-45f74406f6f6
# ╠═7a47ad5b-23ad-4385-abc6-09b92cecf09a
# ╠═03190a64-09cf-4a9b-9b86-e2b50392481c
# ╠═d571712f-981e-4d1e-95ce-b4dbf7de50c6
# ╠═4b379226-1646-4b13-b64a-7b55b66442f5
# ╠═7783216f-6a73-44aa-98f7-e0ca628ba145
# ╠═b51877ef-a8a6-48d4-a521-5b1e01d2be40
# ╠═06b9ba3f-4dff-4fba-afce-923b95fae9c3
# ╠═af30bbcd-d6e4-4bd7-97f0-3243c514194a
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
# ╠═5373eaae-4372-4561-83cf-b62664bf28af
# ╠═0bc7a04d-7534-4e0b-993c-a6ca65f8d363
# ╠═e57ab914-8db4-41ca-9950-1134e35e8c40
# ╠═ab55966a-ce1c-440d-b7d3-d6f902332e10
# ╠═0935e104-bc25-4a7d-a430-18474d6dc332
# ╠═5b0e7e1f-bf59-4500-ada6-fef3edbda03f
# ╠═024984d7-1771-4bbb-8d58-be43c7793710
# ╠═1a92f83f-d67d-461a-92d7-8ffe2cae850d
# ╠═9bf9cf0e-0be9-4d3f-93d7-22da4f1a9d1a
# ╠═1e102fe9-8c3f-4392-9eb7-c9d265876574
# ╠═6e48e199-699c-464f-940e-6b78b67d80a4
# ╠═8c676a07-f179-4a71-ae16-f7791fff054f
# ╠═9d5ea53f-f1a6-4746-9000-caf6c3d8d9d0
# ╠═b19c9703-d57c-4ac8-8b50-7314e0218f2f
# ╠═9bed6c0b-7a48-4a30-ae15-4ce4784570ef
# ╠═4e523bf4-34cf-4a1c-b08b-8bcf8f4119d8
# ╠═df791608-5c03-42d6-b8ef-205882a18553
# ╠═7c5b1955-bb74-430c-a80b-284aecabd3a6
# ╠═92cdb34f-a791-4352-9855-29e7dafa6741
# ╠═fba31c65-0465-458a-8fa5-7df1abdf179c
# ╠═679140af-7598-41fe-aad1-90c50e503218
# ╠═1f92dabf-8597-4ef8-9e05-628e922defbf
# ╠═6efdb3ef-4f6e-4073-bfa1-501732088e5c
# ╠═097f05ad-8723-4693-bfba-c350cd478f54
# ╠═417fa98b-38d1-434a-a44c-2948d71bec42
# ╠═28fbe11a-4132-41e1-8ea8-3b65749f41b1
# ╠═67c99b55-6874-4b3f-adea-caf12504bf46
