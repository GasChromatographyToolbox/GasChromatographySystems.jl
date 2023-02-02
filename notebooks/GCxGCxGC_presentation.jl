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
TP = GasChromatographySystems.TemperatureProgram([0.0, 120.0, 3120.0, 300.0], [40.0, 40.0, 300.0, 300.0])

# ╔═╡ 9bd5cee4-fd2d-44e9-8818-6c4733565c02
(300-40)/(3120/60)

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
	pp[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 120.0, 3120.0, 300.0], [259550.0, 259550.0, 365600.0, 365600.0]) # inlet 
	pp[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 120.0, 3120.0, 300.0], [NaN, NaN, NaN, NaN]) 
	pp[3] = GasChromatographySystems.PressurePoint("p₂", [0.0, 120.0, 3120.0, 300.0], [NaN, NaN, NaN, NaN]) 
	pp[4] = GasChromatographySystems.PressurePoint("p₃", [0.0, 120.0, 3120.0, 300.0], [eps(Float64), eps(Float64), eps(Float64), eps(Float64)]) # outlet
	# modules
	modules = Array{GasChromatographySystems.AbstractModule}(undef, ne(g))
	modules[1] = GasChromatographySystems.ModuleColumn("GC column 1", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", TP)
	modules[2] = GasChromatographySystems.ModuleColumn("GC column 2", 2.0, 0.25e-3, 0.25e-6, "SPB50", TP)
	modules[3] = GasChromatographySystems.ModuleColumn("GC column 3", 0.5, 0.25e-3, 0.25e-6, "Wax", TP)
	# system
	sys = GasChromatographySystems.update_system(GasChromatographySystems.System(g, pp, modules, GasChromatographySystems.Options(ng=true)))
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

# ╔═╡ 9d02024e-10dd-4181-b280-6dcba25b21d4
GasChromatographySystems.plot_graph_with_flow(sys, 750; lay=SquareGrid(cols=4), node_size=90, nlabels_fontsize=20, elabels_fontsize=20, elabels_distance=28)

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
function peaklist_GCxGCxGC(pl_1, pl_2, pl_3)
	CAS1 = GasChromatographySimulator.CAS_identification(pl_1.Name).CAS
	pl_GCxGCxGC = DataFrame(Name = pl_1.Name, CAS = CAS1, tR1 = pl_1.tR, τR1 = pl_1.τR)
	tR2 = Array{Float64}(undef, length(pl_GCxGCxGC.Name))
	τR2 = Array{Float64}(undef, length(pl_GCxGCxGC.Name))
	tR3 = Array{Float64}(undef, length(pl_GCxGCxGC.Name))
	τR3 = Array{Float64}(undef, length(pl_GCxGCxGC.Name))
	Cat = Array{Array{String}}(undef, length(pl_GCxGCxGC.Name))
	#Cat_2 = Array{Union{Missing, String}}(undef, length(pl_GCxGC.Name))
	#Cat_3 = Array{Union{Missing, String}}(undef, length(pl_GCxGC.Name))
	#Cat_4 = Array{Union{Missing, String}}(undef, length(pl_GCxGC.Name))
	CAS2 = GasChromatographySimulator.CAS_identification(pl_2.Name).CAS
	CAS3 = GasChromatographySimulator.CAS_identification(pl_3.Name).CAS
	for i=1:length(pl_GCxGCxGC.Name)
		i2 = findfirst(CAS1[i].==CAS2)
		i3 = findfirst(CAS1[i].==CAS3)
		tR2[i] = pl_2.tR[i2]
		τR2[i] = pl_2.τR[i2]
		tR3[i] = pl_3.tR[i3]
		τR3[i] = pl_3.τR[i3]
		i_db = findfirst(CAS1[i].==db.CAS)
		Cat[i] = unique([string(db.Cat_1[i_db]), string(db.Cat_2[i_db]), string(db.Cat_3[i_db]), string(db.Cat_4[i_db])])
	end
	pl_GCxGCxGC[!, :tR2] = tR2
	pl_GCxGCxGC[!, :τR2] = τR2
	pl_GCxGCxGC[!, :tR3] = tR3
	pl_GCxGCxGC[!, :τR3] = τR3
	pl_GCxGCxGC[!, :Cat] = Cat
	return pl_GCxGCxGC
end

# ╔═╡ 0bc7a04d-7534-4e0b-993c-a6ca65f8d363
pl = peaklist_GCxGCxGC(sim_res[2][1][1], sim_res[2][1][2], sim_res[2][1][3])

# ╔═╡ 024984d7-1771-4bbb-8d58-be43c7793710
begin
	plotly()
	p_t1_t2 = Plots.scatter(pl.tR1, pl.tR2.-pl.tR1, xlabel="t¹ in s", ylabel="t² in s", legend=false)
end

# ╔═╡ cda948c1-db04-4778-bc73-284c9e88ce9f
#savefig(p_t1_t2, "3D_chrom_t1_t2.svg")

# ╔═╡ 222323d5-d909-434e-9b87-921c5308b326
Plots.scatter(pl.tR1, pl.tR3.-pl.tR2, xlabel="t¹ in s", ylabel="t³ in s")

# ╔═╡ e2ca66f9-2271-48f9-8ac0-ee186e456a28
p_t2_t3 = Plots.scatter(pl.tR2.-pl.tR1, pl.tR3.-pl.tR2, xlabel="t² in s", ylabel="t³ in s", legend=false)

# ╔═╡ 5643b1d5-f693-4cc0-8c28-f97865d3f8bb
#savefig(p_t2_t3, "3D_chrom_t2_t3.svg")

# ╔═╡ ede44987-fb1a-458f-b700-f92e5f504aa1
Plots.scatter(pl.tR1, pl.tR2.-pl.tR1, pl.tR3.-pl.tR2, xlabel="t¹ in s", ylabel="t² in s", zlabel="t³ in s", legend=false)

# ╔═╡ 99932526-134e-48fc-8bcd-bf89297fd4f1
begin
	gr()
	p_chrom_3D = Plots.scatter(pl.tR1, pl.tR2.-pl.tR1, pl.tR3.-pl.tR2, xlabel="t¹ in s", ylabel="t² in s", zlabel="t³ in s", legend=false)
end

# ╔═╡ 8cb201e5-141d-4de0-a6e4-7d8e29f226c2
#savefig(p_chrom_3D, "3D_chrom.svg")

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
# ╠═9bd5cee4-fd2d-44e9-8818-6c4733565c02
# ╠═f9041774-42ff-4f86-993a-e3df618cc36e
# ╠═3444976f-7687-4a63-9810-45f74406f6f6
# ╠═7a47ad5b-23ad-4385-abc6-09b92cecf09a
# ╠═03190a64-09cf-4a9b-9b86-e2b50392481c
# ╠═d571712f-981e-4d1e-95ce-b4dbf7de50c6
# ╠═4b379226-1646-4b13-b64a-7b55b66442f5
# ╠═7783216f-6a73-44aa-98f7-e0ca628ba145
# ╠═b51877ef-a8a6-48d4-a521-5b1e01d2be40
# ╠═06b9ba3f-4dff-4fba-afce-923b95fae9c3
# ╠═9d02024e-10dd-4181-b280-6dcba25b21d4
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
# ╠═024984d7-1771-4bbb-8d58-be43c7793710
# ╠═cda948c1-db04-4778-bc73-284c9e88ce9f
# ╠═222323d5-d909-434e-9b87-921c5308b326
# ╠═e2ca66f9-2271-48f9-8ac0-ee186e456a28
# ╠═5643b1d5-f693-4cc0-8c28-f97865d3f8bb
# ╠═ede44987-fb1a-458f-b700-f92e5f504aa1
# ╠═99932526-134e-48fc-8bcd-bf89297fd4f1
# ╠═8cb201e5-141d-4de0-a6e4-7d8e29f226c2
# ╠═67c99b55-6874-4b3f-adea-caf12504bf46
