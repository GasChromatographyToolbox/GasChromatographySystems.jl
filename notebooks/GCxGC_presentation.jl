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

# ╔═╡ af30bbcd-d6e4-4bd7-97f0-3243c514194a
plot_flow_over_time(sys)

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
GasChromatographySystems.plot_graph_with_flow(sys, 0; lay=SquareGrid(cols=3), node_size=90, nlabels_fontsize=20, elabels_fontsize=20, elabels_distance=28)

# ╔═╡ d256dfee-4639-4687-bba0-c1c662e1a12b
function plot_graph_with_flow(sys, t; lay = Stress(), color=:lightblue, node_size=80, arrow_size=20)
	p_func = GasChromatographySystems.pressure_functions(sys)
	F_func = GasChromatographySystems.flow_functions(sys)
	fig, ax, p = GraphMakie.graphplot(sys.g, 
						layout=lay,
						nlabels = [sys.pressurepoints[i].name*"\n$(trunc(Int, p_func[i](t)/1000))kPa" for i in 1:nv(sys.g)],
						nlabels_align=(:center,:center),
						node_size = [node_size for i in 1:nv(sys.g)],
						node_color = [color for i in 1:nv(sys.g)],
						elabels = [sys.modules[i].name*"\n $(round(F_func[i](t)*60*1e6; sigdigits=2)) mL/min" for i in 1:ne(sys.g)],
						elabels_distance = 20,
						arrow_size = arrow_size
					)
	hidedecorations!(ax)
	hidespines!(ax)
	#ax.aspect = DataAspect()
	return fig
end

# ╔═╡ 9c9aed82-3015-419b-ae6a-aa6f97a9d753
plot_pressure_over_time(sys)

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

# ╔═╡ cc85656c-59f9-4134-99d6-3c33103372b8
function change_initial_focussed(par::GasChromatographySimulator.Parameters, pl; τ₀=zeros(length(pl.tR)))
	# copys the parameters `par` and changes the values of par.sub[i].t₀ to pl.tR[]
	CAS_pl = GasChromatographySimulator.CAS_identification(pl.Name).CAS
	newsub = Array{GasChromatographySimulator.Substance}(undef, length(par.sub))
	for i=1:length(par.sub)
		ii = findfirst(par.sub[i].CAS.==CAS_pl)
		newsub[i] = GasChromatographySimulator.Substance(par.sub[i].name, par.sub[i].CAS, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀, par.sub[i].ann, par.sub[i].Cag, pl.tR[ii], 0.0)
	end
	newpar = GasChromatographySimulator.Parameters(par.col, par.prog, newsub, par.opt)
	return newpar
end

# ╔═╡ 53b140dc-d5ba-4808-83c9-af711563be9e
function simulate_along_paths(sys, paths, db_dataframe, selected_solutes, par_sys; t₀=zeros(length(selected_solutes)), τ₀=zeros(length(selected_solutes)), refocus=falses(ne(sys.g)))
	#par_sys = graph_to_parameters(sys, db_dataframe, selected_solutes)
	E = collect(edges(sys.g))
	peaklists = Array{Array{DataFrame,1}}(undef, length(paths))
	solutions = Array{Array{Any,1}}(undef, length(paths))
	path_pos = Array{String}(undef, length(paths))
	new_par_sys = Array{GasChromatographySimulator.Parameters}(undef, length(par_sys))
	visited_E = falses(length(E))
	for i=1:length(paths)
		i_par = GasChromatographySystems.index_parameter(sys.g, paths[i])
		if GasChromatographySystems.path_possible(sys, paths[i]) == true
			path_pos[i] = "path is possible"
			peaklists_ = Array{DataFrame}(undef, length(i_par))
			solutions_ = Array{Any}(undef, length(i_par))
			for j=1:length(i_par)
				if (i>1) && (all(visited_E[1:i_par[j]].==true))
					# look in all previous paths for the correct result -> the simulation correlated to the same edge and where this edge is connected to only previouse visited edges
					i_path = 0
					for k=1:i-1 # previous paths
						i_par_previous = GasChromatographySystems.index_parameter(sys.g, paths[k])[1:j]
						if all(x->x in i_par_previous, i_par[1:j]) == true # all edges up to j are the same between the two paths
							i_path = k
						end
					end
					peaklists_[j] = peaklists[i_path][i_par[j]]
					solutions_[j] = solutions[i_path][i_par[j]]
				else
					if j == 1
						new_par_sys[i_par[j]] = GasChromatographySystems.change_initial(par_sys[i_par[j]], t₀, τ₀)
						peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(new_par_sys[i_par[j]])
					else
						if refocus[i_par[j]] == true
							new_par_sys[i_par[j]] = change_initial_focussed(par_sys[i_par[j]], peaklists_[j-1])
						else
							new_par_sys[i_par[j]] = GasChromatographySystems.change_initial(par_sys[i_par[j]], peaklists_[j-1])
						end
						peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(new_par_sys[i_par[j]])
					end
				end
			end
			visited_E[i_par] .= true
			peaklists[i] = peaklists_
			solutions[i] = solutions_
		else
			neg_flow_modules = sys.modules[findall(paths[i][findall(GasChromatographySystems.positive_flow(sys)[i_par].==false)].==edges(sys.g))]
			str_neg_flow = "path not possible: "
			for j=1:length(neg_flow_modules)
				str_neg_flow = str_neg_flow*"Flow in module >$(neg_flow_modules[j].name)< becomes negative during the program. "
			end
			path_pos[i] = str_neg_flow
		end
	end
	return path_pos, peaklists, solutions, new_par_sys
end

# ╔═╡ 1436c866-9124-4ab4-8ab7-0b04f4de896c
md"""
### 1. GCxGC common temperature program
"""

# ╔═╡ 457df2b2-8bfa-4e04-b4f2-ae1ca8902a4e
paths = GasChromatographySystems.all_paths(sys.g, 1)[2]

# ╔═╡ f106b5fd-e17a-4876-8562-e82670a1ab25
sim_res = GasChromatographySystems.simulate_along_paths(sys, paths, db, selected_solutes, par)

# ╔═╡ d8c30c13-7209-43ea-93c8-ac0310b0a8e5
begin
	par_ = GasChromatographySystems.graph_to_parameters(sys, db, selected_solutes)
	paths_ = GasChromatographySystems.all_paths(sys.g, 1)[2]
	sim_res_ = GasChromatographySystems.simulate_along_paths(sys, paths_, db, selected_solutes, par_)
end

# ╔═╡ 7dac264f-2e91-44bd-b96d-207f072cba9c
sim_res_foc = simulate_along_paths(sys, paths, db, selected_solutes, par; refocus=[true, true])

# ╔═╡ 5373eaae-4372-4561-83cf-b62664bf28af
sim_res[2][1][1]

# ╔═╡ e57ab914-8db4-41ca-9950-1134e35e8c40
function peaklist_GCxGC(pl_1, pl_2)
	CAS1 = GasChromatographySimulator.CAS_identification(pl_1.Name).CAS
	pl_GCxGC = DataFrame(Name = pl_1.Name, CAS = CAS1, tR1 = pl_1.tR, τR1 = pl_1.τR)
	tR2 = Array{Float64}(undef, length(pl_GCxGC.Name))
	τR2 = Array{Float64}(undef, length(pl_GCxGC.Name))
	Cat = Array{Array{String}}(undef, length(pl_GCxGC.Name))
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

# ╔═╡ 5429291b-8efd-4ebb-9d4d-1be143073482
pl_foc = peaklist_GCxGC(sim_res_foc[2][1][1], sim_res_foc[2][1][2])

# ╔═╡ c5888be4-c641-4871-81a5-11d01ace12ed
#savefig(p_GCxGC_cats, "chrom_GCxGC_10.svg")

# ╔═╡ 127069e6-6da1-4742-88b4-a0d398a88b24
function plot_GCxGC(pl_GCxGC, sys; categories = String[])
	t1end = sum(GasChromatographySystems.common_timesteps(sys))
	t2end = ceil(maximum(pl_GCxGC.tR2.-pl_GCxGC.tR1))
	p_gcxgc = Plots.scatter(pl_GCxGC.tR1, pl_GCxGC.tR2.-pl_GCxGC.tR1, xlims=(0.0,t1end), ylims=(0.0, t2end), label="", xlabel=string(sys.modules[1].stationary_phase, ", t¹ in s"), ylabel=string(sys.modules[2].stationary_phase, ", t² in s"))
	if !isempty(categories)
		for i=1:length(categories) 
			pl_f = filter([:Cat] => x -> categories[i] in x, pl_GCxGC)
			Plots.scatter!(p_gcxgc, pl_f.tR1, pl_f.tR2.-pl_f.tR1, label=categories[i], legend=:bottomright)
		end
	end
	return p_gcxgc
end

# ╔═╡ 3986a24a-689c-48b8-8e78-ec5d3842d0d7
plot_GCxGC(pl, sys)

# ╔═╡ 938d8677-7898-43c8-b07c-6e9d7f3c7e2e
begin
	gr()
	p_GCxGC_cats = plot_GCxGC(pl, sys; categories = ["saturated FAME", "unsaturated FAME", "aromatic", "PCB", "alkanes", "alcohols", "alkanones", "phenones"])
	p_GCxGC_cats
end

# ╔═╡ 987efbb5-8a48-4c83-9a20-50e5fbd0a16d
function plot_GCxGC_contour_(pl_GCxGC; split_ratio=ones(length(pl_GCxGC.Name)))
	t¹ = 0.9*minimum(pl_GCxGC.tR1):(1.1*maximum(pl_GCxGC.tR1)-0.9*minimum(pl_GCxGC.tR1))/1000:1.1*maximum(pl_GCxGC.tR1)
	t² = 0.9*minimum(pl_GCxGC.tR2.-pl_GCxGC.tR1):(1.1*maximum(pl_GCxGC.tR2.-pl_GCxGC.tR1)-0.9*minimum(pl_GCxGC.tR2.-pl_GCxGC.tR1))/1000:1.1*maximum(pl_GCxGC.tR2.-pl_GCxGC.tR1)
	chrom2D = Array{Float64}(undef, length(t¹), length(t²), length(pl_GCxGC.Name))
	for j=1:length(t²)
		for i=1:length(t¹)
			for k=1:length(pl_GCxGC.Name)
				chrom2D[i,j,k] = split_ratio[k]*exp(-((t¹[i]-pl_GCxGC.tR1[k])^2/(2*pl_GCxGC.τR1[k]^2)+(t²[j]-(pl_GCxGC.tR2[k]-pl_GCxGC.tR1[k]))^2/(2*pl_GCxGC.τR2[k]^2)))
			end
		end
	end
	#p_2D = Plots.contourf(t¹, t², sum(chrom2D, dims=3)[:,:,1]', fill=true, legend=false, xlabel="tR1 in s", ylabel="tR2 in s", c=:turbo, levels=40)
	return t¹, t², sum(chrom2D, dims=3)[:,:,1]'
end

# ╔═╡ 939e3a45-03cd-4e07-a672-10c32cf2637e
t¹, t², chrom = plot_GCxGC_contour_(pl_foc)

# ╔═╡ 7100eeed-5604-4d81-855e-b08cd40135a6
Plots.contour(t¹, t², chrom, legend=false, xlabel="tR1 in s", ylabel="tR2 in s", c=:turbo, levels=100, xlims=(0.0, 1700.0), ylims=(0.0, 5.0), grid=false)

# ╔═╡ 00fbcd7c-71a1-43f7-b358-abfde851a237
begin
	plotly()
	p_GCxGC = Plots.heatmap(t¹, t², chrom.^(1//3), legend=false, xlabel="tR1 in s", ylabel="tR2 in s", c=:jet1, dpi=500)#, xlims=(0.0, 1700.0), ylims=(0.0, 5.0), grid=false)
	gui()
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

# ╔═╡ e7956916-dccb-46c1-a30b-ef77650e54ce
md"""
### 1b. reverse phases
"""

# ╔═╡ 2c45c22a-842f-4d49-8b0b-d6f9546932c9
par_rev = GasChromatographySystems.graph_to_parameters(sys_rev, db, selected_solutes)

# ╔═╡ eccd7341-eae1-4a37-9cfd-f1817857e572
paths_rev = GasChromatographySystems.all_paths(sys_rev.g, 1)[2]

# ╔═╡ aa0d4fda-3f1b-468c-9568-acadf61b8b95
sim_res_rev = GasChromatographySystems.simulate_along_paths(sys_rev, paths_rev, db, selected_solutes, par_rev)

# ╔═╡ 397a913b-fa6c-4eba-b8b7-7571256a001b
pl_rev = peaklist_GCxGC(sim_res_rev[2][1][1], sim_res_rev[2][1][2])

# ╔═╡ c3f16558-43b8-4822-98af-daf8e06cd1df
p_GCxGC_cats_rev = plot_GCxGC(pl_rev, sys_rev; categories = ["saturated FAME", "unsaturated FAME", "aromatic", "PCB", "alkanes", "alcohols", "alkanones", "phenones"])

# ╔═╡ 834739e9-fa7f-43f0-beac-ae5790575255
#savefig(p_GCxGC_cats_rev, "chrom_GCxGC_8_rev.svg")

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
p_GCxGC_cats_offset = plot_GCxGC(pl_offset, sys_offset; categories = ["saturated FAME", "unsaturated FAME", "aromatic", "PCB", "alkanes", "alcohols", "alkanones", "phenones"])

# ╔═╡ 8d14ae41-a062-4d89-8e5a-e78b62cc0564
#savefig(p_GCxGC_cats_offset, "chrom_GCxGC_8_offset.svg")

# ╔═╡ 67c99b55-6874-4b3f-adea-caf12504bf46
md"""
# End
"""

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
# ╠═6cfda8f9-f918-4451-97c3-084eb1d9d339
# ╠═a8971018-3f53-4276-b004-9245985ca490
# ╠═e4fa29f4-fc83-4aa4-a548-f2683ad5d335
# ╠═38a4081d-49ed-4688-84ba-88218f535911
# ╠═4b379226-1646-4b13-b64a-7b55b66442f5
# ╠═7783216f-6a73-44aa-98f7-e0ca628ba145
# ╠═b51877ef-a8a6-48d4-a521-5b1e01d2be40
# ╠═06b9ba3f-4dff-4fba-afce-923b95fae9c3
# ╠═d256dfee-4639-4687-bba0-c1c662e1a12b
# ╠═9c9aed82-3015-419b-ae6a-aa6f97a9d753
# ╠═a65c868c-fa26-4a7f-be94-2d86ee9518fc
# ╠═b424f5f1-aca8-40e5-aa0d-32fc7c4d4ee6
# ╠═0b0c0577-ec2e-450f-9f19-4a29da213efe
# ╠═3cd3ea64-a45c-49fb-90bb-c807c47eefc0
# ╠═e935ae05-c654-4365-a569-fed418f35fe1
# ╠═d060c9a3-d760-472f-b52b-3c481b39a85c
# ╠═79bc683a-8d51-4d81-a27c-e42663720917
# ╠═cc85656c-59f9-4134-99d6-3c33103372b8
# ╠═53b140dc-d5ba-4808-83c9-af711563be9e
# ╟─1436c866-9124-4ab4-8ab7-0b04f4de896c
# ╠═457df2b2-8bfa-4e04-b4f2-ae1ca8902a4e
# ╠═f106b5fd-e17a-4876-8562-e82670a1ab25
# ╠═d8c30c13-7209-43ea-93c8-ac0310b0a8e5
# ╠═7dac264f-2e91-44bd-b96d-207f072cba9c
# ╠═5373eaae-4372-4561-83cf-b62664bf28af
# ╠═0bc7a04d-7534-4e0b-993c-a6ca65f8d363
# ╠═5429291b-8efd-4ebb-9d4d-1be143073482
# ╠═e57ab914-8db4-41ca-9950-1134e35e8c40
# ╠═3986a24a-689c-48b8-8e78-ec5d3842d0d7
# ╠═938d8677-7898-43c8-b07c-6e9d7f3c7e2e
# ╠═c5888be4-c641-4871-81a5-11d01ace12ed
# ╠═127069e6-6da1-4742-88b4-a0d398a88b24
# ╠═987efbb5-8a48-4c83-9a20-50e5fbd0a16d
# ╠═939e3a45-03cd-4e07-a672-10c32cf2637e
# ╠═7100eeed-5604-4d81-855e-b08cd40135a6
# ╠═00fbcd7c-71a1-43f7-b358-abfde851a237
# ╠═fba31c65-0465-458a-8fa5-7df1abdf179c
# ╠═88bf3e84-6772-453d-8bee-2ccbbcf7cb34
# ╠═0cc7f988-bf98-4e45-93d5-7e3dea48c4ed
# ╠═e7956916-dccb-46c1-a30b-ef77650e54ce
# ╠═2c45c22a-842f-4d49-8b0b-d6f9546932c9
# ╠═eccd7341-eae1-4a37-9cfd-f1817857e572
# ╠═aa0d4fda-3f1b-468c-9568-acadf61b8b95
# ╠═397a913b-fa6c-4eba-b8b7-7571256a001b
# ╠═c3f16558-43b8-4822-98af-daf8e06cd1df
# ╠═834739e9-fa7f-43f0-beac-ae5790575255
# ╠═679140af-7598-41fe-aad1-90c50e503218
# ╠═1f92dabf-8597-4ef8-9e05-628e922defbf
# ╠═6efdb3ef-4f6e-4073-bfa1-501732088e5c
# ╠═097f05ad-8723-4693-bfba-c350cd478f54
# ╠═417fa98b-38d1-434a-a44c-2948d71bec42
# ╠═28fbe11a-4132-41e1-8ea8-3b65749f41b1
# ╠═8d14ae41-a062-4d89-8e5a-e78b62cc0564
# ╠═67c99b55-6874-4b3f-adea-caf12504bf46
