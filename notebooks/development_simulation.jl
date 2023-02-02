### A Pluto.jl notebook ###
# v0.19.17

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
# Development of the simulation of GC systems along the paths
"""

# ╔═╡ f6c6443d-0e39-4e82-b2c4-49f2954e5b49
md"""
## Definition of the system
"""

# ╔═╡ f484559b-95d5-4155-aad6-fdeed3d238bc
default_TP = GasChromatographySystems.TemperatureProgram([0.0, 1800.0], [40.0, 240.0])

# ╔═╡ 0a0fff3d-1744-4a2a-977f-68a784e1ccc9
begin
	g_split = SimpleDiGraph(4)
	add_edge!(g_split, 1, 2) # Inj -> GC column -> Split point
	add_edge!(g_split, 2, 3) # Split point -> TL column -> Det 1
	add_edge!(g_split, 2, 4) # Split point -> TL column -> Det 2
	# pressure points
	pp_split = Array{GasChromatographySystems.PressurePoint}(undef, nv(g_split))
	pp_split[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 1800.0], [200000.0, 250000.0]) # inlet 
	pp_split[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 1800.0], [NaN, NaN]) # 
	pp_split[3] = GasChromatographySystems.PressurePoint("p₃", [0.0, 1800.0], [eps(Float64), eps(Float64)]) # outlet 1 
	pp_split[4] = GasChromatographySystems.PressurePoint("p₄", [0.0, 1800.0], [101300.0, 101300.0]) # outlet 2
	# modules
	modules_split = Array{GasChromatographySystems.AbstractModule}(undef, ne(g_split))
	modules_split[1] = GasChromatographySystems.ModuleColumn("GC column", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	modules_split[2] = GasChromatographySystems.ModuleColumn("TL column 1", 0.5, 0.1e-3, 0.0e-6, "", 300.0)
	modules_split[3] = GasChromatographySystems.ModuleColumn("TL column 2", 1.0, 0.15e-3, 0.15e-6, "SLB5ms", 300.0)
	# system
	sys_split = GasChromatographySystems.update_system(GasChromatographySystems.System(g_split, pp_split, modules_split, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ ea5d611b-6b21-4185-94fc-ce3a6293c0dd
GasChromatographySystems.plot_graph(sys_split; color=:green)

# ╔═╡ 4b379226-1646-4b13-b64a-7b55b66442f5
md"""
## Pressure and flow calculation
"""

# ╔═╡ 63a6a5a0-a4d3-430e-af78-ba446c100db5
GasChromatographySystems.flow_balance(sys_split)

# ╔═╡ 40c2806b-ed48-4948-8e8c-2e7c635b76f9
GasChromatographySystems.solve_balance(sys_split)

# ╔═╡ 73f9a1c4-73c9-4690-a1ad-fcadeb6d575a
p_func = GasChromatographySystems.pressure_functions(sys_split)

# ╔═╡ 46bec177-cdf1-4527-9862-daafe5af6828
F_func = GasChromatographySystems.flow_functions(sys_split)

# ╔═╡ 0ac6c03b-5be1-41ea-ad99-380a670eee62
GasChromatographySystems.plot_graph_with_flow(sys_split, 0)

# ╔═╡ 7783216f-6a73-44aa-98f7-e0ca628ba145
begin
	plotly()
	com_timesteps = GasChromatographySystems.common_timesteps(sys_split)
	trange = 0:sum(com_timesteps)/100:sum(com_timesteps)
	p_flow = Plots.plot(xlabel="time in s", ylabel="flow in mL/min")
	for i=1:ne(sys_split.g)
		Plots.plot!(p_flow, trange, F_func[i].(trange)*60*1e6, label="F_$(sys_split.modules[i].name)")
	end
	p_flow
end

# ╔═╡ a65c868c-fa26-4a7f-be94-2d86ee9518fc
md"""
## Selection of solutes
"""

# ╔═╡ b424f5f1-aca8-40e5-aa0d-32fc7c4d4ee6
db = DataFrame(CSV.File(joinpath(dirname(pwd()), "data", "database_Kcentric.csv"), header=1, silencewarnings=true))

# ╔═╡ 0b0c0577-ec2e-450f-9f19-4a29da213efe
unique(GasChromatographySystems.all_stationary_phases(sys_split))

# ╔═╡ 2eff2573-137d-4a85-9bc9-e43d081b8d77
GasChromatographySystems.common_solutes(db, sys_split)

# ╔═╡ 5e824c58-4d19-47fa-a1f9-27b977cb0761
selected_solutes = ["Pentadecane", "C11 acid methyl ester"]

# ╔═╡ e935ae05-c654-4365-a569-fed418f35fe1
md"""
## Graph to Parameters
"""

# ╔═╡ d060c9a3-d760-472f-b52b-3c481b39a85c
par_split = GasChromatographySystems.graph_to_parameters(sys_split, db, selected_solutes)

# ╔═╡ 79bc683a-8d51-4d81-a27c-e42663720917
md"""
## Simulation
"""

# ╔═╡ 12cec32a-333e-442e-8394-96863c3dc369
function simulate_along_paths_brute_force(sys, paths, db_dataframe, selected_solutes; t₀=zeros(length(selected_solutes)), τ₀=zeros(length(selected_solutes)))
	par_sys = GasChromatographySystems.graph_to_parameters(sys, db_dataframe, selected_solutes)
	E = collect(edges(sys.g))
	peaklists = Array{Array{DataFrame,1}}(undef, length(paths))
	solutions = Array{Array{Any,1}}(undef, length(paths))
	path_pos = Array{String}(undef, length(paths))
	new_par_sys = Array{GasChromatographySimulator.Parameters}(undef, length(par_sys))
	for i=1:length(paths)
		i_par = GasChromatographySystems.index_parameter(sys.g, paths[i])
		if GasChromatographySystems.path_possible(sys, paths[i]) == true
			path_pos[i] = "path is possible"
			peaklists_ = Array{DataFrame}(undef, length(i_par))
			solutions_ = Array{Any}(undef, length(i_par))
			for j=1:length(i_par)
				if j == 1
					new_par_sys[i_par[j]] = GasChromatographySystems.change_initial(par_sys[i_par[j]], t₀, τ₀)
					peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(new_par_sys[i_par[j]])
				else
					new_par_sys[i_par[j]] = GasChromatographySystems.change_initial(par_sys[i_par[j]], peaklists_[j-1].tR, peaklists_[j-1].τR)
					peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(new_par_sys[i_par[j]])
				end
			end
			peaklists[i] = peaklists_
			solutions[i] = solutions_
		else
			neg_flow_modules = sys.modules[findall(paths[i][findall(positive_flow_(sys)[i_par].==false)].==edges(sys.g))]
			str_neg_flow = "path not possible: "
			for j=1:length(neg_flow_modules)
				str_neg_flow = str_neg_flow*"Flow in module >$(neg_flow_modules[j].name)< becomes negative during the program. "
			end
			path_pos[i] = str_neg_flow
		end
	end
	return path_pos, peaklists, solutions, new_par_sys
end

# ╔═╡ cc2cb755-6ecc-4609-920c-51dff152db1f
function simulate_along_paths(sys, paths, db_dataframe, selected_solutes; t₀=zeros(length(selected_solutes)), τ₀=zeros(length(selected_solutes)))
	par_sys = GasChromatographySystems.graph_to_parameters(sys, db_dataframe, selected_solutes)
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
						new_par_sys[i_par[j]] = GasChromatographySystems.change_initial(par_sys[i_par[j]], peaklists_[j-1].tR, peaklists_[j-1].τR)
						peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(new_par_sys[i_par[j]])
					end
				end
			end
			visited_E[i_par] .= true
			peaklists[i] = peaklists_
			solutions[i] = solutions_
		else
			neg_flow_modules = sys.modules[findall(paths[i][findall(positive_flow_(sys)[i_par].==false)].==edges(sys.g))]
			str_neg_flow = "path not possible: "
			for j=1:length(neg_flow_modules)
				str_neg_flow = str_neg_flow*"Flow in module >$(neg_flow_modules[j].name)< becomes negative during the program. "
			end
			path_pos[i] = str_neg_flow
		end
	end
	return path_pos, peaklists, solutions, new_par_sys
end

# ╔═╡ 457df2b2-8bfa-4e04-b4f2-ae1ca8902a4e
paths = GasChromatographySystems.all_paths(sys_split.g, 2)[2]

# ╔═╡ 5d4fb93c-fb08-4157-8690-2e4185aeb120
sim_res = simulate_along_paths_brute_force(sys_split, paths, db, selected_solutes)

# ╔═╡ f106b5fd-e17a-4876-8562-e82670a1ab25
sim_res_ = simulate_along_paths(sys_split, paths, db, selected_solutes)

# ╔═╡ Cell order:
# ╠═e9d64f0e-95b7-11ed-2fee-ab579cd02ebc
# ╠═139bfe64-17ac-455d-84a2-09c9564156cd
# ╠═f6c6443d-0e39-4e82-b2c4-49f2954e5b49
# ╠═f484559b-95d5-4155-aad6-fdeed3d238bc
# ╠═0a0fff3d-1744-4a2a-977f-68a784e1ccc9
# ╠═ea5d611b-6b21-4185-94fc-ce3a6293c0dd
# ╠═4b379226-1646-4b13-b64a-7b55b66442f5
# ╠═63a6a5a0-a4d3-430e-af78-ba446c100db5
# ╠═40c2806b-ed48-4948-8e8c-2e7c635b76f9
# ╠═73f9a1c4-73c9-4690-a1ad-fcadeb6d575a
# ╠═46bec177-cdf1-4527-9862-daafe5af6828
# ╠═0ac6c03b-5be1-41ea-ad99-380a670eee62
# ╠═7783216f-6a73-44aa-98f7-e0ca628ba145
# ╠═a65c868c-fa26-4a7f-be94-2d86ee9518fc
# ╠═b424f5f1-aca8-40e5-aa0d-32fc7c4d4ee6
# ╠═0b0c0577-ec2e-450f-9f19-4a29da213efe
# ╠═2eff2573-137d-4a85-9bc9-e43d081b8d77
# ╠═5e824c58-4d19-47fa-a1f9-27b977cb0761
# ╠═e935ae05-c654-4365-a569-fed418f35fe1
# ╠═d060c9a3-d760-472f-b52b-3c481b39a85c
# ╠═79bc683a-8d51-4d81-a27c-e42663720917
# ╠═12cec32a-333e-442e-8394-96863c3dc369
# ╠═cc2cb755-6ecc-4609-920c-51dff152db1f
# ╠═457df2b2-8bfa-4e04-b4f2-ae1ca8902a4e
# ╠═5d4fb93c-fb08-4157-8690-2e4185aeb120
# ╠═f106b5fd-e17a-4876-8562-e82670a1ab25
