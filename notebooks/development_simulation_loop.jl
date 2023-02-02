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
	# setup of the Graph -> ATTENTION to the right order of the modules
	# look at `collect(edges(sys_loop.g))`
	g_loop = SimpleDiGraph(6)
	add_edge!(g_loop, 1, 2) # Inj->TL1->Split 
	add_edge!(g_loop, 2, 3) # Split -> TL2
	add_edge!(g_loop, 2, 5) # GC2
	add_edge!(g_loop, 3, 4) # GC1
	add_edge!(g_loop, 4, 5) # TL3
	add_edge!(g_loop, 5, 6) # TL4 -> Det
	# pressure points
	pp_loop = Array{GasChromatographySystems.PressurePoint}(undef, nv(g_loop))
	pp_loop[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 1800.0], [600000.0, 600000.0]) # inlet 
	pp_loop[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 1800.0], [NaN, NaN]) # 
	pp_loop[3] = GasChromatographySystems.PressurePoint("p₃", [0.0, 1800.0], [NaN, NaN]) #  
	pp_loop[4] = GasChromatographySystems.PressurePoint("p₄", [0.0, 1800.0], [NaN, NaN]) #
	pp_loop[5] = GasChromatographySystems.PressurePoint("p₅", [0.0, 1800.0], [NaN, NaN]) #
	pp_loop[6] = GasChromatographySystems.PressurePoint("p₆", [0.0, 1800.0], [101300.0, 101300.0]) # outlet
	# modules
	modules_loop = Array{GasChromatographySystems.AbstractModule}(undef, ne(g_loop))
	modules_loop[1] = GasChromatographySystems.ModuleColumn("1 TL1", 0.5, 0.25e-3, 0.25e-6, "Wax", 300.0)
	modules_loop[2] = GasChromatographySystems.ModuleColumn("2 TL2", 1.0, 0.25e-3, 0.25e-6, "Wax", 300.0)
	modules_loop[3] = GasChromatographySystems.ModuleColumn("5 GC2", 2.0, 0.1e-3, 0.1e-6, "SLB5ms", default_TP) 
	modules_loop[4] = GasChromatographySystems.ModuleColumn("3 GC1", 10.0, 0.25e-3, 0.25e-6, "Wax", default_TP)
	modules_loop[5] = GasChromatographySystems.ModuleColumn("4 TL3", 1.0, 0.1e-3, 0.1e-6, "Wax", default_TP)
	
	modules_loop[6] = GasChromatographySystems.ModuleColumn("6 TL4", 1.5, 0.1e-3, 0.1e-6, "", 300.0)
	# system
	sys_loop = GasChromatographySystems.update_system(GasChromatographySystems.System(g_loop, pp_loop, modules_loop, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ ea5d611b-6b21-4185-94fc-ce3a6293c0dd
GasChromatographySystems.plot_graph(sys_loop; color=:green)

# ╔═╡ 4b379226-1646-4b13-b64a-7b55b66442f5
md"""
## Pressure and flow calculation
"""

# ╔═╡ 63a6a5a0-a4d3-430e-af78-ba446c100db5
GasChromatographySystems.flow_balance(sys_loop)

# ╔═╡ 40c2806b-ed48-4948-8e8c-2e7c635b76f9
GasChromatographySystems.solve_balance(sys_loop)

# ╔═╡ 73f9a1c4-73c9-4690-a1ad-fcadeb6d575a
p_func = GasChromatographySystems.pressure_functions(sys_loop)

# ╔═╡ 46bec177-cdf1-4527-9862-daafe5af6828
F_func = GasChromatographySystems.flow_functions(sys_loop)

# ╔═╡ 06d0e7d5-4370-4901-88a5-728001c2bee3
GasChromatographySystems.plot_graph_with_flow(sys_loop, 0)

# ╔═╡ 0ac6c03b-5be1-41ea-ad99-380a670eee62
GasChromatographySystems.plot_graph_with_flow(sys_loop, 1800)

# ╔═╡ 7783216f-6a73-44aa-98f7-e0ca628ba145
begin
	plotly()
	com_timesteps = GasChromatographySystems.common_timesteps(sys_loop)
	trange = 0:sum(com_timesteps)/100:sum(com_timesteps)
	p_flow = Plots.plot(xlabel="time in s", ylabel="flow in mL/min")
	for i=1:ne(sys_loop.g)
		Plots.plot!(p_flow, trange, F_func[i].(trange)*60*1e6, label="F_$(sys_loop.modules[i].name)")
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
unique(GasChromatographySystems.all_stationary_phases(sys_loop))

# ╔═╡ 2eff2573-137d-4a85-9bc9-e43d081b8d77
GasChromatographySystems.common_solutes(db, sys_loop)

# ╔═╡ 5e824c58-4d19-47fa-a1f9-27b977cb0761
selected_solutes = ["Pentadecane", "2,6-Dimethyloctane"]

# ╔═╡ e935ae05-c654-4365-a569-fed418f35fe1
md"""
## Graph to Parameters
"""

# ╔═╡ d060c9a3-d760-472f-b52b-3c481b39a85c
par_split = GasChromatographySystems.graph_to_parameters(sys_loop, db, selected_solutes)

# ╔═╡ 6dce6315-333a-4633-9c35-77a45e041b06
par_split[1].prog.pout_itp(0)

# ╔═╡ bff96e1c-fb9f-48b9-8ea7-2b9295054b8d
par_split[1].prog.pout_itp(1800)

# ╔═╡ 79bc683a-8d51-4d81-a27c-e42663720917
md"""
## Simulation
"""

# ╔═╡ 504dbc06-0c83-4e1c-8003-a1f562d71336
paths = GasChromatographySystems.all_paths(sys_loop.g, 2)[2]

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

# ╔═╡ 347bd70a-967a-4c41-9d52-87a04b2e86ff
simulate_along_paths_brute_force(sys_loop, paths, db, selected_solutes)

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

# ╔═╡ ae0aca4e-2ba4-40bc-9b34-db5ee741557a
simulate_along_paths(sys_loop, paths, db, selected_solutes)

# ╔═╡ 65fee760-73bd-457e-b709-7f59fc7ffe26
md"""
## Manually simulated
"""

# ╔═╡ 98538dc1-c8a9-4cd8-a270-7d50ca450c33
par_sys_ = GasChromatographySystems.graph_to_parameters(sys_loop, db, selected_solutes)

# ╔═╡ ae9f6808-5721-4344-ac8b-4e7022affaa7
GasChromatographySystems.plot_graph_with_flow(sys_loop, 0)

# ╔═╡ c6ead0e2-80fb-483a-b498-986e50915e97


# ╔═╡ bc87a815-0568-41f6-991b-7150191eb861
md"""
## Path 1
"""

# ╔═╡ a6e2998c-cfde-4359-82ff-24d90a4acc8e
paths[1]

# ╔═╡ 084ac0f5-a212-4be8-9461-77d97c52c917
md"""
### 1 TL1
"""

# ╔═╡ ecd3507d-8fca-49c5-8071-cf17608a979b
par_sys_1 = GasChromatographySystems.change_initial(par_sys_[1], zeros(length(selected_solutes)), zeros(length(selected_solutes)))

# ╔═╡ ec6c33ec-2c36-4fb6-9c68-9e871408a798
pl1, sol1 = GasChromatographySimulator.simulate(par_sys_1)

# ╔═╡ e4a2613f-c7e6-4417-8453-a0fda9e9afc3
md"""
### 2 TL2
"""

# ╔═╡ d9df22c8-2245-4c47-95d7-10d1a706a1a3
par_sys_2 = GasChromatographySystems.change_initial(par_sys_[2], pl1.tR, pl1.τR)

# ╔═╡ 5d5afb9c-9333-4d4c-9222-9f27b5c081f9
pl2, sol2 = GasChromatographySimulator.simulate(par_sys_2)

# ╔═╡ 5a30539d-9aa0-46b7-80a8-5c00c9955e25
md"""
### 3 GC1
"""

# ╔═╡ 02aa913f-8c9b-4cbf-bbc8-b0b3dd13aa2c
par_sys_3 = GasChromatographySystems.change_initial(par_sys_[4], pl2.tR, pl2.τR)

# ╔═╡ 364b6f72-b704-49a7-bb86-380be41ea2fc
pl3, sol3 = GasChromatographySimulator.simulate(par_sys_3)

# ╔═╡ 52c21ba0-afa2-40d6-92f4-ed4899c07056
md"""
### 4 TL3
"""

# ╔═╡ aca53f94-3709-4466-a692-3334179e7d31
par_sys_4 = GasChromatographySystems.change_initial(par_sys_[5], pl3.tR, pl3.τR)

# ╔═╡ d15acf8b-cf02-4d27-86c4-fdde9860f749
pl4, sol4 = GasChromatographySimulator.simulate(par_sys_4)

# ╔═╡ bbbea6e1-a1f1-45fc-9253-ee4236f68fd0
md"""
### 6 TL4
"""

# ╔═╡ d7bfce12-440f-47c4-8cba-db98789ce215
par_sys_6 = GasChromatographySystems.change_initial(par_sys_[6], pl4.tR, pl4.τR)

# ╔═╡ d371e4ed-c862-43e1-81d9-225abc057e81
pl6, sol6 = GasChromatographySimulator.simulate(par_sys_6)

# ╔═╡ 8fe53a28-70e5-434e-8032-330f41fb4596
md"""
## Path 2
"""

# ╔═╡ 36ebe09e-3f19-452b-a7ca-740d0992956b
paths[2]

# ╔═╡ e67a7327-5010-4360-8943-f4c0dc32f7af
md"""
## 1 TL1
"""

# ╔═╡ 0c4af679-e705-4c80-a1c9-cfc32a53c6d2
par_sys_1_2 = GasChromatographySystems.change_initial(par_sys_[1], zeros(length(selected_solutes)), zeros(length(selected_solutes)))

# ╔═╡ a385762e-2724-4833-a542-0cc8bab86ad5
pl1_2, sol1_2 = GasChromatographySimulator.simulate(par_sys_1_2)

# ╔═╡ 77bc166f-26f9-47a9-b6a1-36e0e22781ae
md"""
### 5 GC2
"""

# ╔═╡ e855cf55-0aaa-410b-805d-8db172f8ec49
par_sys_5_2 = GasChromatographySystems.change_initial(par_sys_[3], pl1_2.tR, pl1_2.τR)

# ╔═╡ b0018815-7476-43a3-a72c-289d6a7b4316
pl5_2, sol5_2 = GasChromatographySimulator.simulate(par_sys_5_2)

# ╔═╡ 781147c9-b739-4849-a3ea-7b37affc4b8c
md"""
### 6 TL4
"""

# ╔═╡ 257a3374-16e8-44b8-90e8-132ad4cb7e5e
par_sys_6_2 = GasChromatographySystems.change_initial(par_sys_[6], pl5_2.tR, pl5_2.τR)

# ╔═╡ a4583645-9e09-416c-8c54-d86c5d98c1bd
pl6_2, sol6_2 = GasChromatographySimulator.simulate(par_sys_6_2)

# ╔═╡ df855d1b-cd83-45b8-99ec-56861a56cc87
md"""
## Functio step-by-step
"""

# ╔═╡ 3d3d9073-a134-4829-a1d1-299717dc0a57
begin
	sys = sys_loop
	db_dataframe = db
	t₀=zeros(length(selected_solutes))
	τ₀=zeros(length(selected_solutes))
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
	#return path_pos, peaklists, solutions, new_par_sys
end

# ╔═╡ 83062c0b-2621-4d12-a889-599551a1f441
path_pos

# ╔═╡ 39f4c5da-b4c7-4415-bdf0-6cd22eb763da
peaklists

# ╔═╡ 6069b868-a13b-4db9-85f7-1c1621891db8
md"""
# End
"""

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
# ╠═06d0e7d5-4370-4901-88a5-728001c2bee3
# ╠═0ac6c03b-5be1-41ea-ad99-380a670eee62
# ╠═7783216f-6a73-44aa-98f7-e0ca628ba145
# ╠═a65c868c-fa26-4a7f-be94-2d86ee9518fc
# ╠═b424f5f1-aca8-40e5-aa0d-32fc7c4d4ee6
# ╠═0b0c0577-ec2e-450f-9f19-4a29da213efe
# ╠═2eff2573-137d-4a85-9bc9-e43d081b8d77
# ╠═5e824c58-4d19-47fa-a1f9-27b977cb0761
# ╠═e935ae05-c654-4365-a569-fed418f35fe1
# ╠═d060c9a3-d760-472f-b52b-3c481b39a85c
# ╠═6dce6315-333a-4633-9c35-77a45e041b06
# ╠═bff96e1c-fb9f-48b9-8ea7-2b9295054b8d
# ╠═79bc683a-8d51-4d81-a27c-e42663720917
# ╠═504dbc06-0c83-4e1c-8003-a1f562d71336
# ╠═12cec32a-333e-442e-8394-96863c3dc369
# ╠═347bd70a-967a-4c41-9d52-87a04b2e86ff
# ╠═cc2cb755-6ecc-4609-920c-51dff152db1f
# ╠═ae0aca4e-2ba4-40bc-9b34-db5ee741557a
# ╠═65fee760-73bd-457e-b709-7f59fc7ffe26
# ╠═98538dc1-c8a9-4cd8-a270-7d50ca450c33
# ╠═ae9f6808-5721-4344-ac8b-4e7022affaa7
# ╠═c6ead0e2-80fb-483a-b498-986e50915e97
# ╠═bc87a815-0568-41f6-991b-7150191eb861
# ╠═a6e2998c-cfde-4359-82ff-24d90a4acc8e
# ╠═084ac0f5-a212-4be8-9461-77d97c52c917
# ╠═ecd3507d-8fca-49c5-8071-cf17608a979b
# ╠═ec6c33ec-2c36-4fb6-9c68-9e871408a798
# ╠═e4a2613f-c7e6-4417-8453-a0fda9e9afc3
# ╠═d9df22c8-2245-4c47-95d7-10d1a706a1a3
# ╠═5d5afb9c-9333-4d4c-9222-9f27b5c081f9
# ╠═5a30539d-9aa0-46b7-80a8-5c00c9955e25
# ╠═02aa913f-8c9b-4cbf-bbc8-b0b3dd13aa2c
# ╠═364b6f72-b704-49a7-bb86-380be41ea2fc
# ╠═52c21ba0-afa2-40d6-92f4-ed4899c07056
# ╠═aca53f94-3709-4466-a692-3334179e7d31
# ╠═d15acf8b-cf02-4d27-86c4-fdde9860f749
# ╠═bbbea6e1-a1f1-45fc-9253-ee4236f68fd0
# ╠═d7bfce12-440f-47c4-8cba-db98789ce215
# ╠═d371e4ed-c862-43e1-81d9-225abc057e81
# ╠═8fe53a28-70e5-434e-8032-330f41fb4596
# ╠═36ebe09e-3f19-452b-a7ca-740d0992956b
# ╠═e67a7327-5010-4360-8943-f4c0dc32f7af
# ╠═0c4af679-e705-4c80-a1c9-cfc32a53c6d2
# ╠═a385762e-2724-4833-a542-0cc8bab86ad5
# ╠═77bc166f-26f9-47a9-b6a1-36e0e22781ae
# ╠═e855cf55-0aaa-410b-805d-8db172f8ec49
# ╠═b0018815-7476-43a3-a72c-289d6a7b4316
# ╠═781147c9-b739-4849-a3ea-7b37affc4b8c
# ╠═257a3374-16e8-44b8-90e8-132ad4cb7e5e
# ╠═a4583645-9e09-416c-8c54-d86c5d98c1bd
# ╠═df855d1b-cd83-45b8-99ec-56861a56cc87
# ╠═3d3d9073-a134-4829-a1d1-299717dc0a57
# ╠═83062c0b-2621-4d12-a889-599551a1f441
# ╠═39f4c5da-b4c7-4415-bdf0-6cd22eb763da
# ╠═6069b868-a13b-4db9-85f7-1c1621891db8
