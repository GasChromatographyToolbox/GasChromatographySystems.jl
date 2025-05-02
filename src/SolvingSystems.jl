# estimate possible paths in the graphs
function all_paths(g) # brute force method using multiple random walk and only using unique results
	num_paths = number_of_paths(g)
	rand_paths = Any[]
	while length(unique(rand_paths))<num_paths && length(rand_paths)<nv(g)*10
		push!(rand_paths, non_backtracking_randomwalk(g, 1, nv(g)))
	end
	Vp = sort(unique(rand_paths))
	Ep = Array{Array{Graphs.SimpleGraphs.SimpleEdge{Int64}, 1}}(undef, length(Vp))
	for j=1:length(Vp)
		Ep_ = Array{Graphs.SimpleGraphs.SimpleEdge{Int64}, 1}(undef, length(Vp[j])-1)
		for i=1:(length(Vp[j])-1)
			index = findfirst(Vp[j][i].==src.(edges(g)) .&& Vp[j][i+1].==dst.(edges(g)))
			Ep_[i] = collect(edges(g))[index]
		end
		Ep[j] = Ep_
	end
	return Vp, Ep
end

function all_paths(g, num_paths)
	rand_paths = Any[]
	while length(unique(rand_paths))<num_paths && length(rand_paths)<nv(g)*10
		push!(rand_paths, non_backtracking_randomwalk(g, 1, nv(g)))
	end
	Vp = sort(unique(rand_paths))
	Ep = Array{Array{Graphs.SimpleGraphs.SimpleEdge{Int64}, 1}}(undef, length(Vp))
	for j=1:length(Vp)
		Ep_ = Array{Graphs.SimpleGraphs.SimpleEdge{Int64}, 1}(undef, length(Vp[j])-1)
		for i=1:(length(Vp[j])-1)
			index = findfirst(Vp[j][i].==src.(edges(g)) .&& Vp[j][i+1].==dst.(edges(g)))
			Ep_[i] = collect(edges(g))[index]
		end
		Ep[j] = Ep_
	end
	return Vp, Ep
end

# simulate along the paths
function index_parameter(g, path)
	return findall(x->x in path, collect(edges(g)))
end

function common_edges(path1, path2)
	return path2[findall(x->x in path2, path1)]
end

function positive_flow(sys, p2fun; mode="λ")
	F_func = flow_functions(sys, p2fun; mode=mode)
	tend = sum(common_timesteps(sys))
	t = 0:tend/2:tend # !!!perhaps use interval arithmatic to test, if the flow is positive over the interval 0:tend!!!???
	#ipar = index_parameter(sys.g, path)
	pos_Flow = Array{Bool}(undef, length(F_func))
	for i=1:length(F_func)
		if isempty(findall(F_func[i].(t).<=0))
			pos_Flow[i] = true
		else
			pos_Flow[i] = false
		end
	end
	return pos_Flow
end

function path_possible(sys, p2fun, paths; mode="λ")
	#F_func = flow_functions(sys)
	i_E = index_parameter(sys.g, paths)
	if length(i_E) == length(findall(positive_flow(sys, p2fun; mode=mode)[i_E]))
		possible = true
	else
		possible = false
	end
	return possible
end

function change_initial(par::GasChromatographySimulator.Parameters, init_t, init_τ)
	# copys the parameters `par` and changes the values of par.sub[i].τ₀ and par.sub[i].t₀ to init_τ[i] resp. init_t[i]
	newsub = Array{GasChromatographySimulator.Substance}(undef, length(par.sub))
	for i=1:length(par.sub)
		newsub[i] = GasChromatographySimulator.Substance(par.sub[i].name, par.sub[i].CAS, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀, par.sub[i].ann, par.sub[i].Cag, init_t[i], init_τ[i])
	end
	newpar = GasChromatographySimulator.Parameters(par.col, par.prog, newsub, par.opt)
	return newpar
end

function change_initial(par::GasChromatographySimulator.Parameters, prev_pl)
	new_sub = Array{GasChromatographySimulator.Substance}(undef, length(prev_pl.CAS))
	for i=1:length(prev_pl.CAS)
		# filter for correct CAS and annotation (slice number)
		CAS_par = [par.sub[i].CAS for i in 1:length(par.sub)]
		i_sub = findfirst(prev_pl.CAS[i] .== CAS_par)
		#ii_ = common_index(prev_pl, prev_par.sub[i].CAS, join(split(prev_par.sub[i].ann, ", ")[1:end-1], ", "))
		new_sub[i] = GasChromatographySimulator.Substance(par.sub[i_sub].name, par.sub[i_sub].CAS, par.sub[i_sub].Tchar, par.sub[i_sub].θchar, par.sub[i_sub].ΔCp, par.sub[i_sub].φ₀, prev_pl.Annotations[i], par.sub[i_sub].Cag, prev_pl.tR[i], prev_pl.τR[i])
	end
	# here changes of options could be applied
	new_par = GasChromatographySimulator.Parameters(par.col, par.prog, new_sub, par.opt)
	return new_par
end

function change_initial_focussed(par::GasChromatographySimulator.Parameters, pl; τ₀=zeros(length(pl.tR)))
	# copys the parameters `par` and changes the values of par.sub[i].t₀ to pl.tR[]
	CAS_pl = [GasChromatographySimulator.CAS_identification(name).CAS for name in pl.Name] # CAS-numbers of the peaklist entries
	newsub = Array{GasChromatographySimulator.Substance}(undef, length(par.sub))
	for i=1:length(par.sub)
		ii = findfirst(par.sub[i].CAS.==CAS_pl)
		newsub[i] = GasChromatographySimulator.Substance(par.sub[i].name, par.sub[i].CAS, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀, par.sub[i].ann, par.sub[i].Cag, pl.tR[ii], τ₀[ii])
	end
	newpar = GasChromatographySimulator.Parameters(par.col, par.prog, newsub, par.opt)
	return newpar
end

# first column segment
function simulate_ModuleColumn(segment_par, t₀, τ₀)
	new_segment_par = GasChromatographySystems.change_initial(segment_par, t₀, τ₀)
	peaklist, solutions = GasChromatographySimulator.simulate(new_segment_par)
	peaklist[!,:A] = ones(length(peaklist.Name))
	return new_segment_par, peaklist, solutions
end

# later column segments
function simulate_ModuleColumn(segment_par, prev_peaklist)
	new_segment_par = GasChromatographySystems.change_initial(segment_par, prev_peaklist)
	peaklist, solutions = GasChromatographySimulator.simulate(new_segment_par)
	GasChromatographySystems.add_A_to_pl!(peaklist, prev_peaklist)
	return new_segment_par, peaklist, solutions
end

function simulate_along_paths(sys, p2fun, paths, db_dataframe, selected_solutes; t₀=zeros(length(selected_solutes)), τ₀=zeros(length(selected_solutes)), mode="λ")
	par_sys = graph_to_parameters(sys, p2fun, db_dataframe, selected_solutes, mode=mode)
	path_pos, peaklists, solutions, new_par_sys = simulate_along_paths(sys, p2fun, paths, par_sys; t₀=t₀, τ₀=τ₀, mode=mode)
	return path_pos, peaklists, solutions, new_par_sys
end

"""
    simulate_along_paths(sys, p2fun, paths, par_sys; t₀=zeros(length(par_sys[1].sub)), τ₀=zeros(length(par_sys[1].sub)), nτ=6, refocus=falses(ne(sys.g)), τ₀_focus=zeros(length(par_sys[1].sub)), mode="λ", kwargsTM...)

Simulate solute transport along multiple paths in a gas chromatography system.

This function simulates the transport of solutes through a network of GC modules (columns and thermal modulators)
along specified paths. It handles both direct simulation of new segments and reuse of previously simulated
segments when possible.

# Arguments
- `sys`: The GC system structure containing the network of modules
- `p2fun`: Pressure functions for the system
- `paths`: Array of paths to simulate, where each path is a sequence of edges in the system graph
- `par_sys`: Array of parameter sets for each module in the system

# Keyword Arguments
- `t₀`: Initial times for solutes (default: zeros)
- `τ₀`: Initial peak widths for solutes (default: zeros)
- `nτ`: Number of peak widths to consider for slicing in the thermal modulator (default: 6)
- `refocus`: Boolean array indicating which modules should refocus peaks (default: all false)
- `τ₀_focus`: Initial peak widths for refocusing (default: zeros)
- `mode`: Mode for flow calculations ("λ" for permeability or "κ" for restriction) (default: "λ")
- `kwargsTM`: Additional keyword arguments for thermal modulator simulation

# Returns
- `path_pos`: Array of strings indicating if each path is possible and any issues encountered
- `peaklists`: Array of peak lists for each path, containing retention times and peak widths
- `solutions`: Array of solutions for each path, containing detailed simulation results
- `new_par_sys`: Updated parameter sets for the system

# Notes
- Reuses simulation results from previous paths when possible to improve efficiency
- Handles both ModuleColumn and ModuleTM (thermal modulator) types
- Checks for negative flows to determine if paths are possible
- First segment after injection is assumed to be a ModuleColumn
"""
function simulate_along_paths(sys, p2fun, paths, par_sys; t₀=zeros(length(par_sys[1].sub)), τ₀=zeros(length(par_sys[1].sub)), nτ=6, refocus=falses(ne(sys.g)), τ₀_focus=zeros(length(par_sys[1].sub)), mode="λ", kwargsTM...)
	
	E = collect(edges(sys.g))
	peaklists = Array{Array{DataFrame,1}}(undef, length(paths))
	solutions = Array{Array{Any,1}}(undef, length(paths))
	path_pos = Array{String}(undef, length(paths))
	new_par_sys = Array{GasChromatographySimulator.Parameters}(undef, length(par_sys))
	visited_E = falses(length(E))

	# i -> path number
	# j -> segment/module number
	for i=1:length(paths)
		i_par = GasChromatographySystems.index_parameter(sys.g, paths[i])
		if GasChromatographySystems.path_possible(sys, p2fun, paths[i]; mode=mode) == true
			path_pos[i] = "path is possible"
			peaklists_ = Array{DataFrame}(undef, length(i_par))
			solutions_ = Array{Any}(undef, length(i_par))
			#As_ = Array{DataFrame}(undef, length(i_par))
			for j=1:length(i_par)
				if (i>1) && (all(visited_E[1:i_par[j]].==true))
					# was the segment already simulated in a previous simulated path?
					# look in all previous paths for the correct result -> the simulation correlated to the same edge and where this edge is connected to only previouse visited edges
					i_path = 0
					i_edge = 0
					for k=1:i-1 # previous paths
						i_par_previous = GasChromatographySystems.index_parameter(sys.g, paths[k])
						if length(i_par_previous) < j
							j0 = length(i_par_previous)
						else
							j0 = j
						end
						if all(x->x in i_par_previous[1:j0], i_par[1:j]) == true # all edges up to j are the same between the two paths
							i_path = k
							i_edge = findfirst(i_par[j].==i_par_previous)
						end
					end
					# re-use the results
					peaklists_[j] = peaklists[i_path][i_edge]
					solutions_[j] = solutions[i_path][i_edge]
				else # new simulated segments
					if j == 1 # first segment, directly after injection, it is assumed to be a segment of type `ModuleColumn`
						new_par_sys[j], peaklists_[j], solutions_[j] = simulate_ModuleColumn(par_sys[i_par[j]], t₀, τ₀)
					elseif typeof(sys.modules[i_par[j]]) == GasChromatographySystems.ModuleTM
						new_par_sys[j], peaklists_[j], solutions_[j] = simulate_ModuleTM(par_sys[i_par[j]], sys.modules[i_par[j]], peaklists_[j-1]; nτ=nτ, τ₀_focus=τ₀_focus, refocus=refocus, kwargsTM...)
					else # ModuleColumn
						new_par_sys[j], peaklists_[j], solutions_[j] = simulate_ModuleColumn(par_sys[i_par[j]], peaklists_[j-1])
					end
				end
			end
			visited_E[i_par] .= true
			peaklists[i] = peaklists_
			solutions[i] = solutions_
		else # path not possible
			neg_flow_modules = sys.modules[findall(paths[i][findall(positive_flow(sys, p2fun; mode=mode)[i_par].==false)].==edges(sys.g))]
			str_neg_flow = "path not possible: "
			for j=1:length(neg_flow_modules)
				str_neg_flow = str_neg_flow*"Flow in module >$(neg_flow_modules[j].name)< becomes negative during the program. "
			end
			path_pos[i] = str_neg_flow
		end
	end
	return path_pos, peaklists, solutions, new_par_sys
end