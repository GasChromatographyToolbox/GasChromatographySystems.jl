# estimate possible paths in the graphs

"""Map a vertex walk `vp` to the corresponding edge list (edge index order matches `collect(edges(g))`)."""
function _vertex_path_to_edges(g, vp)
	Eg = collect(edges(g))
	Ep_ = Array{Graphs.SimpleGraphs.SimpleEdge{Int64}, 1}(undef, length(vp) - 1)
	for i in 1:(length(vp) - 1)
		index = findfirst(vp[i] .== src.(Eg) .&& vp[i + 1] .== dst.(Eg))
		index === nothing && throw(ArgumentError("no edge from vertex $(vp[i]) to $(vp[i + 1])"))
		Ep_[i] = Eg[index]
	end
	return Ep_
end

"""
    path_is_chromatographic(g, modules, edge_path)

Return `true` if every edge in `edge_path` is a simulation segment (`ModuleColumn` or `ModuleTM`).
Paths containing a `ModuleValve` edge return `false`.
"""
function path_is_chromatographic(g, modules, edge_path)
	Eg = collect(edges(g))
	for e in edge_path
		idx = findfirst(==(e), Eg)
		idx === nothing && return false
		is_simulation_segment(modules[idx]) || return false
	end
	return true
end

function _filter_chromatographic_paths(g, modules, Vp, Ep)
	keep = [path_is_chromatographic(g, modules, Ep[j]) for j in eachindex(Ep)]
	return Vp[keep], Ep[keep]
end

"""
    all_paths(g, modules)

Enumerate vertex paths and edge paths by random walk (brute force, unique results).

Paths that include a `ModuleValve` edge are excluded (hydraulic-only edges; use
`path_is_chromatographic` to test a single path).

# Returns
- `Vp`: Vector of vertex sequences
- `Ep`: Vector of edge sequences (same length as `Vp`)
"""
function all_paths(g, modules::AbstractVector{<:GasChromatographySystems.AbstractModule})
	return all_paths(g, modules, nv(g) * 10)
end

function _collect_random_vertex_paths(g, num_paths::Integer)
	rand_paths = Any[]
	while length(unique(rand_paths)) < num_paths && length(rand_paths) < nv(g) * 10
		push!(rand_paths, non_backtracking_randomwalk(g, 1, nv(g)))
	end
	return sort(unique(rand_paths))
end

"""
    all_paths(g, modules, num_paths)

Like [`all_paths`](@ref)(`g`, `modules`), but stop once `num_paths` unique random walks are collected
(before excluding non-chromatographic paths).
"""
function all_paths(g, modules::AbstractVector{<:GasChromatographySystems.AbstractModule}, num_paths::Integer)
	Vp = _collect_random_vertex_paths(g, num_paths)
	Ep = [_vertex_path_to_edges(g, Vp[j]) for j in eachindex(Vp)]
	return _filter_chromatographic_paths(g, modules, Vp, Ep)
end

"""
    all_paths(sys, num_paths)

Convenience wrapper: [`all_paths`](@ref)(`sys.g`, `sys.modules`, `num_paths`).
"""
function all_paths(sys::GasChromatographySystems.System, num_paths::Integer)
	return all_paths(sys.g, sys.modules, num_paths)
end

"""
    all_paths(sys)

Convenience wrapper: [`all_paths`](@ref)(`sys.g`, `sys.modules`, `nv(sys.g) * 10`).
"""
function all_paths(sys::GasChromatographySystems.System)
	return all_paths(sys.g, sys.modules, nv(sys.g) * 10)
end

# simulate along the paths

"""
    index_parameter(g, path)

Map a path (sequence of graph edges) to module / parameter indices.

Each edge in `collect(edges(g))` corresponds to one entry in `sys.modules` and in the
`graph_to_parameters` output. Returns the indices of all edges that appear in `path`, in
graph edge order (not path order).

# Arguments
- `g`: Graph whose `edges(g)` define the module numbering.
- `path`: Vector of edges along a path (e.g. from [`all_paths`](@ref)).

# Returns
- `Vector{Int}`: Edge indices `i` with `collect(edges(g))[i] ∈ path`.
"""
function index_parameter(
	g::Graphs.AbstractGraph,
	path::AbstractVector{<:Graphs.AbstractEdge},
)::Vector{Int}
	E = collect(edges(g))
	return findall(e -> e in path, E)
end

"""
    common_edges(path1, path2)

Select edges from `path2` at positions where the corresponding edge in `path1` also occurs in `path2`.

For each index `i` in `path1`, if `path1[i]` is contained in `path2`, include `path2[i]` in the result.
(Used when comparing path prefixes; paths are usually aligned so `i` is meaningful for both.)

# Arguments
- `path1`: Reference edge sequence (typically a prefix).
- `path2`: Edge sequence to subselect.

# Returns
- Vector of edges from `path2` (same element type as `path2`).
"""
function common_edges(
	path1::AbstractVector{<:Graphs.AbstractEdge},
	path2::AbstractVector{<:Graphs.AbstractEdge},
)
	indices = findall(x -> x in path2, path1)
	return path2[indices]
end

"""
    positive_flow(sys, p2fun; mode="λ")

Test whether each edge has strictly positive volumetric flow over the program horizon.

Uses [`flow_functions`](@ref) and samples times `t = 0:Δt:t_end` with `Δt = t_end/2` and
`t_end = sum(common_timesteps(sys))`. An edge is marked `true` when `F(t) > 0` at all sample
points (no non-positive values). This is a coarse check, not interval arithmetic over `[0, t_end]`.

# Arguments
- `sys`: Capillary system after flow balance.
- `p2fun`: Squared-pressure solution functions from `build_pressure_squared_functions`.
- `mode`: `"λ"` or `"κ"` (passed to `flow_functions`).

# Returns
- `Vector{Bool}` of length `ne(sys.g)`; `true` if flow on that edge stays positive at the sample times.
"""
function positive_flow(
	sys::GasChromatographySystems.System,
	p2fun,
; mode::AbstractString = "λ",
)::Vector{Bool}
	F_func = flow_functions(sys, p2fun; mode=mode)
	tend = sum(common_timesteps(sys))
	t = 0:tend / 2:tend # !!!perhaps use interval arithmatic to test, if the flow is positive over the interval 0:tend!!!???
	pos_Flow = Vector{Bool}(undef, length(F_func))
	for i in eachindex(F_func)
		pos_Flow[i] = isempty(findall(F_func[i].(t) .<= 0))
	end
	return pos_Flow
end

"""
    path_possible(sys, p2fun, path; mode="λ")

Return whether solute transport along `path` is feasible (forward flow on every edge of the path).

Maps `path` to edge indices with [`index_parameter`](@ref) and requires
[`positive_flow`](@ref) to be `true` on each of those edges. Used by [`simulate_along_paths`](@ref)
to skip paths with backflush or zero-flow segments.

# Arguments
- `sys`: Capillary system.
- `p2fun`: Squared-pressure solution functions from `build_pressure_squared_functions`.
- `path`: Sequence of edges (one outlet path from [`all_paths`](@ref)).
- `mode`: `"λ"` or `"κ"` (passed to `positive_flow`).

# Returns
- `true` if every edge on `path` has positive flow at the sample times; `false` otherwise.
"""
function path_possible(
	sys::GasChromatographySystems.System,
	p2fun,
	path::AbstractVector{<:Graphs.AbstractEdge},
; mode::AbstractString = "λ",
)::Bool
	i_E = index_parameter(sys.g, path)
	pos = positive_flow(sys, p2fun; mode=mode)
	return length(i_E) == count(pos[i_E])
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

"""Index in `par.sub` for a peaklist row (match CAS + annotation when present)."""
function _sub_index_for_peaklist_row(par::GasChromatographySimulator.Parameters, cas, ann)
	ann_str = ann isa AbstractString ? ann : string(ann)
	if !isempty(ann_str)
		for k in eachindex(par.sub)
			if par.sub[k].CAS == cas && par.sub[k].ann == ann_str
				return k
			end
		end
	end
	return findfirst(==(cas), (s.CAS for s in par.sub))
end

function change_initial(par::GasChromatographySimulator.Parameters, prev_pl)
	# Drop non-finite carry-over peaks (failed upstream rows) so Substance construction
	# does not throw on t₀/τ₀ = NaN/Inf. This keeps simulate_along_paths running and
	# surfaces a clearer warning instead of an InexactError/ArgumentError deep in GCSim.
	finite_rows = findall(isfinite.(prev_pl.tR) .&& isfinite.(prev_pl.τR))
	if length(finite_rows) < nrow(prev_pl)
		@warn "change_initial: dropping $(nrow(prev_pl) - length(finite_rows)) non-finite peak rows before downstream simulation."
	end
	isempty(finite_rows) && throw(ArgumentError("change_initial: no finite peaks to pass downstream (all tR/τR are non-finite). Check upstream module simulation and flow/pressure settings."))
	has_ann = :Annotations in propertynames(prev_pl)
	new_sub = GasChromatographySimulator.Substance[]
	for i in finite_rows
		ann_i = has_ann ? prev_pl.Annotations[i] : ""
		i_sub = _sub_index_for_peaklist_row(par, prev_pl.CAS[i], ann_i)
		i_sub === nothing && continue
		# Injection time from peaklist. Default τ₀ = upstream simulated width (column→column, TM→column).
		# Use `par.sub.τ₀` only when this row matches a sliced substance (CAS + annotation),
		# e.g. after valve junction slicing — not on CAS-only fallback (downstream template `par`).
		τ_init = prev_pl.τR[i]
		if !isempty(ann_i) && par.sub[i_sub].ann == ann_i
			τ_init = par.sub[i_sub].τ₀
		end
		push!(new_sub, GasChromatographySimulator.Substance(
			par.sub[i_sub].name,
			par.sub[i_sub].CAS,
			par.sub[i_sub].Tchar,
			par.sub[i_sub].θchar,
			par.sub[i_sub].ΔCp,
			par.sub[i_sub].φ₀,
			prev_pl.Annotations[i],
			par.sub[i_sub].Cag,
			prev_pl.tR[i],
			τ_init,
		))
	end
	isempty(new_sub) && throw(ArgumentError("change_initial: no matching CAS entries found between previous peak list and downstream parameters."))
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

"""
    simulate_ModuleColumn(segment_par, t₀, τ₀)
    simulate_ModuleColumn(segment_par, prev_peaklist)

Simulate solute transport through a column module in a gas chromatography system.

This function simulates the behavior of solutes in a column module, either starting from initial
conditions or continuing from a previous module's results. It calculates retention times, peak
widths, and other chromatographic parameters.

# Arguments
- `segment_par`: Simulation parameters for the column module
- `t₀`: Initial time points (first variant)
- `τ₀`: Initial peak widths (first variant)
- `prev_peaklist`: Peak list from previous module (second variant)

# Returns
- Tuple containing:
  1. Updated simulation parameters
  2. Peak list with chromatographic parameters
  3. Solution trajectories

# Notes
- First variant is used for the initial column segment
- Second variant is used for subsequent column segments
- Peak areas are preserved between segments
- Initial areas are set to 1.0 for the first segment
- For valve-sliced peaks, pass the `Parameters` returned from
  [`apply_valve_junctions_at_vertex`](@ref): `t₀` from peaklist `tR`, `τ₀` from matching
  `par.sub` when CAS and annotation match (not peaklist `τR`). Otherwise upstream `τR` is used.
"""
function simulate_ModuleColumn(segment_par, t₀, τ₀)
	new_segment_par = GasChromatographySystems.change_initial(segment_par, t₀, τ₀)
	peaklist, solutions = GasChromatographySimulator.simulate(new_segment_par)
	peaklist[!,:A] = ones(length(peaklist.Name))
	return new_segment_par, peaklist, solutions
end

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
- Handles `ModuleColumn`, `ModuleTM`, and valve junction slicing at tee vertices (see `apply_valve_junctions_at_vertex`)
- Checks for negative flows to determine if paths are possible
- First segment after injection is assumed to be a `ModuleColumn`
"""
function simulate_along_paths(sys, p2fun, paths, par_sys; t₀=zeros(length(par_sys[1].sub)), τ₀=zeros(length(par_sys[1].sub)), nτ=6, refocus=falses(ne(sys.g)), τ₀_focus=zeros(length(par_sys[1].sub)), mode="λ", kwargsTM...)
	# -------------------------------------------------------------------------
	# Output buffers (one slot per path / per graph edge)
	# -------------------------------------------------------------------------
	E = collect(edges(sys.g))  # fixed edge numbering: index i ↔ sys.modules[i], par_sys[i]
	peaklists = Array{Array{DataFrame, 1}}(undef, length(paths))   # peaklists[i][j] = DataFrame after segment j on path i
	solutions = Array{Array{Any, 1}}(undef, length(paths))         # ODE / modulator solutions, same nesting
	path_pos = Array{String}(undef, length(paths))                 # human-readable status per path
	new_par_sys = Array{GasChromatographySimulator.Parameters}(undef, length(par_sys))  # updated Parameters per edge index

	# visited_E[e] = true once edge e has been simulated on some *feasible* path (enables reuse below)
	visited_E = falses(length(E))

	# Diagnostic helper: report where non-finite tR/τR first appears.
	function _log_nonfinite_peaklist(pl, i_path, j_seg, edge_idx, stage::AbstractString)
		if !(:tR in names(pl)) || !(:τR in names(pl))
			return
		end
		finite_mask = isfinite.(pl.tR) .&& isfinite.(pl.τR)
		n_bad = count(.!finite_mask)
		if n_bad > 0
			mod_name = sys.modules[edge_idx].name
			@warn "Non-finite peaks detected" path_index=i_path segment_index=j_seg edge_index=edge_idx module_name=mod_name stage=stage dropped_or_invalid=n_bad total_rows=nrow(pl)
		end
	end

	# -------------------------------------------------------------------------
	# Outer loop: each candidate outlet path (sequence of graph edges)
	# -------------------------------------------------------------------------
	for i in 1:length(paths)
		# i_par: global edge indices on this path (order follows collect(edges(g)), not necessarily
		#        the order solute travels). Used to index par_sys, sys.modules, and reuse logic.
		i_par = GasChromatographySystems.index_parameter(sys.g, paths[i])

		# path_edge_idx: same edges as i_par but in **path order** (injection → detector).
		# Needed to locate junction vertices between consecutive column segments (DPM tee, etc.).
		path_edge_idx = GasChromatographySystems.edges_along_path_in_order(sys.g, paths[i])

		# ---------------------------------------------------------------------
		# Feasibility: forward flow on every edge of the path over the program horizon
		# ---------------------------------------------------------------------
		if GasChromatographySystems.path_possible(sys, p2fun, paths[i]; mode=mode) == true
			path_pos[i] = "path is possible"

			# Per-path working arrays; j indexes position along this path's module chain
			peaklists_ = Array{DataFrame}(undef, length(i_par))
			solutions_ = Array{Any}(undef, length(i_par))

			# -----------------------------------------------------------------
			# Inner loop: simulate or reuse each segment along the path
			# -----------------------------------------------------------------
			for j in 1:length(i_par)
				# -------------------------------------------------------------
				# Branch A — Reuse: segment already computed on an earlier path
				# -------------------------------------------------------------
				# Conditions:
				#   (1) We are not on the first path (i > 1).
				#   (2) Every graph edge up to and including the current one was already
				#       visited on a previous *successful* path (prefix of the network
				#       is fully covered).
				# Then search paths 1:(i-1) for one whose first j segments share the
				# same edge indices; copy its peak list and solution for segment j.
				if (i > 1) && (all(visited_E[1:i_par[j]].==true))
					i_path = 0   # index of donor path in `paths`
					i_edge = 0   # segment index j on that donor path
					for k in 1:(i - 1)
						i_par_previous = GasChromatographySystems.index_parameter(sys.g, paths[k])
						# Compare only the overlapping prefix if the donor path is shorter
						j0 = min(length(i_par_previous), j)
						if all(x -> x in i_par_previous[1:j0], i_par[1:j])
							i_path = k
							i_edge = findfirst(==(i_par[j]), i_par_previous)
						end
					end
					peaklists_[j] = peaklists[i_path][i_edge]
					solutions_[j] = solutions[i_path][i_edge]

				# -------------------------------------------------------------
				# Branch B — Fresh simulation for this edge / segment
				# -------------------------------------------------------------
				else
					if j == 1
						# First module after injection: initial band conditions (t₀, τ₀).
						# Convention: first edge on a chromatographic path is a ModuleColumn.
						new_par_sys[i_par[j]], peaklists_[j], solutions_[j] =
							simulate_ModuleColumn(par_sys[i_par[j]], t₀, τ₀)
					else
						# Carry the peak list produced by the previous segment on this path.
						pl_in = peaklists_[j - 1]
						_log_nonfinite_peaklist(pl_in, i, j - 1, i_par[j - 1], "input from previous segment")

						# Valve junction (e.g. FastGC×GC DPM tee): path is column-only, but a
						# ModuleValve may be attached at the vertex between two columns. If the
						# valve state varies in time, split peaks by period before the downstream
						# segment (analogous to TM slicing, see ValveJunction.jl).
						p_pos = findfirst(==(i_par[j]), path_edge_idx)
						par_in = par_sys[i_par[j]]
						if p_pos !== nothing && p_pos > 1
							# Junction vertex = destination of upstream path edge
							v_junc = dst(paths[i][p_pos - 1])
							par_in, pl_in, _ = GasChromatographySystems.apply_valve_junctions_at_vertex(
								sys, v_junc, par_in, pl_in; nτ=nτ,
							)
						end

						mod_j = sys.modules[i_par[j]]
						if mod_j isa GasChromatographySystems.ModuleTM
							# Thermal modulator on the path: slice by PM, then modulate (simplifiedTM / ODE).
							new_par_sys[i_par[j]], peaklists_[j], solutions_[j] =
								simulate_ModuleTM(
									par_sys[i_par[j]], mod_j, pl_in;
									nτ=nτ, τ₀_focus=τ₀_focus, refocus=refocus, kwargsTM...,
								)
							_log_nonfinite_peaklist(peaklists_[j], i, j, i_par[j], "output of simulate_ModuleTM")
						else
							# Ordinary column: continue with (possibly valve-sliced) peak list.
							try
								new_par_sys[i_par[j]], peaklists_[j], solutions_[j] =
									simulate_ModuleColumn(par_in, pl_in)
							catch err
								if err isa ArgumentError && occursin("no finite peaks to pass downstream", sprint(showerror, err))
									mod_name = sys.modules[i_par[j]].name
									@error "All peaks non-finite before downstream column simulation" path_index=i segment_index=j edge_index=i_par[j] module_name=mod_name
									@error "Upstream segment context" upstream_segment_index=j-1 upstream_edge_index=i_par[j-1] upstream_module_name=sys.modules[i_par[j-1]].name
									_log_nonfinite_peaklist(pl_in, i, j - 1, i_par[j - 1], "failing input to simulate_ModuleColumn")
								end
								rethrow(err)
							end
							_log_nonfinite_peaklist(peaklists_[j], i, j, i_par[j], "output of simulate_ModuleColumn")
						end
					end
				end
			end

			# Mark all edges on this path as available for reuse by later paths
			visited_E[i_par] .= true
			peaklists[i] = peaklists_
			solutions[i] = solutions_

		# ---------------------------------------------------------------------
		# Path infeasible: backflush or zero flow on at least one path edge
		# ---------------------------------------------------------------------
		else
			# Build a diagnostic message listing modules where F(t) ≤ 0 was detected
			neg_flow_modules = sys.modules[findall(
				paths[i][findall(positive_flow(sys, p2fun; mode=mode)[i_par].==false)].==edges(sys.g),
			)]
			str_neg_flow = "path not possible: "
			for j in 1:length(neg_flow_modules)
				str_neg_flow *= "Flow in module >$(neg_flow_modules[j].name)< becomes negative during the program. "
			end
			path_pos[i] = str_neg_flow
		end
	end
	return path_pos, peaklists, solutions, new_par_sys
end