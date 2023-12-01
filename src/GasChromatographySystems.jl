module GasChromatographySystems

using Graphs
using GraphMakie
using Symbolics
using NetworkLayout
using CairoMakie
using GasChromatographySimulator
using Interpolations
using CSV
using DataFrames
using Plots
using SpecialFunctions
using OrdinaryDiffEq
using Integrals
using LsqFit

# some constants
const Tst = 273.15 # K
const R = 8.31446261815324 # J mol⁻¹ K⁻¹
const Tn = 25.0 + Tst # K
const pn = 101300 # Pa

# Definition of structures and methods
include("./Structures.jl")
include("./Flowcalc.jl")

# functions

# begin - misc, for update_system
# common programs
function common_timesteps(sys)
	com_timesteps = []
	for i=1:nv(sys.g)
		com_timesteps = GasChromatographySimulator.common_time_steps(com_timesteps, sys.pressurepoints[i].time_steps)
	end
	for i=1:ne(sys.g)
		if typeof(sys.modules[i].T) <: TemperatureProgram
			com_timesteps = GasChromatographySimulator.common_time_steps(com_timesteps, sys.modules[i].T.time_steps)
		end
	end
	return com_timesteps
end

function index_modules_with_temperature_program(sys)
	i_tempprog = Int[]
	for i=1:ne(sys.g)
		if typeof(sys.modules[i].T) <: TemperatureProgram
			push!(i_tempprog, i)
		end
	end
	return i_tempprog
end

function match_programs(sys)
	com_times = common_timesteps(sys)
	new_press_steps = Array{Array{Float64,1}}(undef, nv(sys.g))
	for i=1:nv(sys.g)
		new_press_steps[i] = GasChromatographySimulator.new_value_steps(sys.pressurepoints[i].pressure_steps, sys.pressurepoints[i].time_steps, com_times)
	end
	i_tempprog = index_modules_with_temperature_program(sys)
	new_temp_steps = Array{Array{Float64,1}}(undef, length(i_tempprog))
	for i=1:length(i_tempprog)
		new_temp_steps[i] = GasChromatographySimulator.new_value_steps(sys.modules[i_tempprog[i]].T.temp_steps, sys.modules[i_tempprog[i]].T.time_steps, com_times)
	end
	# add for gradient new_a_gf
	return com_times, new_press_steps, new_temp_steps, i_tempprog
end

function update_system(sys)
	new_timesteps, new_pressuresteps, new_temperaturesteps, index_module_tempprog = match_programs(sys)
	new_pp = Array{PressurePoint}(undef, nv(sys.g))
	for i=1:nv(sys.g)
		new_pp[i] = PressurePoint(sys.pressurepoints[i].name, new_timesteps, new_pressuresteps[i])
	end
	new_modules = Array{AbstractModule}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		if typeof(sys.modules[i].T) <: Number
			new_modules[i] = sys.modules[i]
		elseif typeof(sys.modules[i].T) <: TemperatureProgram
			# add/modify for gradient
			ii = findfirst(index_module_tempprog.==i)
			new_tp = TemperatureProgram(new_timesteps, new_temperaturesteps[ii])
			if typeof(sys.modules[i]) == ModuleColumn
				new_modules[i] = ModuleColumn(sys.modules[i].name, sys.modules[i].L, sys.modules[i].d, sys.modules[i].df, sys.modules[i].sp, new_tp, sys.modules[i].F, sys.modules[i].opt)
			elseif typeof(sys.modules[i]) == ModuleTM
				new_modules[i] = ModuleTM(sys.modules[i].name, sys.modules[i].L, sys.modules[i].d, sys.modules[i].df, sys.modules[i].sp, new_tp, sys.modules[i].shift, sys.modules[i].PM, sys.modules[i].ratio, sys.modules[i].Thot, sys.modules[i].Tcold, sys.modules[i].F, sys.modules[i].opt)
			end
		end
	end
	new_sys = System(sys.g, new_pp, new_modules, sys.options)
    return new_sys
end
# end - misc, for update_system


# begin - pressure and flow calculations
# flow balance and pressure calculations
function inlet_vertices(g)
	v = vertices(g)
	inlet = Int[]
	for i=1:length(v)
		if isempty(inneighbors(g, v[i]))
			push!(inlet, i)
		end
	end
	return inlet
end

function outlet_vertices(g)
	v = vertices(g)
	outlet = Int[]
	for i=1:length(v)
		if isempty(outneighbors(g, v[i]))
			push!(outlet, i)
		end
	end
	return outlet
end

function inner_vertices(g)
	v = vertices(g)
	inner = Int[]
	for i=1:length(v)
		if !isempty(outneighbors(g, v[i])) && !isempty(inneighbors(g, v[i]))
			push!(inner, i)
		end
	end
	return inner
end

function module_temperature(module_::ModuleColumn, sys)
	L = if isnan(module_.L)
		1.0
	else
		module_.L
	end
	if typeof(module_.T) <: TemperatureProgram # temperature is a TemperatureProgram
		time_steps = module_.T.time_steps
		temp_steps = module_.T.temp_steps	
		gf = module_.T.gf
		a_gf = module_.T.a_gf
		T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, L)
	elseif typeof(module_.T) <: Number # temperature is a constant value
		time_steps = common_timesteps(sys)
		temp_steps = module_.T.*ones(length(time_steps))
		gf(x) = zero(x).*ones(length(time_steps))
		a_gf = [zeros(length(time_steps)) zeros(length(time_steps)) ones(length(time_steps)) zeros(length(time_steps))]
		T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, L)
	end
	return time_steps, temp_steps, gf, a_gf, T_itp
end

function module_temperature(module_::GasChromatographySystems.ModuleTM, sys)
	if typeof(module_.T) <: GasChromatographySystems.TemperatureProgram # temperature is a TemperatureProgram
		time_steps = module_.T.time_steps
		temp_steps = module_.T.temp_steps
		gf = module_.T.gf
		a_gf = module_.T.a_gf
		T_itp_ = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, module_.L)
	elseif typeof(module_.T) <: Number # temperature is a constant value
		time_steps = GasChromatographySystems.common_timesteps(sys)
		temp_steps = module_.T.*ones(length(time_steps))
		gf(x) = zero(x).*ones(length(time_steps))
		a_gf = [zeros(length(time_steps)) zeros(length(time_steps)) ones(length(time_steps)) zeros(length(time_steps))]
		T_itp_ = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, module_.L)
	end
	T_itp(x,t) = if module_.opt.Tcold_abs == true # cool jet always at Tcold
		therm_mod(t, module_.shift, module_.PM, module_.ratio, module_.Tcold, T_itp_(x, t) .+ module_.Thot .- 273.15; flank=module_.opt.tflank) .+ 273.15 
	else # cool jet always at T_itp_ + Tcold
		therm_mod(t, module_.shift, module_.PM, module_.ratio, T_itp_(x, t) .+ module_.Tcold .- 273.15, T_itp_(x, t) .+ module_.Thot .- 273.15; flank=module_.opt.tflank) .+ 273.15 
	end
	
	spot(x,t) = if module_.opt.ng == false
		GasChromatographySystems.smooth_rectangle(x, 0.0, sys.modules[5].L, T_itp_(x, t), T_itp(x,t); flank=module_.opt.sflank)
	else
		T_itp(x,t)
	end
	return time_steps, temp_steps, gf, a_gf, spot
end

"""
    flow_restrictions(sys)

Calculates the flow restrictions κ of all edges (capliaries) of a system of capillaries.

# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
"""
function flow_restrictions(sys)
	kappas = Array{Function}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		T_itp = module_temperature(sys.modules[i], sys)[5]
		κ(t) = GasChromatographySimulator.flow_restriction(sys.modules[i].L, t, T_itp, sys.modules[i].d, sys.options.gas; ng=sys.modules[i].opt.ng, vis=sys.options.vis)
		kappas[i] = κ
	end
	return kappas
end

"""
    flow_permeabilities(sys)

Calculates the flow permeabilities λ of all edges (capliaries) of a system of capillaries.

# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
"""
function flow_permeabilities(sys)
	lambdas = Array{Function}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		T_itp = module_temperature(sys.modules[i], sys)[5]
		λ(t) = 1/GasChromatographySimulator.flow_restriction(sys.modules[i].L, t, T_itp, sys.modules[i].d, sys.options.gas; ng=sys.modules[i].opt.ng, vis=sys.options.vis)
		lambdas[i] = λ
	end
	return lambdas
end

function pressures_squared(sys)
	#p² = Array{Interpolations.Extrapolation}(undef, nv(sys.g))
	p² = Array{Any}(undef, nv(sys.g))
	for i=1:nv(sys.g)
		g = if isnan(sys.pressurepoints[i].pressure_steps[1])
			GasChromatographySimulator.steps_interpolation(sys.pressurepoints[i].time_steps, NaN.*ones(length(sys.pressurepoints[i].time_steps)))
		else
			GasChromatographySimulator.steps_interpolation(sys.pressurepoints[i].time_steps, identity.(sys.pressurepoints[i].pressure_steps))
		end
		f(t) = g(t).^2 
		p²[i] = f
	end
	return p²
end

# begin - plotting of the graphs and functions
function plot_graph(g, edge_labels, node_labels; lay = Spring(), color=:lightblue, node_size=40, arrow_size=20, arrow_shift=0.8, dataaspect=false, nlabels_fontsize=10, elabels_fontsize=10)
	fig, ax, p = GraphMakie.graphplot(g, 
						layout=lay,
						nlabels=node_labels,
						nlabels_fontsize = nlabels_fontsize, 
						nlabels_align=(:center,:center),
						node_size = [node_size for i in 1:Graphs.nv(g)],
						node_color = [color for i in 1:Graphs.nv(g)],
						elabels = edge_labels,
						elabels_fontsize = elabels_fontsize,
						arrow_size = arrow_size,
						arrow_shift = arrow_shift
					)
	hidedecorations!(ax)
	hidespines!(ax)
	if dataaspect == true
		ax.aspect = DataAspect()
	end
	return fig
end

function plot_graph(sys; lay = Spring(), color=:lightblue, node_size=40, arrow_size=20, arrow_shift=0.8, dataaspect=false, nlabels_fontsize=20, elabels_fontsize=20)
	fig, ax, p = GraphMakie.graphplot(sys.g, 
						layout=lay,
						nlabels = [sys.pressurepoints[i].name for i in 1:Graphs.nv(sys.g)],
						nlabels_align=(:center,:center),
						nlabels_fontsize = nlabels_fontsize,
						node_size = [node_size for i in 1:Graphs.nv(sys.g)],
						node_color = [color for i in 1:Graphs.nv(sys.g)],
						elabels = [sys.modules[i].name for i in 1:Graphs.ne(sys.g)],
						elabels_fontsize = elabels_fontsize,
						arrow_size = arrow_size,
						arrow_shift = arrow_shift
					)
	hidedecorations!(ax)
	hidespines!(ax)
	if dataaspect == true
		ax.aspect = DataAspect()
	end
	return fig
end

function plot_graph_with_flow(sys, t; lay = Spring(), color=:lightblue, node_size=80, arrow_size=20, arrow_shift=0.8, dataaspect=false, nlabels_fontsize=14, elabels_fontsize=14, elabels_distance = 20)
	p_func = pressure_functions(sys)
	F_func = flow_functions(sys)
	fig, ax, p = GraphMakie.graphplot(sys.g, 
						layout=lay,
						nlabels = [sys.pressurepoints[i].name*"\n$(trunc(Int, p_func[i](t)/1000))kPa" for i in 1:nv(sys.g)],
						nlabels_align=(:center,:center),
						nlabels_fontsize = nlabels_fontsize,
						node_size = [node_size for i in 1:nv(sys.g)],
						node_color = [color for i in 1:nv(sys.g)],
						elabels = [sys.modules[i].name*"\n $(round(F_func[i](t)*60*1e6; sigdigits=2)) mL/min" for i in 1:ne(sys.g)],
						elabels_distance = elabels_distance,
						elabels_fontsize = elabels_fontsize,
						arrow_size = arrow_size,
						arrow_shift = arrow_shift
					)
	hidedecorations!(ax)
	hidespines!(ax)
	if dataaspect == true
		ax.aspect = DataAspect()
	end
	return fig
end

function plot_flow_over_time(sys; dt=60.0)
	#plotly()
	F_func = flow_functions(sys)
	com_timesteps = GasChromatographySystems.common_timesteps(sys)
	trange = 0:dt:sum(com_timesteps)
	p_flow = Plots.plot(xlabel="time in s", ylabel="flow in mL/min")
	for i=1:ne(sys.g)
		Plots.plot!(p_flow, trange, F_func[i].(trange).*60e6, label="F_$(sys.modules[i].name)")
	end
	return p_flow
end

function plot_pressure_over_time(sys; dt=60.0)
	#plotly()
	p_func = pressure_functions(sys)
	com_timesteps = GasChromatographySystems.common_timesteps(sys)
	trange = 0:dt:sum(com_timesteps)
	p_pres = Plots.plot(xlabel="time in s", ylabel="pressure in Pa(a)", legend=:topleft)
	for i=1:nv(sys.g)
		Plots.plot!(p_pres, trange, p_func[i].(trange), label="$(sys.pressurepoints[i].name)")
	end
	return p_pres
end

# transform system to GasChromatographySimulator.Parameters
function all_stationary_phases(sys)
	stat_phases = Array{String}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		stat_phases[i] = sys.modules[i].sp
	end
	return stat_phases
end
# end - plotting of the graphs and functions

# begin - system to parameters
function common_solutes(db, sys)
	# gives the soultes with common stationary phases between the database 'db' and the
	# GC-System 'GCsys'
	usp = setdiff(unique(GasChromatographySystems.all_stationary_phases(sys)), [""])
	if length(usp)==0 # no stationary phase
		common_solutes = DataFrame(Name=db.Name, CAS=db.CAS)
	else
		filter_db = Array{DataFrame}(undef, length(usp))
		for i=1:length(usp)
			filter_db[i] = filter([:Phase] => x -> x==usp[i], db)
		end
		if length(usp)==1 # only one stationary phase
			common_db = filter_db[1]
		else
			common_db = innerjoin(filter_db[1], filter_db[2], on=:CAS, makeunique=true)
			if length(usp)>2 # more than two stationary phases
				for i=3:length(usp)
				common_db = innerjoin(common_db, filter_db[i], on=:CAS, makeunique=true)
				end
			end
		end
		CAS = unique(common_db.CAS)
		Name = Array{String}(undef, length(CAS))
		for i=1:length(CAS)
			ii = findfirst(common_db.CAS.==CAS[i])
			Name[i] = common_db.Name[ii]
		end
		common_solutes = DataFrame(Name=Name, CAS=CAS)
	end	
	return common_solutes
end

function graph_to_parameters(sys, db_dataframe, selected_solutes; interp=true, dt=1)
	# ng should be taken from the separat module options 
	E = collect(edges(sys.g))
	srcE = src.(E) # source indices
	dstE = dst.(E) # destination indices
	if interp == true # linear interpolation of pressure functions with step width dt
		p_func = interpolate_pressure_functions(sys; dt=dt)
	else
		p_func = pressure_functions(sys)
	end
	parameters = Array{GasChromatographySimulator.Parameters}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		# column parameters
		col = GasChromatographySimulator.Column(sys.modules[i].L, sys.modules[i].d, [sys.modules[i].d], sys.modules[i].df, [sys.modules[i].df], sys.modules[i].sp, sys.options.gas)

		# program parameters
		pin_steps = sys.pressurepoints[srcE[i]].pressure_steps
		pout_steps = sys.pressurepoints[dstE[i]].pressure_steps
		pin_itp = p_func[srcE[i]]
		pout_itp = p_func[dstE[i]]	
		time_steps, temp_steps, gf, a_gf, T_itp = module_temperature(sys.modules[i], sys)
		prog = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)

		# substance parameters
		sub = GasChromatographySimulator.load_solute_database(db_dataframe, sys.modules[i].sp, sys.options.gas, selected_solutes, NaN.*ones(length(selected_solutes)), NaN.*ones(length(selected_solutes)))

		# option parameters
		#opt = if typeof(sys.modules[i]) == GasChromatographySystems.ModuleTM 
		opt = GasChromatographySimulator.Options(alg=sys.modules[i].opt.alg, abstol=sys.modules[i].opt.abstol, reltol=sys.modules[i].opt.reltol, Tcontrol=sys.modules[i].opt.Tcontrol, odesys=sys.options.odesys, ng=sys.modules[i].opt.ng, vis=sys.options.vis, control=sys.options.control, k_th=sys.options.k_th)
		#else
			# put options for ModuleColumn in separat options structure?
		#	GasChromatographySimulator.Options(alg=sys.modules[i].opt.alg, abstol=sys.modules[i].opt.abstol, reltol=sys.modules[i].opt.reltol, Tcontrol=sys.modules[i].opt.Tcontrol, odesys=sys.options.odesys, ng=sys.modules[i].opt.ng, vis=sys.options.vis, control=sys.options.control, k_th=sys.options.k_th)
		#end

		parameters[i] = GasChromatographySimulator.Parameters(col, prog, sub, opt)
	end
	return parameters
end
# end - system to parameters

# begin - simulation of system
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

function positive_flow(sys)
	F_func = flow_functions(sys)
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

function path_possible(sys, paths)
	#F_func = flow_functions(sys)
	i_E = index_parameter(sys.g, paths)
	if length(i_E) == length(findall(positive_flow(sys)[i_E]))
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

function simulate_along_paths(sys, paths, db_dataframe, selected_solutes; t₀=zeros(length(selected_solutes)), τ₀=zeros(length(selected_solutes)))
	par_sys = graph_to_parameters(sys, db_dataframe, selected_solutes)
	path_pos, peaklists, solutions, new_par_sys = simulate_along_paths(sys, paths, par_sys; t₀=t₀, τ₀=τ₀)
	return path_pos, peaklists, solutions, new_par_sys
end

function simulate_along_paths(sys, paths, par_sys; t₀=zeros(length(par_sys[1].sub)), τ₀=zeros(length(par_sys[1].sub)), nτ=6, refocus=falses(ne(sys.g)), τ₀_focus=zeros(length(par_sys[1].sub)), kwargsTM...)
	
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
		if GasChromatographySystems.path_possible(sys, paths[i]) == true
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
						new_par_sys[i_par[j]] = GasChromatographySystems.change_initial(par_sys[i_par[j]], t₀, τ₀)
						peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(new_par_sys[i_par[j]])
						peaklists_[j][!,:A] = ones(length(peaklists_[j].Name)) # add a relativ area factor, splitting of `A` at split points is not accounted for yet, is only used for the slicing of peaks at modulators
					elseif typeof(sys.modules[i_par[j]]) == GasChromatographySystems.ModuleTM
						# put this in a separate function
						if refocus[i_par[j]] == true
							τ₀=τ₀_focus
						else
							τ₀=peaklists_[j-1].τR # PM?
						end
						new_par_sys[i_par[j]], df_A = GasChromatographySystems.slicing(peaklists_[j-1], sys.modules[i_par[j]].PM, sys.modules[i_par[j]].ratio, sys.modules[i_par[j]].shift, par_sys[i_par[j]]; nτ=nτ, τ₀=τ₀, abstol=sys.modules[i_par[j]].opt.abstol, reltol=sys.modules[i_par[j]].opt.reltol, alg=sys.modules[i_par[j]].opt.alg)
						if sys.modules[i_par[j]].opt.alg == "simplifiedTM"
							peaklists_[j], solutions_[j] = approximate_modulator(sys.modules[i_par[j]].T, new_par_sys[i_par[j]], df_A, sys.modules[i_par[j]].PM, sys.modules[i_par[j]].ratio, sys.modules[i_par[j]].shift, sys.modules[i_par[j]].Thot)
						else
							sol = Array{Any}(undef, length(new_par_sys[i_par[j]].sub))
							for i_sub=1:length(new_par_sys[i_par[j]].sub)
								dt = sys.modules[i_par[j]].opt.dtinit
								sol[i_sub] = GasChromatographySimulator.solving_odesystem_r(new_par_sys[i_par[j]].col, new_par_sys[i_par[j]].prog, new_par_sys[i_par[j]].sub[i_sub], new_par_sys[i_par[j]].opt; kwargsTM..., dt=dt)
								tR = sol[i_sub].u[end][1]
								while (fld(tR + sys.modules[i_par[j]].shift, sys.modules[i_par[j]].PM) > fld(new_par_sys[i_par[j]].sub[i_sub].t₀ + sys.modules[i_par[j]].shift, sys.modules[i_par[j]].PM) && dt > eps()) || (sol[i_sub].retcode != ReturnCode.Success && dt > eps()) # solute elutes not in the same modulation periode or the solving failed
									dt = dt/10 # reduce initial step-width
									dtmax = dt*1000
									#@warn "Retention time $(tR) surpasses modulation period, t₀ = $(new_par_sys[i_par[j]].sub[i_sub].t₀), PM = $(sys.modules[i_par[j]].PM). Initial step-width dt is decreased ($(dt))."
									sol[i_sub] = GasChromatographySimulator.solving_odesystem_r(new_par_sys[i_par[j]].col, new_par_sys[i_par[j]].prog, new_par_sys[i_par[j]].sub[i_sub], new_par_sys[i_par[j]].opt; kwargsTM..., dt=dt, dtmax=dtmax)
									tR = sol[i_sub].u[end][1]
								end
								#if isnan(tR) # does not work
								#	@warn "Simulation result is NaN. Approximate modulator."
								#	sol[i_sub] = approximate_modulator(new_par_sys[i_par[j]], df_A, sys.modules[i_par[j]].PM, sys.modules[i_par[j]].ratio, sys.modules[i_par[j]].shift)[2][i_sub]
								#end
							end
							peaklists_[j] = GasChromatographySimulator.peaklist(sol, new_par_sys[i_par[j]])
							GasChromatographySystems.add_A_to_pl!(peaklists_[j], df_A)
							solutions_[j] = sol
						end
						if maximum(peaklists_[j].τR) > sys.modules[i_par[j]].PM
							return @warn "Peak width of focussed peaks $(maximum(peaklists_[j].τR)) > modulation period $(sys.modules[i_par[j]].PM). Simulation is aborted. alg=$(sys.modules[i_par[j]].opt.alg), T=$(sys.modules[i_par[j]].T), peaklist=$(peaklists_[j])."
						end
					else # ModuleColumn
						new_par_sys[i_par[j]] = GasChromatographySystems.change_initial(par_sys[i_par[j]], peaklists_[j-1])
						peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(new_par_sys[i_par[j]])
						GasChromatographySystems.add_A_to_pl!(peaklists_[j], peaklists_[j-1])
					end
				end
			end
			visited_E[i_par] .= true
			peaklists[i] = peaklists_
			solutions[i] = solutions_
		else # path not possible
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

function simulate_along_one_path(sys, path, par_sys; t₀=zeros(length(par_sys[1].sub)), τ₀=zeros(length(par_sys[1].sub)), nτ=6, refocus=falses(length(par_sys[1].sub)), τ₀_focus=zeros(length(par_sys[1].sub)), pl_thread=true, kwargsTM...)
	
	E = collect(edges(sys.g))
#	peaklists = Array{Array{DataFrame,1}}(undef, length(paths))
#	solutions = Array{Array{Any,1}}(undef, length(paths))
#	path_pos = Array{String}(undef, length(paths))
	new_par_sys = Array{GasChromatographySimulator.Parameters}(undef, length(par_sys))
#	visited_E = falses(length(E))

	# i -> path number
	# j -> segment/module number
#	for i=1:length(paths)
		i_par = index_parameter(sys.g, path)
		if path_possible(sys, path) == true
			path_pos = "path is possible"
			peaklists_ = Array{DataFrame}(undef, length(i_par))
			solutions_ = Array{Any}(undef, length(i_par))
			#As_ = Array{DataFrame}(undef, length(i_par))
			for j=1:length(i_par)
#				if (i>1) && (all(visited_E[1:i_par[j]].==true))
#					# was the segment already simulated in a previous simulated path?
#					# look in all previous paths for the correct result -> the simulation correlated to the same edge and where this edge is connected to only previouse visited edges
#					i_path = 0
#					i_edge = 0
#					for k=1:i-1 # previous paths
#						i_par_previous = GasChromatographySystems.index_parameter(sys.g, paths[k])
#						if length(i_par_previous) < j
#							j0 = length(i_par_previous)
#						else
#							j0 = j
#						end
#						if all(x->x in i_par_previous[1:j0], i_par[1:j]) == true # all edges up to j are the same between the two paths
#							i_path = k
#							i_edge = findfirst(i_par[j].==i_par_previous)
#						end
#					end
#					# re-use the results
#					peaklists_[j] = peaklists[i_path][i_edge]
#					solutions_[j] = solutions[i_path][i_edge]
#				else # new simulated segments
					if j == 1 # first segment, directly after injection, it is assumed to be a segment of type `ModuleColumn`
						new_par_sys[i_par[j]] = change_initial(par_sys[i_par[j]], t₀, τ₀)
						
						
						T_itp = new_par_sys[i_par[j]].prog.T_itp
						Fpin_itp = new_par_sys[i_par[j]].prog.Fpin_itp
						pout_itp = new_par_sys[i_par[j]].prog.pout_itp
						L = new_par_sys[i_par[j]].col.L
						d = new_par_sys[i_par[j]].col.d
						df = new_par_sys[i_par[j]].col.df
						Tchars = [new_par_sys[i_par[j]].sub[x].Tchar for x in 1:length(new_par_sys[i_par[j]].sub)]
						θchars = [new_par_sys[i_par[j]].sub[x].θchar for x in 1:length(new_par_sys[i_par[j]].sub)]
						ΔCps = [new_par_sys[i_par[j]].sub[x].ΔCp for x in 1:length(new_par_sys[i_par[j]].sub)]
						φ₀s = [new_par_sys[i_par[j]].sub[x].φ₀ for x in 1:length(new_par_sys[i_par[j]].sub)]
						Cags = [new_par_sys[i_par[j]].sub[x].Cag for x in 1:length(new_par_sys[i_par[j]].sub)]
						t₀s = [new_par_sys[i_par[j]].sub[x].t₀ for x in 1:length(new_par_sys[i_par[j]].sub)]
						τ₀s = [new_par_sys[i_par[j]].sub[x].τ₀ for x in 1:length(new_par_sys[i_par[j]].sub)]
						gas = new_par_sys[i_par[j]].col.gas
						opt = new_par_sys[i_par[j]].opt
						
						#peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(T_itp, Fpin_itp, pout_itp, L, d, df, Tchars, θchars, ΔCps, φ₀s, Cags, t₀s, τ₀s, gas, opt)
						
						sol_ = Array{Any}(undef, length(new_par_sys[i_par[j]].sub))
						for i=1:length(new_par_sys[i_par[j]].sub)
							sol_[i] = GasChromatographySimulator.solving_odesystem_r(L, d, df, T_itp, Fpin_itp, pout_itp, Tchars[i], θchars[i], ΔCps[i], φ₀s[i], Cags[i], t₀s[i], τ₀s[i], gas, opt)
						end
						peaklists_[j] = GasChromatographySimulator.peaklist(sol_, new_par_sys[i_par[j]]; thread=pl_thread)
						peaklists_[j][!,:A] = ones(length(peaklists_[j].Name)) # add a relativ area factor, splitting of `A` at split points is not accounted for yet, is only used for the slicing of peaks at modulators
						solutions_[j] = sol_
					elseif typeof(sys.modules[i_par[j]]) == ModuleTM
						# put this in a separate function
						if refocus[i_par[j]] == true
							τ₀=τ₀_focus
						else
							τ₀=peaklists_[j-1].τR # PM?
						end
						new_par_sys[i_par[j]], df_A = slicing(peaklists_[j-1], sys.modules[i_par[j]].PM, sys.modules[i_par[j]].ratio, sys.modules[i_par[j]].shift, par_sys[i_par[j]]; nτ=nτ, τ₀=τ₀, abstol=sys.modules[i_par[j]].opt.abstol, reltol=sys.modules[i_par[j]].opt.reltol, alg=sys.modules[i_par[j]].opt.alg)
						if sys.modules[i_par[j]].opt.alg == "simplifiedTM"
							peaklists_[j], solutions_[j] = approximate_modulator(sys.modules[i_par[j]].T, new_par_sys[i_par[j]], df_A, sys.modules[i_par[j]].PM, sys.modules[i_par[j]].ratio, sys.modules[i_par[j]].shift, sys.modules[i_par[j]].Thot)
						else
							sol = Array{Any}(undef, length(new_par_sys[i_par[j]].sub))
							for i_sub=1:length(new_par_sys[i_par[j]].sub)
								dt = sys.modules[i_par[j]].opt.dtinit
								sol[i_sub] = GasChromatographySimulator.solving_odesystem_r(new_par_sys[i_par[j]].col, new_par_sys[i_par[j]].prog, new_par_sys[i_par[j]].sub[i_sub], new_par_sys[i_par[j]].opt; kwargsTM..., dt=dt)
								tR = sol[i_sub].u[end][1]
								while (fld(tR + sys.modules[i_par[j]].shift, sys.modules[i_par[j]].PM) > fld(new_par_sys[i_par[j]].sub[i_sub].t₀ + sys.modules[i_par[j]].shift, sys.modules[i_par[j]].PM) && dt > eps()) || (sol[i_sub].retcode != ReturnCode.Success && dt > eps()) # solute elutes not in the same modulation periode or the solving failed
									dt = dt/10 # reduce initial step-width
									dtmax = dt*1000
									#@warn "Retention time $(tR) surpasses modulation period, t₀ = $(new_par_sys[i_par[j]].sub[i_sub].t₀), PM = $(sys.modules[i_par[j]].PM). Initial step-width dt is decreased ($(dt))."
									sol[i_sub] = GasChromatographySimulator.solving_odesystem_r(new_par_sys[i_par[j]].col, new_par_sys[i_par[j]].prog, new_par_sys[i_par[j]].sub[i_sub], new_par_sys[i_par[j]].opt; kwargsTM..., dt=dt, dtmax=dtmax)
									tR = sol[i_sub].u[end][1]
								end
								#if isnan(tR) # does not work
								#	@warn "Simulation result is NaN. Approximate modulator."
								#	sol[i_sub] = approximate_modulator(new_par_sys[i_par[j]], df_A, sys.modules[i_par[j]].PM, sys.modules[i_par[j]].ratio, sys.modules[i_par[j]].shift)[2][i_sub]
								#end
							end
							peaklists_[j] = GasChromatographySimulator.peaklist(sol, new_par_sys[i_par[j]]; thread=pl_thread)
							add_A_to_pl!(peaklists_[j], df_A)
							solutions_[j] = sol
						end
						if maximum(peaklists_[j].τR) > sys.modules[i_par[j]].PM
							return @warn "Peak width of focussed peaks $(maximum(peaklists_[j].τR)) > modulation period $(sys.modules[i_par[j]].PM). Simulation is aborted. alg=$(sys.modules[i_par[j]].opt.alg), T=$(sys.modules[i_par[j]].T), peaklist=$(peaklists_[j])."
						end
					else # ModuleColumn
						new_par_sys[i_par[j]] = change_initial(par_sys[i_par[j]], peaklists_[j-1])
						peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(new_par_sys[i_par[j]])
						add_A_to_pl!(peaklists_[j], peaklists_[j-1])
					end
#				end
			end
#			visited_E[i_par] .= true
#			peaklists[i] = peaklists_
#			solutions[i] = solutions_
		else # path not possible
			neg_flow_modules = sys.modules[findall(paths[i][findall(positive_flow(sys)[i_par].==false)].==edges(sys.g))]
			str_neg_flow = "path not possible: "
			for j=1:length(neg_flow_modules)
				str_neg_flow = str_neg_flow*"Flow in module >$(neg_flow_modules[j].name)< becomes negative during the program. "
			end
			path_pos = str_neg_flow
		end
#	end
	return path_pos, peaklists_, solutions_, new_par_sys
end
# end - simulation of system

# begin - thermal modulator specific functions
# general smooth rectangle function
function smooth_rectangle(x, a, b, m)
	# a ... mid-position of rising flank
	# b ... mid-position of falling flank
	# m ... width of the flank
	if x.>a-3*m && x.<=a+3*m
		f = 1/2 .*(1 .+erf.((x.-a)./sqrt.(2 .*m.^2)))
	elseif x.>a+3*m && x.<=b-3*m
		f = 1
	elseif x.>b-3*m && x<=b+3*m
		f = 1 .- 1/2 .*(1 .+erf.((x.-b)./sqrt.(2 .*m.^2)))
	else
		f = 0
	end
	return f
end

# smooth rectangle with mitpoint of the rising flank at 'xstart', after 'width' the falling flank, minimum values 'min' and maximum values 'max'. 
function smooth_rectangle(x, xstart, width, min, max; flank=20)
	m = width/flank
	a = xstart + 3*m
	b = xstart + width - 3*m
	val = (max-min)*smooth_rectangle(x, a, b, m) + min
	return val
end

# periodic repeated smoothed rectangle function with period 'PM', a shift by 'shift', 'ratio' of time of Tcold to time of Thot. A small shift is incorporated to move the falling flank from the beginning to the end
function therm_mod(t, shift, PM, ratio, Tcold, Thot; flank=20) 
	# add warning, if flank value is to low -> jumps in the function
	width = (1-ratio)*PM
	tmod = mod(t+shift, PM)
	tstart = ratio*PM
	if flank == Inf # rectangle function
		# ceil() to assure, that the rounding errors of (mod) do not affect the rectangle function at the rising flank (but what is at the falling flank?) 
		return ifelse(ceil(tmod, digits=4) < tstart, Tcold, Thot)
	else # smoothed rectangle
		return smooth_rectangle.(tmod, tstart, width, Tcold, Thot; flank=flank) 
	end
end

# alternative rectangle function
# rect(t, (1+ratio)*PM/2, Tcold, Thot, (1-ratio)*PM)
function rect(x, x0, min, max, width)
	if abs(x-x0) > width/2
		min
	elseif abs(x-x0) == width/2
		(min+max)/2
	elseif abs(x-x0) < width/2
		max
	end
end

# this splicing function is with separating focussed and unfocussed peak segments
# area as an input is needed 'AR'
# A_focussed calculated in relation to it
# A_focussed is the complete area during a modulation period
# not focussed segment is already included in the focussed segment, assuming it will be focussed in the 2nd modulation in a multi stage modulation
function slicing(pl, PM, ratio, shift, par::GasChromatographySimulator.Parameters; nτ=6, τ₀=zeros(length(pl.τR)), abstol=1e-8, reltol=1e-6, alg=OwrenZen5())
	tR = pl.tR
	τR = pl.τR
	AR = pl.A
	tcold = PM*ratio
	thot = PM*(1-ratio)
	init_t_start = (fld.(tR .+ shift .- nτ.*τR, PM)).*PM .- shift # start time of the peaks, rounded down to multiple of PM
	init_t_end = (fld.(tR .+ shift .+ nτ.*τR, PM)).*PM .- shift # end time of the peaks, rounded down to multiple of PM
	n_slice = round.(Int, (init_t_end .- init_t_start)./PM .+ 1) # number of slices for every substance
	sub_TM_focussed = Array{GasChromatographySimulator.Substance}(undef, sum(n_slice))
	A_focussed = Array{Float64}(undef, sum(n_slice))
	#A_unfocussed = Array{Float64}(undef, sum(n_slice))
	g(x,p) = 1/sqrt(2*π*p[2]^2)*exp(-(x-p[1])^2/(2*p[2]^2))
	ii = 1
	Name = Array{String}(undef, sum(n_slice))
	CAS = Array{String}(undef, sum(n_slice))
	Ann_focussed = Array{String}(undef, sum(n_slice))
	t0_foc = Array{Float64}(undef, sum(n_slice))
	for i=1:length(n_slice)
		for j=1:n_slice[i]
			if n_slice[i] == 1 # no slicing, peak fits completly inside the modulation periode
				t₀ = tR[i]
			else
				t₀ = init_t_start[i]+(j-1)*PM # initial start time
			end

			CAS_par = [par.sub[i].CAS for i in 1:length(par.sub)]
			i_sub = findfirst(pl.CAS[i] .== CAS_par)
			sub_TM_focussed[ii] = GasChromatographySimulator.Substance(par.sub[i_sub].name, par.sub[i_sub].CAS, par.sub[i_sub].Tchar, par.sub[i_sub].θchar, par.sub[i_sub].ΔCp, par.sub[i_sub].φ₀, "s$(j)_"*pl.Annotations[i], par.sub[i_sub].Cag, t₀, τ₀[i])

			# Integrals:
			p = [tR[i], τR[i]]
			# approximated integrals
			#prob_focussed = IntegralProblem(g, init_t_start[i]+(j-1)*PM, init_t_start[i]+(j-1)*PM+tcold, p)
			#prob_unfocussed = IntegralProblem(g, init_t_start[i]+(j-1)*PM+tcold, init_t_start[i]+(j-1)*PM+tcold+thot, p)
			prob_focussed = IntegralProblem(g, init_t_start[i]+(j-1)*PM, init_t_start[i]+(j-1)*PM+PM, p)
			A_focussed[ii] = solve(prob_focussed, QuadGKJL(); reltol = 1e-18, abstol = 1e-30).u * AR[i]
			#A_unfocussed[ii] = solve(prob_unfocussed, QuadGKJL(); reltol = 1e-18, abstol = 1e-30).u * AR[i]
			# Areas in the same order as sub_TM_focussed
			Name[ii] = sub_TM_focussed[ii].name
			CAS[ii] = sub_TM_focussed[ii].CAS
			Ann_focussed[ii] = sub_TM_focussed[ii].ann
			t0_foc[ii] = t₀
			ii = ii + 1
		end
	end
	newopt = GasChromatographySimulator.Options(alg, abstol, reltol, par.opt.Tcontrol, par.opt.odesys, par.opt.ng, par.opt.vis, par.opt.control, par.opt.k_th)
	newpar_focussed = GasChromatographySimulator.Parameters(par.col, par.prog, sub_TM_focussed, newopt)

	#df_A_foc = DataFrame(Name=Name, CAS=CAS, Annotations=Ann_focussed, A=A_focussed+A_unfocussed, t0=t0_foc)
	df_A_foc = DataFrame(Name=Name, CAS=CAS, Annotations=Ann_focussed, A=A_focussed, t0=t0_foc)
	
	return newpar_focussed, df_A_foc
end

function add_A_to_pl!(pl, df_A)
	# add the areas in the dataframe df_A to the peaklist, same CAS and Annotations are used to identify the correct solute
#	sort_A_foc = Array{Float64}(undef, length(df_A.A_foc))
#	sort_A_unfoc = Array{Float64}(undef, length(df_A.A_foc))
	sort_A = Array{Float64}(undef, length(df_A.A))
	for i=1:length(df_A.A)
		ii = common_index(pl, df_A.CAS[i], join(split(df_A.Annotations[i], "_")[1:end-1], "_"))
		#ii = common_index(pl, df_A.CAS[i], split(df_A.Annotations[i], ", ")[1])
		sort_A[ii] = df_A.A[i]
	end
	pl[!, :A] = sort_A
	return pl
end

# identfy the common index of a substances CAS number and Annotations in a peaklist 
function common_index(pl, CAS, Annotation)
	ii_CAS = findall(CAS.==pl.CAS)
	ii_ann = findall(occursin.(Annotation, pl.Annotations))
	ii = intersect(ii_CAS, ii_ann)[1]
	return ii
end

function approximate_modulator(T, par, df_A, PM, ratio, shift, Thot;)
	# rectangular function is assumed
	tcold = PM*ratio
	thot = PM*(1-ratio)
	
	sort_df_A = sort(df_A, :t0)

	# start time as multiple of modulation period
	t₀ = fld.(sort_df_A.t0 .+ shift, PM).*PM .- shift
	
	A = sort_df_A.A

	No = [parse(Int, split(sort_df_A.Annotations[x], " ")[end]) for x in 1:length(sort_df_A.Annotations)]

	Name = sort_df_A.Name

	CAS = sort_df_A.CAS

	tR = t₀ .+ tcold

	#TR = par.prog.T_itp.(par.col.L, tR)
	# using par.prog.T_itp can result in wrong temperatures at tR, because of rounding errors for Float64 in `mod()`-function inside the `therm_mod()`-function.
	# take the temperature/temperature program defined for the module and add Thot 
	T_itp = if typeof(T) <: Number
		gf(x) = zero(x).*ones(2)
		GasChromatographySimulator.temperature_interpolation([0.0, 3600.0], [T + Thot, T + Thot], gf, par.col.L)
	else
		GasChromatographySimulator.temperature_interpolation(T.time_steps, T.temp_steps .+ Thot, T.gf, par.col.L)
	end
	TR = T_itp(par.col.L, tR) .- 273.15

	kR = Array{Float64}(undef, length(tR))
	uR = Array{Float64}(undef, length(tR))
	τ₀ = Array{Float64}(undef, length(tR))
	for i=1:length(tR)
		i_sub = findfirst(Name[i] .== df_A.Name)
		kR[i] = GasChromatographySimulator.retention_factor(par.col.L, tR[i], T_itp, par.col.d, par.col.df, par.sub[i_sub].Tchar, par.sub[i_sub].θchar, par.sub[i_sub].ΔCp, par.sub[i_sub].φ₀)

		rM = GasChromatographySimulator.mobile_phase_residency(par.col.L, tR[i], T_itp, par.prog.Fpin_itp, par.prog.pout_itp, par.col.L, par.col.d, par.col.gas; ng=par.opt.ng, vis=par.opt.vis, control=par.opt.control)
		uR[i] = 1/(rM*(1+kR[i]))

		τ₀[i] = par.sub[i_sub].τ₀
	end

	τR = par.col.L./uR

	σR = uR./τR

	Annotations = sort_df_A.Annotations

	pl = sort!(DataFrame(No=No, Name=Name, CAS=CAS, tR=tR, τR=τR, TR=TR, σR=σR, uR=uR, kR=kR, Annotations=Annotations, A=A), [:tR])

	Res = Array{Float64}(undef, length(tR))
	Δs = Array{Float64}(undef, length(tR))
	for i=1:(length(tR)-1)
		Res[i] = (pl.tR[i+1] - pl.tR[i])/(2*(pl.τR[i+1] + pl.τR[i]))
        Δs[i] = (pl.tR[i+1] - pl.tR[i])/(pl.τR[i+1] - pl.τR[i]) * log(pl.τR[i+1]/pl.τR[i])
	end
	pl[!, :Res] = Res
	pl[!, :Δs] = Δs

	sol = Array{NamedTuple{(:t, :u), Tuple{Vector{Float64}, Vector{Tuple{Float64, Float64}}}}}(undef, length(tR))
	for i=1:length(tR)
		sol[i] = (t = [0.0, par.col.L], u = [(t₀[i], τ₀[i]), (tR[i], τR[i]^2)])
	end
	
	return select(pl, [:No, :Name, :CAS, :tR, :τR, :TR, :σR, :uR, :kR, :Res, :Δs, :Annotations, :A]), sol
end

## further functions for evaluation of GCxGC_TM results
function peaklist_GCxGC(pl_end, pl_1D, PM)#, shift)
	# estimated from the gaussian peak fit to apex of projected retention times
	# here the shift is not known
	#fit_D1 = fit_gauss_D1(pl_end, pl_1D, PM, shift)
	#fit_D2 = fit_gauss_D2(pl_end, PM, shift)
	fit_D1 = fit_gauss_D1(pl_end, pl_1D, PM)
	fit_D2 = fit_gauss_D2(pl_end, PM)
	tR1 = Array{Float64}(undef, length(fit_D1.Name))
	tR2 = Array{Float64}(undef, length(fit_D1.Name))
	for i=1:length(fit_D1.Name)
		ii = findfirst(fit_D1.Name[i].==fit_D2.Name)
		tR1[i] = fit_D1.fits[i].param[1]
		tR2[i] = fit_D2.fits[ii].param[1]
	end
	return DataFrame(Name=fit_D1.Name, tR1=tR1, tR2=tR2)
end

function peaklist_GCxGC(pl_end, PM)#, shift)
	# estimated from the weighted (height) means of projected retention times
	# here the shift is not known
	heights = GasChromatographySystems.heights_of_peaks(pl_end)
	#tR1 = fld.(pl_end.tR .+ shift, PM).*PM .- shift
	#tR2 = pl_end.tR .- (fld.(pl_end.tR .+ shift, PM).*PM .- shift)
	tR1 = fld.(pl_end.tR, PM).*PM
	tR2 = pl_end.tR .- (fld.(pl_end.tR, PM).*PM)
	Name = unique(pl_end.Name)
	nsub = length(Name)
	
	sort_heights = Array{Array{Float64}}(undef, nsub)
	sort_tR1 = Array{Array{Float64}}(undef, nsub)
	sort_tR2 = Array{Array{Float64}}(undef, nsub)
	mean_tR1 = Array{Float64}(undef, nsub)
	mean_tR2 = Array{Float64}(undef, nsub)
	for i=1:nsub
		ii_name = findall(Name[i].==pl_end.Name)
		sort_heights[i] = heights[ii_name]
		sort_tR1[i] = tR1[ii_name]
		sort_tR2[i] = tR2[ii_name]
		mean_tR1[i] = sum(sort_tR1[i].*sort_heights[i])/sum(sort_heights[i])
		mean_tR2[i] = sum(sort_tR2[i].*sort_heights[i])/sum(sort_heights[i])
	end
	return DataFrame(Name=Name, tR1=mean_tR1, tR2=mean_tR2, tR1s=sort_tR1, tR2s=sort_tR2, heights=sort_heights)
end

function heights_of_peaks(pl)
	heights = Array{Float64}(undef, length(pl.tR))
 	for i=1:length(pl.tR)
		heights[i] = (GasChromatographySimulator.chromatogram([pl.tR[i]], [pl.tR[i]], [pl.τR[i]])[1].*pl.A[i])
	end
	return heights
end

function fit_envelope(pl_D2, pl_D1)
	heights = heights_of_peaks(pl_D2)
	tR = pl_D2.tR
	Name = unique(pl_D2.Name)
	nsub = length(Name)

	@. model_g(x,p) = p[3]/sqrt(2*π*p[2]^2)*exp(-(x-p[1])^2/(2*p[2]^2))
	
	sort_heights = Array{Array{Float64}}(undef, nsub)
	sort_tR = Array{Array{Float64}}(undef, nsub)
	fits = Array{LsqFit.LsqFitResult}(undef, nsub)
	for i=1:nsub
		ii_name = findall(Name[i].==pl_D2.Name)
		ii = findfirst(Name[i].==pl_D1.Name)
		sort_heights[i] = heights[ii_name]
		sort_tR[i] = tR[ii_name]
		fits[i] = curve_fit(model_g, sort_tR[i], sort_heights[i], [pl_D1.tR[ii], pl_D1.τR[ii], 1.0])
	end
	return DataFrame(Name=Name, tRs=sort_tR, heights=sort_heights, fits=fits)
end

function fit_gauss_D1(pl_D2, pl_D1, PM)#, shift) # shift?
	# here the shift is not known
	heights = heights_of_peaks(pl_D2)
	#tR = fld.(pl_D2.tR .+ shift, PM).*PM .- shift
	tR = fld.(pl_D2.tR, PM).*PM
	Name = unique(pl_D2.Name)
	nsub = length(Name)

	@. model_g(x,p) = p[3]/sqrt(2*π*p[2]^2)*exp(-(x-p[1])^2/(2*p[2]^2))
	
	sort_heights = Array{Array{Float64}}(undef, nsub)
	sort_tR = Array{Array{Float64}}(undef, nsub)
	fits = Array{LsqFit.LsqFitResult}(undef, nsub)
	for i=1:nsub
		ii_name = findall(Name[i].==pl_D2.Name)
		ii = findfirst(Name[i].==pl_D1.Name)
		sort_heights[i] = heights[ii_name]
		sort_tR[i] = tR[ii_name]
		fits[i] = curve_fit(model_g, sort_tR[i], sort_heights[i], [pl_D1.tR[ii], pl_D1.τR[ii], 1.0])
	end
	return DataFrame(Name=Name, tRs=sort_tR, heights=sort_heights, fits=fits)
end

function fit_gauss_D2(pl_D2, PM)#, shift) # shift ?
	# here the shift is not known
	heights = heights_of_peaks(pl_D2)
	#tR = pl_D2.tR .- (fld.(pl_D2.tR .+ shift, PM).*PM .- shift)
	tR = pl_D2.tR .- (fld.(pl_D2.tR, PM).*PM)
	Name = unique(pl_D2.Name)
	nsub = length(Name)

	@. model_g(x,p) = p[3]/sqrt(2*π*p[2]^2)*exp(-(x-p[1])^2/(2*p[2]^2))
	
	sort_heights = Array{Array{Float64}}(undef, nsub)
	sort_tR = Array{Array{Float64}}(undef, nsub)
	fits = Array{LsqFit.LsqFitResult}(undef, nsub)
	for i=1:nsub
		ii_name = findall(Name[i].==pl_D2.Name)
		#ii = findfirst(Name[i].==pl_D1.Name)
		sort_heights[i] = heights[ii_name]
		sort_tR[i] = tR[ii_name]
		mean_tR = sum(sort_tR[i])/length(sort_tR[i])
		fits[i] = curve_fit(model_g, sort_tR[i], sort_heights[i], [mean_tR, mean_tR/10, 1.0])
	end
	return DataFrame(Name=Name, tRs=sort_tR, heights=sort_heights, fits=fits)
end

function chrom(pl; nτ=6)
	tstart = Array{Float64}(undef, length(pl.Name))
	tend = Array{Float64}(undef, length(pl.Name))
	t = Array{Array{Float64}}(undef, length(pl.Name))
	c = Array{Array{Float64}}(undef, length(pl.Name))
	for i=1:length(pl.Name)
		tstart[i] = pl.tR[i] - nτ * pl.τR[i]
		tend[i] = pl.tR[i] + nτ * pl.τR[i]
		t[i] = collect(tstart[i]:(2*nτ*pl.τR[i]/100):tend[i])
		c[i] = GasChromatographySimulator.chromatogram(t[i], [pl.tR[i]], [pl.τR[i]])*pl.A[i]
	end
	t1 = minimum(tstart)
	t2 = maximum(tend)
	dt = 2*nτ*minimum(pl.τR)/100
	t_sum = collect(t1:dt:t2)
	c_sum = fill(0.0, length(t_sum))
	for i=1:length(pl.Name)
		c_ = GasChromatographySimulator.chromatogram(t_sum, [pl.tR[i]], [pl.τR[i]])*pl.A[i]
		c_sum = c_sum .+ c_
	end

	names = unique(pl.Name)
	
	p_chrom = Plots.plot(t_sum, c_sum, xlabel="time in s", label="Chromatogram")
	for i=1:length(pl.Name)
		i_names = findfirst(pl.Name[i].==names)
		#if i > length(names)
		#	lbl = ""
		#else
			lbl = pl.Name[i]
		#end
		Plots.plot!(p_chrom, t[i], c[i], label=lbl, color=i_names+1)
	end
	Plots.plot!(p_chrom, ylims=(-0.02*maximum(c_sum), 1.02*maximum(c_sum)), xlims=(minimum(t_sum), maximum(t_sum)))
	return p_chrom, t_sum, c_sum, t, c 
end

function chrom_marked(pl, PM, ratio, shift; nτ=6)
	p_chrom, t_sum, c_sum, t, c = chrom(pl; nτ=nτ)
	max_y = maximum(c_sum)
	# add modulation period
	n = unique(fld.(t_sum, PM + shift))
	for i=1:length(n)
		Plots.plot!(p_chrom, [n[i]*PM-shift, n[i]*PM-shift], [0.0, 1.1*maximum(c_sum)], c=:black, label="")
		Plots.plot!(p_chrom, [n[i]*PM-shift+PM*ratio, n[i]*PM-shift+PM*ratio], [0.0, 1.1*maximum(c_sum)], c=:black, linestyle=:dash, label="")
	end
	Plots.plot!(p_chrom, ylims=(-0.02*maximum(c_sum), 1.02*maximum(c_sum)), xlims=(minimum(t_sum), maximum(t_sum)))
	return p_chrom, t_sum, c_sum, t, c
end

function collect_chrom(pl_array, sys; markings=true)
	# collect chromatograms for all segments
	# chromatograms of segments, which are followed by a ModuleTM should be marked with the modulation period (coldjet, hotjet)
	# add titles (module name)
	p = Array{Any}(undef, length(pl_array))
	for i=1:length(pl_array)
		if markings == true && i < length(pl_array)
			if typeof(sys.modules[i+1]) == GasChromatographySystems.ModuleTM
				p[i] = chrom_marked(pl_array[i], sys.modules[i+1].PM, sys.modules[i+1].ratio, sys.modules[i+1].shift)[1]
				
			else
				p[i] = chrom(pl_array[i])[1]
			end
		else
			p[i] = chrom(pl_array[i])[1]
		end
		Plots.plot!(p[i], title=sys.modules[i].name)
	end
	return p
end

#=function chrom_slicing(t, c, PM, shift) # correctly account for the shift!!!
	n = Int.(fld.(collect(t).-shift, PM)) # number of the slices
	slices = Array{Array{Float64, 1}, 1}(undef, length(unique(n)))
	t_D2 = Array{Array{Float64, 1}, 1}(undef, length(unique(n)))
	for i=1:length(unique(n))
		i1 = findfirst(unique(n)[i].==n)
		i2 = findlast(unique(n)[i].==n)
		slices[i] = c[i1:i2]
		t_D2[i] = t[i1:i2] .- unique(n)[i] * PM
	end
	t_D1 = (0.0:PM:t[end]).- shift# shift?
	return slices, t_D1, t_D2
end=#

function chrom2d(pl_final, sys, PM)
	# a shift is unknown here
	t = 0.0:0.01:sum(GasChromatographySystems.common_timesteps(sys))
	chrom = Array{Array{Float64}}(undef, length(pl_final.tR))
	for i=1:length(pl_final.tR)
		chrom[i] = GasChromatographySimulator.chromatogram(collect(t), [pl_final.tR[i]], [pl_final.τR[i]]).*pl_final.A[i]
	end
	chrom_sum = chrom[1]
	for i=2:length(chrom)
		chrom_sum = chrom_sum .+ chrom[i]
	end

	t_D1 = (0.0:PM:t[end])
	
	n = Int.(fld.(collect(t), PM))
	slices = Array{Array{Float64, 1}, 1}(undef, length(unique(n)))
	t_D2 = Array{Array{Float64, 1}, 1}(undef, length(unique(n)))
	for i=1:length(unique(n))
		i1 = findfirst(unique(n)[i].==n)
		i2 = findlast(unique(n)[i].==n)
		slices[i] = chrom_sum[i1:i2]
		t_D2[i] = t[i1:i2] .- (unique(n)[i] * PM )
	end
	
	slice_mat = Array{Float64}(undef, length(slices)-1, length(t_D2[1]))
	for j=1:length(t_D2[1])
		for i=1:(length(slices)-1)
			slice_mat[i,j] = slices[i][j]
		end
	end
	return slice_mat, t_D1, t_D2, slices, t, chrom_sum#, chrom_sliced 
end

function comparison_meas_sim(meas, pl_sim)
	Name = meas.Name
	tR1_meas = meas.tR1
	tR2_meas = meas.tR2
	index = [findfirst(meas.Name[x].==pl_sim.Name) for x in 1:length(meas.Name)]
	tR1_sim = Array{Union{Missing,Float64}}(undef, length(meas.Name))
	tR2_sim = Array{Union{Missing,Float64}}(undef, length(meas.Name))
	for i=1:length(meas.Name)
		if isnothing(index[i])
			tR1_sim[i] = missing
			tR2_sim[i] = missing
		else
			tR1_sim[i] = pl_sim.tR1[index[i]]
			tR2_sim[i] = pl_sim.tR2[index[i]]
		end
	end
	comp = DataFrame(Name=Name, tR1_meas=tR1_meas, tR1_sim=tR1_sim, ΔtR1=tR1_meas.-tR1_sim, relΔtR1_percent=(tR1_meas.-tR1_sim)./tR1_meas.*100.0, tR2_meas=tR2_meas, tR2_sim=tR2_sim, ΔtR2=tR2_meas.-tR2_sim, relΔtR2_percent=(tR2_meas.-tR2_sim)./tR2_meas.*100.0)
	#for i=1:length(comp.Name)
	#	ii = findfirst(comp.Name[i].==meas.Name)
	#	if ismissing(meas.tR1[ii])
	#		comp[i, :tR1_meas] = NaN
	#		comp[i, :tR2_meas] = NaN
	#		comp[i, :ΔtR1] = NaN
	#		comp[i, :ΔtR2] = NaN
	#	else
	#		comp[i, :tR1_meas] = meas.tR1[ii]
	#		comp[i, :tR2_meas] = meas.tR2[ii]
	#		comp[i, :ΔtR1] = meas.tR1[ii] - comp[i, :tR1_sim]
	#		comp[i, :ΔtR2] = meas.tR2[ii] - comp[i, :tR2_sim]
	#	end
	#end
	return comp
end

function chrom_slicing(t1, c, PM)
	n = Int.(fld.(collect(t1), PM)) # number of the slices
	slices = Array{Array{Float64, 1}, 1}(undef, length(unique(n)))
	t_D2 = Array{Array{Float64, 1}, 1}(undef, length(unique(n)))
	for i=1:length(unique(n))
		i1 = findfirst(unique(n)[i].==n)
		i2 = findlast(unique(n)[i].==n)
		slices[i] = c[i1:i2]
		t_D2[i] = t1[i1:i2] .- unique(n)[i] * PM
	end
	t_D1 = 0.0:PM:t1[end]
	return slices, t_D1, t_D2
end

function check_area(pl)
	names = unique(pl.Name)
	area = Array{Float64}(undef, length(names))
	for i=1:length(names)
		area[i] = sum(filter([:Name] => x -> x == names[i], pl).A)
	end
	return DataFrame(Name=names, sum_A = area)
end

function check_peakwidths(pl; τ_threshold = 0.5)
	pl_f = filter([:τR] => x -> x > τ_threshold, pl)
	names = unique(pl_f.Name)
	return names
end

function check_duration_modulation(pl_array, par, PM, ratio)
	Δts = duration_in_module(pl_array, par)
	ok_TM1 = Array{Bool}(undef, length(Δts[3].Δt))
	for i=1:length(Δts[3].Δt)
		if (Δts[3].Δt[i] > ratio*PM) && (Δts[3].Δt[i] < PM)
			ok_TM1[i] = true
		else
			ok_TM1[i] = false
		end
	end
	ok_TM2 = Array{Bool}(undef, length(Δts[5].Δt))
	
	for i=1:length(Δts[5].Δt)
		#if (Δts[5].Δt[i] > ratio*PM) && (Δts[5].Δt[i] < PM) # change this condition
		if (Δts[5].Δt[i] < PM)
			ok_TM2[i] = true
		else
			ok_TM2[i] = false
		end
	end
	index_TM1 = findall(ok_TM1.==false)
	index_TM2 = findall(ok_TM2.==false)
	TM1 = DataFrame(Index=index_TM1, Name=Δts[3].Name[index_TM1], CAS=Δts[3].CAS[index_TM1], tR=pl_array[3].tR[index_TM1], τR=pl_array[3].τR[index_TM1], Δt=Δts[3].Δt[index_TM1], Annotations=Δts[3].Annotations[index_TM1])
	TM2 = DataFrame(Index=index_TM2, Name=Δts[5].Name[index_TM2], CAS=Δts[5].CAS[index_TM2], tR=pl_array[5].tR[index_TM2], τR=pl_array[5].τR[index_TM2], Δt=Δts[5].Δt[index_TM2], Annotations=Δts[5].Annotations[index_TM2])
	return TM1, TM2
end

function duration_in_module(pl_array, par)
	Δts = Array{DataFrame}(undef, length(pl_array))
	for i=1:length(pl_array)
		CAS_par = [par[i].sub[x].CAS for x in 1:length(par[i].sub)]
		ann_par = [par[i].sub[x].ann for x in 1:length(par[i].sub)]

		Δt = Array{Float64}(undef, length(CAS_par))
		name = Array{String}(undef, length(CAS_par))
		for j=1:length(CAS_par)
			jj = GasChromatographySystems.common_index(pl_array[i], CAS_par[j], ann_par[j])
			Δt[j] = pl_array[i].tR[jj] - par[i].sub[j].t₀
			name[j] = pl_array[i].Name[jj]
		end
		df_Δt = DataFrame(Name=name, CAS=CAS_par, Δt=Δt, Annotations=ann_par)
		Δts[i] = df_Δt
	end
	return Δts
end
# end - thermal modulator specific functions

# old? 
function plot_GCxGC(pl_GCxGC, sys; categories = String[])
	t1end = sum(GasChromatographySystems.common_timesteps(sys))
	t2end = ceil(maximum(pl_GCxGC.tR2.-pl_GCxGC.tR1))
	p_gcxgc = Plots.scatter(pl_GCxGC.tR1, pl_GCxGC.tR2.-pl_GCxGC.tR1, xlims=(0.0,t1end), ylims=(0.0, t2end), label="", xlabel=string(sys.modules[1].sp, ", t¹ in s"), ylabel=string(sys.modules[2].sp, ", t² in s"))
	if !isempty(categories)
		for i=1:length(categories) 
			pl_f = filter([:Cat] => x -> categories[i] in x, pl_GCxGC)
			Plots.scatter!(p_gcxgc, pl_f.tR1, pl_f.tR2.-pl_f.tR1, label=categories[i], legend=:bottomright)
		end
	end
	return p_gcxgc
end


# begin - specific systems
function SeriesSystem(Ls, ds, dfs, sps, TPs, F, pin, pout; opt=GasChromatographySystems.Options(), kwargs...)
	# add test for correct lengths of input
	# ? make two versions
	# 1. defining flow over the columns (calculate pin)
	# 2. defining inlet pressure (calculate F)
	n = length(Ls)
	g = SimpleDiGraph(n+1)
	for i=1:n
		add_edge!(g, i, i+1) 
	end
	# common time steps
	com_timesteps = []
	for i=1:length(TPs)
		if typeof(TPs[i]) <: GasChromatographySystems.TemperatureProgram
			com_timesteps = GasChromatographySimulator.common_time_steps(com_timesteps, TPs[i].time_steps)
		end
	end
	if isempty(com_timesteps)
		com_timesteps = [0.0, 36000.0]
	end
	# pressure points
	pp = Array{GasChromatographySystems.PressurePoint}(undef, n+1)
	pins = pin*1000.0.*ones(length(com_timesteps))
	nans = NaN.*ones(length(com_timesteps))
	if pout == 0.0
		pouts = eps(Float64).*ones(length(com_timesteps))
	else 
		pouts = pout*1000.0.*ones(length(com_timesteps))
	end
	pp[1] = GasChromatographySystems.PressurePoint("p1", com_timesteps, pins) # inlet
	for i=2:n
		pp[i] = GasChromatographySystems.PressurePoint("p$(i)", com_timesteps, nans) #
	end
	pp[end] = GasChromatographySystems.PressurePoint("p$(n+1)", com_timesteps, pouts) # outlet
	# modules
	modules = Array{GasChromatographySystems.AbstractModule}(undef, n)
	for i=1:n
		if i==1
			modules[i] = GasChromatographySystems.ModuleColumn("$(i) -> $(i+1)", Ls[i], ds[i]*1e-3, dfs[i]*1e-6, sps[i], TPs[i], F/60e6; kwargs...)
		else
			modules[i] = GasChromatographySystems.ModuleColumn("$(i) -> $(i+1)", Ls[i], ds[i]*1e-3, dfs[i]*1e-6, sps[i], TPs[i], NaN; kwargs...)
		end
	end
	# system
	sys = GasChromatographySystems.update_system(GasChromatographySystems.System(g, pp, modules, opt))

	# add test for the defined pressures and flows
	return sys
end

function SeriesSystem(; Ls = [10.0, 5.0, 2.0, 1.0], ds = [0.53, 0.32, 0.25, 0.1], dfs = [0.53, 0.32, 0.25, 0.1], sps = ["Rxi17SilMS", "Rxi17SilMS", "Rxi17SilMS", "Rxi17SilMS"], TPs = [default_TP(), default_TP(), default_TP(), default_TP()], F = NaN, pin = 300.0, pout = 0.0, opt=GasChromatographySystems.Options(), kwargs...)
	sys = SeriesSystem(Ls, ds, dfs, sps, TPs, F, pin, pout; opt=opt, kwargs...)
	return sys
end

#example_SeriesSystem() = SeriesSystem([10.0, 5.0, 2.0, 1.0], [0.53, 0.32, 0.25, 0.1], [0.53, 0.32, 0.25, 0.1], ["Rxi17SilMS", "Rxi17SilMS", "Rxi17SilMS", "Rxi17SilMS"], [default_TP(), default_TP(), default_TP(), default_TP()], NaN, 300.0, 0.0)

function SplitSystem(Ls, ds, dfs, sps, TPs, Fs, pin, pout1, pout2; opt=GasChromatographySystems.Options(), kwargs...)
	g = SimpleDiGraph(4)
	add_edge!(g, 1, 2) # Inj -> GC column -> Split point
	add_edge!(g, 2, 3) # Split point -> TL column -> Det 1
	add_edge!(g, 2, 4) # Split point -> TL column -> Det 2
	# common time steps
	com_timesteps = []
	for i=1:length(TPs)
		if typeof(TPs[i]) <: GasChromatographySystems.TemperatureProgram
			com_timesteps = GasChromatographySimulator.common_time_steps(com_timesteps, TPs[i].time_steps)
		end
	end
	if isempty(com_timesteps)
		com_timesteps = [0.0, 36000.0]
	end
	# pressure points
	pp = Array{GasChromatographySystems.PressurePoint}(undef, nv(g))
	pins = pin*1000.0.*ones(length(com_timesteps))
	nans = NaN.*ones(length(com_timesteps))
	if pout1 == 0.0
		pout1s = eps(Float64).*ones(length(com_timesteps))
	else 
		pout1s = pout1*1000.0.*ones(length(com_timesteps))
	end
	if pout2 == 0.0
		pout2s = eps(Float64).*ones(length(com_timesteps))
	else 
		pout2s = pout2*1000.0.*ones(length(com_timesteps))
	end
	pp[1] = GasChromatographySystems.PressurePoint("p₁", com_timesteps, pins) # inlet 
	pp[2] = GasChromatographySystems.PressurePoint("p₂", com_timesteps, nans) # 
	pp[3] = GasChromatographySystems.PressurePoint("p₃", com_timesteps, pout1s) # outlet 1 
	pp[4] = GasChromatographySystems.PressurePoint("p₄", com_timesteps, pout2s) # outlet 2
	# modules
	modules = Array{GasChromatographySystems.AbstractModule}(undef, ne(g))
	modules[1] = GasChromatographySystems.ModuleColumn("1 -> 2", Ls[1], ds[1]*1e-3, dfs[1]*1e-6, sps[1], TPs[1], Fs[1]/60e6; kwargs...)
	modules[2] = GasChromatographySystems.ModuleColumn("2 -> 3", Ls[2], ds[2]*1e-3, dfs[2]*1e-6, sps[2], TPs[2], Fs[2]/60e6; kwargs...)
	modules[3] = GasChromatographySystems.ModuleColumn("2 -> 4", Ls[3], ds[3]*1e-3, dfs[3]*1e-6, sps[3], TPs[3], Fs[3]/60e6; kwargs...)
	# system
	sys = GasChromatographySystems.update_system(GasChromatographySystems.System(g, pp, modules, opt))

	# add test for the defined pressures and flows
	return sys
end

function SplitSystem(; Ls = [10.0, 1.0, 5.0], ds = [0.25, 0.1, 0.25], dfs = [0.25, 0.0, 0.0], sps = ["Rxi17SilMS", "", ""], TPs = [default_TP(), 300.0, 300.0], Fs = [1.0, NaN, NaN], pin = NaN, pout1 = 0.0, pout2 = 101.3, opt=GasChromatographySystems.Options(), kwargs...)
	sys = SplitSystem(Ls, ds, dfs, sps, TPs, Fs, pin, pout1, pout2; opt=opt, kwargs...)
	return sys
end

# definition GCxGC system with thermal modulator
function GCxGC_TM(L1, d1, df1, sp1, TP1, L2, d2, df2, sp2, TP2, LTL, dTL, dfTL, spTL, TPTL, LM::Array{Float64,1}, dM, dfM, spM, shift, PM, ratioM, HotM, ColdM, TPM, F, pin, pout; opt=GasChromatographySystems.Options(), optTM=ModuleTMOptions(), optCol=ModuleColumnOptions())

	TPs = [TP1, TP2, TPM]
	
	g = SimpleDiGraph(9)
	add_edge!(g, 1, 2) # 1st-D GC
	add_edge!(g, 2, 3) # modulator
	add_edge!(g, 3, 4) # hot/cold 1
	add_edge!(g, 4, 5) # modulator
	add_edge!(g, 5, 6) # hot/cold 2 
	add_edge!(g, 6, 7) # modulator
	add_edge!(g, 7, 8) # 2nd-D GC
	add_edge!(g, 8, 9) # TL
	
	# common time steps
	com_timesteps = []
	for i=1:length(TPs)
		if typeof(TPs[i]) <: TemperatureProgram
			com_timesteps = GasChromatographySimulator.common_time_steps(com_timesteps, TPs[i].time_steps)
		end
	end
	if isempty(com_timesteps)
		com_timesteps = [0.0, 36000.0]
	end
	
	# pressure points
	if length(pin) == 1
		pins = pin*1000.0.*ones(length(com_timesteps))
	else
	end
	nans = NaN.*ones(length(com_timesteps))
	if pout == 0.0
		pouts = eps(Float64).*ones(length(com_timesteps))
	else 
		pouts = pout*1000.0.*ones(length(com_timesteps))
	end
	pp = Array{PressurePoint}(undef, nv(g))
	pp[1] = PressurePoint("p1", com_timesteps, pins) # inlet 
	for i=2:(nv(g)-1)
		pp[i] = PressurePoint("p$(i)", com_timesteps, nans)
	end
	pp[end] = PressurePoint("p$(nv(g))", com_timesteps, pouts) # outlet
	# modules
	modules = Array{AbstractModule}(undef, ne(g))
	modules[1] = ModuleColumn("GC column 1", L1, d1*1e-3, df1*1e-6, sp1, TP1, F/60e6, optCol)
	modules[2] = ModuleColumn("mod in", LM[1], dM*1e-3, dfM*1e-6, spM, TPM, optCol)
	modules[3] = ModuleTM("TM1", LM[2], dM*1e-3, dfM*1e-6, spM, TPM, shift, PM, ratioM, HotM, ColdM, NaN, optTM)
	modules[4] = ModuleColumn("mod loop", LM[3], dM*1e-3, dfM*1e-6, spM, TPM, optCol)
	modules[5] = ModuleTM("TM2", LM[4], dM*1e-3, dfM*1e-6, spM, TPM, shift, PM, ratioM, HotM, ColdM, NaN, optTM)
	modules[6] = ModuleColumn("mod out", LM[5], dM*1e-3, dfM*1e-6, spM, TPM, NaN, optCol)
	modules[7] = ModuleColumn("GC column 2", L2, d2*1e-3, df2*1e-6, sp2, TP2, NaN, optCol)
	modules[8] = ModuleColumn("TL", LTL, dTL*1e-3, dfTL*1e-6, spTL, TPTL, NaN, optCol)
	# system
	sys = update_system(System(g, pp, modules, opt))
	return sys
end

function GCxGC_TM(; L1 = 30.0, d1 = 0.25, df1 = 0.25, sp1 = "ZB1ms", TP1 = default_TP(), L2 = 2.0, d2 = 0.1, df2 = 0.1, sp2 = "Stabilwax", TP2 = default_TP(), LTL = 0.25, dTL = 0.1, dfTL = 0.1, spTL = "Stabilwax", TPTL = 280.0, LM = [0.30, 0.01, 0.90, 0.01, 0.30], dM = 0.1, dfM = 0.1, spM = "Stabilwax", shift = 0.0, PM = 4.0, ratioM = 0.9125, HotM = 30.0, ColdM = -120.0, TPM = default_TP(), F = 0.8, pin = NaN, pout = 0.0, opt=GasChromatographySystems.Options(), optTM=ModuleTMOptions(), optCol=ModuleColumnOptions())
	sys = GCxGC_TM(L1, d1, df1, sp1, TP1, L2, d2, df2, sp2, TP2, LTL, dTL, dfTL, spTL, TPTL, LM::Array{Float64,1}, dM, dfM, spM, shift, PM, ratioM, HotM, ColdM, TPM, F, pin, pout; opt=opt, optTM=optTM, optCol=optCol)
	return sys
end
# end - specific systems


#example_GCxGC_TM() = GCxGC_TM(30.0, 0.25, 0.25, "ZB1ms", default_TP(), 0.1, 0.1, 0.1, "Stabilwax", default_TP(), 0.56, 0.1, 0.1, "Stabilwax", 280.0, [0.30, 0.01, 0.90, 0.01, 0.30], 0.1, 0.1, "Stabilwax", 0.0, 0.0, 4.0, 0.9125, 80.0, -120.0, default_TP(), 0.8, NaN, 0.0)



# begin - Misc
function traces(sol, par, i_select)
	z = sol[i_select].t
	
		tt = Array{Float64}(undef, length(sol[i_select].t))
		ττ = Array{Float64}(undef, length(sol[i_select].t))
		TT = Array{Float64}(undef, length(sol[i_select].t))
		kk = Array{Float64}(undef, length(sol[i_select].t))
		for j=1:length(sol[i_select].t)
			tt[j] = sol[i_select].u[j][1]
			ττ[j] = sol[i_select].u[j][2]
			TT[j] = par.prog.T_itp(sol[i_select].t[j], sol[i_select].u[j][1])
			kk[j] = GasChromatographySimulator.retention_factor(sol[i_select].t[j], sol[i_select].u[j][1], par.col, par.prog, par.sub[i_select], par.opt)
		end
	return trace = DataFrame(z=z, t=tt, τ²=ττ, T=TT, k=kk)
end
# end - Misc

# begin - hold-up times
function holdup_time_functions(sys)
	p_func = GasChromatographySystems.pressure_functions(sys)
	tM_func = Array{Function}(undef, GasChromatographySystems.ne(sys.g))
	E = collect(GasChromatographySystems.edges(sys.g))
	srcE = GasChromatographySystems.src.(E)
	dstE = GasChromatographySystems.dst.(E)
	for i=1:GasChromatographySystems.ne(sys.g)
		pin(t) = p_func[srcE[i]](t)
		pout(t) = p_func[dstE[i]](t)
		T_itp = GasChromatographySystems.module_temperature(sys.modules[i], sys)[5]
		f(t) = GasChromatographySimulator.holdup_time(t, T_itp, pin, pout, sys.modules[i].L, sys.modules[i].d, sys.options.gas; ng=sys.modules[i].opt.ng, vis=sys.options.vis, control=sys.options.control)
		tM_func[i] = f
	end
	return tM_func
end

function holdup_time_path(sys, num_paths)
	tM = holdup_time_functions(sys)
	paths = all_paths(sys.g, num_paths)[2]
	
	tMp = Array{Function}(undef, num_paths)
	for i=1:num_paths
		i_paths = GasChromatographySystems.index_parameter(sys.g, paths[i])
		f(t) = sum([tM[x](t) for x in i_paths])
		tMp[i] = f
	end
	return tMp
end

function holdup_time_functions(sys, p2fun)
	# collecting the hold-up time functions of every edge as function of time t for system `sys` and the squared pressure solution functions `p2fun`. This function should be used, if parameters of the system are to be changes, e.g. column length or diameter, but the structure of the system is the same (same grape, same unknown pressures/flows)
	p_func = pressure_functions(sys, p2fun)
	tM_func = Array{Function}(undef, GasChromatographySystems.ne(sys.g))
	E = collect(GasChromatographySystems.edges(sys.g))
	srcE = GasChromatographySystems.src.(E)
	dstE = GasChromatographySystems.dst.(E)
	for i=1:GasChromatographySystems.ne(sys.g)
		pin(t) = p_func[srcE[i]](t)
		pout(t) = p_func[dstE[i]](t)
		T_itp = GasChromatographySystems.module_temperature(sys.modules[i], sys)[5]
		f(t) = GasChromatographySimulator.holdup_time(t, T_itp, pin, pout, sys.modules[i].L, sys.modules[i].d, sys.options.gas; ng=sys.modules[i].opt.ng, vis=sys.options.vis, control=sys.options.control)
		tM_func[i] = f
	end
	return tM_func
end

function holdup_time_path(sys, p2fun, num_paths)
	# collecting the hold-up time functions of every path as function of time t for system `sys` and the squared pressure solution functions `p2fun`. This function should be used, if parameters of the system are to be changes, e.g. column length or diameter, but the structure of the system is the same (same grape, same unknown pressures/flows)
	tM = holdup_time_functions(sys, p2fun)
	paths = all_paths(sys.g, num_paths)[2]
	
	tMp = Array{Function}(undef, num_paths)
	for i=1:num_paths
		i_paths = GasChromatographySystems.index_parameter(sys.g, paths[i])
		f(t) = sum([tM[x](t) for x in i_paths])
		tMp[i] = f
	end
	return tMp
end

function plot_holdup_time_path_over_time(sys, num_paths; dt=60.0)
	#plotly()
	tMp_func = holdup_time_path(sys, num_paths)
	com_timesteps = GasChromatographySystems.common_timesteps(sys)
	trange = 0:dt:sum(com_timesteps)
	p_tM = Plots.plot(xlabel="time in s", ylabel="hold-up time in s")
	for i=1:num_paths
		Plots.plot!(p_tM, trange, tMp_func[i].(trange), label="path: $(i)")
	end
	return p_tM
end
# end - hold-up times


end # module