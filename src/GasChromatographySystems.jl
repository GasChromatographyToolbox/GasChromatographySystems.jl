module GasChromatographySystems

using Reexport
@reexport using Graphs
using GraphMakie
@reexport using Symbolics
using NetworkLayout
using CairoMakie
using GasChromatographySimulator
using Interpolations
@reexport using CSV
@reexport using DataFrames
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
include("./Systems.jl")
include("./SolvingSystems.jl")
include("./ThermalModulator.jl")
# functions

# begin - misc, for update_system
# common programs
function common_timesteps(sys; default=[0.0, 36000.0])
	com_timesteps = []
	for i=1:nv(sys.g)
		if typeof(sys.pressurepoints[i].P) <: PressureProgram
			com_timesteps = GasChromatographySimulator.common_time_steps(com_timesteps, sys.pressurepoints[i].P.time_steps)
		end
	end
	for i=1:ne(sys.g)
		if typeof(sys.modules[i].T) <: TemperatureProgram
			com_timesteps = GasChromatographySimulator.common_time_steps(com_timesteps, sys.modules[i].T.time_steps)
		end
	end
	if isempty(com_timesteps)
		com_timesteps = default
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

function index_pressurepoints_with_pressure_program(sys)
	i_pressprog = Int[]
	for i=1:nv(sys.g)
		if typeof(sys.pressurepoints[i].P) <: PressureProgram
			push!(i_pressprog, i)
		end
	end
	return i_pressprog
end

function match_programs(sys)
	com_times = GasChromatographySystems.common_timesteps(sys)
	i_pressprog = index_pressurepoints_with_pressure_program(sys)
	new_press_steps = Array{Array{Float64,1}}(undef, length(i_pressprog))
	for i=1:length(i_pressprog)
		new_press_steps[i] = GasChromatographySimulator.new_value_steps(sys.pressurepoints[i_pressprog[i]].P.pressure_steps, sys.pressurepoints[i_pressprog[i]].P.time_steps, com_times)
	end
	i_tempprog = index_modules_with_temperature_program(sys)
	new_temp_steps = Array{Array{Float64,1}}(undef, length(i_tempprog))
	for i=1:length(i_tempprog)
		new_temp_steps[i] = GasChromatographySimulator.new_value_steps(sys.modules[i_tempprog[i]].T.temp_steps, sys.modules[i_tempprog[i]].T.time_steps, com_times)
	end
	# add for gradient new_a_gf
	return com_times, new_press_steps, new_temp_steps, i_pressprog, i_tempprog
end

function update_system(sys)
	new_timesteps, new_pressuresteps, new_temperaturesteps, index_pp_pressprog, index_module_tempprog = match_programs(sys)
	new_pp = Array{PressurePoint}(undef, nv(sys.g))
	for i=1:nv(sys.g)
		if typeof(sys.pressurepoints[i].P) <: Number
			new_pp[i] = sys.pressurepoints[i]
		elseif typeof(sys.pressurepoints[i].P) <: PressureProgram
			ii = findfirst(index_pp_pressprog.==i)
			new_presprog = PressureProgram(new_timesteps, new_pressuresteps[ii])
			new_pp[i] = PressurePoint(sys.pressurepoints[i].name, new_presprog)
		end
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
	new_sys = System(sys.name, sys.g, new_pp, new_modules, sys.options)
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
		T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, L; ng=module_.opt.ng)
	elseif typeof(module_.T) <: Number # temperature is a constant value
		time_steps = if isempty(GasChromatographySystems.common_timesteps(sys))
			[0.0, 36000.0]
		else
			GasChromatographySystems.common_timesteps(sys)
		end
		temp_steps = module_.T.*ones(length(time_steps))
		gf(x) = zero(x).*ones(length(time_steps))
		a_gf = [zeros(length(time_steps)) zeros(length(time_steps)) ones(length(time_steps)) zeros(length(time_steps))]
		T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, L; ng=module_.opt.ng)
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
		time_steps = if isempty(GasChromatographySystems.common_timesteps(sys))
			[0.0, 36000.0]
		else
			GasChromatographySystems.common_timesteps(sys)
		end
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
		g = if typeof(sys.pressurepoints[i].P) <: Number
				GasChromatographySimulator.steps_interpolation([0.0, 36000.0], fill(sys.pressurepoints[i].P, 2))
		elseif typeof(sys.pressurepoints[i].P) <: GasChromatographySystems.PressureProgram
			#if isnan(sys.pressurepoints[i].P.pressure_steps[1])
			#	GasChromatographySimulator.steps_interpolation(sys.pressurepoints[i].P.time_steps, NaN.*ones(length(sys.pressurepoints[i].P.time_steps)))
			#else
				GasChromatographySimulator.steps_interpolation(sys.pressurepoints[i].P.time_steps, identity.(sys.pressurepoints[i].P.pressure_steps))
			#end
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

#=
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
=#

function plot_graph_with_flow(sys, p2fun, t; lay = GasChromatographySystems.Spring(), color=:lightblue, node_size=80, arrow_size=20, arrow_shift=0.8, dataaspect=false, nlabels_fontsize=14, elabels_fontsize=14, elabels_distance = 20)
	p_func = GasChromatographySystems.pressure_functions(sys, p2fun)
	F_func = GasChromatographySystems.flow_functions(sys, p2fun)
	fig, ax, p = GasChromatographySystems.GraphMakie.graphplot(sys.g, 
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
	GasChromatographySystems.hidedecorations!(ax)
	GasChromatographySystems.hidespines!(ax)
	if dataaspect == true
		ax.aspect = DataAspect()
	end
	return fig
end

#=
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
=#

function plot_flow_over_time(sys, p2fun; dt=60.0)
	#plotly()
	F_func = GasChromatographySystems.flow_functions(sys, p2fun)
	com_timesteps = GasChromatographySystems.common_timesteps(sys)
	trange = 0:dt:sum(com_timesteps)
	p_flow = Plots.plot(xlabel="time in s", ylabel="flow in mL/min")
	for i=1:ne(sys.g)
		Plots.plot!(p_flow, trange, F_func[i].(trange).*60e6, label="F_$(sys.modules[i].name)")
	end
	return p_flow
end

#=
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
=#

function plot_pressure_over_time(sys, p2fun; dt=60.0)
	#plotly()
	p_func = GasChromatographySystems.pressure_functions(sys, p2fun)
	com_timesteps = GasChromatographySystems.common_timesteps(sys)
	trange = 0:dt:sum(com_timesteps)
	p_pres = Plots.plot(xlabel="time in s", ylabel="pressure in Pa(a)", legend=:topleft)
	for i=1:nv(sys.g)
		Plots.plot!(p_pres, trange, p_func[i].(trange), label="$(sys.pressurepoints[i].name)")
	end
	return p_pres
end

#=
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
=#

function plot_holdup_time_path_over_time(sys, p2fun, num_paths; dt=60.0)
	#plotly()
	tMp_func = holdup_time_path(sys, p2fun, num_paths)
	com_timesteps = GasChromatographySystems.common_timesteps(sys)
	trange = 0:dt:sum(com_timesteps)
	p_tM = Plots.plot(xlabel="time in s", ylabel="hold-up time in s")
	for i=1:num_paths
		Plots.plot!(p_tM, trange, tMp_func[i].(trange), label="path: $(i)")
	end
	return p_tM
end
# end - plotting of the graphs and functions

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
"""
    therm_mod(t, shift, PM, ratio, Tcold, Thot; flank=20)

Generate a periodic temperature modulation pattern for thermal modulation in gas chromatography.

The function creates a temperature profile that alternates between cold (Tcold) and hot (Thot) phases
with a specified period (PM). The pattern can be either a sharp rectangular function or a smoothed
version with controlled transition flanks.

# Arguments
- `t`: Time point(s) at which to evaluate the temperature
- `shift`: Time shift of the modulation pattern (in seconds)
- `PM`: Modulation period (in seconds)
- `ratio`: Ratio of cold phase duration to total period (0 < ratio < 1)
- `Tcold`: Temperature during the cold phase (in °C)
- `Thot`: Temperature during the hot phase (in °C)

# Keyword Arguments
- `flank`: Controls the smoothness of temperature transitions. Use `Inf` for sharp transitions,
          or a positive number for smoothed transitions (default: 20)

# Returns
- Temperature value(s) at the specified time point(s)

# Notes
- on modulation consists of a cold phase and a hot phase in this order
- The cold phase duration is `ratio * PM`
- The hot phase duration is `(1 - ratio) * PM`
- When `shift = 0`, the pattern starts with the hot phase
- The `flank` parameter controls the width of the temperature transition regions
"""
function therm_mod(t, shift, PM, ratio, Tcold, Thot; flank=20) 
	# add warning, if flank value is to low -> jumps in the function
	
	tcold = ratio*PM
	thot = (1-ratio)*PM
	width = thot
	totalshift = tcold - shift
	tmod = mod(t + totalshift, PM)
	tstart = tcold
	if flank == Inf # rectangle function
		# ceil() with 4 digits precision ensures consistent switching between cold and hot phases
		return ifelse(ceil(tmod, digits=4) < tstart, Tcold, Thot)
	else # smoothed rectangle
		return GasChromatographySystems.smooth_rectangle.(tmod, tstart, width, Tcold, Thot; flank=flank) 
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

# periodic repeated smoothed rectangle function with period 'PM', a shift by 'shift', 'ratio' of time of Tcold to time of Thot. A small shift is incorporated to move the falling flank from the beginning to the end
"""
    mod_number(t, shift, PM, ratio)

Calculate the modulation number for a given time point using the same timing logic as `therm_mod`.

# Arguments
- `t`: Time point at which to calculate the modulation number
- `shift`: Time shift of the modulation pattern (in seconds)
- `PM`: Modulation period (in seconds)
- `ratio`: Ratio of cold phase duration to total period (0 < ratio < 1)

# Returns
- Integer modulation number (1-based)

# Notes
- Uses the same timing logic as `therm_mod`
- Returns 1 for the first modulation period
"""
function mod_number(t, shift, PM, ratio)
    tcold = ratio*PM
    totalshift = tcold - shift
    return cld(t + totalshift, PM)
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
	CAS = pl_end.CAS[[findfirst(Name[i].==pl_end.Name) for i=1:length(Name)]]
	
	sort_heights = Array{Array{Float64}}(undef, nsub)
	sort_tR1 = Array{Array{Float64}}(undef, nsub)
	sort_tR2 = Array{Array{Float64}}(undef, nsub)
	mean_tR1 = Array{Float64}(undef, nsub)
	mean_tR2 = Array{Float64}(undef, nsub)
	sort_τR2 = Array{Array{Float64}}(undef, nsub)
	mean_τR2 = Array{Float64}(undef, nsub)
	for i=1:nsub
		ii_name = findall(Name[i].==pl_end.Name)
		sort_heights[i] = heights[ii_name]
		sort_tR1[i] = tR1[ii_name]
		sort_tR2[i] = tR2[ii_name]
		mean_tR1[i] = sum(sort_tR1[i].*sort_heights[i])/sum(sort_heights[i])
		mean_tR2[i] = sum(sort_tR2[i].*sort_heights[i])/sum(sort_heights[i])
		sort_τR2[i] = pl_end.τR[ii_name]
		mean_τR2[i] = sum(sort_τR2[i].*sort_heights[i])/sum(sort_heights[i])
	end
	return DataFrame(Name=Name, CAS=CAS, tR1=mean_tR1, tR2=mean_tR2, τR2=mean_τR2, tR1s=sort_tR1, tR2s=sort_tR2, τR2s=sort_τR2, heights=sort_heights)
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

end # module