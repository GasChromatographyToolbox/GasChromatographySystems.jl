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

# some constants
const Tst = 273.15 # K
const R = 8.31446261815324 # J mol⁻¹ K⁻¹
const Tn = 25.0 + Tst # K
const pn = 101300 # Pa

# structures

# options structure
struct Options
    mobile_phase::String# gas of the mobile phase
    alg                 # algorithmen for the ODE solver
    abstol              # absolute tolerance for ODE solver
    reltol              # relative tolerance for ODE solver 
    Tcontrol::String    # temperature control at 'inlet' (top) or 'outlet' (bottom) of the column
    odesys::Bool  		# calculate the two ODEs (migration and peak-width) separately (false) or 
                        # combined as a system of ODEs (true)                        
    ng::Bool            # non-gradient calculation, ignores a defined spatial change of d, df or T
    vis::String         # viscosity model 'HP' or 'Blumberg'
    control::String     # control of the 'Flow' or of the inlet 'Pressure' during the program
    k_th                # threshold for the max. possible retention factor
end

function Options(;mobile_phase="He", alg=OwrenZen5(), abstol=1e-6, reltol=1e-3, Tcontrol="inlet", odesys=true, ng=false, vis="Blumberg", control="Pressure", k_th=1e12)
    opt = Options(mobile_phase, alg, abstol, reltol, Tcontrol, odesys, ng, vis, control, k_th)
    return opt
end

# module structure
abstract type AbstractModule end

struct ModuleColumn<:GasChromatographySystems.AbstractModule
	# Module
	# GC column, gradients are possible
	name::String
	length::Float64
	diameter#::Fd # Function
	a_diameter::Array{Float64,1} # Parameters of diameter function, just for information
	film_thickness#::Fdf # Function
	a_film_thickness::Array{Float64,1} # Parameters of film_thickness function, just for information
	stationary_phase::String
	temperature # a number (constant temperature) or a TemperatureProgram structure
	flow # an number (constant flow) or a Function
end

function ModuleColumn(name, L, d, df, sp, tp)
	# function to construct the Column structure
	# for the case of constant diameter and constant film thickness
	# and undefined flow
	col = ModuleColumn(name, L, d, [d], df, [df], sp, tp, NaN)
	return col
end

function ModuleColumn(name, L, d, df, sp, tp, flow)
	# function to construct the Column structure
	# for the case of constant diameter and constant film thickness
	col = ModuleColumn(name, L, d, [d], df, [df], sp, tp, flow)
	return col
end

# add a thermal modulator module

struct ModuleTM<:GasChromatographySystems.AbstractModule
	# Module
	# thermal modulator
	name::String
	length::Float64
	diameter#::Fd # Function
	a_diameter::Array{Float64,1} # Parameters of diameter function, just for information
	film_thickness#::Fdf # Function
	a_film_thickness::Array{Float64,1} # Parameters of film_thickness function, just for information
	stationary_phase::String
	temperature # a number (constant temperature) or a TemperatureProgram structure
	PM::Float64 # a number, modulation periode 
	ratio::Float64 # a number, ratio of the duration between hot and cold jet, approx. as rectangular function
	Thot::Float64 # heating with hot jet
	Tcold::Float64 # cooling with cold jet
	flow # an number (constant flow) or a Function
end

function ModuleTM(name, L, d, df, sp, tp, pm, ratio, Thot, Tcold)
	TM = ModuleTM(name, L, d, [d], df, [df], sp, tp, pm, ratio, Thot, Tcold, NaN)
	return TM
end

function ModuleTM(name, L, d, df, sp, tp, pm, ratio, Thot, Tcold, F)
	TM = ModuleTM(name, L, d, [d], df, [df], sp, tp, pm, ratio, Thot, Tcold, F)
	return TM
end

# add a flow modulator module

# temperature program structure
struct TemperatureProgram{F<:Function}
    timesteps::Array{<:Real,1}
    temperaturesteps::Array{<:Real,1}
    gradient_function::F
    a_gradient_function::Array{<:Real,2} # Parameters of the gradient_function, just for information
    TemperatureProgram(ts,Ts,gf,a) = (length(ts)!=length(Ts) || length(ts)!=length(gf(0.0))) ? error("Mismatch between length(timesteps) = $(length(ts)), length(temperaturesteps) = $(length(Ts)) and length(gradient_function(0.0)) = $(length(gf(0.0)))") : new{typeof(gf)}(ts,Ts,gf,a)
end

function TemperatureProgram(timesteps::Array{<:Real, 1}, temperaturesteps::Array{<:Real, 1})
    # function to construct the TEmperature Program structure
    # without a thermal gradient
    # using as gradient function the exponential model 'gf_exp(x,a_gf,Tcontrol)'
    a_gf = [zeros(length(timesteps)) zeros(length(timesteps)) ones(length(timesteps)) zeros(length(timesteps))]
    gf(x) = GasChromatographySimulator.gradient(x, a_gf)
    prog = TemperatureProgram(timesteps, temperaturesteps, gf, a_gf)
    return prog
end

default_TP() = GasChromatographySystems.TemperatureProgram([0.0, 1800.0], [40.0, 340.0])

# pressure point structure
struct PressurePoint
	# Pressure program, same structure for inlet and outlet
	name::String
	timesteps::Array{Float64,1}
	pressure_steps::Array{Float64,1}
	PressurePoint(n,ts,ps) = length(ts)!=length(ps) ? error("Mismatch between length(timesteps) = $(length(ts)), length(pressure_steps) = $(length(ps))") : new(n,ts,ps)
end

# system structure
struct System
	g::Graphs.SimpleDiGraph{Int}
	pressurepoints::Array{PressurePoint}
	modules::Array{AbstractModule}
	options::Options
	System(g_,pressurepoints_,modules_,options_) = nv(g_)!=length(pressurepoints_) || ne(g_)!=length(modules_) ? error("Mismatch between number of nodes ($(nv(g_))) and number of pressure points ($(length(pressurepoints_))) and/or mismatch between number of edges ($(ne(g_))) and number of modules ($(length(modules_))).") : new(g_,pressurepoints_,modules_,options_)
end

# functions

# rectangular temperature modulation
function therm_mod(t, shift, PM, ratio, Thot, Tcold) 
	return ifelse(mod(t + shift, PM) < ratio*PM, Thot, Tcold)
end

# common programs
function common_timesteps(sys)
	com_timesteps = []
	for i=1:nv(sys.g)
		com_timesteps = GasChromatographySimulator.common_time_steps(com_timesteps, sys.pressurepoints[i].timesteps)
	end
	for i=1:ne(sys.g)
		if typeof(sys.modules[i].temperature) <: TemperatureProgram
			com_timesteps = GasChromatographySimulator.common_time_steps(com_timesteps, sys.modules[i].temperature.timesteps)
		end
	end
	return com_timesteps
end

function index_modules_with_temperature_program(sys)
	i_tempprog = Int[]
	for i=1:ne(sys.g)
		if typeof(sys.modules[i].temperature) <: TemperatureProgram
			push!(i_tempprog, i)
		end
	end
	return i_tempprog
end

function match_programs(sys)
	com_times = common_timesteps(sys)
	new_press_steps = Array{Array{Float64,1}}(undef, nv(sys.g))
	for i=1:nv(sys.g)
		new_press_steps[i] = GasChromatographySimulator.new_value_steps(sys.pressurepoints[i].pressure_steps, sys.pressurepoints[i].timesteps, com_times)
	end
	i_tempprog = index_modules_with_temperature_program(sys)
	new_temp_steps = Array{Array{Float64,1}}(undef, length(i_tempprog))
	for i=1:length(i_tempprog)
		new_temp_steps[i] = GasChromatographySimulator.new_value_steps(sys.modules[i_tempprog[i]].temperature.temperaturesteps, sys.modules[i_tempprog[i]].temperature.timesteps, com_times)
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
		if typeof(sys.modules[i].temperature) <: Number
			new_modules[i] = sys.modules[i]
		elseif typeof(sys.modules[i].temperature) <: TemperatureProgram
			# add/modify for gradient
			ii = findfirst(index_module_tempprog.==i)
			new_tp = TemperatureProgram(new_timesteps, new_temperaturesteps[ii])
			if typeof(sys.modules[i]) == ModuleColumn
				new_modules[i] = ModuleColumn(sys.modules[i].name, sys.modules[i].length, sys.modules[i].diameter, sys.modules[i].film_thickness, sys.modules[i].stationary_phase, new_tp, sys.modules[i].flow)
			elseif typeof(sys.modules[i]) == ModuleTM
				new_modules[i] = ModuleTM(sys.modules[i].name, sys.modules[i].length, sys.modules[i].diameter, sys.modules[i].film_thickness, sys.modules[i].stationary_phase, new_tp, sys.modules[i].PM, sys.modules[i].ratio, sys.modules[i].Thot, sys.modules[i].Tcold, sys.modules[i].flow)
			end
		end
	end
	new_sys = System(sys.g, new_pp, new_modules, sys.options)
    return new_sys
end

# flow balance and pressure calculations
function flow_balance(g, i_n, F)
	#@variables P²[1:nv(g)], κ[1:ne(g)]
	E = collect(edges(g)) # all edges
	srcE = src.(E) # source nodes of the edges
	dstE = dst.(E) # destination nodes of the edges
	# find edges, where node `i_n` is the source
	i_src = findall(srcE.==i_n)
	# find edges, where node `i_n` is the destination
	i_dst = findall(dstE.==i_n)
	balance = 0
	for j=1:length(i_dst) # ingoing flows
		balance = balance + F[i_dst[j]]
	end
	for j=1:length(i_src) # outgoing flows
		balance = balance - F[i_src[j]]
	end
	return balance
end

#=function flow_balance(g, i_n, P², κ)
	#@variables P²[1:nv(g)], κ[1:ne(g)]
	E = collect(edges(g))
	srcE = src.(E) # source nodes of the edges
	dstE = dst.(E) # destination nodes of the edges
	# find edges, where node `i_n` is the source
	i_src = findall(srcE.==i_n)
	# find edges, where node `i_n` is the destination
	i_dst = findall(dstE.==i_n)
	balance = 0
	for j=1:length(i_dst) # ingoing flows
		balance = balance + (P²[srcE[i_dst[j]]]-P²[dstE[i_dst[j]]])/κ[i_dst[j]]
	end
	for j=1:length(i_src) # outgoing flows
		balance = balance - (P²[srcE[i_src[j]]]-P²[dstE[i_src[j]]])/κ[i_src[j]]
	end
	return balance
end=#

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

# first construct the flow balance equations only using the flows over the edges
function flow_balance(sys)
	@variables F[1:ne(sys.g)]
	inner_V = GasChromatographySystems.inner_vertices(sys.g) # one flow balance equation for every inner node
	bal_eq = Array{Symbolics.Equation}(undef, length(inner_V))
	for i=1:length(inner_V)
		bal_eq[i] = flow_balance(sys.g, inner_V[i], F) ~ 0
	end
	return bal_eq
end

function unknown_F(sys)
	i_unknown_flows = Int[]
	for i=1:ne(sys.g)
		if isnan(sys.modules[i].flow)
			push!(i_unknown_flows, i)
		end
	end
	return i_unknown_flows
end

function unknown_p(sys)
	i_unknown_pressures = Int[]
	for i=1:nv(sys.g)
		if isnan(sys.pressurepoints[i].pressure_steps[1])
			push!(i_unknown_pressures, i)
		end
	end
	return i_unknown_pressures
end

# second substitute the unknowns in the flow balance equations and add equations for the known flows
function substitute_unknown_flows(sys)
	@variables A, P²[1:nv(sys.g)], κ[1:ne(sys.g)], F[1:ne(sys.g)]
	E = collect(edges(sys.g)) # all edges
	srcE = src.(E) # source nodes of the edges
	dstE = dst.(E) # destination nodes of the edges
	i_unknown_F = unknown_F(sys) # indices of the modules with an unknown flow
	# create dictionary for the substitution of the unknown flows
	sub_dict = Dict()
	for i=1:length(i_unknown_F)
		j = i_unknown_F[i]
		sub_dict = merge(sub_dict, Dict(F[j] => A*(P²[srcE[j]]-P²[dstE[j]])/κ[j]))
	end
	# index of the known flows
	i_known_F = collect(1:length(edges(sys.g)))[Not(unknown_F(sys))]
	# substitute the unknown flows in all balance equations
	bal_eq = flow_balance(sys)
	sub_bal_eq = Array{Equation}(undef, length(bal_eq)+length(i_known_F))
	for i=1:length(bal_eq)
		sub_bal_eq[i] = substitute(bal_eq[i], sub_dict)
	end
	for i=1:length(i_known_F)
		j = i_known_F[i]
		sub_bal_eq[length(bal_eq)+i] = F[j] - A*(P²[srcE[j]]-P²[dstE[j]])/κ[j] ~ 0
	end
	return sub_bal_eq
end

# not needed???
function unknowns_in_flow_balances(sys)
	@variables P²[1:nv(sys.g)], κ[1:ne(sys.g)]
	i_unknown_p = GasChromatographySystems.unknown_p(sys)
	F_balance = GasChromatographySystems.flow_balance(sys)
	in_equation = Array{Array{Int64,1}}(undef, length(i_unknown_p))
	touched = zeros(Int, length(F_balance))
	for i=1:length(i_unknown_p)
		var_unknown = Symbolics.get_variables(P²[i_unknown_p[i]])
		in_equation_ = Int[]
		for j=1:length(F_balance)
			var_equation = Symbolics.get_variables(F_balance[j])
			counter_touched = touched[j]
			for k=1:length(var_equation)
				if isequal(var_unknown[1], var_equation[k])
					push!(in_equation_, j)
					counter_touched += 1
				end
			end
			touched[j] = counter_touched
		end
		in_equation[i] = in_equation_
	end
	# in_equation ... Array of Arrays, lists the idices of the flow balances in which the unknown is in it.
	# touched ... total number of different unknowns in the flow balance equations
	return in_equation, touched
end

function solve_balance(sys)
	@variables A, P²[1:nv(sys.g)], κ[1:ne(sys.g)], F[1:ne(sys.g)]
	i_unknown_p = GasChromatographySystems.unknown_p(sys) # indices of the nodes with unknown pressure
	#i_unknown_F = unknown_F(sys) # indices of the edges with an unknown flow
	bal_eq = substitute_unknown_flows(sys)
	#num_use_eq = GasChromatographySystems.unknowns_in_flow_balances(sys)[2]
	if length(i_unknown_p) == length(bal_eq)
		sol = Symbolics.solve_for(bal_eq, [P²[i_unknown_p[i]] for i=1:length(i_unknown_p)])
	#elseif length(i_unknown_p) == 1
	#	sol = Symbolics.solve_for(F_balance[i_unknown_p[1]-1], [P²[i_unknown_p[1]]])
	#elseif length(i_unknown_p) == (length(F_balance) - 1)
		#if length(findall(minimum(num_use_eq).==num_use_eq)) == 1
	#		leave_out_eq = findlast(num_use_eq.==minimum(num_use_eq))
		#else
		#	leave_out_eq = 
		#end
	#	sol = Symbolics.solve_for(F_balance[Not(leave_out_eq)], [P²[i_unknown_p[i]] for i=1:length(i_unknown_p)])
	#elseif length(i_unknown_p) < (length(F_balance) - 1)
	#	error("More flow balance equations than unknown pressures. ToDo: leave equations out.")
	elseif length(i_unknown_p) > length(bal_eq)
		error("More unknown pressures than flow balance equations.")
	else # loop of the i_unknown_p -> is this really used???
		# identifie inner nodes which are unknown pressures, there index in the inner_vertices() is the index of the flow balance equations to use 
		# this works only for unknown_p which are inner nodes, unknown_p at outer nodes (e.g. inlet pressure), lead to error 
		inner_V = GasChromatographySystems.inner_vertices(sys.g) 
		bal_eq_i = Array{Int}(undef, length(i_unknown_p))
		for i=1:length(i_unknown_p)
			# if the unknown pressure is not an inner node, than add a equation from the end of the list of balance equation, which should be a flow definition.
			if isnothing(findfirst(i_unknown_p[i].==inner_V))
				bal_eq_i[i] = length(bal_eq)-(i+0)
			else
				bal_eq_i[i] = findfirst(i_unknown_p[i].==inner_V)
			end
		end
		sol = Symbolics.solve_for([bal_eq[x] for x in bal_eq_i], [P²[i_unknown_p[i]] for i=1:length(i_unknown_p)])
	end
	return sol
end

function flow_restriction(t, module_, options)
	if typeof(module_.temperature) <: TemperatureProgram
		T_itp = GasChromatographySimulator.temperature_interpolation(module_.temperature.timesteps, module_.temperature.temperaturesteps, module_.temperature.gradient_function, module_.length)
	elseif typeof(module_.temperature) <: Number
		gf(x) = [zero(x), zero(x)]
		T_itp = GasChromatographySimulator.temperature_interpolation([0.0, 3600.0], [module_.temperature, module_.temperature], gf, module_.length)
	end
	κ = GasChromatographySimulator.flow_restriction(module_.length, t, T_itp, module_.diameter, options.mobile_phase; ng=options.ng, vis=options.vis)
	return κ
end

function flow_restrictions(sys)
	kappas = Array{Function}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		f(t) = flow_restriction(t, sys.modules[i], sys.options)
		kappas[i] = f
	end
	return kappas
end

function pressures_squared(sys)
	#p² = Array{Interpolations.Extrapolation}(undef, nv(sys.g))
	p² = Array{Any}(undef, nv(sys.g))
	for i=1:nv(sys.g)
		g = if isnan(sys.pressurepoints[i].pressure_steps[1])
			GasChromatographySimulator.steps_interpolation(sys.pressurepoints[i].timesteps, NaN.*ones(length(sys.pressurepoints[i].timesteps)))
		else
			GasChromatographySimulator.steps_interpolation(sys.pressurepoints[i].timesteps, identity.(sys.pressurepoints[i].pressure_steps))
		end
		f(t) = g(t).^2 
		p²[i] = f
	end
	return p²
end

function solve_pressure(sys)
	@variables A, P²[1:nv(sys.g)], κ[1:ne(sys.g)], F[1:ne(sys.g)]
	#balance = flow_balance(sys)
	i_unknown_p = GasChromatographySystems.unknown_p(sys)
	i_known_F = collect(1:length(edges(sys.g)))[Not(unknown_F(sys))]
	solutions = solve_balance(sys)
	a = π/256 * GasChromatographySystems.Tn/GasChromatographySystems.pn
	κs = GasChromatographySystems.flow_restrictions(sys)
	p²s = GasChromatographySystems.pressures_squared(sys)
	p_solution = Array{Function}(undef, length(i_unknown_p))
	for i=1:length(i_unknown_p)
		sub_dict(t) = merge(Dict((P²[j] => p²s[j](t) for j=setdiff(1:nv(sys.g), i_unknown_p))), Dict((κ[j] => κs[j](t) for j=1:ne(sys.g))), Dict(A => a), Dict(F[j] => sys.modules[j].flow for j in i_known_F))
		f(t) = sqrt(substitute(solutions[i], sub_dict(t)))
		p_solution[i] = f
	end
	return p_solution, i_unknown_p
end

function pressure_functions(sys)
	pres, unk = solve_pressure(sys)
	#p²s = pressures_squared(sys)
	p_func = Array{Any}(undef, nv(sys.g))
	for i=1:nv(sys.g)
		f = if i in unk
			pres[findfirst(unk.==i)]
		else
			GasChromatographySimulator.steps_interpolation(sys.pressurepoints[i].timesteps, identity.(sys.pressurepoints[i].pressure_steps))
		end
	p_func[i] = f
	end
	return p_func
end

function interpolate_pressure_functions(sys; dt=1)
	tsteps = GasChromatographySystems.common_timesteps(sys)
	tend = sum(tsteps)
	if all(degree(sys.g).<3) # no split/merge nodes, straight graph -> pressure at nodes is linear
		trange = cumsum(tsteps)
	else
		trange = 0:dt:tend
	end
	p_func = GasChromatographySystems.pressure_functions(sys)
	p_itp = Array{Any}(undef, length(p_func))
	for i=1:length(p_func)
		if all(isnan.(sys.pressurepoints[i].pressure_steps))
			p_itp[i] = LinearInterpolation((trange, ), p_func[i].(trange), extrapolation_bc=Flat())
		else
			p_itp[i] = p_func[i]
		end
	end
	return p_itp
end

function flow_functions(sys)
	p_func = pressure_functions(sys)
	F_func = Array{Function}(undef, ne(sys.g))
	E = collect(edges(sys.g))
	srcE = src.(E)
	dstE = dst.(E)
	for i=1:ne(sys.g)
		pin(t) = p_func[srcE[i]](t)
		pout(t) = p_func[dstE[i]](t)
		if typeof(sys.modules[i].temperature) <: GasChromatographySystems.TemperatureProgram
			T_itp = GasChromatographySimulator.temperature_interpolation(sys.modules[i].temperature.timesteps, sys.modules[i].temperature.temperaturesteps, sys.modules[i].temperature.gradient_function, sys.modules[i].length)
		elseif typeof(sys.modules[i].temperature) <: Number
			gf(x) = [zero(x), zero(x)]
			T_itp = GasChromatographySimulator.temperature_interpolation([0.0, 3600.0], [sys.modules[i].temperature, sys.modules[i].temperature], gf, sys.modules[i].length)
		end
		f(t) = GasChromatographySimulator.flow(t, T_itp, pin, pout, sys.modules[i].length, sys.modules[i].diameter, sys.options.mobile_phase; ng=sys.options.ng, vis=sys.options.vis, control=sys.options.control)
		F_func[i] = f
	end
	return F_func
end

# plotting of the graphs
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
		stat_phases[i] = sys.modules[i].stationary_phase
	end
	return stat_phases
end

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
	E = collect(edges(sys.g))
	srcE = src.(E)
	dstE = dst.(E)
	if interp == true
		p_func = interpolate_pressure_functions(sys; dt=dt)
	else
		p_func = pressure_functions(sys)
	end
	parameters = Array{GasChromatographySimulator.Parameters}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		col = GasChromatographySimulator.Column(sys.modules[i].length, sys.modules[i].diameter, [sys.modules[i].diameter], sys.modules[i].film_thickness, [sys.modules[i].film_thickness], sys.modules[i].stationary_phase, sys.options.mobile_phase)

		pin_steps = sys.pressurepoints[srcE[i]].pressure_steps
		pout_steps = sys.pressurepoints[dstE[i]].pressure_steps
		pin_itp = p_func[srcE[i]]
		pout_itp = p_func[dstE[i]]	
		if typeof(sys.modules[i].temperature) <: TemperatureProgram
			time_steps = sys.modules[i].temperature.timesteps
			temp_steps = sys.modules[i].temperature.temperaturesteps	
			gf = sys.modules[i].temperature.gradient_function
			a_gf = sys.modules[i].temperature.a_gradient_function
			T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, sys.modules[i].length)		
		elseif typeof(sys.modules[i].temperature) <: Number
			time_steps = common_timesteps(sys)
			temp_steps = sys.modules[i].temperature.*ones(length(time_steps))
			gf(x) = zero(x).*ones(length(time_steps))
			a_gf = [zeros(length(time_steps)) zeros(length(time_steps)) ones(length(time_steps)) zeros(length(time_steps))]
			T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, sys.modules[i].length)
		end
		prog = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)

		sub = GasChromatographySimulator.load_solute_database(db_dataframe, sys.modules[i].stationary_phase, sys.options.mobile_phase, selected_solutes, NaN.*ones(length(selected_solutes)), NaN.*ones(length(selected_solutes)))

		opt = GasChromatographySimulator.Options(sys.options.alg, sys.options.abstol, sys.options.reltol, sys.options.Tcontrol, sys.options.odesys, sys.options.ng, sys.options.vis, sys.options.control, sys.options.k_th)

		parameters[i] = GasChromatographySimulator.Parameters(col, prog, sub, opt)
	end
	return parameters
end

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
	#i_par = Array{Int}(undef, length(path))
	#for i=1:length(path)
	#	i_par[i] = findall(src.(path)[i].==src.(edges(g)))[findfirst(findall(src.(path)[i].==src.(edges(g))).==findall(dst.(path)[i].==dst.(edges(g))))]
	#end
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

function change_initial(par::GasChromatographySimulator.Parameters, pl)
	# copys the parameters `par` and changes the values of par.sub[i].τ₀ and par.sub[i].t₀ to pl.tR[], pl.τR[]
	CAS_pl = GasChromatographySimulator.CAS_identification(pl.Name).CAS
	newsub = Array{GasChromatographySimulator.Substance}(undef, length(par.sub))
	for i=1:length(par.sub)
		ii = findfirst(par.sub[i].CAS.==CAS_pl)
		newsub[i] = GasChromatographySimulator.Substance(par.sub[i].name, par.sub[i].CAS, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀, par.sub[i].ann, par.sub[i].Cag, pl.tR[ii], pl.τR[ii])
	end
	newpar = GasChromatographySimulator.Parameters(par.col, par.prog, newsub, par.opt)
	return newpar
end

function change_initial_focussed(par::GasChromatographySimulator.Parameters, pl; τ₀=zeros(length(pl.tR)))
	# copys the parameters `par` and changes the values of par.sub[i].t₀ to pl.tR[]
	CAS_pl = GasChromatographySimulator.CAS_identification(pl.Name).CAS
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

function simulate_along_paths(sys, paths, par_sys; t₀=zeros(length(par_sys[1].sub)), τ₀=zeros(length(par_sys[1].sub)), refocus=falses(ne(sys.g)), τ₀_focus=zeros(length(par_sys[1].sub)))
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
					peaklists_[j] = peaklists[i_path][i_edge]
					solutions_[j] = solutions[i_path][i_edge]
				else
					if j == 1
						new_par_sys[i_par[j]] = GasChromatographySystems.change_initial(par_sys[i_par[j]], t₀, τ₀)
						peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(new_par_sys[i_par[j]])
					else
						if refocus[i_par[j]] == true
							new_par_sys[i_par[j]] = change_initial_focussed(par_sys[i_par[j]], peaklists_[j-1]; τ₀=τ₀_focus)
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

function peaklist_GCxGC(pl_1, pl_2)
	pl_GCxGC = DataFrame(Name = pl_1.Name, tR1 = pl_1.tR, τR1 = pl_1.τR)
	tR2 = Array{Float64}(undef, length(pl_GCxGC.Name))
	τR2 = Array{Float64}(undef, length(pl_GCxGC.Name))
	CAS1 = GasChromatographySimulator.CAS_identification(pl_1.Name).CAS
	CAS2 = GasChromatographySimulator.CAS_identification(pl_2.Name).CAS
	for i=1:length(pl_GCxGC.Name)
		ii = findfirst(CAS1[i].==CAS2)
		tR2[i] = pl_2.tR[ii]
		τR2[i] = pl_2.τR[ii]
	end
	pl_GCxGC[!, :tR2] = tR2
	pl_GCxGC[!, :τR2] = τR2
	return pl_GCxGC
end

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

function SeriesSystem(Ls, ds, dfs, sps, TPs, F, pin, pout; opt=GasChromatographySystems.Options(ng=true))
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
			com_timesteps = GasChromatographySimulator.common_time_steps(com_timesteps, TPs[i].timesteps)
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
		modules[i] = GasChromatographySystems.ModuleColumn("$(i) -> $(i+1)", Ls[i], ds[i]*1e-3, dfs[i]*1e-6, sps[i], TPs[i], F/60e6)
	end
	# system
	sys = GasChromatographySystems.update_system(GasChromatographySystems.System(g, pp, modules, opt))

	# add test for the defined pressures and flows
	return sys
end

example_SeriesSystem() = SeriesSystem([10.0, 5.0, 2.0, 1.0], [0.53, 0.32, 0.25, 0.1], [0.53, 0.32, 0.25, 0.1], ["Rxi17SilMS", "Rxi17SilMS", "Rxi17SilMS", "Rxi17SilMS"], [default_TP(), default_TP(), default_TP(), default_TP()], NaN, 300.0, 0.0)

function SplitSystem(Ls, ds, dfs, sps, TPs, Fs, pin, pout1, pout2; opt=GasChromatographySystems.Options(ng=true))
	g = SimpleDiGraph(4)
	add_edge!(g, 1, 2) # Inj -> GC column -> Split point
	add_edge!(g, 2, 3) # Split point -> TL column -> Det 1
	add_edge!(g, 2, 4) # Split point -> TL column -> Det 2
	# common time steps
	com_timesteps = []
	for i=1:length(TPs)
		if typeof(TPs[i]) <: GasChromatographySystems.TemperatureProgram
			com_timesteps = GasChromatographySimulator.common_time_steps(com_timesteps, TPs[i].timesteps)
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
	modules[1] = GasChromatographySystems.ModuleColumn("1 -> 2", Ls[1], ds[1]*1e-3, dfs[1]*1e-6, sps[1], TPs[1], Fs[1]/60e6)
	modules[2] = GasChromatographySystems.ModuleColumn("2 -> 3", Ls[2], ds[2]*1e-3, dfs[2]*1e-6, sps[2], TPs[2], Fs[2]/60e6)
	modules[3] = GasChromatographySystems.ModuleColumn("2 -> 4", Ls[3], ds[3]*1e-3, dfs[3]*1e-6, sps[3], TPs[3], Fs[3]/60e6)
	# system
	sys = GasChromatographySystems.update_system(GasChromatographySystems.System(g, pp, modules, opt))

	# add test for the defined pressures and flows
	return sys
end

example_SplitSystem() = SplitSystem([10.0, 1.0, 5.0], [0.25, 0.1, 0.25], [0.25, 0.0, 0.0], ["Rxi17SilMS", "", ""], [default_TP(), 300.0, 300.0], [1.0, NaN, NaN], NaN, 0.0, 101.3)

function GCxGC_TM_simp(L1, d1, df1, sp1, TP1, L2, d2, df2, sp2, TP2, F, pin, pout; opt=GasChromatographySystems.Options(ng=true))
	# ? make two versions
	# 1. defining flow over the columns (calculate pin)
	# 2. defining inlet pressure (calculate F)
	Ls = [L1, L2]
	ds = [d1, d2]
	dfs = [df1, df2]
	sps = [sp1, sp2]
	TPs = [TP1, TP2]
	n = 2
	g = SimpleDiGraph(n+1)
	for i=1:n
		add_edge!(g, i, i+1) 
	end
	# common time steps
	com_timesteps = []
	for i=1:length(TPs)
		if typeof(TPs[i]) <: GasChromatographySystems.TemperatureProgram
			com_timesteps = GasChromatographySimulator.common_time_steps(com_timesteps, TPs[i].timesteps)
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
	pp[2] = GasChromatographySystems.PressurePoint("p2", com_timesteps, nans) #
	pp[3] = GasChromatographySystems.PressurePoint("p3", com_timesteps, pouts) # outlet
	# modules
	modules = Array{GasChromatographySystems.AbstractModule}(undef, n)
	modules[1] = GasChromatographySystems.ModuleColumn("GC1", Ls[1], ds[1]*1e-3, dfs[1]*1e-6, sps[1], TPs[1], F/60e6)
	modules[2] = GasChromatographySystems.ModuleColumn("GC2", Ls[2], ds[2]*1e-3, dfs[2]*1e-6, sps[2], TPs[2], NaN/60e6)
	# system
	sys = GasChromatographySystems.update_system(GasChromatographySystems.System(g, pp, modules, opt))
	# add test for the defined pressures and flows
	return sys
end

example_GCxGC_TM_simp() = GCxGC_TM_simp(30.0, 0.25, 0.25, "Rxi17SilMS", default_TP(), 2.0, 0.25, 0.25, "Wax", default_TP(), NaN, 200.0, 0.0)

function GCxGC_FM_simp(L1, d1, df1, sp1, TP1, F1, L2, d2, df2, sp2, TP2, F2, pin, pout)
	# ? make two versions
	# 1. defining flow over the columns (calculate pin)
	# 2. defining inlet pressure (calculate F)
	g = SimpleDiGraph(4)
	add_edge!(g, 1, 3) # Inj -> GC1 -> FM
	add_edge!(g, 2, 3) # Modpress -> TL -> FM
	add_edge!(g, 3, 4) # FM -> GC2 -> Det
	# common time steps
	com_timesteps = []
	TPs = [TP1, TP2]
	for i=1:length(TPs)
		if typeof(TPs[i]) <: GasChromatographySystems.TemperatureProgram
			com_timesteps = GasChromatographySimulator.common_time_steps(com_timesteps, TPs[i].timesteps)
		end
	end
	if isempty(com_timesteps)
		com_timesteps = [0.0, 36000.0]
	end
	# pressure points
	pp = Array{GasChromatographySystems.PressurePoint}(undef, nv(g))
	pins = pin*1000.0.*ones(length(com_timesteps))
	nans = NaN.*ones(length(com_timesteps))
	if pout == 0.0
		pouts = eps(Float64).*ones(length(com_timesteps))
	else 
		pouts = pout*1000.0.*ones(length(com_timesteps))
	end
	pp[1] = GasChromatographySystems.PressurePoint("p₁", com_timesteps, pins) # inlet 
	pp[2] = GasChromatographySystems.PressurePoint("p₂", com_timesteps, nans) # FM inlet
	pp[3] = GasChromatographySystems.PressurePoint("p₃", com_timesteps, nans) # FM
	pp[4] = GasChromatographySystems.PressurePoint("p₄", com_timesteps, pouts) # outlet
	# modules
	modules = Array{GasChromatographySystems.AbstractModule}(undef, ne(g))
	modules[1] = GasChromatographySystems.ModuleColumn("GC1", L1, d1*1e-3, df1*1e-6, sp1, TP1, F1/60e6)
	modules[2] = GasChromatographySystems.ModuleColumn("FM inlet", 0.1, 0.5*1e-3, 0.0, "", TP1, NaN/60e6)
	modules[3] = GasChromatographySystems.ModuleColumn("GC2", L2, d2*1e-3, df2*1e-6, sp2, TP2, F2/60e6)
	# system
	sys = GasChromatographySystems.update_system(GasChromatographySystems.System(g, pp, modules, GasChromatographySystems.Options(ng=true)))

	# add test for the defined pressures and flows
	return sys
end

example_GCxGC_FM_simp() = GCxGC_FM_simp(30.0, 0.25, 0.25, "Rxi17SilMS", default_TP(), 1.0, 2.0, 0.25, 0.25, "Wax", default_TP(), 2.0, NaN, 0.0)

end # module