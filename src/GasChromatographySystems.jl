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

# structures

# options structure
struct Options
    mobile_phase::String# gas of the mobile phase
    odesys::Bool  		# calculate the two ODEs (migration and peak-width) separately (false) or 
                        # combined as a system of ODEs (true)                        
    vis::String         # viscosity model 'HP' or 'Blumberg'
    control::String     # control of the 'Flow' or of the inlet 'Pressure' during the program
    k_th                # threshold for the max. possible retention factor
end

function Options(;mobile_phase="He", odesys=true, vis="Blumberg", control="Pressure", k_th=1e12)
    opt = Options(mobile_phase, odesys, vis, control, k_th)
    return opt
end

# module structure
abstract type AbstractModule end

struct ModuleColumnOpt
	alg                 # algorithmen for the ODE solver
    abstol              # absolute tolerance for ODE solver
    reltol              # relative tolerance for ODE solver 
	ng::Bool            # non-gradient calculation, ignores a defined spatial change of d, df or T
	Tcontrol::String    # temperature control at 'inlet' (top) or 'outlet' (bottom) of the column
end

function ModuleColumnOpt(; alg=OwrenZen5(), abstol=1e-8, reltol=1e-6, ng=false, Tcontrol="inlet")
	ModuleColumnOpt(alg, abstol, reltol, ng, Tcontrol)
end

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
	opt::ModuleColumnOpt
end

function ModuleColumn(name, L, d, df, sp, tp, opt::ModuleColumnOpt)
	# function to construct the Column structure
	# for the case of constant diameter and constant film thickness
	# and undefined flow
	col = ModuleColumn(name, L, d, [d], df, [df], sp, tp, NaN, opt)
	return col
end

function ModuleColumn(name, L, d, df, sp, tp, flow, opt::ModuleColumnOpt)
	# function to construct the Column structure
	# for the case of constant diameter and constant film thickness
	col = ModuleColumn(name, L, d, [d], df, [df], sp, tp, flow, opt)
	return col
end

function ModuleColumn(name, L, d, df, sp, tp; alg=OwrenZen5(), abstol=1e-8, reltol=1e-6, ng=false, Tcontrol="inlet")
	opt = ModuleColumnOpt(; alg=alg, abstol=abstol, reltol=reltol, ng=ng, Tcontrol=Tcontrol)
	TM = ModuleColumn(name, L, d, [d], df, [df], sp, tp, NaN, opt)
	return TM
end

function ModuleColumn(name, L, d, df, sp, tp, flow; alg=OwrenZen5(), abstol=1e-8, reltol=1e-6, ng=false, Tcontrol="inlet")
	opt = ModuleColumnOpt(; alg=alg, abstol=abstol, reltol=reltol, ng=ng, Tcontrol=Tcontrol)
	TM = ModuleColumn(name, L, d, [d], df, [df], sp, tp, flow, opt)
	return TM
end

# add a thermal modulator module
struct ModuleTMopt
	Tcold_abs::Bool 	# Tcold as absolute value (`true`) or as relative difference to the program temperature (`false`) 
	sflank::Float64 	# flank factor for the spatial smoothed rectangle, values between 12 and Inf? 
	tflank::Float64 	# flank factor for the temporal smoothed rectangle, values between 12 and Inf?
	alg 				# solver algorithm for the simulation, typical `Vern9()`, `OwrenZen5()`. With the option "simplifiedTM" no solving of ODEs is used for the modulator spot but an approximation is used, assuming a rectangle modulation.
	abstol::Float64 	# absolute tolerance for the solving of the ODEs, typical value 1e-10
	reltol::Float64 	# relativ tolerance for solving of the ODEs, typical value 1e-8
	dtinit::Float64 	# initial step width for the solving of the ODEs, typical value `module length * 1e-6
	ng::Bool 			# model modulation spot also as smoothed rectangle over the length (`false`) or as uniform (`true`)
	Tcontrol::String    # temperature control at 'inlet' (top) or 'outlet' (bottom) of the column
end

function ModuleTMopt(; Tcold_abs=true, sflank=40, tflank=20, alg=Vern9(), abstol=1e-10, reltol=1e-8, dtinit=1e-6, ng=true, Tcontrol="inlet")
	ModuleTMopt(Tcold_abs, sflank, tflank, alg, abstol, reltol, dtinit, ng, Tcontrol)
end

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
	shift::Float64
	PM::Float64 # a number, modulation periode 
	ratio::Float64 # a number, ratio of the duration between cold and hot jet, approx. as rectangular function
	Thot::Float64 # heating with hot jet
	Tcold::Float64 # cooling with cold jet
	flow # an number (constant flow) or a Function
	opt::ModuleTMopt
end

function ModuleTM(name, L, d, df, sp, tp, shift, pm, ratio, Thot, Tcold, opt::ModuleTMopt)
	TM = ModuleTM(name, L, d, [d], df, [df], sp, tp, shift, pm, ratio, Thot, Tcold, NaN, opt)
	return TM
end

function ModuleTM(name, L, d, df, sp, tp, shift, pm, ratio, Thot, Tcold, F, opt::ModuleTMopt)
	TM = ModuleTM(name, L, d, [d], df, [df], sp, tp, shift, pm, ratio, Thot, Tcold, F, opt)
	return TM
end

function ModuleTM(name, L, d, df, sp, tp, shift, pm, ratio, Thot, Tcold; Tcold_abs=true, sflank=40, tflank=20, alg=Vern9(), abstol=1e-10, reltol=1e-8, dtinit=1e-6, ng=false, Tcontrol="inlet")
	opt = ModuleTMopt(; Tcold_abs=Tcold_abs, sflank=sflank, tflank=tflank, alg=alg, abstol=abstol, reltol=reltol, dtinit=dtinit, ng=ng, Tcontrol=Tcontrol)
	TM = ModuleTM(name, L, d, df, sp, tp, shift, pm, ratio, Thot, Tcold, NaN, opt)
	return TM
end

function ModuleTM(name, L, d, df, sp, tp, shift, pm, ratio, Thot, Tcold, F; Tcold_abs=true, sflank=40, tflank=20, alg=Vern9(), abstol=1e-10, reltol=1e-8, dtinit=1e-6, ng=false, Tcontrol="inlet")
	opt = ModuleTMopt(; Tcold_abs=Tcold_abs, sflank=sflank, tflank=tflank, alg=alg, abstol=abstol, reltol=reltol, dtinit=dtinit, ng=ng, Tcontrol=Tcontrol)
	TM = ModuleTM(name, L, d, df, sp, tp, shift, pm, ratio, Thot, Tcold, F, opt)
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
#function therm_mod(t, shift, PM, ratio, Thot, Tcold) 
#	return ifelse(mod(t + shift, PM) < ratio*PM, Thot, Tcold)
#end

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
				new_modules[i] = ModuleColumn(sys.modules[i].name, sys.modules[i].length, sys.modules[i].diameter, sys.modules[i].film_thickness, sys.modules[i].stationary_phase, new_tp, sys.modules[i].flow, sys.modules[i].opt)
			elseif typeof(sys.modules[i]) == ModuleTM
				new_modules[i] = ModuleTM(sys.modules[i].name, sys.modules[i].length, sys.modules[i].diameter, sys.modules[i].film_thickness, sys.modules[i].stationary_phase, new_tp, sys.modules[i].shift, sys.modules[i].PM, sys.modules[i].ratio, sys.modules[i].Thot, sys.modules[i].Tcold, sys.modules[i].flow, sys.modules[i].opt)
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

function flow_restrictions(sys)
	kappas = Array{Function}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		T_itp = module_temperature(sys.modules[i], sys)[5]
		κ(t) = GasChromatographySimulator.flow_restriction(sys.modules[i].length, t, T_itp, sys.modules[i].diameter, sys.options.mobile_phase; ng=sys.modules[i].opt.ng, vis=sys.options.vis)
		kappas[i] = κ
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
		T_itp = GasChromatographySystems.module_temperature(sys.modules[i], sys)[5]
		f(t) = GasChromatographySimulator.flow(t, T_itp, pin, pout, sys.modules[i].length, sys.modules[i].diameter, sys.options.mobile_phase; ng=sys.modules[i].opt.ng, vis=sys.options.vis, control=sys.options.control)
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

function module_temperature(module_::ModuleColumn, sys)
	if typeof(module_.temperature) <: TemperatureProgram # temperature is a TemperatureProgram
		time_steps = module_.temperature.timesteps
		temp_steps = module_.temperature.temperaturesteps	
		gf = module_.temperature.gradient_function
		a_gf = module_.temperature.a_gradient_function
		T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, module_.length)
	elseif typeof(module_.temperature) <: Number # temperature is a constant value
		time_steps = common_timesteps(sys)
		temp_steps = module_.temperature.*ones(length(time_steps))
		gf(x) = zero(x).*ones(length(time_steps))
		a_gf = [zeros(length(time_steps)) zeros(length(time_steps)) ones(length(time_steps)) zeros(length(time_steps))]
		T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, module_.length)
	end
	return time_steps, temp_steps, gf, a_gf, T_itp
end

function module_temperature(module_::GasChromatographySystems.ModuleTM, sys)
	if typeof(module_.temperature) <: GasChromatographySystems.TemperatureProgram # temperature is a TemperatureProgram
		time_steps = module_.temperature.timesteps
		temp_steps = module_.temperature.temperaturesteps
		gf = module_.temperature.gradient_function
		a_gf = module_.temperature.a_gradient_function
		T_itp_ = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, module_.length)
	elseif typeof(module_.temperature) <: Number # temperature is a constant value
		time_steps = GasChromatographySystems.common_timesteps(sys)
		temp_steps = module_.temperature.*ones(length(time_steps))
		gf(x) = zero(x).*ones(length(time_steps))
		a_gf = [zeros(length(time_steps)) zeros(length(time_steps)) ones(length(time_steps)) zeros(length(time_steps))]
		T_itp_ = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, module_.length)
	end
	T_itp(x,t) = if module_.opt.Tcold_abs == true # cool jet always at Tcold
		therm_mod(t, module_.shift, module_.PM, module_.ratio, module_.Tcold, T_itp_(x, t) .+ module_.Thot .- 273.15; flank=module_.opt.tflank) .+ 273.15 
	else # cool jet always at T_itp_ + Tcold
		therm_mod(t, module_.shift, module_.PM, module_.ratio, T_itp_(x, t) .+ module_.Tcold .- 273.15, T_itp_(x, t) .+ module_.Thot .- 273.15; flank=module_.opt.tflank) .+ 273.15 
	end
	
	spot(x,t) = if module_.opt.ng == false
		GasChromatographySystems.smooth_rectangle(x, 0.0, sys.modules[5].length, T_itp_(x, t), T_itp(x,t); flank=module_.opt.sflank)
	else
		T_itp(x,t)
	end
	return time_steps, temp_steps, gf, a_gf, spot
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
		col = GasChromatographySimulator.Column(sys.modules[i].length, sys.modules[i].diameter, [sys.modules[i].diameter], sys.modules[i].film_thickness, [sys.modules[i].film_thickness], sys.modules[i].stationary_phase, sys.options.mobile_phase)

		# program parameters
		pin_steps = sys.pressurepoints[srcE[i]].pressure_steps, optCol=ModuleColumnOpt()
		pout_steps = sys.pressurepoints[dstE[i]].pressure_steps
		pin_itp = p_func[srcE[i]]
		pout_itp = p_func[dstE[i]]	
		time_steps, temp_steps, gf, a_gf, T_itp = module_temperature(sys.modules[i], sys)
		prog = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)

		# substance parameters
		sub = GasChromatographySimulator.load_solute_database(db_dataframe, sys.modules[i].stationary_phase, sys.options.mobile_phase, selected_solutes, NaN.*ones(length(selected_solutes)), NaN.*ones(length(selected_solutes)))

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
							peaklists_[j], solutions_[j] = approximate_modulator(new_par_sys[i_par[j]], df_A, sys.modules[i_par[j]].PM, sys.modules[i_par[j]].ratio, sys.modules[i_par[j]].shift)
						else
							sol = Array{Any}(undef, length(new_par_sys[i_par[j]].sub))
							for i_sub=1:length(new_par_sys[i_par[j]].sub)
								dt = sys.modules[i_par[j]].opt.dtinit
								sol[i_sub] = GasChromatographySimulator.solving_odesystem_r(new_par_sys[i_par[j]].col, new_par_sys[i_par[j]].prog, new_par_sys[i_par[j]].sub[i_sub], new_par_sys[i_par[j]].opt; kwargsTM..., dt=dt)
								tR = sol[i_sub].u[end][1]
								while (fld(tR, sys.modules[i_par[j]].PM) > fld(new_par_sys[i_par[j]].sub[i_sub].t₀, sys.modules[i_par[j]].PM) && dt > eps()) || (sol[i_sub].retcode != ReturnCode.Success && dt > eps()) # solute elutes not in the same modulation periode or the solving failed
									dt = dt/10 # reduce initial step-width
									#@warn "Retention time $(tR) surpasses modulation period, t₀ = $(new_par_sys[i_par[j]].sub[i_sub].t₀), PM = $(sys.modules[i_par[j]].PM). Initial step-width dt is decreased ($(dt))."
									sol[i_sub] = GasChromatographySimulator.solving_odesystem_r(new_par_sys[i_par[j]].col, new_par_sys[i_par[j]].prog, new_par_sys[i_par[j]].sub[i_sub], new_par_sys[i_par[j]].opt; kwargsTM..., dt=dt)
									tR = sol[i_sub].u[end][1]
								end
								if isnan(tR) # does not work
									@warn "Simulation result is NaN. Approximate modulator."
									sol[i_sub] = approximate_modulator(new_par_sys[i_par[j]], df_A, sys.modules[i_par[j]].PM, sys.modules[i_par[j]].ratio, sys.modules[i_par[j]].shift)[2][i_sub]
								end
							end
							peaklists_[j] = GasChromatographySimulator.peaklist(sol, new_par_sys[i_par[j]])
							GasChromatographySystems.add_A_to_pl!(peaklists_[j], df_A)
							solutions_[j] = sol
						end
						if maximum(peaklists_[j].τR) > sys.modules[i_par[j]].PM
							return @warn "Peak width of focussed peaks > modulation period. Simulation is aborted."
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

function approximate_modulator(par, df_A, PM, ratio, shift)
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

	TR = par.prog.T_itp.(par.col.L, tR)

	kR = Array{Float64}(undef, length(tR))
	uR = Array{Float64}(undef, length(tR))
	τ₀ = Array{Float64}(undef, length(tR))
	for i=1:length(tR)
		i_sub = findfirst(Name[i] .== df_A.Name)
		kR[i] = GasChromatographySimulator.retention_factor(par.col.L, tR[i], par.prog.T_itp, par.col.d, par.col.df, par.sub[i_sub].Tchar, par.sub[i_sub].θchar, par.sub[i_sub].ΔCp, par.sub[i_sub].φ₀)

		rM = GasChromatographySimulator.mobile_phase_residency(par.col.L, tR[i], par.prog.T_itp, par.prog.Fpin_itp, par.prog.pout_itp, par.col.L, par.col.d, par.col.gas; ng=par.opt.ng, vis=par.opt.vis, control=par.opt.control)
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

#=function peaklist_GCxGC(pl_1, pl_2)
	pl_GCxGC = DataFrame(Name = pl_1.Name, tR1 = pl_1.tR, τR1 = pl_1.τR)
	tR2 = Array{Float64}(undef, length(pl_GCxGC.Name))
	τR2 = Array{Float64}(undef, length(pl_GCxGC.Name))
	CAS1 = [GasChromatographySimulator.CAS_identification(name).CAS for name in pl_1.Name] 
	CAS2 = [GasChromatographySimulator.CAS_identification(name).CAS for name in pl_2.Name]
	for i=1:length(pl_GCxGC.Name)
		ii = findfirst(CAS1[i].==CAS2)
		tR2[i] = pl_2.tR[ii]
		τR2[i] = pl_2.τR[ii]
	end
	pl_GCxGC[!, :tR2] = tR2
	pl_GCxGC[!, :τR2] = τR2
	return pl_GCxGC
end=#

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
		modules[i] = GasChromatographySystems.ModuleColumn("$(i) -> $(i+1)", Ls[i], ds[i]*1e-3, dfs[i]*1e-6, sps[i], TPs[i], F/60e6; kwargs...)
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

#example_SplitSystem() = SplitSystem([10.0, 1.0, 5.0], [0.25, 0.1, 0.25], [0.25, 0.0, 0.0], ["Rxi17SilMS", "", ""], [default_TP(), 300.0, 300.0], [1.0, NaN, NaN], NaN, 0.0, 101.3)
#=
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

function GCxGC_TM_simp(; L1 = 30.0, d1 = 0.25, df1 = 0.25, sp1 = "Rxi17SilMS", TP1 = default_TP(), L2 = 2.0, d2 = 0.25, df2 = 0.25, sp2 = "Wax", TP2 = default_TP(), F = NaN, pin = 200.0, pout = 0.0, opt=GasChromatographySystems.Options(ng=true))
	sys = GCxGC_TM_simp(L1, d1, df1, sp1, TP1, L2, d2, df2, sp2, TP2, F, pin, pout; opt=opt)
	return sys
end

#example_GCxGC_TM_simp() = GCxGC_TM_simp(30.0, 0.25, 0.25, "Rxi17SilMS", default_TP(), 2.0, 0.25, 0.25, "Wax", default_TP(), NaN, 200.0, 0.0)

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

function GCxGC_FM_simp(; L1 = 30.0, d1 = 0.25, df1 = 0.25, sp1 = "Rxi17SilMS", TP1 = default_TP(), F1 = 1.0, L2 = 2.0, d2 = 0.25, df2 = 0.25, sp2 = "Wax", TP2 = default_TP(), F2 = 2.0, pin = NaN, pout = 0.0)
	sys = GCxGC_FM_simp(L1, d1, df1, sp1, TP1, F1, L2, d2, df2, sp2, TP2, F2, pin, pout)
	return sys
end

#example_GCxGC_FM_simp() = GCxGC_FM_simp(30.0, 0.25, 0.25, "Rxi17SilMS", default_TP(), 1.0, 2.0, 0.25, 0.25, "Wax", default_TP(), 2.0, NaN, 0.0)
=#
# thermal modulation

# definitions for the periodic smoothed rectangle function
# attention to the used values

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
		return ifelse(tmod < tstart, Tcold, Thot)
	else # smoothed rectangle
		return smooth_rectangle.(tmod, tstart, width, Tcold, Thot; flank=flank) 
	end
end

# definition GCxGC system with thermal modulator
function GCxGC_TM(L1, d1, df1, sp1, TP1, L2, d2, df2, sp2, TP2, LTL, dTL, dfTL, spTL, TPTL, LM::Array{Float64,1}, dM, dfM, spM, shift, PM, ratioM, HotM, ColdM, TPM, F, pin, pout; opt=GasChromatographySystems.Options(), optTM=ModuleTMopt(), optCol=ModuleColumnOpt())

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
			com_timesteps = GasChromatographySimulator.common_time_steps(com_timesteps, TPs[i].timesteps)
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

function GCxGC_TM(; L1 = 30.0, d1 = 0.25, df1 = 0.25, sp1 = "ZB1ms", TP1 = default_TP(), L2 = 2.0, d2 = 0.1, df2 = 0.1, sp2 = "Stabilwax", TP2 = default_TP(), LTL = 0.25, dTL = 0.1, dfTL = 0.1, spTL = "Stabilwax", TPTL = 280.0, LM = [0.30, 0.01, 0.90, 0.01, 0.30], dM = 0.1, dfM = 0.1, spM = "Stabilwax", shift = 0.0, PM = 4.0, ratioM = 0.9125, HotM = 30.0, ColdM = -120.0, TPM = default_TP(), F = 0.8, pin = NaN, pout = 0.0, opt=GasChromatographySystems.Options(), optTM=ModuleTMopt(), optCol=ModuleColumnOpt())
	sys = GCxGC_TM(L1, d1, df1, sp1, TP1, L2, d2, df2, sp2, TP2, LTL, dTL, dfTL, spTL, TPTL, LM::Array{Float64,1}, dM, dfM, spM, shift, PM, ratioM, HotM, ColdM, TPM, F, pin, pout; opt=opt, optTM=optTM, optCol=optCol)
	return sys
end

#example_GCxGC_TM() = GCxGC_TM(30.0, 0.25, 0.25, "ZB1ms", default_TP(), 0.1, 0.1, 0.1, "Stabilwax", default_TP(), 0.56, 0.1, 0.1, "Stabilwax", 280.0, [0.30, 0.01, 0.90, 0.01, 0.30], 0.1, 0.1, "Stabilwax", 0.0, 0.0, 4.0, 0.9125, 80.0, -120.0, default_TP(), 0.8, NaN, 0.0)

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
	n_slice = Int.((init_t_end .- init_t_start)./PM .+ 1) # number of slices for every substance
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

## further functions for evaluation of GCxGC_TM results
function peaklist_GCxGC(pl_end, pl_1D, PM, shift)
	fit_D1 = fit_gauss_D1(pl_end, pl_1D, PM, shift)
	fit_D2 = fit_gauss_D2(pl_end, PM, shift)
	tR1 = Array{Float64}(undef, length(fit_D1.Name))
	tR2 = Array{Float64}(undef, length(fit_D1.Name))
	for i=1:length(fit_D1.Name)
		ii = findfirst(fit_D1.Name[i].==fit_D2.Name)
		tR1[i] = fit_D1.fits[i].param[1]
		tR2[i] = fit_D2.fits[ii].param[1]
	end
	return DataFrame(Name=fit_D1.Name, tR1=tR1, tR2=tR2)
end

function peaklist_GCxGC_weighted_mean(pl_end, PM, shift)
	heights = GasChromatographySystems.heights_of_peaks(pl_end)
	tR1 = fld.(pl_end.tR .+ shift, PM).*PM .- shift
	tR2 = pl_end.tR .- (fld.(pl_end.tR .+ shift, PM).*PM .- shift)
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

function fit_gauss_D1(pl_D2, pl_D1, PM, shift) # shift?
	heights = heights_of_peaks(pl_D2)
	tR = fld.(pl_D2.tR .+ shift, PM).*PM .- shift
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

function fit_gauss_D2(pl_D2, PM, shift) # shift ?
	heights = heights_of_peaks(pl_D2)
	tR = pl_D2.tR .- (fld.(pl_D2.tR .+ shift, PM).*PM .- shift)
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

function chrom_slicing(t, c, PM, shift) # correctly account for the shift!!!
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
end

function chrom2d(pl_final, sys)
	t_ = 0.0:0.01:sum(sys.modules[1].temperature.timesteps)
	chrom_sliced = Array{Array{Float64}}(undef, length(pl_final.tR))
	for i=1:length(pl_final.tR)
		chrom_sliced[i] = GasChromatographySimulator.chromatogram(collect(t_), [pl_final.tR[i]], [pl_final.τR[i]]).*pl_final.A[i]
	end
	chrom_sliced_sum = chrom_sliced[1]
	for i=2:length(chrom_sliced)
		chrom_sliced_sum = chrom_sliced_sum .+ chrom_sliced[i]
	end
	#Plots.plot(t_, chrom_sliced_sum)
	# determin the index of the ModuleTM
	c_slices, t_D1, t_D2 = chrom_slicing(t_, chrom_sliced_sum, sys.modules[5].PM, sys.modules[5].shift)
	slice_mat = Array{Float64}(undef, length(c_slices)-1, length(t_D2[1]))
	for j=1:length(t_D2[1])
		for i=1:(length(c_slices)-1)
			slice_mat[i,j] = c_slices[i][j]
		end
	end
	return slice_mat, t_D1, t_D2, c_slices, t_, chrom_sliced_sum, chrom_sliced 
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

end # module