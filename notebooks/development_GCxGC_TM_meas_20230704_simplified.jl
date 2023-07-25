### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 23c71f14-efdf-11ed-33f2-95304516114d
begin
	import Pkg
    # activate the shared project environment
	Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
	#Pkg.upgrade_manifest()
	#Pkg.precompile()
    Pkg.instantiate()

	using CSV, DataFrames
	using Plots, CairoMakie, GraphMakie
	using Graphs, NetworkLayout, Symbolics
	using GasChromatographySimulator
	using PlutoUI
	include(joinpath(dirname(pwd()), "src", "GasChromatographySystems.jl"))
	TableOfContents(depth=4)
end

# ╔═╡ d7f045e5-47c7-44b6-8bce-ce9dc3154fff
using Integrals

# ╔═╡ d8d6880e-13c3-4cc9-9dcc-cf616b3d81b1
using OrdinaryDiffEq

# ╔═╡ 08d9257e-6d49-4cf6-97a9-f8b546d9e933
using LsqFit

# ╔═╡ d419921b-1472-4b66-bae4-469537259814
md"""
# Development of a thermal modulated GCxGC system **simplified**

Simulation will be made for the `ModuleColumn` segments in the normal way. For `ModuleTM` no simulation will be made. The elution time from the modulator point will be at t₀ + tcold. The peak width will be approximated as (length modulator point)/(solute velocity(t₀ + tcold + ϵ)). The area of the not focussed part of the sliced peak at the first modulator point will be added to the previous focussed part, assuming this part would be focussed at the second modulator point. The time needed to migrate between the modulator points and a potential time shift between them is used to check, if a breaktrough is possible and a warning will be given. 
"""

# ╔═╡ dff75c18-eee5-4f70-9509-0e8228932819
md"""
## Periodic rectangle temperature function

for the simplified approace, the realization of the temperature modulation is not so much important.
"""

# ╔═╡ 326d2249-dfcc-40b6-b8c0-d9ab5ff88627
function therm_mod(t, shift, PM, ratio, Tcold, Thot) 
	return ifelse(mod(t+shift, PM) < ratio*PM, Tcold, Thot)
end

# ╔═╡ 1f5c1ca8-7b61-4dc4-91e6-9fb187073313
begin
	_t = 0.0:0.001:20.0
	_shift = 0.0
	_PM = 4.0
	_ratio = (_PM-0.35)/_PM
	_Thot = +80.0
	_Tcold = -80.0
	_T = Array{Float64}(undef, length(_t))
	for i=1:length(_t)
		_T[i] = therm_mod(_t[i], _shift, _PM, _ratio, _Tcold, _Thot)
	end
	Plots.plot(_t, _T)
end

# ╔═╡ 9a6dea14-f8d4-4b79-9c02-268c0c799515
# differentiation of this function is zero everywehre

# ╔═╡ 60d86bdd-5527-49da-8c35-ecd508171e6c
md"""
## Definition System
"""

# ╔═╡ b63ea500-c2bb-49c3-9915-032fd3e33e14
begin
	tsteps_, Tsteps_ = GasChromatographySimulator.conventional_program([50.0, 1.0, 3.0, 200.0, 0.0, 10.0, 225.0, 10.0]) 
	GCxGC_TP = GasChromatographySystems.TemperatureProgram(tsteps_, Tsteps_)
end

# ╔═╡ 522a7582-0ff1-42d6-960a-adb952a225d2
function update_system(sys)
	new_timesteps, new_pressuresteps, new_temperaturesteps, index_module_tempprog = GasChromatographySystems.match_programs(sys)
	new_pp = Array{GasChromatographySystems.PressurePoint}(undef, nv(sys.g))
	for i=1:nv(sys.g)
		new_pp[i] = GasChromatographySystems.PressurePoint(sys.pressurepoints[i].name, new_timesteps, new_pressuresteps[i])
	end
	new_modules = Array{GasChromatographySystems.AbstractModule}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		if typeof(sys.modules[i].temperature) <: Number
			new_modules[i] = sys.modules[i]
		elseif typeof(sys.modules[i].temperature) <: GasChromatographySystems.TemperatureProgram
			# add/modify for gradient
			ii = findfirst(index_module_tempprog.==i)
			new_tp = GasChromatographySystems.TemperatureProgram(new_timesteps, new_temperaturesteps[ii])
			if typeof(sys.modules[i]) == GasChromatographySystems.ModuleColumn
				new_modules[i] = GasChromatographySystems.ModuleColumn(sys.modules[i].name, sys.modules[i].length, sys.modules[i].diameter, sys.modules[i].film_thickness, sys.modules[i].stationary_phase, new_tp, sys.modules[i].flow)
			elseif typeof(sys.modules[i]) == GasChromatographySystems.ModuleTM
				new_modules[i] = GasChromatographySystems.ModuleTM(sys.modules[i].name, sys.modules[i].length, sys.modules[i].diameter, sys.modules[i].film_thickness, sys.modules[i].stationary_phase, new_tp, sys.modules[i].shift, sys.modules[i].PM, sys.modules[i].ratio, sys.modules[i].Thot, sys.modules[i].Tcold, sys.modules[i].flow)
			end
		end
	end
	new_sys = GasChromatographySystems.System(sys.g, new_pp, new_modules, sys.options)
    return new_sys
end

# ╔═╡ 80150c58-f9b9-43d9-bd7b-e38388f29fd8
function GCxGC_TM(L1, d1, df1, sp1, TP1, L2, d2, df2, sp2, TP2, LTL, dTL, dfTL, spTL, TPTL, LM::Array{Float64,1}, dM, dfM, spM, shiftM, PM, ratioM, HotM, ColdM, TPM, F, pin, pout; opt=GasChromatographySystems.Options(ng=true))

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
		if typeof(TPs[i]) <: GasChromatographySystems.TemperatureProgram
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
	pp = Array{GasChromatographySystems.PressurePoint}(undef, nv(g))
	pp[1] = GasChromatographySystems.PressurePoint("p₁", com_timesteps, pins) # inlet 
	for i=2:(nv(g)-1)
		pp[i] = GasChromatographySystems.PressurePoint("p$(i)", com_timesteps, nans)
	end
	pp[end] = GasChromatographySystems.PressurePoint("p$(nv(g))", com_timesteps, pouts) # outlet
	# modules
	modules = Array{GasChromatographySystems.AbstractModule}(undef, ne(g))
	modules[1] = GasChromatographySystems.ModuleColumn("GC column 1", L1, d1*1e-3, df1*1e-6, sp1, TP1, F/60e6)
	modules[2] = GasChromatographySystems.ModuleColumn("mod in", LM[1], dM*1e-3, dfM*1e-6, spM, TPM)
	modules[3] = GasChromatographySystems.ModuleTM("TM1", LM[2], dM*1e-3, dfM*1e-6, spM, TPM, shiftM, PM, ratioM, HotM, ColdM, NaN)
	modules[4] = GasChromatographySystems.ModuleColumn("mod loop", LM[3], dM*1e-3, dfM*1e-6, spM, TPM)
	modules[5] = GasChromatographySystems.ModuleTM("TM2", LM[4], dM*1e-3, dfM*1e-6, spM, TPM, shiftM, PM, ratioM, HotM, ColdM, NaN)
	#modules[5] = GasChromatographySystems.ModuleColumn("TM2c", LM[4], dM*1e-3, dfM*1e-6, spM, TPM, NaN) # simulate 2nd Modulator point with the oven temperature -> single stage modulator
	modules[6] = GasChromatographySystems.ModuleColumn("mod out", LM[5], dM*1e-3, dfM*1e-6, spM, TPM, NaN)
	modules[7] = GasChromatographySystems.ModuleColumn("GC column 2", L2, d2*1e-3, df2*1e-6, sp2, TP2, NaN)
	modules[8] = GasChromatographySystems.ModuleColumn("TL", LTL, dTL*1e-3, dfTL*1e-6, spTL, TPTL, NaN)
	# system
	sys = update_system(GasChromatographySystems.System(g, pp, modules, opt))
	return sys
end

# ╔═╡ 7bb430cd-6bc5-4c48-85f8-3ffce3220d50
begin
	PM = 4.0
	hotjet = 0.35
	coldjet = PM - hotjet
	
#	sys = GCxGC_TM(30.0, 0.25, 0.25, "ZB1ms", GCxGC_TP, 0.1, 0.1, 0.1, "Stabilwax", GCxGC_TP, 0.56, 0.1, 0.1, "Stabilwax", 280.0, [0.30, 0.005, 0.90, 0.005, 0.30], 0.1, 0.1, "Stabilwax", 0.0, PM, coldjet/PM, 80.0, -120.0, GCxGC_TP, 0.8, NaN, 0.0)
	# old:
#	sys = GCxGC_TM(31.5, 0.25, 0.25, "ZB1ms", GCxGC_TP, 0.1, 0.1, 0.1, "Stabilwax", GCxGC_TP, 0.25, 0.1, 0.1, "Stabilwax", 280.0, [0.60, 0.01, 0.30, 0.01, 0.50], 0.1, 0.1, "Stabilwax", 0.0, PM, coldjet/PM, 80.0, -120.0, GCxGC_TP, 0.65, NaN, 0.0)
	sys = GCxGC_TM(31.5, 0.25, 0.25, "ZB1ms", GCxGC_TP, 0.56, 0.1, 0.1, "Stabilwax", GCxGC_TP, 0.1, 0.1, 0.1, "Stabilwax", 280.0, [0.30, 0.01, 0.90, 0.01, 0.30], 0.1, 0.1, "Stabilwax", 0.0, PM, coldjet/PM, 25.0, -120.0, GCxGC_TP, 0.8, NaN, 0.0)
end

# ╔═╡ ca604ac6-fdb1-41d7-83d2-a0c9b2440dad
GasChromatographySystems.plot_graph(sys)

# ╔═╡ b9e3acb8-a7ec-40ef-ba58-932f79b08fb0
#=function flow_restriction(t, module_::GasChromatographySystems.ModuleColumn, options)
	if typeof(module_.temperature) <: GasChromatographySystems.TemperatureProgram
		T_itp = GasChromatographySimulator.temperature_interpolation(module_.temperature.timesteps, module_.temperature.temperaturesteps, module_.temperature.gradient_function, module_.length)
	elseif typeof(module_.temperature) <: Number
		gf(x) = [zero(x), zero(x)]
		T_itp = GasChromatographySimulator.temperature_interpolation([0.0, 3600.0], [module_.temperature, module_.temperature], gf, module_.length)
	end
	κ = GasChromatographySimulator.flow_restriction(module_.length, t, T_itp, module_.diameter, options.mobile_phase; ng=options.ng, vis=options.vis)
	return κ
end=#

# ╔═╡ 9605b6ad-7cdc-47ee-9edf-c3ce531ae9a2
#=function flow_restriction(t, module_::ModuleTM, options)
	T_itp(x, t) = modulator_temperature(x, t, module_)
	κ = GasChromatographySimulator.flow_restriction(module_.length, t, T_itp, module_.diameter, options.mobile_phase; ng=options.ng, vis=options.vis)
	return κ
end=#

# ╔═╡ 06487ba0-8752-40a8-b5d8-07fb08514ffe
#=function flow_restrictions(sys)
	kappas = Array{Function}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		f(t) = flow_restriction(t, sys.modules[i], sys.options)
		kappas[i] = f
	end
	return kappas
end=#

# ╔═╡ a5c2b339-4879-4367-9d03-ada1e7e6dd46
sys

# ╔═╡ bdbd8cc9-29fb-43ce-850b-8d18980808d3
md"""
## Simulation
"""

# ╔═╡ 4adbf7e9-232d-4889-92c6-75e8bde0d88d
begin
	db = DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/Estimator/GCsim_d_renamed.csv", header=1, silencewarnings=true))
	db.Tchar = db.Tchar .- 273.15
	db[!, :No] = collect(1:length(db.Name))
	db
end

# ╔═╡ fc36ffe3-b1df-42c8-ad14-e3b90f27a272
selected_solutes_ = GasChromatographySystems.common_solutes(db, sys).Name

# ╔═╡ 43f27a81-dca0-4188-baf4-3f8140d9e57c
selected_solutes = selected_solutes_[[14,16,17]]

# ╔═╡ 2ec6f339-c7c4-41c1-869e-8cd94d993952
#=function modulator_temperature(x, t, module_::ModuleTM)
	if typeof(module_.temperature) <: GasChromatographySystems.TemperatureProgram
		T_itp_ = GasChromatographySimulator.temperature_interpolation(module_.temperature.timesteps, module_.temperature.temperaturesteps, module_.temperature.gradient_function, module_.length)
	elseif typeof(module_.temperature) <: Number
		gf(x) = [zero(x), zero(x)]
		T_itp_ = GasChromatographySimulator.temperature_interpolation([0.0, 3600.0], [module_.temperature, module_.temperature], gf, module_.length)
	end
	T_itp(x, t) = therm_mod(t, module_.shift, module_.PM, module_.ratio, module_.Tcold, T_itp_(x, t) .+ module_.Thot .- 273.15) .+ 273.15 # cool jet always at Tcold
	#T_itp(x, t) = therm_mod(t, module_.shift, module_.PM, module_.ratio, T_itp_(x, t) .+ module_.Tcold .- 273.15, T_itp_(x, t) .+ module_.Thot .- 273.15) .+ 273.15 # cool jet always at T_itp_ + Tcold
	return T_itp(x, t)
end=#

# ╔═╡ ac60e9eb-0dfb-4a8d-a6bb-5289110b40cd
begin

	function module_temperature(module_::GasChromatographySystems.ModuleColumn, sys; Tcold_abs=true)
		if typeof(module_.temperature) <: GasChromatographySystems.TemperatureProgram # temperature is a TemperatureProgram
			time_steps = module_.temperature.timesteps
			temp_steps = module_.temperature.temperaturesteps	
			gf = module_.temperature.gradient_function
			a_gf = module_.temperature.a_gradient_function
			T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, module_.length)
		elseif typeof(module_.temperature) <: Number # temperature is a constant value
			time_steps = GasChromatographySystems.common_timesteps(sys)
			temp_steps = module_.temperature.*ones(length(time_steps))
			gf(x) = zero(x).*ones(length(time_steps))
			a_gf = [zeros(length(time_steps)) zeros(length(time_steps)) ones(length(time_steps)) zeros(length(time_steps))]
			T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, module_.length)
		end
		return time_steps, temp_steps, gf, a_gf, T_itp
	end

	function module_temperature(module_::GasChromatographySystems.ModuleTM, sys; Tcold_abs=true)
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
		T_itp(x,t) = if Tcold_abs == true # cool jet always at Tcold
			therm_mod(t, module_.shift, module_.PM, module_.ratio, module_.Tcold, T_itp_(x, t) .+ module_.Thot .- 273.15) .+ 273.15 
		else # cool jet always at T_itp_ + Tcold
			therm_mod(t, module_.shift, module_.PM, module_.ratio, T_itp_(x, t) .+ module_.Tcold .- 273.15, T_itp_(x, t) .+ module_.Thot .- 273.15) .+ 273.15 
		end
		return time_steps, temp_steps, gf, a_gf, T_itp
	end
end

# ╔═╡ 54c85e34-941d-41a1-9780-cf156f82d765
# not essential function
function flow_restrictions(sys)
	kappas = Array{Function}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		T_itp = module_temperature(sys.modules[i], sys)[5]
		κ(t) = GasChromatographySimulator.flow_restriction(sys.modules[i].length, t, T_itp, sys.modules[i].diameter, sys.options.mobile_phase; ng=sys.options.ng, vis=sys.options.vis)
		kappas[i] = κ
	end
	return kappas
end

# ╔═╡ c366df65-59a4-427a-beee-bf0bef5fb409
begin
	κ = flow_restrictions(sys)
	Plots.plot(0.0:0.01:100.0, κ[3].(0.0:0.01:100.0), title="flow restriction at modulator point", xlabel="time in s", ylabel="flow restriction in Pa s K m⁻³", label="")
end

# ╔═╡ 607265de-c784-4ed0-8e89-ebb206785b0b
begin
	TM_itp = module_temperature(sys.modules[3], sys)[5]
	Plots.plot(500.0:0.01:600.0, TM_itp.(0.0, 500.0:0.01:600.0).-273.115, title="temperature at modulator point", xlabel="time in s", ylabel="temperature in °C", label="")
end

# ╔═╡ 7277a646-fd70-4032-9691-569eddeafec6
function graph_to_parameters(sys, db_dataframe, selected_solutes; interp=true, dt=1, Tcold_abs=true)
	E = collect(edges(sys.g))
	srcE = src.(E) # source indices
	dstE = dst.(E) # destination indices
	if interp == true # linear interpolation of pressure functions with step width dt
		p_func = GasChromatographySystems.interpolate_pressure_functions(sys; dt=dt)
	else
		p_func = GasChromatographySystems.pressure_functions(sys)
	end
	parameters = Array{GasChromatographySimulator.Parameters}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		# column parameters
		col = GasChromatographySimulator.Column(sys.modules[i].length, sys.modules[i].diameter, [sys.modules[i].diameter], sys.modules[i].film_thickness, [sys.modules[i].film_thickness], sys.modules[i].stationary_phase, sys.options.mobile_phase)

		# program parameters
		pin_steps = sys.pressurepoints[srcE[i]].pressure_steps
		pout_steps = sys.pressurepoints[dstE[i]].pressure_steps
		pin_itp = p_func[srcE[i]]
		pout_itp = p_func[dstE[i]]	
#		if typeof(sys.modules[i].temperature) <: GasChromatographySystems.TemperatureProgram # temperature is a TemperatureProgram
#			time_steps = sys.modules[i].temperature.timesteps
#			temp_steps = sys.modules[i].temperature.temperaturesteps	
#			gf = sys.modules[i].temperature.gradient_function
#			a_gf = sys.modules[i].temperature.a_gradient_function
#			T_itp_ = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, sys.modules[i].length)
#		elseif typeof(sys.modules[i].temperature) <: Number # temperature is a constant value
#			time_steps = GasChromatographySystems.common_timesteps(sys)
#			temp_steps = sys.modules[i].temperature.*ones(length(time_steps))
#			gf(x) = zero(x).*ones(length(time_steps))
#			a_gf = [zeros(length(time_steps)) zeros(length(time_steps)) ones(length(time_steps)) zeros(length(time_steps))]
#			T_itp_ = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, sys.modules[i].length)
#		end
#		T_itp(x, t) = if typeof(sys.modules[i]) == GasChromatographySystems.ModuleColumn
#			T_itp_(x,t) 
#		elseif typeof(sys.modules[i]) == ModuleTM
#			# cool jet always at Tcold
#			therm_mod(t, sys.modules[i].shift, sys.modules[i].PM, sys.modules[i].ratio, sys.modules[i].Tcold, T_itp_(x, t) .+ sys.modules[i].Thot .- 273.15) .+ 273.15 
#			# cool jet always at T_itp_ + Tcold
#			#therm_mod(t, sys.modules[i].shift, sys.modules[i].PM, sys.modules[i].ratio, T_itp_(x, t) .+ sys.modules[i].Tcold .- 273.15, T_itp_(x, t) .+ sys.modules[i].Thot .- 273.15) .+ 273.15 
#			##modulator_temperature(x, t, sys.modules[i])
#		end
		time_steps, temp_steps, gf, a_gf, T_itp = module_temperature(sys.modules[i], sys; Tcold_abs=Tcold_abs)
		prog = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)

		# substance parameters
		sub = GasChromatographySimulator.load_solute_database(db_dataframe, sys.modules[i].stationary_phase, sys.options.mobile_phase, selected_solutes, NaN.*ones(length(selected_solutes)), NaN.*ones(length(selected_solutes)))

		# option parameters
		opt = GasChromatographySimulator.Options(sys.options.alg, sys.options.abstol, sys.options.reltol, sys.options.Tcontrol, sys.options.odesys, sys.options.ng, sys.options.vis, sys.options.control, sys.options.k_th)

		parameters[i] = GasChromatographySimulator.Parameters(col, prog, sub, opt)
	end
	return parameters
end

# ╔═╡ 793560f7-d364-4c68-81ec-994441a41059
par = graph_to_parameters(sys, db, selected_solutes; Tcold_abs=true)

# ╔═╡ 6993afed-4519-4c8a-9cfc-87c72723c444
begin
	gr()#plotly()
	p_TM_Prog = Plots.plot(500.0:0.01:800.0, par[3].prog.T_itp.(0.0, 500.0:0.01:800.0).-273.15)
	Plots.plot!(p_TM_Prog, 500.0:0.01:800.0, par[4].prog.T_itp.(0.0, 500.0:0.01:800.0).-273.15)
	#Plots.plot!(p_TM_Prog, 500.0:0.01:600.0, par[5].prog.T_itp.(0.0, 500.0:0.01:600.0).-273.15)
	Plots.plot!(p_TM_Prog, xlabel="time in s", ylabel="temperature in °C", legend=false, title="temperature at modulator point")
end

# ╔═╡ a14d1ca2-46ef-4296-aa1d-06590eca97cf
paths = GasChromatographySystems.all_paths(sys.g, 1)[2]

# ╔═╡ 69598506-cf92-4626-b6e1-d8f57bf8388f
md"""
### Simulation with the old function

Only results up to the modulator point are assumed to be correct. These correct results are used as starting point for a step-by-step simulation through the different modules.
"""

# ╔═╡ 1a0695c5-0a67-4158-a3bd-c135be31e408
sim = GasChromatographySystems.simulate_along_paths(sys, [paths[1][1:2]], par)

# ╔═╡ 6556a858-65ef-49a6-bdaa-562a2b568a70
md"""
### (deactivated) Simulation using the new function

The new function will reproduce the step-by-step simulation.
"""

# ╔═╡ 67ba928d-5af0-422e-bd58-c6302c3e41af
[paths[1][1:1]]

# ╔═╡ 3e8d2e43-1fa9-40e4-a132-7228664b7546
sys.modules[3]

# ╔═╡ e3e34deb-9249-4811-b74a-44a4c2d06ac2
md"""
## Modulation step-by-step

Using the results of `sim` (Simulation with the old function) just before the 1st modulator spot
"""

# ╔═╡ 54e20666-b217-415d-b2fd-305980588668
md"""
### Peak before 1st Modulator spot

Results after segment 2.
"""

# ╔═╡ ac3db64b-82c8-455a-97aa-6bbddcc03830
sim[2][1][2]

# ╔═╡ 7e088401-a8d4-4648-9380-241777534fe2
#savefig(p_chrom_before, "peak_before_modulator_new.svg")

# ╔═╡ 82c69944-c4e0-4699-87c0-29500b33ba78
# not essential function
function chrom_before_modulation(pl; nτ=6)
	tstart = Array{Float64}(undef, length(pl.Name))
	tend = Array{Float64}(undef, length(pl.Name))
	t = Array{Array{Float64}}(undef, length(pl.Name))
	c = Array{Array{Float64}}(undef, length(pl.Name))
	for i=1:length(pl.Name)
		tstart[i] = pl.tR[i] - nτ * pl.τR[i]
		tend[i] = pl.tR[i] + nτ * pl.τR[i]
		t[i] = collect(tstart[i]:(tend[i]-tstart[i])/100:tend[i])
		c[i] = GasChromatographySimulator.chromatogram(t[i], [pl.tR[i]], [pl.τR[i]])
	end
	t1 = minimum(tstart)
	t2 = maximum(tend)
	t_sum = collect(t1:(t2-t1)/1000:t2)
	c_sum = GasChromatographySimulator.chromatogram(t_sum, pl.tR, pl.τR)
	return t, c, t_sum, c_sum
end

# ╔═╡ ee2e3313-9d05-4815-97a2-5b292dc47b68
begin
	t_befores, c_befores, t_before, c_before = chrom_before_modulation(sim[2][1][2])
	p_chrom_before = Plots.plot(t_before, c_before, xlabel="time in s", label="")
	for i=1:length(c_befores)
		Plots.plot!(p_chrom_before, t_befores[i], c_befores[i], label=sim[2][1][2].Name[i])
	end
	Plots.plot!(p_chrom_before, xticks=0:4:3000)
	p_chrom_before
end

# ╔═╡ ce958d12-3dad-4e3c-ab57-f84136610870
#=begin

	function change_initial(par::GasChromatographySimulator.Parameters, init_t, init_τ)
		# copys the parameters `par` and changes the values of par.sub[i].τ₀ and par.sub[i].t₀ to init_τ[i] resp. init_t[i]
		newsub = Array{GasChromatographySimulator.Substance}(undef, length(par.sub))
		for i=1:length(par.sub)
			newsub[i] = GasChromatographySimulator.Substance(par.sub[i].name, par.sub[i].CAS, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀, par.sub[i].ann, par.sub[i].Cag, init_t[i], init_τ[i])
		end
		newpar = GasChromatographySimulator.Parameters(par.col, par.prog, newsub, par.opt)
		return newpar
	end
	
	function change_initial(par::GasChromatographySimulator.Parameters, prev_par::GasChromatographySimulator.Parameters, init_t, init_τ)
		# copys the parameters `par` and changes the values of par.sub[i].τ₀ and par.sub[i].t₀ to init_τ[i] resp. init_t[i]
		newsub = Array{GasChromatographySimulator.Substance}(undef, length(prev_par.sub))
		for i=1:length(prev_par.sub)
			newsub[i] = GasChromatographySimulator.Substance(prev_par.sub[i].name, prev_par.sub[i].CAS, prev_par.sub[i].Tchar, prev_par.sub[i].θchar, prev_par.sub[i].ΔCp, prev_par.sub[i].φ₀, prev_par.sub[i].ann, prev_par.sub[i].Cag, init_t[i], init_τ[i])
		end
		newpar = GasChromatographySimulator.Parameters(par.col, par.prog, newsub, par.opt)
		return newpar
	end
	
	function change_initial(par::GasChromatographySimulator.Parameters, prev_par::GasChromatographySimulator.Parameters, pl)
		# copys the parameters `par` and changes the values of par.sub[i].τ₀ and par.sub[i].t₀ to pl.tR[], pl.τR[]
		#ann_pl = GasChromatographySimulator.CAS_identification(pl.Name).CAS
		newsub = Array{GasChromatographySimulator.Substance}(undef, length(prev_par.sub))
		for i=1:length(prev_par.sub)
			ii = intersect(findall(prev_par.sub[i].ann.==pl.ann), findall(prev_par.sub[i].CAS.==pl.CAS))[1] # annotations (slice number) and CAS have to be the same
			newsub[i] = GasChromatographySimulator.Substance(prev_par.sub[i].name, prev_par.sub[i].CAS, prev_par.sub[i].Tchar, prev_par.sub[i].θchar, prev_par.sub[i].ΔCp, prev_par.sub[i].φ₀, prev_par.sub[i].ann, prev_par.sub[i].Cag, pl.tR[ii], pl.τR[ii])
		end
		newpar = GasChromatographySimulator.Parameters(par.col, par.prog, newsub, par.opt)
		return newpar
	end
	
	function change_initial_focussed(par::GasChromatographySimulator.Parameters, prev_par::GasChromatographySimulator.Parameters, pl; τ₀=zeros(length(pl.tR)))
		# copys the parameters `par` and changes the values of par.sub[i].t₀ to pl.tR[]
		#CAS_pl = GasChromatographySimulator.CAS_identification(pl.Name).CAS
		newsub = Array{GasChromatographySimulator.Substance}(undef, length(prev_par.sub))
		for i=1:length(prev_par.sub)
			ii = intersect(findall(prev_par.sub[i].ann.==pl.ann), findall(prev_par.sub[i].CAS.==pl.CAS))[1] # annotations (slice number) and CAS have to be the same
			newsub[i] = GasChromatographySimulator.Substance(prev_par.sub[i].name, prev_par.sub[i].CAS, prev_par.sub[i].Tchar, prev_par.sub[i].θchar, prev_par.sub[i].ΔCp, prev_par.sub[i].φ₀, prev_par.sub[i].ann, prev_par.sub[i].Cag, pl.tR[ii], τ₀[ii])
		end
		newpar = GasChromatographySimulator.Parameters(par.col, par.prog, newsub, par.opt)
		return newpar
	end
end=#

# ╔═╡ 70597ee8-760d-4c21-af6c-981922208d10
md"""
### 1st Modulator spot
"""

# ╔═╡ b83330dd-4df3-433a-ab29-273e39f8b32c
begin
	# increase the peak width befor the modulator by 50%
	#sim[2][1][2][!,:τR] = 1.5.*sim[2][1][2][!,:τR]
	sim[2][1][2]
end

# ╔═╡ fc7f8b92-f0ee-4a89-ab7c-0c8f5ba1bb0f
begin
	p_chrom_before_mark = Plots.plot(t_before, c_before, xlabel="time in s", label="")
	for i=1:length(c_befores)
		Plots.plot!(p_chrom_before_mark, t_befores[i], c_befores[i], label=sim[2][1][2].Name[i])
	end
	Plots.plot!(p_chrom_before_mark, xticks=0:4:3000)
	#for j=1:length(c_befores)
	#	for i=1:n_slice[j]
	#		Plots.plot!(p_chrom_before_mark, [init_t_start[j]+(i-1)*sys.modules[3].PM, init_t_start[j]+(i-1)*sys.modules[3].PM], [0.0, 0.4], c=j+1, label="")
	#		Plots.plot!(p_chrom_before_mark, [init_t_start[j]+(i-1)*sys.modules[3].PM, init_t_start[j]+(i-1)*sys.modules[3].PM].+sys.modules[3].PM.*sys.modules[3].ratio, [0.0, 0.4], c=:gray, linestyle=:dash, label="")
			#Plots.plot!(p_chrom_before_mark, [1892, 1892], [0.0, 0.4], c=j+1)
	#	end
	#end

	#add marker for PM and tcold
	n_ = unique(fld.(t_before, sys.modules[3].PM))
	for i=1:length(n_)
		Plots.plot!(p_chrom_before_mark, [n_[i]*sys.modules[3].PM, n_[i]*sys.modules[3].PM], [0.0, 0.4], c=:black, label="")
		Plots.plot!(p_chrom_before_mark, [n_[i]*sys.modules[3].PM+sys.modules[3].PM*sys.modules[3].ratio, n_[i]*sys.modules[3].PM+sys.modules[3].PM*sys.modules[3].ratio], [0.0, 0.4], c=:black, linestyle=:dash, label="")
	end
	p_chrom_before_mark
end

# ╔═╡ 00886bb1-39ff-48f2-9329-2a195f705a30
md"""
#### Sliced focussed peaks

- use the `slicing()` function to get the approximated areas 
- set elution time from the modulator point as ``t_0 + t_{cold}``
- set peak width as ``τ ≈ L_M / u_R``
"""

# ╔═╡ 49f4e7f4-2520-4a97-93c8-495de159e7c4
par_slice_TM1, df_A_TM1 = GasChromatographySystems.slicing(sim[2][1][2], sys.modules[3].PM, sys.modules[3].ratio, sys.modules[3].shift, par[3])

# ╔═╡ 76ca88a4-24a1-476e-8316-916593552523
function approximate_modulator(par, df_A, PM, ratio, shift)
	# rectangular function is assumed
	tcold = PM*ratio
	thot = PM*(1-ratio)
	
	sort_df_A = sort(df_A, :t0)

	# start time as multiple of modulation period
	t₀ = fld.(sort_df_A.t0, PM).*PM .+ shift # shift?
	
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

# ╔═╡ 0b2c8f8b-4648-4af3-a45d-900c5c382860
function simulate_along_paths(sys, paths, par_sys; t₀=zeros(length(par_sys[1].sub)), τ₀=zeros(length(par_sys[1].sub)), nτ=6, abstolTM=1e-10, reltolTM=1e-8, algTM=Vern9(), refocus=falses(ne(sys.g)), τ₀_focus=zeros(length(par_sys[1].sub)))
	
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
						if refocus[i_par[j]] == true
							τ₀=τ₀_focus
						else
							τ₀=peaklists_[j-1].τR # PM?
						end
						
						new_par_sys[i_par[j]], df_A = GasChromatographySystems.slicing(peaklists_[j-1], sys.modules[i_par[j]].PM, sys.modules[i_par[j]].ratio, sys.modules[i_par[j]].shift, par_sys[i_par[j]]; nτ=nτ, τ₀=τ₀, abstol=abstolTM, reltol=reltolTM, alg=algTM)
						if algTM == "simplifiedTM"
							peaklists_[j], solutions_[j] = approximate_modulator(new_par_sys[i_par[j]], df_A, sys.modules[i_par[j]].PM, sys.modules[i_par[j]].ratio, sys.modules[i_par[j]].shift)
						else
							peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(new_par_sys[i_par[j]])
							GasChromatographySystems.add_A_to_pl!(peaklists_[j], df_A)
						end
						if maximum(peaklists_[j].τR) > sys.modules[i_par[j]].PM
							return @warn "Peak width of focussed peaks > modulation period. Simulation is aborted."
							end
					else # ModuleColumn
						new_par_sys[i_par[j]] = GasChromatographySystems.change_initial(par_sys[i_par[j]], peaklists_[j-1])
						#else
						#	new_par_sys[i_par[j]] = GasChromatographySystems.change_initial(par_sys[i_par[j]], peaklists_[j-1])
						#end
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

# ╔═╡ fcdc5015-31d0-4808-8218-f33901437c0b
sim_ = simulate_along_paths(sys, paths, par, algTM="simplifiedTM")

# ╔═╡ 1f62b350-332e-4acf-b907-75362157e4da
sim_[2][1][3]

# ╔═╡ f1bf21dd-07ee-4e4e-beb5-7d9225b3fc6b
sim_[2][1][5]

# ╔═╡ ceb6ba3f-aece-4ee3-b757-5904fa04c7ac
sim_[2][1][8]

# ╔═╡ 88fb12fa-3a88-4577-8709-e9f8ec130b2b
sim__ = simulate_along_paths(sys, [paths[1][1:3]], par)

# ╔═╡ 4c385876-bb3b-48fb-aebe-c5554d21df29
sol_TM1 = approximate_modulator(par_slice_TM1, df_A_TM1, sys.modules[3].PM, sys.modules[3].ratio, sys.modules[3].shift)

# ╔═╡ 3c22e8aa-376d-4eec-9650-3bb102d4a0b8
# identfy the common index of a substances CAS number and Annotations in a peaklist 
function common_index(pl, CAS, Annotation)
	ii_CAS = findall(CAS.==pl.CAS)
	ii_ann = findall(occursin.(Annotation, pl.Annotations))
	ii = intersect(ii_CAS, ii_ann)[1]
	return ii
end

# ╔═╡ 846c3cab-90d4-4e73-9ce5-c1b41f937cdf
function add_A_to_pl!(pl, df_A)
	# add the areas in the dataframe df_A to the peaklist, same CAS and Annotations are used to identify the correct solute
	sort_A = Array{Float64}(undef, length(df_A.A))
	for i=1:length(df_A.A)
		ii = common_index(pl, df_A.CAS[i], join(split(df_A.Annotations[i], ", ")[1:end-1], ", "))
		#ii = common_index(pl, df_A.CAS[i], split(df_A.Annotations[i], ", ")[1])
		sort_A[ii] = df_A.A[i]
	end
	pl[!, :A] = sort_A
	return pl
end

# ╔═╡ 06c3c74b-4ced-405e-89be-a05092fab027
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

# ╔═╡ 67a4dff7-2634-4ae2-adb0-515807988b17
md"""
#### Chromatogram after 1st TM
"""

# ╔═╡ 84111543-21e2-4891-9050-0c102cd46a52
#=begin
	chrom_sliced_1st = Array{Array{Float64}}(undef, length(sol_1st_TM[1].tR))
	height_1st = Array{Float64}(undef, length(sol_1st_TM[1].tR))
 	for i=1:length(sol_1st_TM[1].tR)
		chrom_sliced_1st[i] = GasChromatographySimulator.chromatogram(collect(t_), [sol_1st_TM[1].tR[i]], [newpar.col.L ./ sol_1st_TM[1].uR[i]]).*sol_1st_TM[1].A[i]
		height_1st[i] = (GasChromatographySimulator.chromatogram([sol_1st_TM[1].tR[i]], [sol_1st_TM[1].tR[i]], [newpar.col.L ./ sol_1st_TM[1].uR[i]])[1].*sol_1st_TM[1].A[i])
	end
	#Plots.plot(t_, chrom_sliced)
	height_1st
end=#

# ╔═╡ 9c6573e8-bc0f-4c40-9eea-9733726fbd76
function chrom_after_modulation(pl, par; nτ=6, width_sim=false)
	tstart = Array{Float64}(undef, length(pl.Name))
	tend = Array{Float64}(undef, length(pl.Name))
	t = Array{Array{Float64}}(undef, length(pl.Name))
	c = Array{Array{Float64}}(undef, length(pl.Name))
	τR = Array{Float64}(undef, length(pl.Name))
	for i=1:length(pl.Name)
		# Ableitung der Rechteckfuntion (Temperatur modulation) ist immer NULL -> der schlagartige Temperaturanstieg wird bei der Peakbreitensimulation nicht berücksichtigt, um dies zu berücksichtigen wird für width_sim=false die Peakbreite abgeschätz von der Länge des Modulator Spots und der Geschwindigkeit der Substanz bei Elution (im Hot Jet)
		τR[i] = if width_sim == true
			pl.τR[i]
		else
			par.col.L/pl.uR[i] 
		end
		tstart[i] = pl.tR[i] - nτ * τR[i]
		tend[i] = pl.tR[i] + nτ * τR[i]	
	end
	t1 = minimum(tstart[isnan.(tstart).==false])
	t2 = maximum(tend[isnan.(tstart).==false])
	t_sum = collect(t1:(t2-t1)/1000:t2)
	c_sum = fill(0.0, 1001)
	for i=1:length(pl.Name)
		if isnan(tstart[i])
			c[i] = fill(NaN,1001)
		else
			c[i] = GasChromatographySimulator.chromatogram(t_sum, [pl.tR[i]], [τR[i]])*pl.A[i]
			c_sum = c_sum .+ c[i]
		end
	end
	return t_sum, c, c_sum
end

# ╔═╡ 095ee9ba-d974-4625-a47c-3634d91c6a5d
begin
	gr()
	# focussed peaks
	t_TM1, c_TM1s, c_TM1 = chrom_after_modulation(sol_TM1[1], par_slice_TM1, width_sim=false)
	p_chrom_TM1 = Plots.plot(t_TM1, c_TM1, xlabel="time in s")
	
	Plots.plot!(p_chrom_TM1, xticks=0:4:3000, legend=false)
	Plots.plot!(p_chrom_TM1)#, xlims=(1899.5, 1900.2), ylims=(-0.01, 0.5))
	p_chrom_TM1
end

# ╔═╡ f6093dce-30a7-4628-8362-0c5ed8c55ebc
md"""
### Loop between Modulator spots
"""

# ╔═╡ beb15214-b1bb-4b39-a2aa-86f0727bc2c4
md"""
#### Sliced focused peaks
"""

# ╔═╡ 2c1fc71c-5f5d-4e18-84f5-42980d461753
newpar_loop = GasChromatographySystems.change_initial(par[4], sol_TM1[1])

# ╔═╡ 74121702-9f4a-410b-961b-7865ab3af941
begin
	sol_loop = GasChromatographySimulator.simulate(newpar_loop)
	add_A_to_pl!(sol_loop[1], sol_TM1[1])
end

# ╔═╡ 5d0cb969-8a98-4d28-aa5c-233723755f41
function chrom_after_loop(pl; nτ=6)
	tstart = Array{Float64}(undef, length(pl.Name))
	tend = Array{Float64}(undef, length(pl.Name))
	t = Array{Array{Float64}}(undef, length(pl.Name))
	c = Array{Array{Float64}}(undef, length(pl.Name))
	for i=1:length(pl.Name)
		tstart[i] = pl.tR[i] - nτ * pl.τR[i]
		tend[i] = pl.tR[i] + nτ * pl.τR[i]	
	end
	t1 = minimum(tstart[isnan.(tstart).==false])
	t2 = maximum(tend[isnan.(tstart).==false])
	t_sum = collect(t1:(t2-t1)/10000:t2)
	c_sum = fill(0.0, 10001)
	for i=1:length(pl.Name)
		if isnan(tstart[i])
			c[i] = fill(NaN,10001)
		else
			c[i] = GasChromatographySimulator.chromatogram(t_sum, [pl.tR[i]], [pl.τR[i]])*pl.A[i]
			c_sum = c_sum .+ c[i]
		end
	end
	return t_sum, c, c_sum
end

# ╔═╡ 96200092-8afd-4ce3-8287-79806a331656
# time in loop
ΔtR_loop = sol_loop[1].tR .- sol_TM1[1].tR

# ╔═╡ 66eaa639-416c-4167-9e00-14dcbf583a4c
md"""
#### Chromatogram after loop
"""

# ╔═╡ bb483a0e-72d7-405e-bf47-c5fa878e0994
begin
	plotly()
	
	# focussed peaks
	t_loop, c_loops, c_loop = chrom_after_loop(sol_loop[1])

	
	p_chrom_loop = Plots.plot(t_loop, c_loop, xlabel="time in s")
	for i=1:length(c_loops)
		if sol_loop[1].Name[i] == selected_solutes[1]
			color = :blue
		elseif sol_loop[1].Name[i] == selected_solutes[2]
			color = :red
		else
			color = :orange
		end
		Plots.plot!(p_chrom_loop, t_loop, c_loops[i], label=sol_loop[1].Name[i], c=color)
	end
	Plots.plot!(p_chrom_loop, xticks=0:4:3000, legend=false)

	p_chrom_loop
end

# ╔═╡ 53aff8ef-4585-4d03-921b-7f19827dc04e
#savefig(p_chrom_loop, "Chrom_after_loop_new.svg")

# ╔═╡ b31fa530-922f-4862-b1cf-04d2ec34b1c1
md"""
Slight separation of the different solute on the loop.
"""

# ╔═╡ 86274690-dfe8-4f71-8084-66fe5f03f932
md"""
### 2nd Modulator spot
"""

# ╔═╡ 6ace764f-1691-4f9a-a2bb-54635468ce8c
begin	
	p_chrom_loop_mark = Plots.plot(t_loop, c_loop, xlabel="time in s")
	for i=1:length(c_loops)
		if sol_loop[1].Name[i] == selected_solutes[1]
			color = :blue
		elseif sol_loop[1].Name[i] == selected_solutes[2]
			color = :red
		else
			color = :orange
		end
		Plots.plot!(p_chrom_loop_mark, t_loop, c_loops[i], label=sol_loop[1].Name[i], c=color)
	end
	Plots.plot!(p_chrom_loop_mark, xticks=0:4:3000, legend=false)
	
	
	#add marker for PM and tcold
	n = unique(fld.(t_loop, sys.modules[5].PM))
	for i=1:length(n)
		Plots.plot!(p_chrom_loop_mark, [n[i]*sys.modules[5].PM, n[i]*sys.modules[5].PM], [0.0, 10.0], c=:black)
		Plots.plot!(p_chrom_loop_mark, [n[i]*sys.modules[5].PM+sys.modules[5].PM*sys.modules[5].ratio, n[i]*sys.modules[5].PM+sys.modules[5].PM*sys.modules[5].ratio], [0.0, 10.0], c=:black, linestyle=:dash)
	end
	p_chrom_loop_mark
end

# ╔═╡ 2d71af52-1457-4bd3-a73e-424944ab34b9
par_slice_TM2, df_A_TM2 = GasChromatographySystems.slicing(sol_loop[1], sys.modules[5].PM, sys.modules[5].ratio, sys.modules[5].shift, par[5])

# ╔═╡ d0b3d8b5-3012-43e4-9c09-bdaa97d0e3f6
sol_TM2 = approximate_modulator(par_slice_TM2, df_A_TM2, sys.modules[5].PM, sys.modules[5].ratio, sys.modules[5].shift)

# ╔═╡ db04ef60-3490-49c3-a5c1-f198c4015b06
md"""
#### Chromatogram after 2nd TM
"""

# ╔═╡ ebc1b905-eef4-4c66-abb6-b070d19f17d9
p_chrom_TM1

# ╔═╡ 672e5876-0f0d-483a-9936-604b3fc7f97a
begin
	plotly()
	# focussed peaks
	t_TM2, c_TM2s, c_TM2 = chrom_after_modulation(sol_TM2[1], par_slice_TM2, width_sim=false)
	p_chrom_TM2 = Plots.plot(t_TM2, c_TM2, xlabel="time in s")
	# unfocussed peaks 
	#t_unfoc, c_unfoc, h_unfoc, t_unfoc_sum, c_unfoc_sum = rec_chrom_after_modulation(sol_1st_unfoc_TM[1])
	#for i=1:length(c_TM1s)
	#	Plots.plot!(p_chrom_TM1, t_TM1, c_TM1s[i], label=sol_1st_TM[1].Name[i])
	#end
	Plots.plot!(p_chrom_TM2, xticks=0:4:3000, legend=false)
	#Plots.plot!(p_chrom_TM1, t_unfoc_sum, c_unfoc_sum)
	Plots.plot!(p_chrom_TM2)#, xlims=(1899.5, 1900.2), ylims=(-0.01, 0.5))
	p_chrom_TM2
end

# ╔═╡ 231edaa6-8ef3-4a85-b48e-e0e1829d124d
#savefig(p_chrom_TM2, "Chrom_after_2nd_TM_new.svg")

# ╔═╡ 64e635f8-f5fd-4935-b44e-b02ebdd96b64
md"""
### Connection to GC2
"""

# ╔═╡ 78dbaf44-0468-417a-8c69-90c139e2e2ad
newpar_mod_out = GasChromatographySystems.change_initial(par[6], sol_TM2[1])

# ╔═╡ bb01e67d-fbca-461e-b592-508bc5480b94
begin
	sol_mod_out = GasChromatographySimulator.simulate(newpar_mod_out)
	add_A_to_pl!(sol_mod_out[1], sol_TM2[1])
end

# ╔═╡ 195f8556-e4e8-4b85-86d1-f70b5a6ce0bf
md"""
#### Chromatogram after Connection to GC2
"""

# ╔═╡ 088c0c0f-0722-4498-8306-8286069fc24b
begin
	plotly()
	t_mod_out, c_mod_outs, c_mod_out = chrom_after_loop(sol_mod_out[1])
	p_chrom_mod_out = Plots.plot(t_mod_out, c_mod_out, xlabel="time in s")
	for i=1:length(c_mod_outs)
		if sol_mod_out[1].Name[i] == selected_solutes[1]
			color = :blue
		elseif sol_mod_out[1].Name[i] == selected_solutes[2]
			color = :red
		else
			color = :orange
		end
		Plots.plot!(p_chrom_mod_out, t_mod_out, c_mod_outs[i], label=sol_mod_out[1].Name[i], c=color)
	end
	Plots.plot!(p_chrom_mod_out, xticks=0:4:3000, legend=false)
	p_chrom_mod_out
end

# ╔═╡ 29cbe2b2-34fb-45dc-849b-2f678b3f8797
#savefig(p_chrom_mod_out, "Chrom_after_mod_out_new.svg")

# ╔═╡ 14bf7364-ce16-467e-894b-43d5ad1d91dc
md"""
### GC2
"""

# ╔═╡ 98713f9d-123a-4eda-8a8c-72c6ada1a3a4
newpar_gc2 = GasChromatographySystems.change_initial(par[7], sol_mod_out[1])

# ╔═╡ caef921f-f06c-471b-a56a-8415ac2d4fa6
begin
	sol_gc2 = GasChromatographySimulator.simulate(newpar_gc2)
	add_A_to_pl!(sol_gc2[1], sol_mod_out[1])
end 

# ╔═╡ 05a0e4f4-8a1b-4f3c-84b0-c2491840d20c
# time in gc2
sol_gc2[1].tR .- sol_mod_out[1].tR

# ╔═╡ 2465e19f-36f8-4a7d-93ed-d6b872532c40
md"""
#### Chromatogram after GC2
"""

# ╔═╡ 1bc4ae12-6f10-499a-b20f-35d5bd297378
begin
	gr()
	t_gc2, c_gc2s, c_gc2 = chrom_after_loop(sol_gc2[1])
	p_chrom_gc2 = Plots.plot(t_gc2, c_gc2, xlabel="time in s")
	for i=1:length(c_gc2s)
		if sol_gc2[1].Name[i] == selected_solutes[1]
			color = :blue
		elseif sol_gc2[1].Name[i] == selected_solutes[2]
			color = :red
		else
			color = :orange
		end
		Plots.plot!(p_chrom_gc2, t_gc2, c_gc2s[i], label=sol_gc2[1].Name[i], c=color)
	end
	Plots.plot!(p_chrom_gc2, xticks=0:4:3000, legend=false)
	p_chrom_gc2
end

# ╔═╡ 935485a4-947e-42ad-ae92-9146c3f551ff
#savefig(p_chrom_gc2, "Chrom_after_gc2_new.svg")

# ╔═╡ 16e6bf4c-05ee-466d-95d6-e679e4410bc1
md"""
### TL
"""

# ╔═╡ 6accd799-2058-4756-9e31-c0b8b6bc6130
newpar_tl = GasChromatographySystems.change_initial(par[8], sol_gc2[1])

# ╔═╡ 32e40a19-ebab-47a6-b4ca-9300c2b04dab
begin
	sol_tl = GasChromatographySimulator.simulate(newpar_tl)
	add_A_to_pl!(sol_tl[1], sol_gc2[1])
end 

# ╔═╡ 93ee192a-7374-43f4-a1a1-a1265fda42a4
# time in tl
sol_tl[1].tR .- sol_gc2[1].tR

# ╔═╡ 07922786-3664-4308-a939-6d9cd9cffb33
md"""
#### Chromatogram after TL - Detector
"""

# ╔═╡ 0205e873-f3b5-4100-a49b-5824fc5eecf3
sol_tl[1].tR[10].-fld.(sol_tl[1].tR[10], PM).*PM

# ╔═╡ ec41da66-1aba-4e38-8c69-1cc236157279
hh = (GasChromatographySimulator.chromatogram([sol_tl[1].tR[10]], [sol_tl[1].tR[10]], [sol_tl[1].τR[10]])[1].*sol_tl[1].A[10])

# ╔═╡ 8b500dc4-98de-46af-9a83-f5ff7e8190ee
sim[2][1][2].Name, sim[2][1][2].tR, sim[2][1][2].τR

# ╔═╡ 9d6944b7-7794-4912-b8ab-ce91672cf0de
function dot_color(name, selected_solutes)
		if name == selected_solutes[1]
			color = :blue
		elseif name == selected_solutes[2]
			color = :red
		else
			color = :orange
		end
	return color
end

# ╔═╡ a0b9918a-1ea2-463a-8587-2b4ca0f1230c
#savefig(p_chrom_tl, "Chrom_after_tl_fit.svg")

# ╔═╡ 6c4c9478-1392-442c-9a26-b654b64c1b14
function heights_of_peaks(pl)
	heights = Array{Float64}(undef, length(pl.tR))
 	for i=1:length(pl.tR)
		heights[i] = (GasChromatographySimulator.chromatogram([pl.tR[i]], [pl.tR[i]], [pl.τR[i]])[1].*pl.A[i])
	end
	#Plots.plot(t_, chrom_sliced)
	heights
end

# ╔═╡ 50792b68-8acf-4349-afd5-66267fa9153b
heights_tl = heights_of_peaks(sol_tl[1])

# ╔═╡ df1c6430-8ea8-4bdd-a9be-078a0357fe05
@. model_g(x,p) = p[3]/sqrt(2*π*p[2]^2)*exp(-(x-p[1])^2/(2*p[2]^2))

# ╔═╡ c2074cee-177c-4ece-bf67-0262e504536a
function fit_gauss_D1(pl_D2, pl_D1)
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

# ╔═╡ f1504c87-b24f-4550-9a89-c7759b635218
begin
	gr()
	t_tl, c_tls, c_tl = chrom_after_loop(sol_tl[1])
	p_chrom_tl_proj1 = Plots.plot(t_tl, c_tl, xlabel="time in s")
	for i=1:length(c_tls)
		#if sol_tl[1].Name[i] == selected_solutes[1]
		#	color = :blue
		#elseif sol_tl[1].Name[i] == selected_solutes[2]
		#	color = :red
		#else
		#	color = :orange
		#end
		Plots.plot!(p_chrom_tl_proj1, t_tl, c_tls[i], label=sol_tl[1].Name[i], c=dot_color(sol_tl[1].Name[i], selected_solutes))
		#Plots.scatter!(p_chrom_tl, (sol_tl[1].tR[i], heights_tl[i]), msize=2, c=color)
	end
	fit = fit_gauss_D1(sol_tl[1], sim[2][1][2])
	for i=1:length(fit.Name)
		Plots.scatter!(p_chrom_tl_proj1, fit.tRs[i], fit.heights[i], msize=2, c=dot_color(fit.Name[i], selected_solutes))
		Plots.plot!(p_chrom_tl_proj1, t_tl, model_g(t_tl, fit.fits[i].param), c=dot_color(fit.Name[i], selected_solutes), linestyle=:dash)
	end
	Plots.plot!(p_chrom_tl_proj1, xticks=0:4:3000, legend=false)
	p_chrom_tl_proj1
end

# ╔═╡ d0f14e55-10ad-4ddb-a743-54a3ca28aa85
fit.Name

# ╔═╡ 95897aaf-d2bf-4edd-9cf3-8bf7353fd089
fit.fits

# ╔═╡ fdc6cf54-526f-4535-810d-4fba304cefb3
function fit_gauss_D2(pl_D2, PM)
	heights = heights_of_peaks(pl_D2)
	tR = pl_D2.tR .- fld.(pl_D2.tR, PM).*PM
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

# ╔═╡ aecb6572-cf17-428d-845f-ba3cf2f0a5d4
begin
	plotly()
p_chrom_tl_proj2 = Plots.plot(legend=false, xlabel="time in s")
	for i=1:length(c_tls)
		Plots.plot!(p_chrom_tl_proj2, t_tl.-fld.(t_tl, PM).*PM, c_tls[i], c=dot_color(sol_tl[1].Name[i], selected_solutes))
	end
		fit_2D = fit_gauss_D2(sol_tl[1], PM)
		
		for i=1:length(fit_2D.Name)
			Plots.scatter!(p_chrom_tl_proj2, fit_2D.tRs[i], fit_2D.heights[i], msize=2, c=dot_color(fit_2D.Name[i], selected_solutes))
		end
	p_chrom_tl_proj2
end

# ╔═╡ c8bc0ebc-feee-4bb9-8336-41c24eaf71d6
begin
	plotly()
	p_test = Plots.plot(legend=false, xlabel="time in s")
	for i=1:length(fit_2D.Name)
			Plots.scatter!(p_test, fit_2D.tRs[i], fit_2D.heights[i], msize=2, c=dot_color(fit_2D.Name[i], selected_solutes))
			Plots.plot!(p_test, t_tl.-fld.(t_tl, PM).*PM, model_g(t_tl.-fld.(t_tl, PM).*PM, fit_2D.fits[i].param), c=dot_color(fit_2D.Name[i], selected_solutes), linestyle=:dash)
	end
	p_test
end

# ╔═╡ bf4f1950-4f14-43e6-88cc-635fb079ab13
fit_gauss_D2(sol_tl[1], PM)

# ╔═╡ 78a325ce-a934-4afb-bb31-e1eb98e5349a
sol_final = sol_tl;

# ╔═╡ a76b6b3e-a41c-413e-9b4d-920823e9edda
md"""
## Chromatograms
"""

# ╔═╡ c784c310-3945-4911-8742-932e836f0db0


# ╔═╡ ac2e34ed-a519-48fc-a87e-fb8dde6f3ff5
md"""
### Result from step-by-step
"""

# ╔═╡ 63ec7d6b-df9f-40b8-a988-bbcc2adc52bf
md"""
Before 1st thermal modulator (after GC1 and mod in):
$(embed_display(p_chrom_before_mark))

After 1st thermal modulator:
$(embed_display(p_chrom_TM1))

After loop, before 2nd thermal modulator:
$(embed_display(p_chrom_loop_mark))

After 2nd thermal modulator:
$(embed_display(p_chrom_TM2))

After mod outlet:
$(embed_display(p_chrom_mod_out))

After GC2:

$(embed_display(p_chrom_gc2))

After TL, at detector:
$(embed_display(p_chrom_tl_proj1))

Projection on 2nd dimension:
$(embed_display(p_chrom_tl_proj2))
"""

# ╔═╡ e2eaa21a-fbac-4d00-8663-c69ee4e0a8a9
md"""
## Measurement
"""

# ╔═╡ 17c4fbd1-4a9d-4f7a-922a-a72eb97db33a
md"""
# End
"""

# ╔═╡ Cell order:
# ╠═23c71f14-efdf-11ed-33f2-95304516114d
# ╠═d7f045e5-47c7-44b6-8bce-ce9dc3154fff
# ╠═d419921b-1472-4b66-bae4-469537259814
# ╠═dff75c18-eee5-4f70-9509-0e8228932819
# ╠═326d2249-dfcc-40b6-b8c0-d9ab5ff88627
# ╠═1f5c1ca8-7b61-4dc4-91e6-9fb187073313
# ╠═9a6dea14-f8d4-4b79-9c02-268c0c799515
# ╠═60d86bdd-5527-49da-8c35-ecd508171e6c
# ╠═b63ea500-c2bb-49c3-9915-032fd3e33e14
# ╠═7bb430cd-6bc5-4c48-85f8-3ffce3220d50
# ╠═80150c58-f9b9-43d9-bd7b-e38388f29fd8
# ╠═522a7582-0ff1-42d6-960a-adb952a225d2
# ╠═ca604ac6-fdb1-41d7-83d2-a0c9b2440dad
# ╠═b9e3acb8-a7ec-40ef-ba58-932f79b08fb0
# ╠═9605b6ad-7cdc-47ee-9edf-c3ce531ae9a2
# ╠═06487ba0-8752-40a8-b5d8-07fb08514ffe
# ╠═54c85e34-941d-41a1-9780-cf156f82d765
# ╠═a5c2b339-4879-4367-9d03-ada1e7e6dd46
# ╠═607265de-c784-4ed0-8e89-ebb206785b0b
# ╠═c366df65-59a4-427a-beee-bf0bef5fb409
# ╠═bdbd8cc9-29fb-43ce-850b-8d18980808d3
# ╠═4adbf7e9-232d-4889-92c6-75e8bde0d88d
# ╠═fc36ffe3-b1df-42c8-ad14-e3b90f27a272
# ╠═43f27a81-dca0-4188-baf4-3f8140d9e57c
# ╠═2ec6f339-c7c4-41c1-869e-8cd94d993952
# ╠═ac60e9eb-0dfb-4a8d-a6bb-5289110b40cd
# ╠═7277a646-fd70-4032-9691-569eddeafec6
# ╠═793560f7-d364-4c68-81ec-994441a41059
# ╠═6993afed-4519-4c8a-9cfc-87c72723c444
# ╠═a14d1ca2-46ef-4296-aa1d-06590eca97cf
# ╟─69598506-cf92-4626-b6e1-d8f57bf8388f
# ╠═1a0695c5-0a67-4158-a3bd-c135be31e408
# ╠═6556a858-65ef-49a6-bdaa-562a2b568a70
# ╠═fcdc5015-31d0-4808-8218-f33901437c0b
# ╠═d8d6880e-13c3-4cc9-9dcc-cf616b3d81b1
# ╠═67ba928d-5af0-422e-bd58-c6302c3e41af
# ╠═88fb12fa-3a88-4577-8709-e9f8ec130b2b
# ╠═0b2c8f8b-4648-4af3-a45d-900c5c382860
# ╠═1f62b350-332e-4acf-b907-75362157e4da
# ╠═f1bf21dd-07ee-4e4e-beb5-7d9225b3fc6b
# ╠═ceb6ba3f-aece-4ee3-b757-5904fa04c7ac
# ╠═3e8d2e43-1fa9-40e4-a132-7228664b7546
# ╠═e3e34deb-9249-4811-b74a-44a4c2d06ac2
# ╠═54e20666-b217-415d-b2fd-305980588668
# ╠═ac3db64b-82c8-455a-97aa-6bbddcc03830
# ╠═ee2e3313-9d05-4815-97a2-5b292dc47b68
# ╠═7e088401-a8d4-4648-9380-241777534fe2
# ╠═82c69944-c4e0-4699-87c0-29500b33ba78
# ╠═ce958d12-3dad-4e3c-ab57-f84136610870
# ╠═70597ee8-760d-4c21-af6c-981922208d10
# ╠═b83330dd-4df3-433a-ab29-273e39f8b32c
# ╠═fc7f8b92-f0ee-4a89-ab7c-0c8f5ba1bb0f
# ╠═00886bb1-39ff-48f2-9329-2a195f705a30
# ╠═49f4e7f4-2520-4a97-93c8-495de159e7c4
# ╠═4c385876-bb3b-48fb-aebe-c5554d21df29
# ╠═76ca88a4-24a1-476e-8316-916593552523
# ╠═846c3cab-90d4-4e73-9ce5-c1b41f937cdf
# ╠═3c22e8aa-376d-4eec-9650-3bb102d4a0b8
# ╠═06c3c74b-4ced-405e-89be-a05092fab027
# ╠═67a4dff7-2634-4ae2-adb0-515807988b17
# ╠═84111543-21e2-4891-9050-0c102cd46a52
# ╠═9c6573e8-bc0f-4c40-9eea-9733726fbd76
# ╠═095ee9ba-d974-4625-a47c-3634d91c6a5d
# ╠═f6093dce-30a7-4628-8362-0c5ed8c55ebc
# ╠═beb15214-b1bb-4b39-a2aa-86f0727bc2c4
# ╠═2c1fc71c-5f5d-4e18-84f5-42980d461753
# ╠═74121702-9f4a-410b-961b-7865ab3af941
# ╠═5d0cb969-8a98-4d28-aa5c-233723755f41
# ╠═96200092-8afd-4ce3-8287-79806a331656
# ╠═66eaa639-416c-4167-9e00-14dcbf583a4c
# ╠═bb483a0e-72d7-405e-bf47-c5fa878e0994
# ╠═53aff8ef-4585-4d03-921b-7f19827dc04e
# ╠═b31fa530-922f-4862-b1cf-04d2ec34b1c1
# ╠═86274690-dfe8-4f71-8084-66fe5f03f932
# ╟─6ace764f-1691-4f9a-a2bb-54635468ce8c
# ╠═2d71af52-1457-4bd3-a73e-424944ab34b9
# ╠═d0b3d8b5-3012-43e4-9c09-bdaa97d0e3f6
# ╠═db04ef60-3490-49c3-a5c1-f198c4015b06
# ╠═ebc1b905-eef4-4c66-abb6-b070d19f17d9
# ╠═672e5876-0f0d-483a-9936-604b3fc7f97a
# ╠═231edaa6-8ef3-4a85-b48e-e0e1829d124d
# ╠═64e635f8-f5fd-4935-b44e-b02ebdd96b64
# ╠═78dbaf44-0468-417a-8c69-90c139e2e2ad
# ╠═bb01e67d-fbca-461e-b592-508bc5480b94
# ╠═195f8556-e4e8-4b85-86d1-f70b5a6ce0bf
# ╠═088c0c0f-0722-4498-8306-8286069fc24b
# ╠═29cbe2b2-34fb-45dc-849b-2f678b3f8797
# ╠═14bf7364-ce16-467e-894b-43d5ad1d91dc
# ╠═98713f9d-123a-4eda-8a8c-72c6ada1a3a4
# ╠═caef921f-f06c-471b-a56a-8415ac2d4fa6
# ╠═05a0e4f4-8a1b-4f3c-84b0-c2491840d20c
# ╠═2465e19f-36f8-4a7d-93ed-d6b872532c40
# ╠═1bc4ae12-6f10-499a-b20f-35d5bd297378
# ╠═935485a4-947e-42ad-ae92-9146c3f551ff
# ╠═16e6bf4c-05ee-466d-95d6-e679e4410bc1
# ╠═6accd799-2058-4756-9e31-c0b8b6bc6130
# ╠═32e40a19-ebab-47a6-b4ca-9300c2b04dab
# ╠═93ee192a-7374-43f4-a1a1-a1265fda42a4
# ╠═07922786-3664-4308-a939-6d9cd9cffb33
# ╠═f1504c87-b24f-4550-9a89-c7759b635218
# ╠═0205e873-f3b5-4100-a49b-5824fc5eecf3
# ╠═ec41da66-1aba-4e38-8c69-1cc236157279
# ╠═aecb6572-cf17-428d-845f-ba3cf2f0a5d4
# ╠═c8bc0ebc-feee-4bb9-8336-41c24eaf71d6
# ╠═bf4f1950-4f14-43e6-88cc-635fb079ab13
# ╠═8b500dc4-98de-46af-9a83-f5ff7e8190ee
# ╠═d0f14e55-10ad-4ddb-a743-54a3ca28aa85
# ╠═95897aaf-d2bf-4edd-9cf3-8bf7353fd089
# ╠═9d6944b7-7794-4912-b8ab-ce91672cf0de
# ╠═a0b9918a-1ea2-463a-8587-2b4ca0f1230c
# ╠═6c4c9478-1392-442c-9a26-b654b64c1b14
# ╠═50792b68-8acf-4349-afd5-66267fa9153b
# ╠═08d9257e-6d49-4cf6-97a9-f8b546d9e933
# ╠═df1c6430-8ea8-4bdd-a9be-078a0357fe05
# ╠═c2074cee-177c-4ece-bf67-0262e504536a
# ╠═fdc6cf54-526f-4535-810d-4fba304cefb3
# ╠═78a325ce-a934-4afb-bb31-e1eb98e5349a
# ╠═a76b6b3e-a41c-413e-9b4d-920823e9edda
# ╠═c784c310-3945-4911-8742-932e836f0db0
# ╠═ac2e34ed-a519-48fc-a87e-fb8dde6f3ff5
# ╠═63ec7d6b-df9f-40b8-a988-bbcc2adc52bf
# ╠═e2eaa21a-fbac-4d00-8663-c69ee4e0a8a9
# ╠═17c4fbd1-4a9d-4f7a-922a-a72eb97db33a
