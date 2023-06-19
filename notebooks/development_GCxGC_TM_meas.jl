### A Pluto.jl notebook ###
# v0.19.19

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
	TableOfContents()
end

# ╔═╡ d7f045e5-47c7-44b6-8bce-ce9dc3154fff
using Integrals

# ╔═╡ 08d9257e-6d49-4cf6-97a9-f8b546d9e933
using LsqFit

# ╔═╡ 4dd08d9d-7068-4198-b5d4-5377a787dc8a
using Interpolations

# ╔═╡ d419921b-1472-4b66-bae4-469537259814
md"""
# Development of a thermal modulated GCxGC system in detail
"""

# ╔═╡ dff75c18-eee5-4f70-9509-0e8228932819
md"""
## Periodic rectangle temperature function
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

# ╔═╡ 3a6b6c3c-c839-47ed-bf86-dd1f3149c7f9
# chang this definiton, so that the hot jet comes to the end of the modulation period, not the beginning -> revers Thot and Tcold and ratio is the ratio of PM for which Tcold is applied

# ╔═╡ d0adda7c-8980-4d74-8b79-c3c5edd8131f
md"""
## `ModuleTM`
"""

# ╔═╡ ab33d65f-f0ba-4e5d-b660-2ecd106d8d21
begin
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
		ratio::Float64 # a number, ratio of the duration between hot and cold jet, approx. as rectangular function
		Thot::Float64 # heating with hot jet
		Tcold::Float64 # cooling with cold jet
		flow # an number (constant flow) or a Function
	end

	function ModuleTM(name, L, d, df, sp, tp, shift, pm, ratio, Thot, Tcold)
		TM = ModuleTM(name, L, d, [d], df, [df], sp, tp, shift, pm, ratio, Thot, Tcold, NaN)
		return TM
	end

	function ModuleTM(name, L, d, df, sp, tp, shift, pm, ratio, Thot, Tcold, F)
		TM = ModuleTM(name, L, d, [d], df, [df], sp, tp, shift, pm, ratio, Thot, Tcold, F)
		return TM
	end
end

# ╔═╡ 60d86bdd-5527-49da-8c35-ecd508171e6c
md"""
## Definition System
"""

# ╔═╡ b63ea500-c2bb-49c3-9915-032fd3e33e14
begin
	tsteps_, Tsteps_ = GasChromatographySimulator.conventional_program([50.0, 1.0, 3.0, 200.0, 0.0, 10.0, 225.0, 10.0]) 
	GCxGC_TP = GasChromatographySystems.TemperatureProgram(tsteps_, Tsteps_)
end

# ╔═╡ 3964e72a-4a96-4740-bd66-a224a0be16eb
#savefig(p_meas, "2D_chrom_ketones_new.svg")

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
			elseif typeof(sys.modules[i]) == ModuleTM
				new_modules[i] = ModuleTM(sys.modules[i].name, sys.modules[i].length, sys.modules[i].diameter, sys.modules[i].film_thickness, sys.modules[i].stationary_phase, new_tp, sys.modules[i].shift, sys.modules[i].PM, sys.modules[i].ratio, sys.modules[i].Thot, sys.modules[i].Tcold, sys.modules[i].flow)
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
	modules[3] = ModuleTM("TM1", LM[2], dM*1e-3, dfM*1e-6, spM, TPM, shiftM, PM, ratioM, HotM, ColdM, NaN)
	modules[4] = GasChromatographySystems.ModuleColumn("mod loop", LM[3], dM*1e-3, dfM*1e-6, spM, TPM)
	#modules[5] = ModuleTM("TM2", LM[4], dM*1e-3, dfM*1e-6, spM, TPM, shiftM, PM, ratioM, HotM, ColdM, NaN)
	modules[5] = GasChromatographySystems.ModuleColumn("TM2c", LM[4], dM*1e-3, dfM*1e-6, spM, TPM, NaN) # simulate 2nd Modulator point with the oven temperature -> single stage modulator
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
	# settings for first saved 2D chroms
	#sys = GCxGC_TM(30.0, 0.25, 0.25, "ZB1ms", GCxGC_TP, 0.82, 0.1, 0.1, "Stabilwax", GCxGC_TP, 0.25, 0.1, 0.1, "Stabilwax", 280.0, [0.10, 0.01, 0.80, 0.01, 0.01], 0.1, 0.1, "Stabilwax", 0.0, PM, coldjet/PM, 80.0, -80.0, GCxGC_TP, 0.8, NaN, 0.0)
	
	#sys = GCxGC_TM(31.5, 0.25, 0.25, "ZB1ms", GCxGC_TP, 0.1, 0.1, 0.1, "Stabilwax", GCxGC_TP, 0.25, 0.1, 0.1, "Stabilwax", 280.0, [0.50, 0.01, 0.30, 0.01, 0.38], 0.1, 0.1, "Stabilwax", 0.0, PM, coldjet/PM, 80.0, -80.0, GCxGC_TP, 0.62, NaN, 0.0)

	# new adapted settings
	#sys = GCxGC_TM(31.5, 0.25, 0.25, "ZB1ms", GCxGC_TP, 0.1, 0.1, 0.1, "Stabilwax", GCxGC_TP, 0.25, 0.1, 0.1, "Stabilwax", 280.0, [0.60, 0.01, 0.30, 0.01, 0.50], 0.1, 0.1, "Stabilwax", 0.0, PM, coldjet/PM, 80.0, -80.0, GCxGC_TP, 0.65, NaN, 0.0)
	sys = GCxGC_TM(31.5, 0.25, 0.25, "ZB1ms", GCxGC_TP, 0.1, 0.1, 0.1, "Stabilwax", GCxGC_TP, 0.25, 0.1, 0.1, "Stabilwax", 280.0, [0.60, 0.01, 0.30, 0.01, 0.50], 0.1, 0.1, "Stabilwax", 0.0, PM, coldjet/PM, 80.0, -80.0, GCxGC_TP, 0.65, NaN, 0.0)
end

# ╔═╡ ca604ac6-fdb1-41d7-83d2-a0c9b2440dad
GasChromatographySystems.plot_graph(sys)

# ╔═╡ b9e3acb8-a7ec-40ef-ba58-932f79b08fb0
function flow_restriction(t, module_::GasChromatographySystems.ModuleColumn, options)
	if typeof(module_.temperature) <: GasChromatographySystems.TemperatureProgram
		T_itp = GasChromatographySimulator.temperature_interpolation(module_.temperature.timesteps, module_.temperature.temperaturesteps, module_.temperature.gradient_function, module_.length)
	elseif typeof(module_.temperature) <: Number
		gf(x) = [zero(x), zero(x)]
		T_itp = GasChromatographySimulator.temperature_interpolation([0.0, 3600.0], [module_.temperature, module_.temperature], gf, module_.length)
	end
	κ = GasChromatographySimulator.flow_restriction(module_.length, t, T_itp, module_.diameter, options.mobile_phase; ng=options.ng, vis=options.vis)
	return κ
end

# ╔═╡ 95f57b5b-652c-40c2-b66d-137bb7b2d878
function modulator_temperature(x, t, module_::ModuleTM)
	if typeof(module_.temperature) <: GasChromatographySystems.TemperatureProgram
		T_itp_ = GasChromatographySimulator.temperature_interpolation(module_.temperature.timesteps, module_.temperature.temperaturesteps, module_.temperature.gradient_function, module_.length)
	elseif typeof(module_.temperature) <: Number
		gf(x) = [zero(x), zero(x)]
		T_itp_ = GasChromatographySimulator.temperature_interpolation([0.0, 3600.0], [module_.temperature, module_.temperature], gf, module_.length)
	end
	T_itp(x, t) = therm_mod(t, module_.shift, module_.PM, module_.ratio, module_.Tcold, T_itp_(x, t) .+ module_.Thot .- 273.15) .+ 273.15
	return T_itp(x, t)
end

# ╔═╡ 9605b6ad-7cdc-47ee-9edf-c3ce531ae9a2
function flow_restriction(t, module_::ModuleTM, options)
	T_itp(x, t) = modulator_temperature(x, t, module_)
	κ = GasChromatographySimulator.flow_restriction(module_.length, t, T_itp, module_.diameter, options.mobile_phase; ng=options.ng, vis=options.vis)
	return κ
end

# ╔═╡ 06487ba0-8752-40a8-b5d8-07fb08514ffe
function flow_restrictions(sys)
	kappas = Array{Function}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		f(t) = flow_restriction(t, sys.modules[i], sys.options)
		kappas[i] = f
	end
	return kappas
end

# ╔═╡ a5c2b339-4879-4367-9d03-ada1e7e6dd46
sys

# ╔═╡ 7fc154cd-d1a4-4f14-8750-e92a58da35f0
TM_itp(x, t) = modulator_temperature(x, t, sys.modules[3])

# ╔═╡ daf5e89f-343c-4076-80fd-2dd52ba14f29
TM_itp(0.0,0.0)

# ╔═╡ 051e9ad6-51de-40eb-a3cf-2b929d594400
TM_itp(0.0,3.0)

# ╔═╡ cf03f9b5-62e6-46af-924e-3a16638c768d
TM_itp(0.0,3.1)

# ╔═╡ 607265de-c784-4ed0-8e89-ebb206785b0b
begin
	Plots.plot(500.0:0.01:600.0, TM_itp.(0.0, 500.0:0.01:600.0).-273.115)
end

# ╔═╡ 2e4612e7-ce0f-4ae7-97b0-7076d9bf5651
κ = flow_restrictions(sys)

# ╔═╡ c366df65-59a4-427a-beee-bf0bef5fb409
Plots.plot(0.0:0.01:100.0, κ[3].(0.0:0.01:100.0))

# ╔═╡ c1cb420c-7869-4da8-a240-1583545bde7c
#GasChromatographySystems.plot_pressure_over_time(sys)

# ╔═╡ 93f19f34-5766-4863-bfba-fde045afcdb9
pf = GasChromatographySystems.pressure_functions(sys)

# ╔═╡ 5428a416-48e4-4ebd-a493-54fd30cd2556
#Plots.plot(10.0:0.01:20.0, pf[3].(10.0:0.01:20.0))

# ╔═╡ 6378cb76-57fd-4e59-a183-fd2c745d0f69
Ff = GasChromatographySystems.flow_functions(sys)

# ╔═╡ bdbd8cc9-29fb-43ce-850b-8d18980808d3
md"""
## Simulation
"""

# ╔═╡ 4adbf7e9-232d-4889-92c6-75e8bde0d88d
begin
	db = DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/Estimator/GCsim_d_renamed.csv", header=1, silencewarnings=true))
	db.Tchar = db.Tchar .- 273.15
	db
end

# ╔═╡ fc36ffe3-b1df-42c8-ad14-e3b90f27a272
selected_solutes_ = GasChromatographySystems.common_solutes(db, sys).Name

# ╔═╡ 43f27a81-dca0-4188-baf4-3f8140d9e57c
selected_solutes = selected_solutes_[[14,16,17]]

# ╔═╡ d4bbf7ce-a804-469d-a648-b52908daae51
#selected_solutes = selected_solutes_[29:35]

# ╔═╡ 7277a646-fd70-4032-9691-569eddeafec6
function graph_to_parameters(sys, db_dataframe, selected_solutes; interp=true, dt=1)
	E = collect(edges(sys.g))
	srcE = src.(E)
	dstE = dst.(E)
	if interp == true
		p_func = GasChromatographySystems.interpolate_pressure_functions(sys; dt=dt)
	else
		p_func = GasChromatographySystems.pressure_functions(sys)
	end
	parameters = Array{GasChromatographySimulator.Parameters}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		col = GasChromatographySimulator.Column(sys.modules[i].length, sys.modules[i].diameter, [sys.modules[i].diameter], sys.modules[i].film_thickness, [sys.modules[i].film_thickness], sys.modules[i].stationary_phase, sys.options.mobile_phase)

		pin_steps = sys.pressurepoints[srcE[i]].pressure_steps
		pout_steps = sys.pressurepoints[dstE[i]].pressure_steps
		pin_itp = p_func[srcE[i]]
		pout_itp = p_func[dstE[i]]	
		if typeof(sys.modules[i].temperature) <: GasChromatographySystems.TemperatureProgram
			time_steps = sys.modules[i].temperature.timesteps
			temp_steps = sys.modules[i].temperature.temperaturesteps	
			gf = sys.modules[i].temperature.gradient_function
			a_gf = sys.modules[i].temperature.a_gradient_function
			T_itp_ = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, sys.modules[i].length)
		elseif typeof(sys.modules[i].temperature) <: Number
			time_steps = GasChromatographySystems.common_timesteps(sys)
			temp_steps = sys.modules[i].temperature.*ones(length(time_steps))
			gf(x) = zero(x).*ones(length(time_steps))
			a_gf = [zeros(length(time_steps)) zeros(length(time_steps)) ones(length(time_steps)) zeros(length(time_steps))]
			T_itp_ = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, sys.modules[i].length)
		end
		
		T_itp(x, t) = if typeof(sys.modules[i]) == GasChromatographySystems.ModuleColumn
			T_itp_(x,t) 
		elseif typeof(sys.modules[i]) == ModuleTM
			modulator_temperature(x, t, sys.modules[i])
		end
		
		prog = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)

		sub = GasChromatographySimulator.load_solute_database(db_dataframe, sys.modules[i].stationary_phase, sys.options.mobile_phase, selected_solutes, NaN.*ones(length(selected_solutes)), NaN.*ones(length(selected_solutes)))

		opt = GasChromatographySimulator.Options(sys.options.alg, sys.options.abstol, sys.options.reltol, sys.options.Tcontrol, sys.options.odesys, sys.options.ng, sys.options.vis, sys.options.control, sys.options.k_th)

		parameters[i] = GasChromatographySimulator.Parameters(col, prog, sub, opt)
	end
	return parameters
end

# ╔═╡ 793560f7-d364-4c68-81ec-994441a41059
par = graph_to_parameters(sys, db, selected_solutes)

# ╔═╡ d40361fe-9e0c-4df8-a09c-0c5bf143dbf6
par[3].prog.T_itp(0.0, 0.0)

# ╔═╡ 5d63fd50-208e-4db1-9609-4f08c493d90c
par[1].prog.T_itp(0.0, 0.0)

# ╔═╡ 6993afed-4519-4c8a-9cfc-87c72723c444
begin
	gr()#plotly()
	p_TM_Prog = Plots.plot(500.0:0.01:800.0, par[3].prog.T_itp.(0.0, 500.0:0.01:800.0).-273.15)
	Plots.plot!(p_TM_Prog, 500.0:0.01:800.0, par[4].prog.T_itp.(0.0, 500.0:0.01:800.0).-273.15)
	#Plots.plot!(p_TM_Prog, 500.0:0.01:600.0, par[5].prog.T_itp.(0.0, 500.0:0.01:600.0).-273.15)
	Plots.plot!(p_TM_Prog, xlabel="time in s", ylabel="temperature in °C", legend=false)
end

# ╔═╡ e1785bfc-dfc7-4ae4-b524-b93608845384
Plots.plot(2204.0:0.01:2224.0, par[3].prog.T_itp.(0.0, 2204.0:0.01:2224.0).-273.15)

# ╔═╡ b0469761-2537-47ec-b308-cc7be689ca36
#savefig(p_TM_Prog, "detail_temperatureprogram_TM_and_oven.svg")

# ╔═╡ a14d1ca2-46ef-4296-aa1d-06590eca97cf
paths = GasChromatographySystems.all_paths(sys.g, 1)[2]

# ╔═╡ 69598506-cf92-4626-b6e1-d8f57bf8388f
md"""
### Simulation with the old function
"""

# ╔═╡ 1a0695c5-0a67-4158-a3bd-c135be31e408
sim = GasChromatographySystems.simulate_along_paths(sys, paths, par)

# ╔═╡ 6556a858-65ef-49a6-bdaa-562a2b568a70
md"""
### (deactivated) Simulation using the new function
"""

# ╔═╡ fcdc5015-31d0-4808-8218-f33901437c0b
#sim_ = simulate_along_paths(sys, paths, par, abstol=1e-6, reltol=1e-3)

# ╔═╡ 1f62b350-332e-4acf-b907-75362157e4da
sim_[2][1][3]

# ╔═╡ f1bf21dd-07ee-4e4e-beb5-7d9225b3fc6b
sim_[2][1][5]

# ╔═╡ 6e8f2101-8245-483c-9c05-e4b16c260607
function slicing(tR, τR, PM, par::GasChromatographySimulator.Parameters, prev_par::GasChromatographySimulator.Parameters; nτ=6, τ₀=zeros(length(τR)), abstol=1e-12, reltol=1e-9)
	# slicing of the incoming peaks at a thermal modulator, assuming that everything from a peak moving into the thermal modulator during the time period of PM is captured on the modulator (in reality, the part which comes to the modulator during the hot jet is not captured, but should be captured at the second modulator point)
	init_t_start = (fld.(tR.-nτ.*τR, PM)).*PM # start time of the peaks
	init_t_end = (fld.(tR.+nτ.*τR, PM)).*PM # end time of the peaks
	n_slice = Int.((init_t_end .- init_t_start)./PM .+ 1) # number of slices for every substance
	sub_TM = Array{GasChromatographySimulator.Substance}(undef, sum(n_slice))
	A = Array{Float64}(undef, sum(n_slice))
	g(x,p) = 1/sqrt(2*π*p[2]^2)*exp(-(x-p[1])^2/(2*p[2]^2))
	ii = 1
	for i=1:length(n_slice)
		for j=1:n_slice[i]
			sub_TM[ii] = GasChromatographySimulator.Substance(prev_par.sub[i].name, prev_par.sub[i].CAS, prev_par.sub[i].Tchar, prev_par.sub[i].θchar, prev_par.sub[i].ΔCp, prev_par.sub[i].φ₀, prev_par.sub[i].ann*", slice$(j)", prev_par.sub[i].Cag, init_t_start[i]+(j-1)*PM, τ₀[i])
			p = [tR[i], τR[i]]
			prob = IntegralProblem(g, init_t_start[i]+(j-1)*PM, init_t_start[i]+j*PM, p)
			A[ii] = solve(prob, QuadGKJL()).u
			ii = ii + 1
		end
	end
	newopt = GasChromatographySimulator.Options(par.opt.alg, abstol, reltol, par.opt.Tcontrol, par.opt.odesys, par.opt.ng, par.opt.vis, par.opt.control, par.opt.k_th)
	newpar = GasChromatographySimulator.Parameters(par.col, par.prog, sub_TM, newopt)
	
	return newpar, A
end

# ╔═╡ e3e34deb-9249-4811-b74a-44a4c2d06ac2
md"""
## Modulation step-by-step

Using the results of `sim` (Simulation with the old function) just before the 1st modulator spot
"""

# ╔═╡ ce317588-2f42-477e-a0a6-cc40047e1506
# reconstruction of function simulate_along_paths(sys, paths, par_sys; t₀=zeros(length(par_sys[1].sub)), τ₀=zeros(length(par_sys[1].sub)), refocus=falses(ne(sys.g)), τ₀_focus=zeros(length(par_sys[1].sub)), nτ=6):

# ╔═╡ 54e20666-b217-415d-b2fd-305980588668
md"""
### Peak before 1st Modulator spot
"""

# ╔═╡ ac3db64b-82c8-455a-97aa-6bbddcc03830
sim[2][1][2]

# ╔═╡ 7e088401-a8d4-4648-9380-241777534fe2
#savefig(p_chrom_before, "peak_before_modulator_new.svg")

# ╔═╡ 82c69944-c4e0-4699-87c0-29500b33ba78
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
begin

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
end

# ╔═╡ 70597ee8-760d-4c21-af6c-981922208d10
md"""
### 1st Modulator spot
"""

# ╔═╡ b83330dd-4df3-433a-ab29-273e39f8b32c
begin
	# increase the peak width befor the modulator by 50%
	sim[2][1][2][!,:τR] = 1.5.*sim[2][1][2][!,:τR]
	sim[2][1][2]
end

# ╔═╡ 325a8e2e-ae2f-4dd9-a32b-00c7adf47461
begin
	# if next module == moduleTM
	par_prev = par[2]
	# initial peak width?
	τ0 = 0.0 # or PM?
#	init_t_start = (fld.(sim[2][1][2].tR.-6*sim[2][1][2].τR, PM).-1).*PM # is the -1 here correct?????
#	init_t_end = (fld.(sim[2][1][2].tR.+6*sim[2][1][2].τR, PM).-1).*PM 
	nτ = 6
	init_t_start = (fld.(sim[2][1][2].tR.-nτ.*sim[2][1][2].τR, PM)).*PM # start time of the peaks
	init_t_end = (fld.(sim[2][1][2].tR.+nτ.*sim[2][1][2].τR, PM)).*PM # end time of the peaks
	n_slice = Int.((init_t_end .- init_t_start)./PM .+ 1) # number of slices for every substance
	sub_TM = Array{GasChromatographySimulator.Substance}(undef, sum(n_slice))
	ii = 1
	for i=1:length(n_slice)
		for j=1:n_slice[i]
			sub_TM[ii] = GasChromatographySimulator.Substance(par_prev.sub[i].name, par_prev.sub[i].CAS, par_prev.sub[i].Tchar, par_prev.sub[i].θchar, par_prev.sub[i].ΔCp, par_prev.sub[i].φ₀, par_prev.sub[i].ann*", slice$(j)", par_prev.sub[i].Cag, init_t_start[i]+(j-1)*PM, τ0)
			ii = ii + 1
		end
	end
	newopt = GasChromatographySimulator.Options(par[3].opt.alg, 1e-18, 1e-14, par[3].opt.Tcontrol, par[3].opt.odesys, par[3].opt.ng, par[3].opt.vis, par[3].opt.control, par[3].opt.k_th)
	#newopt = GasChromatographySimulator.Options(BS3(), 1e-6, 1e-3, par[3].opt.Tcontrol, par[3].opt.odesys, par[3].opt.ng, par[3].opt.vis, par[3].opt.control, par[3].opt.k_th)
	newpar = GasChromatographySimulator.Parameters(par[3].col, par[3].prog, sub_TM, newopt)
	
end

# ╔═╡ c626e4de-3c65-4230-bc23-2fa6ea9bbbd9
par[3]

# ╔═╡ 03aeeaac-0a83-4d34-944f-d3457904d69d
newpar

# ╔═╡ e4172084-2dd7-4daf-b2b9-cc84891c06f2
par_slice_1st_TM, A_slice_1st_TM = slicing(sim[2][1][2].tR, sim[2][1][2].τR, PM, par[3], par[2]; nτ=6)

# ╔═╡ 846c3cab-90d4-4e73-9ce5-c1b41f937cdf
function add_A_to_pl!(pl, A, par)
	sort_A = Array{Float64}(undef, length(par.sub))
	for i=1:length(par.sub)
		ii_name = findall(par.sub[i].CAS.==pl.CAS)
		ii_ann = findall(par.sub[i].ann.==pl.ann)
		ii = intersect(ii_name, ii_ann)[1]
		sort_A[ii] = A[i]
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

# ╔═╡ 7708f828-0ea3-4568-8614-ff0a724b9eea
# from GasChromatographySimulator
# modified
function peaklist(sol, par)
	n = length(par.sub)
    # sol is solution from ODE system
    Name = Array{String}(undef, n)
	CAS = Array{String}(undef, n) # mod
	ann = Array{String}(undef, n) # mod
    tR = Array{Float64}(undef, n)
    TR = Array{Float64}(undef, n)
    σR = Array{Float64}(undef, n)
    uR = Array{Float64}(undef, n)
    τR = Array{Float64}(undef, n)
    kR = Array{Float64}(undef, n)
    Res = fill(NaN, n)
    Δs = fill(NaN, n)
    Threads.@threads for i=1:n
        Name[i] = par.sub[i].name
		CAS[i] = par.sub[i].CAS # mod 
		ann[i] = par.sub[i].ann # mod
        if sol[i].t[end]==par.col.L
            tR[i] = sol[i].u[end][1]
            TR[i] = par.prog.T_itp(par.col.L, tR[i]) - 273.15 
            uR[i] = 1/GasChromatographySimulator.residency(par.col.L, tR[i], par.prog.T_itp, par.prog.Fpin_itp, par.prog.pout_itp, par.col.L, par.col.d, par.col.df, par.col.gas, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀; vis=par.opt.vis, control=par.opt.control, k_th=par.opt.k_th)
            τR[i] = sqrt(sol[i].u[end][2])
            σR[i] = τR[i]*uR[i]
            kR[i] = GasChromatographySimulator.retention_factor(par.col.L, tR[i], par.prog.T_itp, par.col.d, par.col.df, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀; k_th=par.opt.k_th)
        else
            tR[i] = NaN
            TR[i] = NaN
            uR[i] = NaN
            τR[i] = NaN
            σR[i] = NaN
            kR[i] = NaN
        end
    end  
    df = sort!(DataFrame(Name = Name, CAS = CAS, tR = tR, τR = τR, TR=TR, σR = σR, uR = uR, kR = kR, ann=ann), [:tR])
    Threads.@threads for i=1:n-1
        Res[i] = (df.tR[i+1] - df.tR[i])/(2*(df.τR[i+1] + df.τR[i]))
        Δs[i] = (df.tR[i+1] - df.tR[i])/(df.τR[i+1] - df.τR[i]) * log(df.τR[i+1]/df.τR[i])
    end
    df[!, :Res] = Res
    df[!, :Δs] = Δs 
	
    return df
end

# ╔═╡ f48d790f-6ef0-4e31-a80d-8acaae8e454f
# from GasChromatographySimulator
function simulate(par)#mod
    if par.opt.odesys==true
        #sol = solve_system_multithreads(par, dt=dt)#mod
		sol = GasChromatographySimulator.solve_system_multithreads(par)
    	pl = peaklist(sol, par)
        return pl, sol
	else
		sol, peak = GasChromatographySimulator.solve_multithreads(par)
    	pl = peaklist(sol, peak, par)
        return pl, sol, peak
	end
end

# ╔═╡ bac5d026-0b84-412f-98bb-8ddedb2b92d9
function simulate_along_paths(sys, paths, par_sys; t₀=zeros(length(par_sys[1].sub)), τ₀=zeros(length(par_sys[1].sub)), refocus=falses(ne(sys.g)), τ₀_focus=zeros(length(par_sys[1].sub)), nτ=6, abstol=1e-14, reltol=1e-9)
	#par_sys = graph_to_parameters(sys, db_dataframe, selected_solutes)
	E = collect(edges(sys.g))
	peaklists = Array{Array{DataFrame,1}}(undef, length(paths))
	solutions = Array{Array{Any,1}}(undef, length(paths))
	path_pos = Array{String}(undef, length(paths))
	new_par_sys = Array{GasChromatographySimulator.Parameters}(undef, length(par_sys))
	visited_E = falses(length(E))
	#As = Array{Array{DataFrame,1}}(undef, length(paths))
	for i=1:length(paths)
		i_par = GasChromatographySystems.index_parameter(sys.g, paths[i])
		if GasChromatographySystems.path_possible(sys, paths[i]) == true
			path_pos[i] = "path is possible"
			peaklists_ = Array{DataFrame}(undef, length(i_par))
			solutions_ = Array{Any}(undef, length(i_par))
			#As_ = Array{DataFrame}(undef, length(i_par))
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
						peaklists_[j], solutions_[j] = simulate(new_par_sys[i_par[j]])
						peaklists_[j][!,:A] = ones(length(peaklists_[j].Name))
						# increase peakwidth of 1st dimension by 100%
						peaklists_[j][!,:τR] = 2.0.*peaklists_[j][!,:τR]
					elseif typeof(sys.modules[i_par[j]]) == ModuleTM
						if refocus[i_par[j]] == true
							τ₀=τ₀_focus
						else
							τ₀=peaklists_[j-1].τR # PM?
						end
						new_par_sys[i_par[j]], A = slicing(peaklists_[j-1].tR, peaklists_[j-1].τR, sys.modules[i_par[j]].PM, par_sys[i_par[j]], new_par_sys[i_par[j-1]]; nτ=nτ, τ₀=τ₀, abstol=abstol, reltol=reltol)
						peaklists_[j], solutions_[j] = simulate(new_par_sys[i_par[j]])
						# par.sub for all future modules has to be filled with the slices!
						peaklists_[j][!,:A] = A
					else
						if refocus[i_par[j]] == true
							new_par_sys[i_par[j]] = change_initial_focussed(par_sys[i_par[j]], new_par_sys[i_par[j-1]], peaklists_[j-1]; τ₀=τ₀_focus)
						elseif typeof(sys.modules[i_par[j-1]]) == ModuleTM
							# if previous module is a ModuleTM, than the initial peak width is set to equivalent of the spatial width of the peak equal to the length of the ModuleTM
							# spacial width of the previous modulator spot * solute residency
							τ₀_TM = new_par_sys[i_par[j-1]].col.L ./ peaklists_[j-1].uR #alternative could be the distance migrated until the hot jet is activated
							new_par_sys[i_par[j]] = change_initial(par_sys[i_par[j]], new_par_sys[i_par[j-1]], peaklists_[j-1].tR, τ₀_TM)
						else
							new_par_sys[i_par[j]] = change_initial(par_sys[i_par[j]], new_par_sys[i_par[j-1]], peaklists_[j-1].tR, peaklists_[j-1].τR)
						end
						peaklists_[j], solutions_[j] = simulate(new_par_sys[i_par[j]])
						peaklists_[j][!,:A] = ones(length(peaklists_[j].Name))
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

# ╔═╡ f436afe9-0440-4149-be61-b426e243b34d
begin
	sol_1st_TM = simulate(newpar)
	
	add_A_to_pl!(sol_1st_TM[1], A_slice_1st_TM, newpar)
end

# ╔═╡ f087274e-5147-48f9-8200-4d6ed90397c9
par[3].col.L ./ sol_1st_TM[1].uR

# ╔═╡ 90d32247-3e00-41ee-bd2e-b8e698c2613c
sol_1st_TM[1].ann

# ╔═╡ d98548cc-30bd-48c5-b408-0227caaa5fd1
sol_1st_TM[1]

# ╔═╡ 402e9e25-3d5b-46c2-8b5e-4a472fb357f7
traces(sol_1st_TM[2], newpar, 37)

# ╔═╡ 9c5ec167-212d-413a-a2bc-fb00d015af51
begin
	p_Tt = Plots.plot(xlabel="time in s", ylabel="temperature in °C", legend=false)
	for i=33:41
		trace = traces(sol_1st_TM[2], newpar, i)
		Plots.plot!(p_Tt, trace.t, trace.T.-273.15)
	end
	p_Tt
end

# ╔═╡ 67a4dff7-2634-4ae2-adb0-515807988b17
md"""
### Chromatogram after 1st TM
"""

# ╔═╡ 1e057e8c-66f3-43d1-9674-6d030b094a74
sol_1st_TM[1]

# ╔═╡ 6fa14ed4-f9ee-4279-a020-6abe60f857d7
sol_1st_TM[1]

# ╔═╡ 9c6573e8-bc0f-4c40-9eea-9733726fbd76
function chrom_after_modulation(pl, par; nτ=6)
	tstart = Array{Float64}(undef, length(pl.Name))
	tend = Array{Float64}(undef, length(pl.Name))
	t = Array{Array{Float64}}(undef, length(pl.Name))
	c = Array{Array{Float64}}(undef, length(pl.Name))
	τR = Array{Float64}(undef, length(pl.Name))
	for i=1:length(pl.Name)
		τR[i] = par.col.L/pl.uR[i] 
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

# ╔═╡ d5b4ac08-43e0-4d2c-ab7c-9fea2a926ff3
chrom_after_modulation(sol_1st_TM[1], newpar )

# ╔═╡ 095ee9ba-d974-4625-a47c-3634d91c6a5d
begin
	t_TM1, c_TM1s, c_TM1 = chrom_after_modulation(sol_1st_TM[1], newpar)
	p_chrom_TM1 = Plots.plot(t_TM1, c_TM1, xlabel="time in s")
	#for i=1:length(c_TM1s)
	#	Plots.plot!(p_chrom_TM1, t_TM1, c_TM1s[i], label=sol_1st_TM[1].Name[i])
	#end
	Plots.plot!(p_chrom_TM1, xticks=0:4:3000, legend=false)
	p_chrom_TM1
end

# ╔═╡ 471878ef-b115-411c-9d9d-764d6db08069
#savefig(p_chrom_TM1, "Chrom_after_first_TM_new.svg")

# ╔═╡ f6093dce-30a7-4628-8362-0c5ed8c55ebc
md"""
### Loop between Modulator spots
"""

# ╔═╡ beb15214-b1bb-4b39-a2aa-86f0727bc2c4
#=begin
	# loop between the modulator spots
	par_prev_ = newpar
	# initial peak width?
	# spacial width of the previous modulator spot * solute residency
	τ0_ = par_prev_.col.L ./ sol_1st_TM[1].uR
	sub_loop = Array{GasChromatographySimulator.Substance}(undef, length(sub_TM))
	for i=1:length(sub_TM)
		# should check with CAS/name and ann for the same solute !!!!
		sub_loop[i] = GasChromatographySimulator.Substance(par_prev_.sub[i].name, par_prev_.sub[i].CAS, par_prev_.sub[i].Tchar, par_prev_.sub[i].θchar, par_prev_.sub[i].ΔCp, par_prev_.sub[i].φ₀, par_prev_.sub[i].ann, par_prev_.sub[i].Cag, sol_1st_TM[1].tR[i], τ0_[i])
	end
	newpar_loop = GasChromatographySimulator.Parameters(par[4].col, par[4].prog, sub_loop, par[4].opt)
end=#

# ╔═╡ 6a3f4043-63ba-44ef-b9ce-a391bc697fe6
function change_parameters(par_, prev_par, prev_pl; prev_TM=true)
	if prev_TM == true # if the previous module is a thermal modulator, than the new initial peak width is approximated as length of modulator spot divided by the velocity of the solute at the end of the modulator spot 
		τ0_ = prev_par.col.L./prev_pl.uR
	else
		τ0_ = prev_pl.τR
	end
	new_sub = Array{GasChromatographySimulator.Substance}(undef, length(prev_par.sub))
	for i=1:length(prev_par.sub)
		# filter for correct CAS and annotation (slice number)
		ii_name = findall(prev_par.sub[i].CAS.==prev_pl.CAS)
		ii_ann = findall(prev_par.sub[i].ann.==prev_pl.ann)
		ii_ = intersect(ii_name, ii_ann)[1]
		new_sub[ii_] = GasChromatographySimulator.Substance(prev_par.sub[i].name, prev_par.sub[i].CAS, prev_par.sub[i].Tchar, prev_par.sub[i].θchar, prev_par.sub[i].ΔCp, prev_par.sub[i].φ₀, prev_par.sub[i].ann, prev_par.sub[i].Cag, prev_pl.tR[ii_], τ0_[ii_])
	end
	# here changes of options could be applied
	new_par = GasChromatographySimulator.Parameters(par_.col, par_.prog, new_sub, par_.opt)
	return new_par
end

# ╔═╡ 2c1fc71c-5f5d-4e18-84f5-42980d461753
newpar_loop = change_parameters(par[4], newpar, sol_1st_TM[1], prev_TM=true)

# ╔═╡ 74121702-9f4a-410b-961b-7865ab3af941
begin
	sol_loop = simulate(newpar_loop)
	add_A_to_pl!(sol_loop[1], sol_1st_TM[1].A, newpar_loop)
end

# ╔═╡ 975b24c5-6557-4d7c-a3be-4d6784fcf66d
sol_loop[1].ann

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

# ╔═╡ a848d367-8a2e-4eb1-9fa9-bb94e83d34ff
sol_1st_TM[1]

# ╔═╡ c5f01fec-7bb6-4157-a125-f7c929dd69e6
sol_1st_TM[1].Name[5], sol_1st_TM[1].ann[5], sol_1st_TM[1].A[5]

# ╔═╡ eefc168f-7de7-45fa-bd25-ea3e7df4b685
sol_loop[1]

# ╔═╡ c2b343a8-f794-42f3-a237-477f00d2d6b3
sol_loop[1].Name[5], sol_loop[1].ann[5], sol_loop[1].A[5]

# ╔═╡ 96200092-8afd-4ce3-8287-79806a331656
# time in loop
sol_loop[1].tR .- sol_1st_TM[1].tR

# ╔═╡ e662f5a1-79ee-4adb-8d98-ff0b743c95da
begin
	gr()
	t_loop, c_loops, c_loop = chrom_after_loop(sol_loop[1])
	p_chrom_loop = Plots.plot(t_loop, c_loop, xlabel="time in s")
	for i=1:length(c_TM1s)
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

# ╔═╡ 0baa9993-b50b-45ce-a1b7-ed6e45826187
#savefig(p_chrom_loop, "Chrom_after_loop_new.svg")

# ╔═╡ 86274690-dfe8-4f71-8084-66fe5f03f932
md"""
### 2nd Modulator spot
"""

# ╔═╡ e1ce8764-cec4-4c43-9354-ccceea375f3f
# for now, this 2nd modulator spot is deactivated and works with the oven temperature

# ╔═╡ 12c5ce26-4a91-4b1b-b86e-0465c1c155a3
#=begin
	# if next module == moduleTM
	par_prev__ = newpar_loop
	# initial peak width?
	τ0__ = sol_loop[1].τR
	init_t_start__ = (fld.(sol_loop[1].tR.-6*sol_loop[1].τR, PM)).*PM
	init_t_end__ = (fld.(sol_loop[1].tR.+6*sol_loop[1].τR, PM)).*PM 
	n_slice__ = Int.((init_t_end__ .- init_t_start__)./PM .+ 1)
	sub_TM__ = Array{GasChromatographySimulator.Substance}(undef, sum(n_slice__))
	if sum(n_slice__) == length(par_prev__.sub)
		for i=1:sum(n_slice__)
			sub_TM__[i] = GasChromatographySimulator.Substance(par_prev__.sub[i].name, par_prev__.sub[i].CAS, par_prev__.sub[i].Tchar, par_prev__.sub[i].θchar, par_prev__.sub[i].ΔCp, par_prev__.sub[i].φ₀, par_prev__.sub[i].ann, par_prev__.sub[i].Cag, init_t_start__[i], τ0__[i])
		end
	else
		ii__ = 1
		for i=1:length(n_slice__)
			for j=1:n_slice__[i] # additional slices can appear
				#sub_TM__[ii] = GasChromatographySimulator.Substance(par_prev__.sub[i].name, par_prev__.sub[i].CAS, par_prev__.sub[i].Tchar, par_prev__.sub[i].θchar, par_prev__.sub[i].ΔCp, par_prev__.sub[i].φ₀, par_prev__.sub[i].ann*", slice$(j)", par_prev__.sub[i].Cag, init_t_start__[i]+(j-1)*PM, τ0__[ii__])
				ii__ = ii__ + 1
			end
		end
	end
	newopt__ = GasChromatographySimulator.Options(par[5].opt.alg, 1e-14, 1e-6, par[5].opt.Tcontrol, par[5].opt.odesys, par[5].opt.ng, par[5].opt.vis, par[5].opt.control, par[5].opt.k_th)
	newpar__ = GasChromatographySimulator.Parameters(par[5].col, par[5].prog, sub_TM__, newopt__)
	
end=#

# ╔═╡ ece4ada0-9944-40e1-b90d-8c1cd34209f6
newpar_2nd_TM = change_parameters(par[5], newpar_loop, sol_loop[1], prev_TM=false)

# ╔═╡ 31c843d2-e4ad-47b0-b90e-409b9bdb0f4c
par_slice_2nd_TM, A_slice_2nd_TM = slicing(sol_loop[1].tR, τ0___, PM, par[5], newpar_loop; nτ=6)

# ╔═╡ 4c85cfd1-cb1b-449d-a101-e537abe25e86
begin
	sol_2nd_TM = simulate(newpar_2nd_TM)
	add_A_to_pl!(sol_2nd_TM[1], sol_loop[1].A, newpar_2nd_TM)
end

# ╔═╡ dca2c405-93dc-460c-8a06-9892154d503a
# time in 2nd TM (at oven temperature)
sol_2nd_TM[1].tR .- sol_loop[1].tR

# ╔═╡ 6010a66c-d23c-4127-afae-62b6fddcb541
begin
	gr()
	t_2nd_TM, c_2nd_TMs, c_2nd_TM = chrom_after_loop(sol_2nd_TM[1])
	p_chrom_2nd_TM = Plots.plot(t_2nd_TM, c_2nd_TM, xlabel="time in s")
	for i=1:length(c_2nd_TMs)
		if sol_2nd_TM[1].Name[i] == selected_solutes[1]
			color = :blue
		elseif sol_2nd_TM[1].Name[i] == selected_solutes[2]
			color = :red
		else
			color = :orange
		end
		Plots.plot!(p_chrom_2nd_TM, t_2nd_TM, c_2nd_TMs[i], label=sol_2nd_TM[1].Name[i], c=color)
	end
	Plots.plot!(p_chrom_2nd_TM, xticks=0:4:3000, legend=false)
	p_chrom_2nd_TM
end

# ╔═╡ 231edaa6-8ef3-4a85-b48e-e0e1829d124d
#savefig(p_chrom_2nd_TM, "Chrom_after_2nd_TM_new.svg")

# ╔═╡ 64e635f8-f5fd-4935-b44e-b02ebdd96b64
md"""
### connection to GC2
"""

# ╔═╡ f388ca37-0167-4459-858b-619f1d2ce2d7
newpar_mod_out = change_parameters(par[6], newpar_2nd_TM, sol_2nd_TM[1], prev_TM=false)

# ╔═╡ bb01e67d-fbca-461e-b592-508bc5480b94
begin
	sol_mod_out = simulate(newpar_mod_out)
	add_A_to_pl!(sol_mod_out[1], sol_2nd_TM[1].A, newpar_mod_out)
end

# ╔═╡ c959c1ae-2c46-412d-b23b-fadb87d4d3be
# time in mod out
sol_mod_out[1].tR .- sol_2nd_TM[1].tR

# ╔═╡ 088c0c0f-0722-4498-8306-8286069fc24b
begin
	gr()
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
newpar_gc2 = change_parameters(par[7], newpar_mod_out, sol_mod_out[1], prev_TM=false)

# ╔═╡ caef921f-f06c-471b-a56a-8415ac2d4fa6
begin
	sol_gc2 = simulate(newpar_gc2)
	add_A_to_pl!(sol_gc2[1], sol_mod_out[1].A, newpar_gc2)
end 

# ╔═╡ 05a0e4f4-8a1b-4f3c-84b0-c2491840d20c
# time in gc2
sol_gc2[1].tR .- sol_mod_out[1].tR

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
newpar_tl = change_parameters(par[8], newpar_gc2, sol_gc2[1], prev_TM=false)

# ╔═╡ 32e40a19-ebab-47a6-b4ca-9300c2b04dab
begin
	sol_tl = simulate(newpar_tl)
	add_A_to_pl!(sol_tl[1], sol_gc2[1].A, newpar_tl)
end 

# ╔═╡ 93ee192a-7374-43f4-a1a1-a1265fda42a4
# time in tl
sol_tl[1].tR .- sol_gc2[1].tR

# ╔═╡ 8b500dc4-98de-46af-9a83-f5ff7e8190ee
sim[2][1][2].Name, sim[2][1][2].tR, sim[2][1][2].τR

# ╔═╡ b0490a91-8a34-4fee-b819-de4f6ecb7a27
1897.52 - 1901.32

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
	p_chrom_tl = Plots.plot(t_tl, c_tl, xlabel="time in s")
	for i=1:length(c_tls)
		#if sol_tl[1].Name[i] == selected_solutes[1]
		#	color = :blue
		#elseif sol_tl[1].Name[i] == selected_solutes[2]
		#	color = :red
		#else
		#	color = :orange
		#end
		Plots.plot!(p_chrom_tl, t_tl, c_tls[i], label=sol_tl[1].Name[i], c=dot_color(sol_tl[1].Name[i], selected_solutes))
		#Plots.scatter!(p_chrom_tl, (sol_tl[1].tR[i], heights_tl[i]), msize=2, c=color)
	end
	fit = fit_gauss_D1(sol_tl[1], sim[2][1][2])
	for i=1:length(fit.Name)
		Plots.scatter!(p_chrom_tl, fit.tRs[i], fit.heights[i], msize=2, c=dot_color(fit.Name[i], selected_solutes))
		Plots.plot!(p_chrom_tl, t_tl, model_g(t_tl, fit.fits[i].param), c=dot_color(fit.Name[i], selected_solutes), linestyle=:dash)
	end
	Plots.plot!(p_chrom_tl, xticks=0:4:3000, legend=false)
	p_chrom_tl
end

# ╔═╡ d0f14e55-10ad-4ddb-a743-54a3ca28aa85
fit.Name

# ╔═╡ 95897aaf-d2bf-4edd-9cf3-8bf7353fd089
fit.fits

# ╔═╡ 78a325ce-a934-4afb-bb31-e1eb98e5349a
sol_final = sol_tl;

# ╔═╡ a76b6b3e-a41c-413e-9b4d-920823e9edda
md"""
## Chromatograms
"""

# ╔═╡ ac2e34ed-a519-48fc-a87e-fb8dde6f3ff5
md"""
### Result from step-by-step
"""

# ╔═╡ 915d8972-5423-4dde-8039-9a912bd0ae61
slice_tR1 = fld.(sol_final[1].tR, PM).*PM

# ╔═╡ 556c157c-321d-4eb8-b2ef-ee8c4fe0f038
slice_tR2 = sol_final[1].tR .- fld.(sol_final[1].tR, PM).*PM

# ╔═╡ a5eb5c69-de38-4dbb-83f7-4b751d26562d
sliced_results = DataFrame(Name=sol_final[1].Name, tR1=slice_tR1, tR2=slice_tR2, A=sol_final[1].A)

# ╔═╡ c78ea8fb-d145-4862-973b-06bd13c03f52
findall(selected_solutes[3].==sliced_results.Name)

# ╔═╡ 80a012dd-738e-4332-857a-90203867a03f
begin
	approx_max = DataFrame(Name=selected_solutes)
	tR1_max = Array{Float64}(undef, length(selected_solutes))
	tR2_max = Array{Float64}(undef, length(selected_solutes))
	A_max = Array{Float64}(undef, length(selected_solutes))
	for i=1:length(selected_solutes)
		ii_solutes = findall(selected_solutes[i].==sliced_results.Name)
		As = sliced_results.A[ii_solutes]
		ii = ii_solutes[findfirst(maximum(As).==As)]
		tR1_max[i] = sliced_results.tR1[ii]
		tR2_max[i] = sliced_results.tR2[ii]
		A_max[i] = sliced_results.A[ii]
	end
	approx_max[!,:tR1] = tR1_max
	approx_max[!,:tR2] = tR2_max
	approx_max[!,:A] = A_max
	approx_max
end

# ╔═╡ ae9bb1ec-fc67-46de-91fe-6e3f47844ce0
begin
	t_ = 0.0:0.01:sum(GCxGC_TP.timesteps)
	c_ = GasChromatographySimulator.chromatogram(collect(t_), sol_final[1].tR, sol_final[1].τR)
	p_chrom_1 = Plots.plot(t_, c_)
end

# ╔═╡ 84111543-21e2-4891-9050-0c102cd46a52
begin
	chrom_sliced_1st = Array{Array{Float64}}(undef, length(sol_1st_TM[1].tR))
	height_1st = Array{Float64}(undef, length(sol_1st_TM[1].tR))
 	for i=1:length(sol_1st_TM[1].tR)
		chrom_sliced_1st[i] = GasChromatographySimulator.chromatogram(collect(t_), [sol_1st_TM[1].tR[i]], [newpar.col.L ./ sol_1st_TM[1].uR[i]]).*sol_1st_TM[1].A[i]
		height_1st[i] = (GasChromatographySimulator.chromatogram([sol_1st_TM[1].tR[i]], [sol_1st_TM[1].tR[i]], [newpar.col.L ./ sol_1st_TM[1].uR[i]])[1].*sol_1st_TM[1].A[i])
	end
	#Plots.plot(t_, chrom_sliced)
	height_1st
end

# ╔═╡ c1529873-e916-4c31-afc4-7d08fd1a18f1
begin
	chrom_sliced = Array{Array{Float64}}(undef, length(sol_final[1].tR))
	for i=1:length(sol_final[1].tR)
		chrom_sliced[i] = GasChromatographySimulator.chromatogram(collect(t_), [sol_final[1].tR[i]], [sol_final[1].τR[i]]).*sol_final[1].A[i]
	end
	#Plots.plot(t_, chrom_sliced)
end

# ╔═╡ b2b45284-bc7b-4c40-ad54-a78f9afc9ff5
begin
	gr()
	p_sliced_1D = Plots.plot(t_, chrom_sliced[1], legend=false)
	for i=2:length(chrom_sliced)
		Plots.plot!(p_sliced_1D, t_, chrom_sliced[i])
	end
	Plots.plot!(p_sliced_1D, xlims=(1880, 1930))
	p_sliced_1D
end

# ╔═╡ 85ecccaa-f38f-4d54-b63a-33abf5390dc5
begin
	chrom_sliced_sum = chrom_sliced[1]
	for i=2:length(chrom_sliced)
		chrom_sliced_sum = chrom_sliced_sum .+ chrom_sliced[i]
	end
	Plots.plot(t_, chrom_sliced_sum)
end

# ╔═╡ 6075d217-9072-46e6-9e59-47c8dad10a55


# ╔═╡ 6021015b-107d-41e0-b258-43219028e292


# ╔═╡ 5af9ebbe-b126-4ad4-9a44-c934d0d18cb3
md"""
### (deactivated) Result from new function
"""

# ╔═╡ f2522fe2-e5e4-44bc-8e5a-528b86dfc290
begin
	t__ = 0.0:0.001:sum(GCxGC_TP.timesteps)
	c_s = Array{Array{Float64}}(undef, length(sim_[2][1][end].tR))
	for i=1:length(sim_[2][1][end].tR)
		c_s[i] = GasChromatographySimulator.chromatogram(collect(t__), [sim_[2][1][end].tR[i]], [sim_[2][1][end].τR[i]])
	end
	c__ = sum(sim_[2][1][3].A .* c_s)
	p_chrom_1_ = Plots.plot(t__, c__)
end

# ╔═╡ 20d112f4-f5cc-4e6b-89c5-a0b3ad4ff4c5
sim_[2][1][3].A

# ╔═╡ 622a22d0-ecac-402f-96a5-79c6dc1683d3
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

# ╔═╡ 5a510b13-491d-48d3-84b1-1335aaf6cd06
begin
	c_slices, t_D1, t_D2 = chrom_slicing(t_, chrom_sliced_sum, PM)
	#p_chrom_2 = Plots.plot(t_D2[1], c_slices[1])
	#for i=2:length(c_slices)
	#	Plots.plot!(p_chrom_2, t_D2[i], c_slices[i].+i*10)
	#end
	#p_chrom_2
end

# ╔═╡ e077613b-6082-4d97-a1d1-bfd90614b00c
c_slices

# ╔═╡ ef63df5b-da93-44cb-809f-ba88c91230ce
t_D2

# ╔═╡ 8b490f44-8c18-44b8-aa37-c53b88000e07
length(c_slices), length(t_D2[1])

# ╔═╡ f1c65fa8-714d-48f7-a898-76e3feb61435
begin 
	lengths_c_slices = Array{Int}(undef, length(c_slices))
	for i=1:length(c_slices)
		lengths_c_slices[i] = length(c_slices[i])
	end
	unique(lengths_c_slices)
end

# ╔═╡ a7be08ba-b312-497d-90c1-f7253802276b
begin
	slice_mat = Array{Float64}(undef, length(c_slices)-1, length(t_D2[1]))
	for j=1:length(t_D2[1])
		for i=1:(length(c_slices)-1)
			slice_mat[i,j] = c_slices[i][j]
		end
	end
	slice_mat
end

# ╔═╡ 5a61f329-39cf-4554-88c1-3277832c0a26
Plots.heatmap(slice_mat', c=:jet1)

# ╔═╡ 4b7f4614-39c6-482d-a1ff-1da58b75a52f
Plots.contour(slice_mat')#, c=:jet1)

# ╔═╡ 8f1103d3-38ec-42a9-99fa-faa348f543c5
size(slice_mat)

# ╔═╡ 9cfdf1e7-b150-4680-a832-e03d76d4172a
length(t_D1)

# ╔═╡ 976f9f2b-2d1b-4f5f-8d3c-d8132bd7fd8e
length(t_D2[1])

# ╔═╡ 39a9ed84-c33b-4a0a-b489-4e0ac4fa5f39
t_D2

# ╔═╡ af6e107a-6a36-4a23-8a1b-9c58dff9d22f
begin
	gr()
	Plots.contour(t_D1[1:end-1], t_D2[1], slice_mat', levels=10, xlims=(1884, 1912), ylims=(0,4), xlabel="tR1 in s", ylabel="tR2 in s", size=(400,300), xticks=1884:4:1912)#, c=:jet1)
	#savefig("2D_chrom_3terpenes.svg")
end

# ╔═╡ 0cadf6e0-c0a7-4632-9d4f-e1754123b75a
itp_slice_mat = LinearInterpolation((t_D1[1:end-1], t_D2[1],), slice_mat, extrapolation_bc=Flat())

# ╔═╡ fd4026bf-7085-4a63-a074-c1bbf49bde2a
begin
	c_slices_, t_D1_, t_D2_ = chrom_slicing(t__, c__, 6.0)
	p_chrom_2_ = Plots.plot(t_D2_[1], c_slices_[1])
	for i=2:length(c_slices_)
		Plots.plot!(p_chrom_2_, t_D2_[i], c_slices_[i].+i*10)
	end
	p_chrom_2_
end

# ╔═╡ 5558e8e6-bb35-497b-9d61-ba0f3518e75e
md"""
## Simplified GCxGC
"""

# ╔═╡ 0ccf11db-c934-4592-95ae-d28bd856f97f
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

# ╔═╡ 21d1659a-e50f-4908-8b7b-bcee33f7acae
sys_simp = GCxGC_TM_simp(30.0, 0.25, 0.25, "ZB1ms", GCxGC_TP, 3.02, 0.1, 0.1, "Stabilwax", GCxGC_TP, 1.0, NaN, 0.0)

# ╔═╡ 9bb3302a-16c7-442c-b6a1-9ee1b9d20587
par_simp = graph_to_parameters(sys_simp, db, selected_solutes)

# ╔═╡ 04ce95ea-affb-45e1-b09c-dfd412a7291b
paths_simp = GasChromatographySystems.all_paths(sys_simp.g, 1)[2]

# ╔═╡ b360f18a-fb8a-4e66-b13f-67e14d008238
sim_simp = simulate_along_paths(sys_simp, paths_simp, par_simp)

# ╔═╡ 3e0f5a36-4db1-4ded-863b-2f0bfe01108a
pl_simp = GasChromatographySystems.peaklist_GCxGC(sim_simp[2][1][1], sim_simp[2][1][end])

# ╔═╡ 690988a3-e9f0-4c01-adb1-399ffbe2f093
pl = GasChromatographySystems.peaklist_GCxGC(sim[2][1][1], sim[2][1][end])

# ╔═╡ e41444ed-8012-4263-aa49-1e15399538a4
begin
	Plots.scatter(pl_simp.tR1, pl_simp.tR2.-pl_simp.tR1, xlabel="1st dim time in s", ylabel="2nd dim time in s")
	Plots.scatter!(pl_simp.tR2, pl_simp.tR2.-pl_simp.tR1)
end

# ╔═╡ 1f1fe7e1-d27d-485c-ad65-022b0862cc16
md"""
### Result from old function
"""

# ╔═╡ b5598366-a6d2-49b3-b71f-1aefe5b9ef26
p_gcxgc = GasChromatographySystems.plot_GCxGC(pl, sys)

# ╔═╡ 5c0103b4-5d18-4100-8fa6-a7e8fdbc7ef1
function plot_GCxGC(pl_GCxGC, sys, PM; categories = String[])
	t1end = sum(GasChromatographySystems.common_timesteps(sys))
	t2end = PM
	p_gcxgc = Plots.scatter(pl_GCxGC.tR1 .+ (pl.tR2.-pl.tR1) .- mod.(pl.tR2.-pl.tR1, PM), mod.(pl.tR2.-pl.tR1, PM), xlims=(0.0,t1end), ylims=(0.0, t2end), label="", xlabel=string(sys.modules[1].stationary_phase, ", t¹ in s"), ylabel=string(sys.modules[2].stationary_phase, ", t² in s"))
	if !isempty(categories)
		for i=1:length(categories) 
			pl_f = filter([:Cat] => x -> categories[i] in x, pl_GCxGC)
			Plots.scatter!(p_gcxgc, pl_f.tR1 .+ (pl_f.tR2.-pl_f.tR1) .- mod.(pl_f.tR2.-pl_f.tR1, PM), mod.(pl_f.tR2.-pl_f.tR1, PM), label=categories[i], legend=:bottomright)
		end
	end
	return p_gcxgc
end

# ╔═╡ f74033df-d0d0-4641-9109-f0ff7b07abd9
plot_GCxGC(pl, sys, 6.0)

# ╔═╡ 6e970b15-ff40-427c-9dd2-517e7ede00f8
function plot_GCxGC!(p_gcxgc, pl_GCxGC, sys; categories = String[])
	t1end = sum(GasChromatographySystems.common_timesteps(sys))
	t2end = ceil(maximum(pl_GCxGC.tR2.-pl_GCxGC.tR1))
	Plots.scatter!(p_gcxgc, pl_GCxGC.tR1, pl_GCxGC.tR2.-pl_GCxGC.tR1, xlims=(0.0,t1end), ylims=(0.0, t2end), label="", xlabel=string(sys.modules[1].stationary_phase, ", t¹ in s"), ylabel=string(sys.modules[2].stationary_phase, ", t² in s"))
	if !isempty(categories)
		for i=1:length(categories) 
			pl_f = filter([:Cat] => x -> categories[i] in x, pl_GCxGC)
			Plots.scatter!(p_gcxgc, pl_f.tR1, pl_f.tR2.-pl_f.tR1, label=categories[i], legend=:bottomright)
		end
	end
	return p_gcxgc
end

# ╔═╡ 9c3223ab-1dcc-40b6-8996-c8dbc22a2c6a
begin
	gr()
	_offset_D2 = 0.0
	p_sliced = Plots.scatter(slice_tR1, slice_tR2.+_offset_D2, ylims=(0.0, 6.0))
	Plots.scatter!(p_sliced, sim[2][1][1].tR, 2.0.*ones(length(sim[2][1][1].tR)))
	Plots.scatter!(p_sliced, approx_max.tR1, approx_max.tR2, c=:black, mshape=:cross)
	plot_GCxGC!(p_sliced, pl_simp, sys_simp)
end

# ╔═╡ 18cc22ce-0f7c-4683-ab2a-cae5a25ad9c4
begin
	slice_tR1_ = fld.(sim_[2][1][end].tR, PM).*PM
	slice_tR2_ = sim_[2][1][end].tR .- fld.(sim_[2][1][end].tR, PM).*PM
	p_sliced_ = Plots.scatter(slice_tR1_, slice_tR2_, ylims=(0.0, 6.0))
	Plots.scatter!(p_sliced_, sim_[2][1][1].tR, 2.0.*ones(length(sim_[2][1][1].tR)))
	plot_GCxGC!(p_sliced_, pl_simp, sys_simp)
end

# ╔═╡ 6ad69424-5cd5-4ac2-ae15-c3220b6ce449
plot_GCxGC!(p_gcxgc, pl_simp, sys_simp)

# ╔═╡ e5d6f0da-f249-4c09-b259-427fe240af34
mod.(pl.tR2.-pl.tR1, PM)

# ╔═╡ 2ca1ac9f-a4f7-48d3-ad4a-43dc01d97b18
(pl.tR2.-pl.tR1) .- mod.(pl.tR2.-pl.tR1, PM)

# ╔═╡ e2eaa21a-fbac-4d00-8663-c69ee4e0a8a9
md"""
## Measurement
"""

# ╔═╡ e9f8df08-7538-4dd4-8bbb-4a6649ae5b1a
begin
	meas = DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/Estimator/ZB1xWax/meas_GCxGC.csv", header=1, silencewarnings=true))

end

# ╔═╡ f55d7b20-2966-4a52-95fd-2f11a10d387a
function comparison_meas_sim(meas, pl_sim)
	comp = DataFrame(Name=pl_sim.Name, tR1_meas=pl_sim.tR1, tR1_sim=pl_sim.tR1, ΔtR1=pl_sim.tR1, tR2_meas=pl_sim.tR2.-pl_sim.tR1, tR2_sim=pl_sim.tR2.-pl_sim.tR1, ΔtR2=pl_sim.tR2.-pl_sim.tR1)
	for i=1:length(comp.Name)
		ii = findfirst(comp.Name[i].==meas.Name)
		if ismissing(meas.tR1[ii])
			comp[i, :tR1_meas] = NaN
			comp[i, :tR2_meas] = NaN
			comp[i, :ΔtR1] = NaN
			comp[i, :ΔtR2] = NaN
		else
			comp[i, :tR1_meas] = meas.tR1[ii]
			comp[i, :tR2_meas] = meas.tR2[ii]
			comp[i, :ΔtR1] = meas.tR1[ii] - comp[i, :tR1_sim]
			comp[i, :ΔtR2] = meas.tR2[ii] - comp[i, :tR2_sim]
		end
	end
	return comp
end

# ╔═╡ 6244a106-9527-4088-92d1-98ce746edd6f
comp = comparison_meas_sim(meas, pl_sim)

# ╔═╡ 301cf349-817d-41e2-bae4-4a4a285e9fd3
begin
	chrom = sort!(DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/Estimator/ZB1xWax/CSV/KetAlkPhenHR3Mod1.csv", header=1, silencewarnings=true)), :RT)

	#chrom = sort!(DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/Estimator/ZB1xWax/CSV/PhenoneHR3Mod1.csv", header=1, silencewarnings=true)), :RT)
end

# ╔═╡ 2aea4d11-3e6b-4c1e-9f0c-e6427e4257cf
names(chrom)

# ╔═╡ ffa85cd8-3fea-4aa6-8ea5-2ac37f0cce39
begin
	gr()
	Plots.plot(chrom.RT, chrom.TIC)
end

# ╔═╡ 27fd2299-ecf4-4ac0-a1d5-40d4a3d8d4f7
chrom.RT.*60

# ╔═╡ 7b3a6e3b-0700-49fb-b6c8-af58446e1d00
(60.01383-59.94718)*60

# ╔═╡ a65c3220-579d-4395-962e-67440ed5a7cd
c_slices_m, t_D1_m, t_D2_m = chrom_slicing(chrom.RT.*60, chrom.TIC, PM)

# ╔═╡ 8b2d1722-b25d-4083-9afd-b097574b7fd0
length(t_D1_m)

# ╔═╡ cbd5fccd-3ece-4b42-b911-5da849759490
begin
	itps = Array{Interpolations.Extrapolation}(undef, length(c_slices_m))
	for i=1:length(c_slices_m)
		itps[i] = LinearInterpolation((t_D2_m[i],), c_slices_m[i], extrapolation_bc=Flat())
	end
end

# ╔═╡ cb289770-50aa-44bc-b2f8-114be1924f82
itps[100]

# ╔═╡ 71868e87-b5f4-4c0c-ae5a-aa6c86e3331a
t_D2_m[100]

# ╔═╡ be644292-4b04-4ac5-ab50-acd7613e5c09
c_slices_m[100]

# ╔═╡ 6ccbc3c3-a54f-4266-8653-684ec0e398f5
begin
	Plots.plot(t_D2_m[100], c_slices_m[100])
	Plots.plot!(0.0:0.001:4.0, itps[100].(0.0:0.001:4.0))
end

# ╔═╡ f8624fbf-e1c7-4c03-9c66-130507a6dbd9
begin
	chrom_mat = Array{Float64}(undef, length(c_slices_m), length(0.0:0.001:4.0))
	for i=1:length(c_slices_m)
		chrom_mat[i,:] = itps[i].(0.0:0.001:4.0)
	end
	chrom_mat
end

# ╔═╡ c70c267f-7867-43ff-8810-b66d9219d7e9
size(chrom_mat)

# ╔═╡ ab513cdd-4f71-4860-a7d5-fd6744d68b1a
meas

# ╔═╡ 456081bc-dcb9-4819-9370-28d8194237a8
begin
	gr()
	chrom_mat_mod = Array{Float64}(undef, size(chrom_mat)[1], size(chrom_mat)[2])
	for j=1:size(chrom_mat)[2]
		for i=1:size(chrom_mat)[1]
			if chrom_mat[i,j] > 25000.0
				chrom_mat_mod[i,j] = 25000.0
			else
				chrom_mat_mod[i,j] = chrom_mat[i,j]
			end
		end
	end
	#gr()
	p_meas = Plots.heatmap(0.0:4.0:(size(chrom_mat)[1]-1)*4.0, 0.0:0.001:4.0, chrom_mat_mod', c=:jet1, colorbar=false);
	# add markers for the measured RTs
	# shift by 15*PM in the 1st dimension
	# why?
	# - measured retention times from different single measurement, where the retention times match the 2D chromatogram
	# - these retention times are shifted by 15 PM respectively to the 2D chromatogram of the measurement of all samples together
	offset_D1 = 15.0
	Plots.scatter!(meas.tR1[1:5].-offset_D1*PM, meas.tR2[1:5], markersize=4, c=:red, label="measured");# alcohols
	Plots.scatter!(meas.tR1[6:22].-offset_D1*PM, meas.tR2[6:22], markersize=4, c=:orange, label="measured") # terpenes
	Plots.scatter!(meas.tR1[23:28].-offset_D1*PM, meas.tR2[23:28], markersize=4, c=:yellow, label="measured"); # phenones
	Plots.scatter!(meas.tR1[29:35].-offset_D1*PM, meas.tR2[29:35], markersize=4, c=:lawngreen, label="measured"); # ketones

	# simulation
	offset_D2 = 0.0
	Plots.scatter!(approx_max.tR1.-offset_D1*PM, approx_max.tR2.+offset_D2, c=:lightblue, m=:diamond, markeralpha=1, msize=3, label="simulated", xlims=(0.0, 3500.0), ylims=(0.0, 4.0), xlabel="tR1 in s", ylabel="tR2 in s");
	#Plots.plot!(p_meas, xlims=(0.0, 3500.0), ylims=(0.0, 4.0), xlabel="tR1 in s", ylabel="tR2 in s")
	p_meas
end

# ╔═╡ 2e66e3a0-90fa-4c6a-9db6-7bbc7edc364a
#savefig(p_meas, "2D_chrom_terpenes.svg")

# ╔═╡ 8628510e-28a1-4823-800e-8bf12f563d74
meas

# ╔═╡ 947e6168-ba5d-4a4d-a780-96b9d7ac0760
0.0:4.0:size(chrom_mat)[1]*4.0

# ╔═╡ 847a7180-78cb-409a-8434-8e225e986367
length(0.0:4.0:(size(chrom_mat)[1]-1)*4.0)

# ╔═╡ 6aa00136-693c-417f-8d90-a5c1c5ffde79
begin
	#Plots.scatter!(p_meas, slice_tR1.-2*PM, slice_tR2)
	#p_meas
end

# ╔═╡ 26e07009-33e7-4198-99ae-3a7c5fce2172
#=begin
	plotly()
	Plots.heatmap(chrom_mat_mod', c=:jet1)
end=#

# ╔═╡ 37ca35d1-50b9-43f8-8ca1-0f388a3f391f
begin
	# tR2 of heptanol depending on L of GC2
	df_tR2_heptanol = DataFrame(L=[2, 1, 3, 0.8, 1.5, 1.9, 1.3, 1.7, 2.2, 1.1, 1.05], tR2=[1.55, 2.37, 1, 1.77, 3.91, 1.21, 3.28, 0.48, 2.21, 2.67, 2.5])
	df_tR2_heptanol_sort = sort(df_tR2_heptanol, :L)
	Plots.scatter(df_tR2_heptanol_sort.L, df_tR2_heptanol_sort.tR2)
	Plots.plot!(df_tR2_heptanol_sort.L, df_tR2_heptanol_sort.tR2)
end

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
# ╠═3a6b6c3c-c839-47ed-bf86-dd1f3149c7f9
# ╟─d0adda7c-8980-4d74-8b79-c3c5edd8131f
# ╠═ab33d65f-f0ba-4e5d-b660-2ecd106d8d21
# ╠═60d86bdd-5527-49da-8c35-ecd508171e6c
# ╠═b63ea500-c2bb-49c3-9915-032fd3e33e14
# ╠═7bb430cd-6bc5-4c48-85f8-3ffce3220d50
# ╠═3964e72a-4a96-4740-bd66-a224a0be16eb
# ╠═80150c58-f9b9-43d9-bd7b-e38388f29fd8
# ╠═522a7582-0ff1-42d6-960a-adb952a225d2
# ╠═ca604ac6-fdb1-41d7-83d2-a0c9b2440dad
# ╠═b9e3acb8-a7ec-40ef-ba58-932f79b08fb0
# ╠═95f57b5b-652c-40c2-b66d-137bb7b2d878
# ╠═9605b6ad-7cdc-47ee-9edf-c3ce531ae9a2
# ╠═06487ba0-8752-40a8-b5d8-07fb08514ffe
# ╠═a5c2b339-4879-4367-9d03-ada1e7e6dd46
# ╠═7fc154cd-d1a4-4f14-8750-e92a58da35f0
# ╠═daf5e89f-343c-4076-80fd-2dd52ba14f29
# ╠═051e9ad6-51de-40eb-a3cf-2b929d594400
# ╠═cf03f9b5-62e6-46af-924e-3a16638c768d
# ╠═607265de-c784-4ed0-8e89-ebb206785b0b
# ╠═2e4612e7-ce0f-4ae7-97b0-7076d9bf5651
# ╠═c366df65-59a4-427a-beee-bf0bef5fb409
# ╠═c1cb420c-7869-4da8-a240-1583545bde7c
# ╠═93f19f34-5766-4863-bfba-fde045afcdb9
# ╠═5428a416-48e4-4ebd-a493-54fd30cd2556
# ╠═6378cb76-57fd-4e59-a183-fd2c745d0f69
# ╠═bdbd8cc9-29fb-43ce-850b-8d18980808d3
# ╠═4adbf7e9-232d-4889-92c6-75e8bde0d88d
# ╠═fc36ffe3-b1df-42c8-ad14-e3b90f27a272
# ╠═43f27a81-dca0-4188-baf4-3f8140d9e57c
# ╠═d4bbf7ce-a804-469d-a648-b52908daae51
# ╠═7277a646-fd70-4032-9691-569eddeafec6
# ╠═793560f7-d364-4c68-81ec-994441a41059
# ╠═d40361fe-9e0c-4df8-a09c-0c5bf143dbf6
# ╠═5d63fd50-208e-4db1-9609-4f08c493d90c
# ╠═6993afed-4519-4c8a-9cfc-87c72723c444
# ╠═e1785bfc-dfc7-4ae4-b524-b93608845384
# ╠═b0469761-2537-47ec-b308-cc7be689ca36
# ╠═a14d1ca2-46ef-4296-aa1d-06590eca97cf
# ╠═69598506-cf92-4626-b6e1-d8f57bf8388f
# ╠═1a0695c5-0a67-4158-a3bd-c135be31e408
# ╠═6556a858-65ef-49a6-bdaa-562a2b568a70
# ╠═fcdc5015-31d0-4808-8218-f33901437c0b
# ╠═1f62b350-332e-4acf-b907-75362157e4da
# ╠═f1bf21dd-07ee-4e4e-beb5-7d9225b3fc6b
# ╠═6e8f2101-8245-483c-9c05-e4b16c260607
# ╠═bac5d026-0b84-412f-98bb-8ddedb2b92d9
# ╠═e3e34deb-9249-4811-b74a-44a4c2d06ac2
# ╠═ce317588-2f42-477e-a0a6-cc40047e1506
# ╠═54e20666-b217-415d-b2fd-305980588668
# ╠═ac3db64b-82c8-455a-97aa-6bbddcc03830
# ╠═ee2e3313-9d05-4815-97a2-5b292dc47b68
# ╠═7e088401-a8d4-4648-9380-241777534fe2
# ╠═82c69944-c4e0-4699-87c0-29500b33ba78
# ╟─ce958d12-3dad-4e3c-ab57-f84136610870
# ╠═70597ee8-760d-4c21-af6c-981922208d10
# ╠═b83330dd-4df3-433a-ab29-273e39f8b32c
# ╠═325a8e2e-ae2f-4dd9-a32b-00c7adf47461
# ╠═c626e4de-3c65-4230-bc23-2fa6ea9bbbd9
# ╠═03aeeaac-0a83-4d34-944f-d3457904d69d
# ╠═e4172084-2dd7-4daf-b2b9-cc84891c06f2
# ╠═f436afe9-0440-4149-be61-b426e243b34d
# ╠═846c3cab-90d4-4e73-9ce5-c1b41f937cdf
# ╠═f087274e-5147-48f9-8200-4d6ed90397c9
# ╠═90d32247-3e00-41ee-bd2e-b8e698c2613c
# ╠═d98548cc-30bd-48c5-b408-0227caaa5fd1
# ╠═402e9e25-3d5b-46c2-8b5e-4a472fb357f7
# ╠═9c5ec167-212d-413a-a2bc-fb00d015af51
# ╠═06c3c74b-4ced-405e-89be-a05092fab027
# ╠═f48d790f-6ef0-4e31-a80d-8acaae8e454f
# ╠═7708f828-0ea3-4568-8614-ff0a724b9eea
# ╠═67a4dff7-2634-4ae2-adb0-515807988b17
# ╠═1e057e8c-66f3-43d1-9674-6d030b094a74
# ╠═84111543-21e2-4891-9050-0c102cd46a52
# ╠═6fa14ed4-f9ee-4279-a020-6abe60f857d7
# ╠═9c6573e8-bc0f-4c40-9eea-9733726fbd76
# ╠═d5b4ac08-43e0-4d2c-ab7c-9fea2a926ff3
# ╠═095ee9ba-d974-4625-a47c-3634d91c6a5d
# ╠═471878ef-b115-411c-9d9d-764d6db08069
# ╠═f6093dce-30a7-4628-8362-0c5ed8c55ebc
# ╠═beb15214-b1bb-4b39-a2aa-86f0727bc2c4
# ╠═2c1fc71c-5f5d-4e18-84f5-42980d461753
# ╠═6a3f4043-63ba-44ef-b9ce-a391bc697fe6
# ╠═74121702-9f4a-410b-961b-7865ab3af941
# ╠═975b24c5-6557-4d7c-a3be-4d6784fcf66d
# ╠═5d0cb969-8a98-4d28-aa5c-233723755f41
# ╠═a848d367-8a2e-4eb1-9fa9-bb94e83d34ff
# ╠═c5f01fec-7bb6-4157-a125-f7c929dd69e6
# ╠═eefc168f-7de7-45fa-bd25-ea3e7df4b685
# ╠═c2b343a8-f794-42f3-a237-477f00d2d6b3
# ╠═96200092-8afd-4ce3-8287-79806a331656
# ╠═e662f5a1-79ee-4adb-8d98-ff0b743c95da
# ╠═0baa9993-b50b-45ce-a1b7-ed6e45826187
# ╠═86274690-dfe8-4f71-8084-66fe5f03f932
# ╠═e1ce8764-cec4-4c43-9354-ccceea375f3f
# ╠═12c5ce26-4a91-4b1b-b86e-0465c1c155a3
# ╠═ece4ada0-9944-40e1-b90d-8c1cd34209f6
# ╠═31c843d2-e4ad-47b0-b90e-409b9bdb0f4c
# ╠═4c85cfd1-cb1b-449d-a101-e537abe25e86
# ╠═dca2c405-93dc-460c-8a06-9892154d503a
# ╠═6010a66c-d23c-4127-afae-62b6fddcb541
# ╠═231edaa6-8ef3-4a85-b48e-e0e1829d124d
# ╠═64e635f8-f5fd-4935-b44e-b02ebdd96b64
# ╠═f388ca37-0167-4459-858b-619f1d2ce2d7
# ╠═bb01e67d-fbca-461e-b592-508bc5480b94
# ╠═c959c1ae-2c46-412d-b23b-fadb87d4d3be
# ╠═088c0c0f-0722-4498-8306-8286069fc24b
# ╠═29cbe2b2-34fb-45dc-849b-2f678b3f8797
# ╠═14bf7364-ce16-467e-894b-43d5ad1d91dc
# ╠═98713f9d-123a-4eda-8a8c-72c6ada1a3a4
# ╠═caef921f-f06c-471b-a56a-8415ac2d4fa6
# ╠═05a0e4f4-8a1b-4f3c-84b0-c2491840d20c
# ╠═1bc4ae12-6f10-499a-b20f-35d5bd297378
# ╠═935485a4-947e-42ad-ae92-9146c3f551ff
# ╠═16e6bf4c-05ee-466d-95d6-e679e4410bc1
# ╠═6accd799-2058-4756-9e31-c0b8b6bc6130
# ╠═32e40a19-ebab-47a6-b4ca-9300c2b04dab
# ╠═93ee192a-7374-43f4-a1a1-a1265fda42a4
# ╠═f1504c87-b24f-4550-9a89-c7759b635218
# ╠═8b500dc4-98de-46af-9a83-f5ff7e8190ee
# ╠═d0f14e55-10ad-4ddb-a743-54a3ca28aa85
# ╠═95897aaf-d2bf-4edd-9cf3-8bf7353fd089
# ╠═b0490a91-8a34-4fee-b819-de4f6ecb7a27
# ╠═9d6944b7-7794-4912-b8ab-ce91672cf0de
# ╠═a0b9918a-1ea2-463a-8587-2b4ca0f1230c
# ╠═6c4c9478-1392-442c-9a26-b654b64c1b14
# ╠═50792b68-8acf-4349-afd5-66267fa9153b
# ╠═08d9257e-6d49-4cf6-97a9-f8b546d9e933
# ╠═df1c6430-8ea8-4bdd-a9be-078a0357fe05
# ╠═c2074cee-177c-4ece-bf67-0262e504536a
# ╠═78a325ce-a934-4afb-bb31-e1eb98e5349a
# ╠═a76b6b3e-a41c-413e-9b4d-920823e9edda
# ╠═ac2e34ed-a519-48fc-a87e-fb8dde6f3ff5
# ╠═915d8972-5423-4dde-8039-9a912bd0ae61
# ╠═556c157c-321d-4eb8-b2ef-ee8c4fe0f038
# ╠═a5eb5c69-de38-4dbb-83f7-4b751d26562d
# ╠═c78ea8fb-d145-4862-973b-06bd13c03f52
# ╠═80a012dd-738e-4332-857a-90203867a03f
# ╠═9c3223ab-1dcc-40b6-8996-c8dbc22a2c6a
# ╠═ae9bb1ec-fc67-46de-91fe-6e3f47844ce0
# ╠═c1529873-e916-4c31-afc4-7d08fd1a18f1
# ╠═b2b45284-bc7b-4c40-ad54-a78f9afc9ff5
# ╠═85ecccaa-f38f-4d54-b63a-33abf5390dc5
# ╠═5a510b13-491d-48d3-84b1-1335aaf6cd06
# ╠═e077613b-6082-4d97-a1d1-bfd90614b00c
# ╠═ef63df5b-da93-44cb-809f-ba88c91230ce
# ╠═8b490f44-8c18-44b8-aa37-c53b88000e07
# ╠═f1c65fa8-714d-48f7-a898-76e3feb61435
# ╠═a7be08ba-b312-497d-90c1-f7253802276b
# ╠═5a61f329-39cf-4554-88c1-3277832c0a26
# ╠═4b7f4614-39c6-482d-a1ff-1da58b75a52f
# ╠═8f1103d3-38ec-42a9-99fa-faa348f543c5
# ╠═9cfdf1e7-b150-4680-a832-e03d76d4172a
# ╠═976f9f2b-2d1b-4f5f-8d3c-d8132bd7fd8e
# ╠═39a9ed84-c33b-4a0a-b489-4e0ac4fa5f39
# ╠═af6e107a-6a36-4a23-8a1b-9c58dff9d22f
# ╠═6075d217-9072-46e6-9e59-47c8dad10a55
# ╠═0cadf6e0-c0a7-4632-9d4f-e1754123b75a
# ╠═6021015b-107d-41e0-b258-43219028e292
# ╠═5af9ebbe-b126-4ad4-9a44-c934d0d18cb3
# ╠═18cc22ce-0f7c-4683-ab2a-cae5a25ad9c4
# ╠═f2522fe2-e5e4-44bc-8e5a-528b86dfc290
# ╠═20d112f4-f5cc-4e6b-89c5-a0b3ad4ff4c5
# ╠═fd4026bf-7085-4a63-a074-c1bbf49bde2a
# ╠═622a22d0-ecac-402f-96a5-79c6dc1683d3
# ╠═5558e8e6-bb35-497b-9d61-ba0f3518e75e
# ╠═21d1659a-e50f-4908-8b7b-bcee33f7acae
# ╠═0ccf11db-c934-4592-95ae-d28bd856f97f
# ╠═9bb3302a-16c7-442c-b6a1-9ee1b9d20587
# ╠═04ce95ea-affb-45e1-b09c-dfd412a7291b
# ╠═b360f18a-fb8a-4e66-b13f-67e14d008238
# ╠═3e0f5a36-4db1-4ded-863b-2f0bfe01108a
# ╠═690988a3-e9f0-4c01-adb1-399ffbe2f093
# ╠═e41444ed-8012-4263-aa49-1e15399538a4
# ╠═1f1fe7e1-d27d-485c-ad65-022b0862cc16
# ╠═b5598366-a6d2-49b3-b71f-1aefe5b9ef26
# ╠═6ad69424-5cd5-4ac2-ae15-c3220b6ce449
# ╠═5c0103b4-5d18-4100-8fa6-a7e8fdbc7ef1
# ╠═f74033df-d0d0-4641-9109-f0ff7b07abd9
# ╠═6e970b15-ff40-427c-9dd2-517e7ede00f8
# ╠═e5d6f0da-f249-4c09-b259-427fe240af34
# ╠═2ca1ac9f-a4f7-48d3-ad4a-43dc01d97b18
# ╠═e2eaa21a-fbac-4d00-8663-c69ee4e0a8a9
# ╠═e9f8df08-7538-4dd4-8bbb-4a6649ae5b1a
# ╠═f55d7b20-2966-4a52-95fd-2f11a10d387a
# ╠═6244a106-9527-4088-92d1-98ce746edd6f
# ╠═301cf349-817d-41e2-bae4-4a4a285e9fd3
# ╠═2aea4d11-3e6b-4c1e-9f0c-e6427e4257cf
# ╠═ffa85cd8-3fea-4aa6-8ea5-2ac37f0cce39
# ╠═27fd2299-ecf4-4ac0-a1d5-40d4a3d8d4f7
# ╠═7b3a6e3b-0700-49fb-b6c8-af58446e1d00
# ╠═a65c3220-579d-4395-962e-67440ed5a7cd
# ╠═8b2d1722-b25d-4083-9afd-b097574b7fd0
# ╠═4dd08d9d-7068-4198-b5d4-5377a787dc8a
# ╠═cbd5fccd-3ece-4b42-b911-5da849759490
# ╠═cb289770-50aa-44bc-b2f8-114be1924f82
# ╠═71868e87-b5f4-4c0c-ae5a-aa6c86e3331a
# ╠═be644292-4b04-4ac5-ab50-acd7613e5c09
# ╠═6ccbc3c3-a54f-4266-8653-684ec0e398f5
# ╠═f8624fbf-e1c7-4c03-9c66-130507a6dbd9
# ╠═c70c267f-7867-43ff-8810-b66d9219d7e9
# ╠═ab513cdd-4f71-4860-a7d5-fd6744d68b1a
# ╠═456081bc-dcb9-4819-9370-28d8194237a8
# ╠═2e66e3a0-90fa-4c6a-9db6-7bbc7edc364a
# ╠═8628510e-28a1-4823-800e-8bf12f563d74
# ╠═947e6168-ba5d-4a4d-a780-96b9d7ac0760
# ╠═847a7180-78cb-409a-8434-8e225e986367
# ╠═6aa00136-693c-417f-8d90-a5c1c5ffde79
# ╠═26e07009-33e7-4198-99ae-3a7c5fce2172
# ╠═37ca35d1-50b9-43f8-8ca1-0f388a3f391f
# ╠═17c4fbd1-4a9d-4f7a-922a-a72eb97db33a
