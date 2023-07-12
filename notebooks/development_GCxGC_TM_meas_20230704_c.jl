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

# ╔═╡ 95cc1f3d-eff3-4cc8-8ba1-7ee04e887b84
using SpecialFunctions

# ╔═╡ 4a20d35c-2156-4e5d-a94c-dc37df64545a
using OrdinaryDiffEq

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
## Periodic smoothed rectangle temperature function
"""

# ╔═╡ 326d2249-dfcc-40b6-b8c0-d9ab5ff88627
#=function therm_mod(t, shift, PM, ratio, Tcold, Thot) 
	return ifelse(mod(t+shift, PM) < ratio*PM, Tcold, Thot)
end=#

# ╔═╡ 9a6dea14-f8d4-4b79-9c02-268c0c799515
# differentiation of this function is zero everywehre

# ╔═╡ e44bcb5f-353e-4ef1-bf4d-d226a864ad48
function smooth_rectangle(x, a, b, m)
	# a ... mid-position of rising flank
	# b ... mid-position of falling flank
	# m ... width of the flank
	if x.>a-6*m && x.<=a+6*m
		f = 1/2 .*(1 .+erf.((x.-a)./sqrt.(2 .*m.^2)))
	elseif x.>a+6*m && x.<=b-6*m
		f = 1
	elseif x.>b-6*m && x<=b+6*m
		f = 1 .- 1/2 .*(1 .+erf.((x.-b)./sqrt.(2 .*m.^2)))
	else
		f = 0
	end
	return f
end

# ╔═╡ 7ccb0031-e7e0-43e9-952e-4f80501716ec
function smooth_rectangle_(x, xstart, width, min, max; flank=20)
	a = xstart #- Δx₁/2
	b = xstart + width#/2
	m = width/flank
	val = (max-min)*smooth_rectangle(x, a, b, m) + min
	return val
end

# ╔═╡ 9efd0ffb-99f7-4500-9ea2-762b91c82448
x₁ = 10.0

# ╔═╡ 5cf107e3-be39-44ae-a7a4-f9463828376a
Δx₁ = 1.0

# ╔═╡ 54058244-e33a-42b6-badb-72cd3fc9696e
begin
	plotly()
	x = 0.0:0.001:20.0
	a = x₁ #- Δx₁/2
	b = x₁ + Δx₁#/2
	m = Δx₁/200

	Plots.plot(x, smooth_rectangle_.(x, x₁, Δx₁, 0.0, 1.0))
	Plots.plot!(x, smooth_rectangle_.(x, x₁-4, 2*Δx₁, 0.0, 1.0))
	Plots.plot!(x, smooth_rectangle_.(x, 4-0.35, 0.35, -10.0, 20.0))
end

# ╔═╡ ab143f29-0a96-4ae2-a178-6bf773e6cd26
a, a-6*m, a+6*m

# ╔═╡ 1a081250-c6a6-42a2-9839-f2fd6c6ea0a2
x₁, a+6*m, b-6*m

# ╔═╡ 0c1e27ab-3b5e-4f42-8352-73756042efa4
b, b-6*m, b+6*m

# ╔═╡ 196c80bb-55c9-4e71-8776-b99397930a4e
6*m

# ╔═╡ 06fa229a-d4e5-4662-a8b9-33c3874ac3f0
function therm_mod(t, shift, PM, ratio, Tcold, Thot) 
	width = (1-ratio)*PM
	tmod = mod(t+shift-width/2, PM)
	tstart = ratio*PM
	#return smooth_rectangle_.(tmod+width/2, tstart, width, Tcold, Thot)
	return smooth_rectangle_.(tmod+2/3*width, tstart, width, Tcold, Thot)
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

# ╔═╡ 8d4391c3-3712-4aeb-9868-3d1d0fd09745
function periodic_smooth_rectangle_(x, xstart, width, period, shift, Tcold, Thot)
	#xmod = mod(x+shift-width/2, period)
	xmod = mod(x+shift-width/2, period)
	#val = smooth_rectangle_.(xmod+width/2, xstart, width, Tcold, Thot)
	val = smooth_rectangle_.(xmod+4*width/6, xstart, width, Tcold, Thot)
	return val
end

# ╔═╡ 1f4ca07c-bac7-4d16-8570-e1b11222640d
_Tcold, _Thot

# ╔═╡ 402ed0ed-2126-47fd-b044-9542f997a050
80*80

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
		ratio::Float64 # a number, ratio of the duration between cold and hot jet, approx. as rectangular function
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
	modules[5] = ModuleTM("TM2", LM[4], dM*1e-3, dfM*1e-6, spM, TPM, shiftM, PM, ratioM, HotM, ColdM, NaN)
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
	
#	sys = GCxGC_TM(30.0, 0.25, 0.25, "ZB1ms", GCxGC_TP, 0.1, 0.1, 0.1, "Stabilwax", GCxGC_TP, 2.0, 0.1, 0.1, "Stabilwax", 280.0, [0.30, 0.01, 0.90, 0.01, 0.01], 0.1, 0.1, "Stabilwax", 0.0, PM, coldjet/PM, 80.0, -120.0, GCxGC_TP, 0.8, NaN, 0.0)
	# old:
	sys = GCxGC_TM(31.5, 0.25, 0.25, "ZB1ms", GCxGC_TP, 0.1, 0.1, 0.1, "Stabilwax", GCxGC_TP, 0.25, 0.1, 0.1, "Stabilwax", 280.0, [0.60, 0.01, 0.30, 0.01, 0.50], 0.1, 0.1, "Stabilwax", 0.0, PM, coldjet/PM, 80.0, -120.0, GCxGC_TP, 0.65, NaN, 0.0)
end

# ╔═╡ f6ce117e-e382-439d-a943-862282714b90
begin
	plotly()
	#_t = 0.0:0.001:20.0
	#_shift = 0.0
	#_PM = 4.0
	#_ratio = (_PM-0.35)/_PM
	#_Thot = +80.0
	#_Tcold = -80.0
	_T_ = Array{Float64}(undef, length(_t))
	for i=1:length(_t)
		_T_[i] = periodic_smooth_rectangle_(_t[i], _ratio*_PM, (1-_ratio)*PM, _PM, _shift.+1, _Tcold, _Thot)
	end
	Plots.plot(_t, _T_)
	Plots.plot!(_t, periodic_smooth_rectangle_.(_t, _ratio*_PM, (1-_ratio)*PM, _PM, _shift, 0.0, 30))
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

	function module_temperature(module_::ModuleTM, sys; Tcold_abs=true)
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
par = graph_to_parameters(sys, db, selected_solutes; Tcold_abs=false)

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
sim = GasChromatographySystems.simulate_along_paths(sys, paths, par)

# ╔═╡ 6556a858-65ef-49a6-bdaa-562a2b568a70
md"""
### (deactivated) Simulation using the new function

The new function will reproduce the step-by-step simulation.
"""

# ╔═╡ fcdc5015-31d0-4808-8218-f33901437c0b
#sim_ = simulate_along_paths(sys, paths, par, abstol=1e-6, reltol=1e-3)

# ╔═╡ 1f62b350-332e-4acf-b907-75362157e4da
sim_[2][1][3]

# ╔═╡ f1bf21dd-07ee-4e4e-beb5-7d9225b3fc6b
sim_[2][1][5]

# ╔═╡ ceb6ba3f-aece-4ee3-b757-5904fa04c7ac
sim_[2][1][8]

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
	sim[2][1][2][!,:τR] = 1.5.*sim[2][1][2][!,:τR]
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
"""

# ╔═╡ 325a8e2e-ae2f-4dd9-a32b-00c7adf47461
begin
	# if next module == moduleTM
	par_prev = par[2]
	# initial peak width?
	τ0 = 0.0 # or PM?
#	init_t_start = (fld.(sim[2][1][2].tR.-6*sim[2][1][2].τR, PM).-1).*PM # is the -1 here correct?????
#	init_t_end = (fld.(sim[2][1][2].tR.+6*sim[2][1][2].τR, PM).-1).*PM 
	nτ = 6
	shift_ = 4.0
	init_t_start = (fld.(sim[2][1][2].tR.-nτ.*sim[2][1][2].τR, sys.modules[3].PM)).*sys.modules[3].PM .+ sys.modules[3].shift # start time of the peaks
	init_t_end = (fld.(sim[2][1][2].tR.+nτ.*sim[2][1][2].τR, sys.modules[3].PM)).*sys.modules[3].PM .+ sys.modules[3].shift# end time of the peaks
	n_slice = Int.((init_t_end .- init_t_start)./sys.modules[3].PM .+ 1) # number of slices for every substance
	sub_TM = Array{GasChromatographySimulator.Substance}(undef, sum(n_slice))
	sub_unsliced = Array{GasChromatographySimulator.Substance}(undef, length(n_slice))
	ii = 1
	for i=1:length(n_slice)
		for j=1:n_slice[i]
			sub_TM[ii] = GasChromatographySimulator.Substance(par_prev.sub[i].name, par_prev.sub[i].CAS, par_prev.sub[i].Tchar, par_prev.sub[i].θchar, par_prev.sub[i].ΔCp, par_prev.sub[i].φ₀, "slice$(j), "*par_prev.sub[i].ann, par_prev.sub[i].Cag, init_t_start[i]+(j-1)*sys.modules[3].PM, τ0)
			ii = ii + 1
		end
		sub_unsliced[i] =  GasChromatographySimulator.Substance(par_prev.sub[i].name, par_prev.sub[i].CAS, par_prev.sub[i].Tchar, par_prev.sub[i].θchar, par_prev.sub[i].ΔCp, par_prev.sub[i].φ₀, par_prev.sub[i].ann, par_prev.sub[i].Cag, sim[2][1][2].tR[i], sim[2][1][2].τR[i])
	end
	newopt = GasChromatographySimulator.Options(par[3].opt.alg, 1e-18, 1e-14, par[3].opt.Tcontrol, par[3].opt.odesys, par[3].opt.ng, par[3].opt.vis, par[3].opt.control, par[3].opt.k_th)
	#newopt = GasChromatographySimulator.Options(BS3(), 1e-6, 1e-3, par[3].opt.Tcontrol, par[3].opt.odesys, par[3].opt.ng, par[3].opt.vis, par[3].opt.control, par[3].opt.k_th)
	newpar = GasChromatographySimulator.Parameters(par[3].col, par[3].prog, sub_TM, newopt)
	newpar_unsliced = GasChromatographySimulator.Parameters(par[3].col, par[3].prog, sub_unsliced, newopt)

	newpar
end

# ╔═╡ 3388fdf6-dff9-4d44-9db3-5bf77a57980d
# this splicing function is with separating focussed and unfocussed peak segments
function slicing(tR, τR, PM, ratio, shift, par::GasChromatographySimulator.Parameters, prev_par::GasChromatographySimulator.Parameters; nτ=6, τ₀=zeros(length(τR)), abstol=1e-8, reltol=1e-6, alg=OwrenZen5())
	tcold = PM*ratio
	thot = PM*(1-ratio)
	init_t_start = (fld.(tR.-nτ.*τR, PM)).*PM .+ shift # start time of the peaks, rounded to multiple of PM
	init_t_end = (fld.(tR.+nτ.*τR, PM)).*PM .+ shift # end time of the peaks, rounded to multiple of PM
	n_slice = Int.((init_t_end .- init_t_start)./PM .+ 1) # number of slices for every substance
	sub_TM_focussed = Array{GasChromatographySimulator.Substance}(undef, sum(n_slice))
	sub_TM_unfocussed = Array{GasChromatographySimulator.Substance}(undef, sum(n_slice))
	A_focussed = Array{Float64}(undef, sum(n_slice))
	A_unfocussed = Array{Float64}(undef, sum(n_slice))
	g(x,p) = 1/sqrt(2*π*p[2]^2)*exp(-(x-p[1])^2/(2*p[2]^2))
	ii = 1
	Name = Array{String}(undef, sum(n_slice))
	CAS = Array{String}(undef, sum(n_slice))
	Ann_focussed = Array{String}(undef, sum(n_slice))
	Ann_unfocussed = Array{String}(undef, sum(n_slice))
	t0_foc = Array{Float64}(undef, sum(n_slice))
	t0_unfoc = Array{Float64}(undef, sum(n_slice))
	for i=1:length(n_slice)
		for j=1:n_slice[i]
			t₀ = init_t_start[i]+(j-1)*PM # initial start time
			
			sub_TM_focussed[ii] = GasChromatographySimulator.Substance(prev_par.sub[i].name, prev_par.sub[i].CAS, prev_par.sub[i].Tchar, prev_par.sub[i].θchar, prev_par.sub[i].ΔCp, prev_par.sub[i].φ₀, "TM1_f$(j), "*prev_par.sub[i].ann, prev_par.sub[i].Cag, t₀, τ₀[i])

			sub_TM_unfocussed[ii] = GasChromatographySimulator.Substance(prev_par.sub[i].name, prev_par.sub[i].CAS, prev_par.sub[i].Tchar, prev_par.sub[i].θchar, prev_par.sub[i].ΔCp, prev_par.sub[i].φ₀, "TM1_uf$(j), "*prev_par.sub[i].ann, prev_par.sub[i].Cag, t₀+tcold, thot)
			# Integrals:
			p = [tR[i], τR[i]]
		
			prob_focussed = IntegralProblem(g, t₀, t₀+tcold, p)
			prob_unfocussed = IntegralProblem(g, t₀+tcold, t₀+tcold+thot, p)
			A_focussed[ii] = solve(prob_focussed, QuadGKJL()).u
			A_unfocussed[ii] = solve(prob_unfocussed, QuadGKJL()).u
			# Areas in the same order as sub_TM_focussed rep. sub_TM_unfocussed
			Name[ii] = sub_TM_focussed[ii].name
			CAS[ii] = sub_TM_focussed[ii].CAS
			Ann_focussed[ii] = sub_TM_focussed[ii].ann
			Ann_unfocussed[ii] = sub_TM_unfocussed[ii].ann
			t0_foc[ii] = t₀
			t0_unfoc[ii] = t₀ + tcold
			ii = ii + 1
		end
	end
	#alg = OrdinaryDiffEq.ROCK4()#par.opt.alg
	newopt = GasChromatographySimulator.Options(alg, abstol, reltol, par.opt.Tcontrol, par.opt.odesys, par.opt.ng, par.opt.vis, par.opt.control, par.opt.k_th)
	newpar_focussed = GasChromatographySimulator.Parameters(par.col, par.prog, sub_TM_focussed, newopt)
	newpar_unfocussed = GasChromatographySimulator.Parameters(par.col, par.prog, sub_TM_unfocussed, newopt)

	df_A_foc = DataFrame(Name=Name, CAS=CAS, Annotations=Ann_focussed, A=A_focussed, t0=t0_foc)
	df_A_unfoc = DataFrame(Name=Name, CAS=CAS, Annotations=Ann_unfocussed, A=A_unfocussed, t0=t0_unfoc)
	
	return newpar_focussed, newpar_unfocussed, df_A_foc, df_A_unfoc
end

# ╔═╡ e4172084-2dd7-4daf-b2b9-cc84891c06f2
par_slice_1st_TM, par_slice_unfoc_1st_TM, A_slice_1st_TM, A_slice_unfoc_1st_TM = slicing(sim[2][1][2].tR, sim[2][1][2].τR, sys.modules[3].PM, sys.modules[3].ratio, sys.modules[3].shift, par[3], par[2]; nτ=6, abstol=1e-10, reltol=1e-8, alg=Vern9()) # competly slicing up the peaks at the modulator

# ╔═╡ 65c2d80b-152c-411d-b4e2-eff77e29acba
md"""
Smoothed rectangle function leads for once to a more accurate simulation of peak width, but it is more costly to simulate and the rate of failing seems higher.
"""

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

# ╔═╡ 4418401d-c418-422c-a341-7f01e7a8a3dd
function add_t0_to_pl!(pl, df_t0)
	# add the initial time in the dataframe df_t0 to the peaklist, same CAS and Annotations are used to identify the correct solute
	sort_t0 = Array{Float64}(undef, length(df_t0.t0))
	for i=1:length(df_t0.t0)
		ii = common_index(pl, df_t0.CAS[i], join(split(df_t0.Annotations[i], ", ")[1:end-1], ", "))
		#ii = common_index(pl, df_A.CAS[i], split(df_A.Annotations[i], ", ")[1])
		sort_t0[ii] = df_t0.t0[i]
	end
	pl[!, :t0] = sort_t0
	return pl
end

# ╔═╡ f436afe9-0440-4149-be61-b426e243b34d
begin
	sol_1st_TM = GasChromatographySimulator.simulate(par_slice_1st_TM)
	
	add_A_to_pl!(sol_1st_TM[1], A_slice_1st_TM)
	add_t0_to_pl!(sol_1st_TM[1], A_slice_1st_TM)
end

# ╔═╡ 629d4ad1-9961-4d2c-bb41-3c3832bc20fd
sol_1st_TM[1].tR .- sol_1st_TM[1].t0

# ╔═╡ 91cf4b42-140d-490f-949a-b655c2254c83
# indices, where modulation period is jumped
findall(sol_1st_TM[1].tR .- sol_1st_TM[1].t0 .> sys.modules[3].PM)

# ╔═╡ ca0bbf4a-5cef-4827-b040-cdad9d1dacb8
md"""
#### Test different solvers
"""

# ╔═╡ 7a81f5de-a2bc-4828-b5df-92c9b7950991
md"""
- OwrenZen3
"""

# ╔═╡ 9c30bd27-c649-4ed9-9660-23de44b8ba06
begin
	par_slice_OwrenZen3 = slicing(sim[2][1][2].tR, sim[2][1][2].τR, sys.modules[3].PM, sys.modules[3].ratio, sys.modules[3].shift, par[3], par[2]; nτ=6, abstol=1e-10, reltol=1e-8, alg=OwrenZen3()) # competly slicing up the peaks at the modulator
	sol_OwrenZen3 = GasChromatographySimulator.simulate(par_slice_OwrenZen3[1])
	
	add_A_to_pl!(sol_OwrenZen3[1], par_slice_OwrenZen3[3])
	add_t0_to_pl!(sol_OwrenZen3[1], par_slice_OwrenZen3[3])
end

# ╔═╡ ccd40fbd-ee6b-4571-a4ad-d2bc6d07a987
# indices, where modulation period is jumped
findall(sol_OwrenZen3[1].tR .- sol_OwrenZen3[1].t0 .> sys.modules[3].PM)

# ╔═╡ 5ae945d5-831a-46b4-b0bb-c461b3301ace
md"""
- OwrenZen4
"""

# ╔═╡ ed8d11a2-8583-4e42-8576-e3bdcb29b8d5
begin
	par_slice_OwrenZen4 = slicing(sim[2][1][2].tR, sim[2][1][2].τR, sys.modules[3].PM, sys.modules[3].ratio, sys.modules[3].shift, par[3], par[2]; nτ=6, abstol=1e-10, reltol=1e-8, alg=OwrenZen4()) # competly slicing up the peaks at the modulator
	sol_OwrenZen4 = GasChromatographySimulator.simulate(par_slice_OwrenZen4[1])
	
	add_A_to_pl!(sol_OwrenZen4[1], par_slice_OwrenZen4[3])
	add_t0_to_pl!(sol_OwrenZen4[1], par_slice_OwrenZen4[3])
end

# ╔═╡ 4f9e4d5e-be9a-432d-8032-0176033975f6
# indices, where modulation period is jumped
findall(sol_OwrenZen4[1].tR .- sol_OwrenZen4[1].t0 .> sys.modules[3].PM)

# ╔═╡ 0fac473b-42c1-41ca-a673-4e3e48e825d8
md"""
- OwrenZen5
"""

# ╔═╡ d8949e4a-77aa-485e-8911-e9faf7bbffc5
begin
	par_slice_OwrenZen5 = slicing(sim[2][1][2].tR, sim[2][1][2].τR, sys.modules[3].PM, sys.modules[3].ratio, sys.modules[3].shift, par[3], par[2]; nτ=6, abstol=1e-10, reltol=1e-8, alg=OwrenZen5()) # competly slicing up the peaks at the modulator
	sol_OwrenZen5 = GasChromatographySimulator.simulate(par_slice_OwrenZen5[1])
	
	add_A_to_pl!(sol_OwrenZen5[1], par_slice_OwrenZen5[3])
	add_t0_to_pl!(sol_OwrenZen5[1], par_slice_OwrenZen5[3])
end

# ╔═╡ 40bd41a7-cf33-4da1-9869-332c16ba597c
# indices, where modulation period is jumped
findall(sol_OwrenZen5[1].tR .- sol_OwrenZen5[1].t0 .> sys.modules[3].PM)

# ╔═╡ 406b56cb-a215-4e89-a1ab-37233156ad9a
md"""
- RadauIIA3
"""

# ╔═╡ c34d2c93-9a01-4674-91f4-37cfccb399e4
begin
	par_slice_RadauIIA3 = slicing(sim[2][1][2].tR, sim[2][1][2].τR, sys.modules[3].PM, sys.modules[3].ratio, sys.modules[3].shift, par[3], par[2]; nτ=6, abstol=1e-10, reltol=1e-8, alg=RadauIIA3()) # competly slicing up the peaks at the modulator
	sol_RadauIIA3 = GasChromatographySimulator.simulate(par_slice_RadauIIA3[1])
	
	add_A_to_pl!(sol_RadauIIA3[1], par_slice_RadauIIA3[3])
	add_t0_to_pl!(sol_RadauIIA3[1], par_slice_RadauIIA3[3])
end

# ╔═╡ bca11846-b79f-48b6-86e0-3fdb1b4c1fac
# indices, where modulation period is jumped
findall(sol_RadauIIA3[1].tR .- sol_RadauIIA3[1].t0 .> sys.modules[3].PM)

# ╔═╡ 4721202c-0b9c-43ac-9561-d2b040c5967e
# solutions are backwards
sol_RadauIIA3[1].tR .- sol_RadauIIA3[1].t0

# ╔═╡ 2b7e5fca-28d5-4801-beea-d60c0533e82f
md"""
- RadauIIA5
"""

# ╔═╡ daf6224d-b893-43e0-9197-3b24232fccdb
begin
	par_slice_RadauIIA5 = slicing(sim[2][1][2].tR, sim[2][1][2].τR, sys.modules[3].PM, sys.modules[3].ratio, sys.modules[3].shift, par[3], par[2]; nτ=6, abstol=1e-10, reltol=1e-8, alg=RadauIIA5()) # competly slicing up the peaks at the modulator
	sol_RadauIIA5 = GasChromatographySimulator.simulate(par_slice_RadauIIA5[1])
	
	add_A_to_pl!(sol_RadauIIA5[1], par_slice_RadauIIA5[3])
	add_t0_to_pl!(sol_RadauIIA5[1], par_slice_RadauIIA5[3])
end

# ╔═╡ f19311a3-2719-41b7-8970-6f519213adcc
# indices, where modulation period is jumped
findall(sol_RadauIIA5[1].tR .- sol_RadauIIA5[1].t0 .> sys.modules[3].PM)

# ╔═╡ 014a923d-95ce-4dd7-bcb2-97f1aee40e71
md"""
- Rodas5
"""

# ╔═╡ 2d4c7225-2e85-4a68-a9fa-5fd405d76ca4
begin
	par_slice_Rodas5 = slicing(sim[2][1][2].tR, sim[2][1][2].τR, sys.modules[3].PM, sys.modules[3].ratio, sys.modules[3].shift, par[3], par[2]; nτ=6, abstol=1e-10, reltol=1e-8, alg=Rodas5()) # competly slicing up the peaks at the modulator
	sol_Rodas5 = GasChromatographySimulator.simulate(par_slice_Rodas5[1])
	
	add_A_to_pl!(sol_Rodas5[1], par_slice_Rodas5[3])
	add_t0_to_pl!(sol_Rodas5[1], par_slice_Rodas5[3])
end

# ╔═╡ 330def8f-12ac-4573-ad4b-373a6f8a6ad0
# indices, where modulation period is jumped
findall(sol_Rodas5[1].tR .- sol_Rodas5[1].t0 .> sys.modules[3].PM)

# ╔═╡ 4f7975da-85a8-4e42-a14d-086219952863
md"""
- Rodas5P
"""

# ╔═╡ 4521feaf-80fc-41dd-8d92-ca96a3200278
begin
	par_slice_Rodas5P = slicing(sim[2][1][2].tR, sim[2][1][2].τR, sys.modules[3].PM, sys.modules[3].ratio, sys.modules[3].shift, par[3], par[2]; nτ=6, abstol=1e-10, reltol=1e-8, alg=Rodas5P()) # competly slicing up the peaks at the modulator
	sol_Rodas5P = GasChromatographySimulator.simulate(par_slice_Rodas5P[1])
	
	add_A_to_pl!(sol_Rodas5P[1], par_slice_Rodas5P[3])
	add_t0_to_pl!(sol_Rodas5P[1], par_slice_Rodas5P[3])
end

# ╔═╡ 81957543-aede-40ee-aea0-3720dcfcf7f7
# indices, where modulation period is jumped
findall(sol_Rodas5P[1].tR .- sol_Rodas5P[1].t0 .> sys.modules[3].PM)

# ╔═╡ 30236719-29f5-48fa-8670-d753aebfc549
md"""
- Tsit5
"""

# ╔═╡ c8f4396f-03bb-45e1-84bb-df4eea13c22d
begin
	par_slice_Tsit5 = slicing(sim[2][1][2].tR, sim[2][1][2].τR, sys.modules[3].PM, sys.modules[3].ratio, sys.modules[3].shift, par[3], par[2]; nτ=6, abstol=1e-10, reltol=1e-8, alg=Tsit5()) # competly slicing up the peaks at the modulator
	sol_Tsit5 = GasChromatographySimulator.simulate(par_slice_Tsit5[1])
	
	add_A_to_pl!(sol_Tsit5[1], par_slice_Tsit5[3])
	add_t0_to_pl!(sol_Tsit5[1], par_slice_Tsit5[3])
end

# ╔═╡ 0bcfb830-7127-4eee-a696-e54596453ed1
# indices, where modulation period is jumped
findall(sol_Tsit5[1].tR .- sol_Tsit5[1].t0 .> sys.modules[3].PM)

# ╔═╡ 9e6a5920-6a48-4531-8fd6-1d2358c574bc
md"""
- DP8
"""

# ╔═╡ 8a41c6e8-49a0-4450-a977-27d5e658b20d
begin
	par_slice_DP8 = slicing(sim[2][1][2].tR, sim[2][1][2].τR, sys.modules[3].PM, sys.modules[3].ratio, sys.modules[3].shift, par[3], par[2]; nτ=6, abstol=1e-10, reltol=1e-8, alg=DP8()) # competly slicing up the peaks at the modulator
	sol_DP8 = GasChromatographySimulator.simulate(par_slice_DP8[1])
	
	add_A_to_pl!(sol_DP8[1], par_slice_DP8[3])
	add_t0_to_pl!(sol_DP8[1], par_slice_DP8[3])
end

# ╔═╡ 4d4eb359-0396-4e16-8268-1f7ae6a02fc5
# indices, where modulation period is jumped
findall(sol_DP8[1].tR .- sol_DP8[1].t0 .> sys.modules[3].PM)

# ╔═╡ 296d9824-27e2-4d70-88f5-74f9b25a8348
md"""
- Vern8
"""

# ╔═╡ 73e76fd5-a799-4340-b89f-1a452f96aa48
begin
	par_slice_Vern8 = slicing(sim[2][1][2].tR, sim[2][1][2].τR, sys.modules[3].PM, sys.modules[3].ratio, sys.modules[3].shift, par[3], par[2]; nτ=6, abstol=1e-10, reltol=1e-8, alg=Vern8()) # competly slicing up the peaks at the modulator
	sol_Vern8 = GasChromatographySimulator.simulate(par_slice_Vern8[1])
	
	add_A_to_pl!(sol_Vern8[1], par_slice_Vern8[3])
	add_t0_to_pl!(sol_Vern8[1], par_slice_Vern8[3])
end

# ╔═╡ d4b3911c-7ea8-4e02-8600-ae8b4c475e49
# indices, where modulation period is jumped
findall(sol_Vern8[1].tR .- sol_Vern8[1].t0 .> sys.modules[3].PM)

# ╔═╡ c994a411-c4ad-4294-92d9-1806292f90d2
md"""
- Vern9
"""

# ╔═╡ db5c5f61-556a-4e7e-bdb1-751a74c7d370
begin
	par_slice_Vern9 = slicing(sim[2][1][2].tR, sim[2][1][2].τR, sys.modules[3].PM, sys.modules[3].ratio, sys.modules[3].shift, par[3], par[2]; nτ=6, abstol=1e-10, reltol=1e-8, alg=Vern9()) # competly slicing up the peaks at the modulator
	sol_Vern9 = GasChromatographySimulator.simulate(par_slice_Vern9[1])
	
	add_A_to_pl!(sol_Vern9[1], par_slice_Vern9[3])
	add_t0_to_pl!(sol_Vern9[1], par_slice_Vern9[3])
end

# ╔═╡ dbac14bf-2143-4de8-9bd1-915319690f3d
# indices, where modulation period is jumped
findall(sol_Vern9[1].tR .- sol_Vern9[1].t0 .> sys.modules[3].PM)

# ╔═╡ 47ddd70c-681c-44dc-aee2-4b1a38eb64f6
md"""
#### Sliced, not focussed rectangular peaks
"""

# ╔═╡ f3a793ea-833e-4155-8d5a-ba73db5cb9ee


# ╔═╡ 3f384f8f-4b97-4cf0-b28f-e567402dcf93
par_slice_unfoc_1st_TM.sub

# ╔═╡ aa8dca86-a3d0-449d-88ef-a598ed469bc2
A_slice_unfoc_1st_TM

# ╔═╡ 6ce76586-68e3-4058-aebd-d721af6fe59d
begin
	sol_1st_unfoc_TM = GasChromatographySimulator.simulate(par_slice_unfoc_1st_TM)
	
	add_A_to_pl!(sol_1st_unfoc_TM[1], A_slice_unfoc_1st_TM)
end

# ╔═╡ 3c20b2b4-ca9a-4b90-abb4-4cc6fd25a733
md"""
#### Misc
"""

# ╔═╡ f087274e-5147-48f9-8200-4d6ed90397c9
par[3].col.L ./ sol_1st_TM[1].uR # new initial peak width for the next module

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

# ╔═╡ 9e84e9e1-76f2-4598-ad62-950513d1b07b
begin
	plotly()
	p_Tt_RadauIIA3 = Plots.plot(xlabel="time in s", ylabel="temperature in °C", legend=false, title="RadauIIA3")
	for i=1:length(sol_RadauIIA3[2])
		trace = traces(sol_RadauIIA3[2], par_slice_RadauIIA3[1], i)
		Plots.plot!(p_Tt_RadauIIA3, trace.t, trace.T.-273.15)
	end
	p_Tt_RadauIIA3
end

# ╔═╡ a08277ae-7d69-4da0-8088-5f52c6ca6b73
begin
	plotly()
	p_Tt_RadauIIA5 = Plots.plot(xlabel="time in s", ylabel="temperature in °C", legend=false, title="RadauIIA5")
	for i=1:length(sol_RadauIIA5[2])
		trace = traces(sol_RadauIIA5[2], par_slice_RadauIIA5[1], i)
		Plots.plot!(p_Tt_RadauIIA5, trace.t, trace.T.-273.15)
	end
	p_Tt_RadauIIA5
end

# ╔═╡ fe255b07-c5cd-4778-bb95-080d29b80d2e
begin
	plotly()
	p_Tt_Vern9 = Plots.plot(xlabel="time in s", ylabel="temperature in °C", legend=false, title="Vern9")
	p_τt_Vern9 = Plots.plot(xlabel="time in s", ylabel="τ in s", legend=false, title="Vern9")
	for i=1:length(sol_Vern9[2])
		trace = traces(sol_Vern9[2], par_slice_Vern9[1], i)
		Plots.plot!(p_Tt_Vern9, trace.t, trace.T.-273.15)
		Plots.plot!(p_τt_Vern9, trace.t, sqrt.(trace.τ²))
	end
	Plots.plot(p_Tt_Vern9, p_τt_Vern9) 
end

# ╔═╡ 402e9e25-3d5b-46c2-8b5e-4a472fb357f7
traces(sol_1st_TM[2], par_slice_1st_TM, 1)

# ╔═╡ 9c5ec167-212d-413a-a2bc-fb00d015af51
begin
	p_Tt = Plots.plot(xlabel="time in s", ylabel="temperature in °C", legend=false)
	for i=1:length(sol_1st_TM[2])
		trace = traces(sol_1st_TM[2], newpar, i)
		Plots.plot!(p_Tt, trace.t, trace.T.-273.15)
	end
	p_Tt
end

# ╔═╡ 740e977e-4cf8-496d-bb1a-000871bc4211
begin
	p_τt = Plots.plot(xlabel="time in s", ylabel="τ in s", legend=false)
	for i=3:3#length(sol_1st_TM[2])
		trace = traces(sol_1st_TM[2], newpar, i)
		Plots.plot!(p_τt, trace.t, sqrt.(trace.τ²))
	end
	p_τt
	# Ableitung der Rechteckfuntion (Temperatur modulation) ist immer NULL -> der schlagartige Temperaturanstieg wird bei der Peakbreitensimulation nicht berücksichtigt
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
	t_sum = collect(t1:minimum(τR)/10:t2)
	c_sum = fill(0.0, length(t_sum))
	for i=1:length(pl.Name)
		if isnan(tstart[i])
			c[i] = fill(NaN,length(t_sum))
		else
			c[i] = GasChromatographySimulator.chromatogram(t_sum, [pl.tR[i]], [τR[i]])*pl.A[i]
			c_sum = c_sum .+ c[i]
		end
	end
	return t_sum, c, c_sum
end

# ╔═╡ d28d487a-10f0-45e3-8df3-e8c07692813c
A_foc = sum(sol_1st_TM[1].A)

# ╔═╡ 3fd0d9d9-740e-4c8f-8cc3-6c0a41242611
A_unfoc = sum(sol_1st_unfoc_TM[1].A)

# ╔═╡ 0ea91a94-2e95-4e7b-814c-b6058894c608
A_total = A_foc + A_unfoc

# ╔═╡ 30feef76-872e-4f5c-b303-c661846e7fcc
begin
	percentunfoc = round(A_unfoc/A_total*100; digits=2)
	percentfoc = round(A_foc/A_total*100; digits=2)
	rat = sys.modules[3].ratio
	md"""
	Area of the unfocussed peaks is much smaller, $percentunfoc %, than the area of the focussed peaks, $percentfoc %. But it is not neglectable. This distribution corresponds to the ratio of time of cold jet to time of hot jet, which is $rat. 
	"""
end

# ╔═╡ 471878ef-b115-411c-9d9d-764d6db08069
#savefig(p_chrom_TM1, "Chrom_after_first_TM_new.svg")

# ╔═╡ 4c2c3d26-951a-412f-833d-83b153b4354e
function rect(x, w, h)
	# rect starting at x=0
	# w ... width 
	# h ... height
	rec = if abs(x-w/2) > w/2
		0
	elseif abs(x-w/2) == w/2
		h/2
	elseif abs(x-w/2) < w/2
		h
	end
	return rec
end

# ╔═╡ 47bd51af-7f68-4967-aac2-c80e4ee73ff7
# plot of rectancular peaks to approximate the unfocussed slices
# width = τR
# height = A/width
function rec_chrom_after_modulation(pl)
	#tstart = Array{Float64}(undef, length(pl.Name))
	#tend = Array{Float64}(undef, length(pl.Name))
	t = Array{Array{Float64}}(undef, length(pl.Name))
	c = Array{Array{Float64}}(undef, length(pl.Name))
	h = Array{Float64}(undef, length(pl.Name))
	ts = Array{Array{Float64}}(undef, 4*length(pl.Name))
	cs = Array{Array{Float64}}(undef, 4*length(pl.Name))
	for i=1:length(pl.Name)
		h[i] = pl.A[i]/pl.τR[i]
		t[i] = (pl.tR[i]-pl.τR[i]/100):pl.τR[i]/100:(pl.tR[i]+pl.τR[i]+pl.τR[i]/100)
		c[i] = rect.(t[i].-pl.tR[i], pl.τR[i], h[i])	
	end
	t1 = t[1][1]
	t2 = t[end][end]
	t_sum = collect(t1:(t2-t1)/10000:t2)
	c_sum = fill(0.0, 10001)
	for i=1:length(pl.Name)
#		if isnan(tstart[i])
#			c[i] = fill(NaN,1001)
#		else
			c_ = rect.(t_sum.-pl.tR[i], pl.τR[i], h[i])
			c_sum = c_sum .+ c_
#		end
	end
	return t, c, h, t_sum, c_sum
end

# ╔═╡ 095ee9ba-d974-4625-a47c-3634d91c6a5d
begin
	gr()
	# focussed peaks
	t_TM1, c_TM1s, c_TM1 = chrom_after_modulation(sol_1st_TM[1], par_slice_1st_TM, width_sim=true)
	p_chrom_TM1 = Plots.plot(t_TM1, c_TM1, xlabel="time in s")
	# unfocussed peaks 
	t_unfoc, c_unfoc, h_unfoc, t_unfoc_sum, c_unfoc_sum = rec_chrom_after_modulation(sol_1st_unfoc_TM[1])
	#for i=1:length(c_TM1s)
	#	Plots.plot!(p_chrom_TM1, t_TM1, c_TM1s[i], label=sol_1st_TM[1].Name[i])
	#end
	Plots.plot!(p_chrom_TM1, xticks=0:4:3000, legend=false)
	Plots.plot!(p_chrom_TM1, t_unfoc_sum, c_unfoc_sum)
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

# ╔═╡ 6a3f4043-63ba-44ef-b9ce-a391bc697fe6
function change_initial(par_, prev_par, prev_pl; prev_TM=false)
	if prev_TM == true 
		# Ableitung der Rechteckfuntion (Temperatur modulation) ist immer NULL -> der schlagartige Temperaturanstieg wird bei der Peakbreitensimulation nicht berücksichtigt
		# if the previous module is a thermal modulator, than the new initial peak width is approximated as length of modulator spot divided by the velocity of the solute at the end of the modulator spot
		τ0_ = prev_par.col.L./prev_pl.uR
	else
		τ0_ = prev_pl.τR
	end
	new_sub = Array{GasChromatographySimulator.Substance}(undef, length(prev_par.sub))
	for i=1:length(prev_par.sub)
		# filter for correct CAS and annotation (slice number)
		ii_ = common_index(prev_pl, prev_par.sub[i].CAS, join(split(prev_par.sub[i].ann, ", ")[1:end-1], ", "))
		new_sub[ii_] = GasChromatographySimulator.Substance(prev_par.sub[i].name, prev_par.sub[i].CAS, prev_par.sub[i].Tchar, prev_par.sub[i].θchar, prev_par.sub[i].ΔCp, prev_par.sub[i].φ₀, prev_par.sub[i].ann, prev_par.sub[i].Cag, prev_pl.tR[ii_], τ0_[ii_])
	end
	# here changes of options could be applied
	new_par = GasChromatographySimulator.Parameters(par_.col, par_.prog, new_sub, par_.opt)
	return new_par
end

# ╔═╡ bac5d026-0b84-412f-98bb-8ddedb2b92d9
# combining all the step-by-step simulations 
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

# ╔═╡ 2c1fc71c-5f5d-4e18-84f5-42980d461753
newpar_loop = change_initial(par[4], par_slice_1st_TM, sol_1st_TM[1], prev_TM=false)

# ╔═╡ 74121702-9f4a-410b-961b-7865ab3af941
begin
	sol_loop = GasChromatographySimulator.simulate(newpar_loop)
	add_A_to_pl!(sol_loop[1], sol_1st_TM[1])
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
ΔtR_loop = sol_loop[1].tR .- sol_1st_TM[1].tR

# ╔═╡ 486273df-608e-4f21-ab97-d5466f38f6e8
md"""
#### Sliced not focused rectangular peaks
"""

# ╔═╡ 6b6e1ada-aa27-4156-967a-49e5f076727f
newpar_loop_unfoc = change_initial(par[4], par_slice_unfoc_1st_TM, sol_1st_unfoc_TM[1], prev_TM=false)

# ╔═╡ 093e0c79-5f78-4616-b765-48b46d64eec5
begin
	sol_loop_unfoc = GasChromatographySimulator.simulate(newpar_loop_unfoc)
	add_A_to_pl!(sol_loop_unfoc[1], sol_1st_unfoc_TM[1])
end

# ╔═╡ 66eaa639-416c-4167-9e00-14dcbf583a4c
md"""
#### Chromatogram after loop
"""

# ╔═╡ bb483a0e-72d7-405e-bf47-c5fa878e0994
begin
	plotly()
	
	# focussed peaks
	t_loop, c_loops, c_loop = chrom_after_loop(sol_loop[1])
	# unfocussed peaks 
	t_loop_unfoc, c_loop_unfoc, h_loop_unfoc, t_loop_unfoc_sum, c_loop_unfoc_sum = rec_chrom_after_modulation(sol_loop_unfoc[1])
	
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
	
	Plots.plot!(p_chrom_loop, t_loop_unfoc_sum, c_loop_unfoc_sum, linestyle=:dash, c=1)
	for i=1:length(c_loop_unfoc)
		if sol_loop_unfoc[1].Name[i] == selected_solutes[1]
			color = :blue
		elseif sol_loop_unfoc[1].Name[i] == selected_solutes[2]
			color = :red
		else
			color = :orange
		end
		Plots.plot!(p_chrom_loop, t_loop_unfoc[i], c_loop_unfoc[i], label=sol_loop_unfoc[1].Name[i], c=color, linestyle=:dash)
	end
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

# ╔═╡ 0d51d9c5-8a70-4d4d-b247-4c5fe07b7d1c
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
	
	Plots.plot!(p_chrom_loop_mark, t_loop_unfoc_sum, c_loop_unfoc_sum, linestyle=:dash, c=1)
	for i=1:length(c_loop_unfoc)
		if sol_loop_unfoc[1].Name[i] == selected_solutes[1]
			color = :blue
		elseif sol_loop_unfoc[1].Name[i] == selected_solutes[2]
			color = :red
		else
			color = :orange
		end
		Plots.plot!(p_chrom_loop_mark, t_loop_unfoc[i], c_loop_unfoc[i], label=sol_loop_unfoc[1].Name[i], c=color, linestyle=:dash)
	end

	#add marker for PM and tcold
	n = unique(fld.(t_loop, sys.modules[5].PM))
	for i=1:length(n)
		Plots.plot!(p_chrom_loop_mark, [n[i]*sys.modules[5].PM, n[i]*sys.modules[5].PM], [0.0, 10.0], c=:black)
		Plots.plot!(p_chrom_loop_mark, [n[i]*sys.modules[5].PM+sys.modules[5].PM*sys.modules[5].ratio, n[i]*sys.modules[5].PM+sys.modules[5].PM*sys.modules[5].ratio], [0.0, 10.0], c=:black, linestyle=:dash)
	end
	p_chrom_loop_mark
end

# ╔═╡ e1ce8764-cec4-4c43-9354-ccceea375f3f
# for now, this 2nd modulator spot is deactivated and works with the oven temperature

# ╔═╡ 269dde5f-9611-458e-a731-ce37cd5bdb97
tR_foc = sol_loop[1].tR

# ╔═╡ 086c2c55-bd1f-4fb3-ad05-e7f6abbf87b8
τR_foc = sol_loop[1].τR

# ╔═╡ a05a4081-db0d-406c-bee6-f5399c49b123
tR_unfoc = sol_loop_unfoc[1].tR

# ╔═╡ 191670bd-952d-4521-b8eb-88202db5f98f
τR_unfoc = sol_loop_unfoc[1].τR

# ╔═╡ 6f727703-05c3-4c1b-afad-084bb59dfb9b
shift = sys.modules[5].shift

# ╔═╡ cd1a2bed-0d5a-4530-8b47-51ca8e46709a
(1-sys.modules[5].ratio)*sys.modules[5].PM

# ╔═╡ becced0a-364f-4414-bcec-4b2b009d46b5
# warning
# 1. time in loop should be bigger than time of hot jet
all(ΔtR_loop .- nτ .* τR_foc .> (1-sys.modules[5].ratio)*sys.modules[5].PM)

# ╔═╡ 0924d96b-c598-45a3-ba08-9d0ae95e4a3a
# warning
# 2. time in loop + time of hot jet (roughly the width of the unfocussed peak) should be smaller than time of cold jet
ΔtR_loop .+ (1-sys.modules[5].ratio)*sys.modules[5].PM .< sys.modules[5].ratio*sys.modules[5].PM

# ╔═╡ 2390b6a8-8927-466d-a392-f76c4cd2c7f2
function warning_2nd_TM(ΔtR_loop, τR_foc_loop, PM, ratio)
	if all(ΔtR_loop .- nτ .* τR_foc .> (1-sys.modules[5].ratio)*sys.modules[5].PM) == false
		@warn "Time in loop is shorter than time of hot jet. Breaktrough in the modulator possible."
	elseif all(ΔtR_loop .+ 1.05*(1-sys.modules[5].ratio)*sys.modules[5].PM .< sys.modules[5].ratio*sys.modules[5].PM) == false
		@warn "Time in loop (plus width of unfocussed part) is longer than time of cold jet. Breaktrough in the modulator possible."
	end
end

# ╔═╡ f071a5bc-922d-4a17-81c1-2e809ff380a5
warning_2nd_TM(ΔtR_loop, τR_foc, sys.modules[5].PM, sys.modules[5].ratio)

# ╔═╡ 5ea43b31-b41b-4e0c-9f4e-dafa8560505e
ΔtR_loop .+ 1.0*(1-sys.modules[5].ratio)*sys.modules[5].PM .< sys.modules[5].ratio*sys.modules[5].PM

# ╔═╡ 30c43181-3ccb-4ec1-bab9-31f35f2727ce
ΔtR_loop .+ 1.0*(1-sys.modules[5].ratio)*sys.modules[5].PM

# ╔═╡ 854277d4-9a74-4abe-b827-b37ff532c6db
sys.modules[5].ratio*sys.modules[5].PM

# ╔═╡ bb24116c-5828-40ce-9e8e-673dceb19a66
init_t_start_foc = (fld.(tR_foc.-nτ.*τR_foc, PM)).*PM .+ shift

# ╔═╡ 723b771c-d727-4438-8092-06a9ad060acc
# warning
tR_foc.+(1-sys.modules[5].ratio)*sys.modules[5].PM .- init_t_start_foc

# ╔═╡ 3699039a-62fe-43e5-be1d-b1b166e9afbe
init_t_end_foc = (fld.(tR_foc.+nτ.*τR_foc, PM)).*PM .+ shift

# ╔═╡ d5f9d884-0274-4f71-b932-0d39fc457b14
n_slice_foc = Int.((init_t_end_foc .- init_t_start_foc)./PM .+ 1)

# ╔═╡ 71bc9925-13fb-4823-9f7c-4d9585114df5
init_t_start_unfoc = (fld.(tR_unfoc.-1.0.*τR_unfoc, PM)).*PM .+ shift 

# ╔═╡ 06771585-9592-4542-979e-e9a271028094
init_t_end_unfoc = (fld.(tR_unfoc.+1.0.*τR_unfoc, PM)).*PM .+ shift

# ╔═╡ 3b8841e3-60dd-4cc2-8e93-9cc274a77225
n_slice_unfoc = Int.((init_t_end_unfoc .- init_t_start_unfoc)./PM .+ 1)

# ╔═╡ f7b718fb-488a-4eea-9e5d-b6d7b4cb29da
A_foc_ = sol_loop[1].A

# ╔═╡ ce48d436-4ebe-4887-b513-fb2562cd76d8
A_unfoc_ = sol_loop_unfoc[1].A

# ╔═╡ 914aa8c6-6ae1-4f56-a978-60496e3fbbc1
n_slice_foc == n_slice_unfoc

# ╔═╡ 0f213e62-ed46-4633-bdc6-3d3b729ef12a
init_t_start_foc == init_t_start_unfoc

# ╔═╡ 43d79322-59f5-4d07-9a42-2abb964cc0a9
sub_TM2 = Array{GasChromatographySimulator.Substance}(undef, sum(n_slice_foc))

# ╔═╡ 5355575b-d3a6-436d-9c68-93a969cad4aa
A_TM2 = Array{Float64}(undef, sum(n_slice))

# ╔═╡ 3e95800f-69c3-445d-936c-1da337f77dbd
τ₀=zeros(length(τR_foc))

# ╔═╡ 9b3cd294-d6b3-42cd-81a1-07ad6119e012
for i=1:length(n_slice_foc)
		#for j=1:n_slice_foc[i] # n_slice_foc[i] should always be = 1
			t₀ = init_t_start_foc[i]#+(j-1)*PM # initial start time
			
			sub_TM2[i] = GasChromatographySimulator.Substance(newpar_loop.sub[i].name, newpar_loop.sub[i].CAS, newpar_loop.sub[i].Tchar, newpar_loop.sub[i].θchar, newpar_loop.sub[i].ΔCp, newpar_loop.sub[i].φ₀, "TM2_f1_"*newpar_loop.sub[i].ann, newpar_loop.sub[i].Cag, t₀, τ₀[i])

			A_TM2[i] = A_foc_[i] + A_unfoc_[i]
		#end
	end

# ╔═╡ f845cfe5-bdde-4bc8-851f-3b05d14d75ad
sub_TM2

# ╔═╡ 16630a1c-d7e5-405d-b3fc-c802b2f48c85
A_TM2

# ╔═╡ 44c7cc8e-a02b-437a-8da0-89f94cac878f
par

# ╔═╡ 96a0bd4f-cfe1-4dc3-8059-ee270e979ddb
par_2nd_TM_ = GasChromatographySimulator.Parameters(par[5].col, par[5].prog, sub_TM2, GasChromatographySimulator.Options(par[5].opt.alg, 1e-8, 1e-5, par[5].opt.Tcontrol, par[5].opt.odesys, par[5].opt.ng, par[5].opt.vis, par[5].opt.control, par[5].opt.k_th))

# ╔═╡ a205623b-bd86-4cec-923f-7a87ede6d408
function slicing_2nd(pl_foc, pl_unfoc, PM, ratio, shift, par::GasChromatographySimulator.Parameters, prev_par::GasChromatographySimulator.Parameters, prev_par_unfoc::GasChromatographySimulator.Parameters; nτ=6, τ₀=zeros(length(τR_foc)), abstol=1e-8, reltol=1e-5, alg=OwrenZen5())
	
	tcold = PM*ratio
	thot = PM*(1-ratio)
	tR_foc = pl_foc.tR
	tR_unfoc = pl_unfoc.tR
	τR_foc = pl_foc.τR
	τR_unfoc = pl_unfoc.τR
	A_foc = pl_foc.A
	A_unfoc = pl_unfoc.A
	# focussed peaks
	init_t_start_foc = (fld.(tR_foc.-nτ.*τR_foc, PM)).*PM .+ shift # start time of the peaks, rounded down to multiple of PM
	init_t_end_foc = (fld.(tR_foc.+nτ.*τR_foc, PM)).*PM .+ shift # end time of the peaks, rounded down to multiple of PM
	n_slice_foc = Int.((init_t_end_foc .- init_t_start_foc)./PM .+ 1) # number of slices for every substance, should always be one, otherwise a warning should be given

	# unfocussed peaks
	init_t_start_unfoc = (fld.(tR_unfoc.-1.0.*τR_unfoc, PM)).*PM .+ shift # start time of the peaks, rounded down to multiple of PM
	init_t_end_unfoc = (fld.(tR_unfoc.+1.0.*τR_unfoc, PM)).*PM .+ shift # end time of the peaks, rounded down to multiple of PM
	n_slice_unfoc = Int.((init_t_end_unfoc .- init_t_start_unfoc)./PM .+ 1) # number of slices for every substance, should always be one, otherwise a warning should be given

	# focussed and unfocussed peaks of the same substance in the same modulation periode have to be combined, there areas added up

	#if init_t_start_foc == init_t_start_unfoc # in every slice one focussed and one unfocussed peak of the same substance should be present
		# the order in both should be the same
	#---------
	sub_TM2 = Array{GasChromatographySimulator.Substance}(undef, sum(n_slice_foc))
	#sub_TM_unfocussed = Array{GasChromatographySimulator.Substance}(undef, sum(n_slice))
	A_TM2 = Array{Float64}(undef, sum(n_slice_foc))
	#A_unfocussed = Array{Float64}(undef, sum(n_slice))
	#g(x,p) = 1/sqrt(2*π*p[2]^2)*exp(-(x-p[1])^2/(2*p[2]^2))
	#ii = 1
	Name = Array{String}(undef, sum(n_slice_foc))
	CAS = Array{String}(undef, sum(n_slice_foc))
	Ann = Array{String}(undef, sum(n_slice_foc))
	t0 = Array{Float64}(undef, sum(n_slice_foc))
	for i=1:length(n_slice_foc)
		#for j=1:n_slice_foc[i] # n_slice_foc[i] should always be = 1
			t₀ = init_t_start_foc[i]#+(j-1)*PM # initial start time
			
			sub_TM2[i] = GasChromatographySimulator.Substance(prev_par.sub[i].name, prev_par.sub[i].CAS, prev_par.sub[i].Tchar, prev_par.sub[i].θchar, prev_par.sub[i].ΔCp, prev_par.sub[i].φ₀, "TM2_f1, "*prev_par.sub[i].ann, prev_par.sub[i].Cag, t₀, τ₀[i])

			ii_foc = common_index(pl_foc, prev_par.sub[i].CAS, join(split(prev_par.sub[i].ann, ", ")[1:end-1], ", ")) # find the correct index of the areas  
			ii_unfoc = common_index(pl_unfoc, prev_par_unfoc.sub[i].CAS, join(split(prev_par_unfoc.sub[i].ann, ", ")[1:end-1], ", ")) 
			A_TM2[i] = A_foc[ii_foc] + A_unfoc[ii_unfoc]
			Name[i] = sub_TM2[i].name
			CAS[i] = sub_TM2[i].CAS
			Ann[i] = sub_TM2[i].ann
			t0[i] = t₀
		#end
	end
	newopt = GasChromatographySimulator.Options(alg, abstol, reltol, par.opt.Tcontrol, par.opt.odesys, par.opt.ng, par.opt.vis, par.opt.control, par.opt.k_th)
	newpar_TM2 = GasChromatographySimulator.Parameters(par.col, par.prog, sub_TM2, newopt)

	df_A_TM2 = DataFrame(Name = Name, CAS = CAS, Annotations = Ann, A = A_TM2, t0 = t0)
	
	return newpar_TM2, df_A_TM2
end

# ╔═╡ 768be212-504b-4530-8456-8a719ff2707f
par_2nd_TM, df_A_2nd_TM = slicing_2nd(sol_loop[1], sol_loop_unfoc[1], sys.modules[5].PM, sys.modules[5].ratio, sys.modules[5].shift, par[5], newpar_loop, newpar_loop_unfoc, abstol=1e-10, reltol=1e-8, alg=Vern9())

# ╔═╡ e78ffe7d-3f9f-471c-925c-4feae33231b8
begin
	sol_2nd_TM = GasChromatographySimulator.simulate(par_2nd_TM)
	add_A_to_pl!(sol_2nd_TM[1], df_A_2nd_TM)
	add_t0_to_pl!(sol_2nd_TM[1], df_A_2nd_TM)
	sol_2nd_TM[1]
end

# ╔═╡ f17a9a69-a9bf-4a7f-b6a6-39ed6350135c
sol_2nd_TM[1].tR .- sol_2nd_TM[1].t0

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
	t_TM2, c_TM2s, c_TM2 = chrom_after_modulation(sol_2nd_TM[1], par_2nd_TM, width_sim=true)
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
newpar_mod_out = change_initial(par[6], par_2nd_TM, sol_2nd_TM[1], prev_TM=false)

# ╔═╡ 0400ee95-e5f2-48f3-bc57-8f3eb24a2e21
sol_2nd_TM[1]

# ╔═╡ bb01e67d-fbca-461e-b592-508bc5480b94
begin
	sol_mod_out = GasChromatographySimulator.simulate(newpar_mod_out)
	add_A_to_pl!(sol_mod_out[1], sol_2nd_TM[1])
end

# ╔═╡ 71875c85-58b3-46bc-ad5a-787bfb6efcc3
filter([:Name] => x -> x == "(S)-(+)-Carvone", sol_mod_out[1])

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
newpar_gc2 = change_initial(par[7], newpar_mod_out, sol_mod_out[1], prev_TM=false)

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
newpar_tl = change_initial(par[8], newpar_gc2, sol_gc2[1], prev_TM=false)

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
	fit_1D = fit_gauss_D1(sol_tl[1], sim[2][1][2])
	for i=1:length(fit_1D.Name)
		Plots.scatter!(p_chrom_tl_proj1, fit_1D.tRs[i], fit_1D.heights[i], msize=2, c=dot_color(fit_1D.Name[i], selected_solutes))
		Plots.plot!(p_chrom_tl_proj1, t_tl, model_g(t_tl, fit_1D.fits[i].param), c=dot_color(fit_1D.Name[i], selected_solutes), linestyle=:dash)
	end
	Plots.plot!(p_chrom_tl_proj1, xticks=0:4:3000, legend=false)
	p_chrom_tl_proj1
end

# ╔═╡ 3a938dd0-1799-450c-9341-6373af448192
fit_1D

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
	mins = minimum.(fit_2D.tRs)
	maxs = maximum.(fit_2D.tRs)
	for i=1:length(fit_2D.Name)
		Plots.scatter!(p_chrom_tl_proj2, fit_2D.tRs[i], fit_2D.heights[i], msize=2, c=dot_color(fit_2D.Name[i], selected_solutes))

		t = 0.9*mins[i]:(1.1*maxs[i]-0.9*mins[i])/1000.0:1.1*maxs[i]
		Plots.plot!(p_chrom_tl_proj2, t, model_g(t, fit_2D.fits[i].param), c=dot_color(fit_2D.Name[i], selected_solutes), linestyle=:dash)
	end
	p_chrom_tl_proj2
end

# ╔═╡ 3361d83c-8beb-465e-8af6-25c0d78dc981
fit_2D

# ╔═╡ 0205e873-f3b5-4100-a49b-5824fc5eecf3
begin
	name = Array{String}(undef, length(selected_solutes))
	tR1s = Array{Float64}(undef, length(selected_solutes))
	tR2s = Array{Float64}(undef, length(selected_solutes))
	for i=1:length(selected_solutes)
		name[i] = fit_1D.Name[i]
		tR1s[i] = fit_1D.fits[i].param[1]
		ii = findfirst(name[i].==fit_2D.Name)
		tR2s[i] = fit_2D.fits[ii].param[1]
	end
	pl_GCxGC = DataFrame(Name=name, tR1=tR1s, tR2=tR2s)
end

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

# ╔═╡ 8361678e-60f9-4a0c-b027-71304fe52c4a
pl_GCxGC

# ╔═╡ 4a7bf73a-7087-4f35-b4be-54ad38355a64
fld(minimum(pl_GCxGC.tR1),4)*4

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

# ╔═╡ 85ecccaa-f38f-4d54-b63a-33abf5390dc5
begin
	t_ = 0.0:0.01:sum(GCxGC_TP.timesteps)
	chrom_sliced = Array{Array{Float64}}(undef, length(sol_final[1].tR))
	for i=1:length(sol_final[1].tR)
		chrom_sliced[i] = GasChromatographySimulator.chromatogram(collect(t_), [sol_final[1].tR[i]], [sol_final[1].τR[i]]).*sol_final[1].A[i]
	end
	chrom_sliced_sum = chrom_sliced[1]
	for i=2:length(chrom_sliced)
		chrom_sliced_sum = chrom_sliced_sum .+ chrom_sliced[i]
	end
	Plots.plot(t_, chrom_sliced_sum)
	c_slices, t_D1, t_D2 = chrom_slicing(t_, chrom_sliced_sum, PM)
	slice_mat = Array{Float64}(undef, length(c_slices)-1, length(t_D2[1]))
	for j=1:length(t_D2[1])
		for i=1:(length(c_slices)-1)
			slice_mat[i,j] = c_slices[i][j]
		end
	end
	slice_mat
end

# ╔═╡ af6e107a-6a36-4a23-8a1b-9c58dff9d22f
begin
	gr()
	Plots.contour(t_D1[1:end-1], t_D2[1], slice_mat', levels=10, xlabel="tR1 in s", ylabel="tR2 in s")#, c=:jet1)
	Plots.scatter!(pl_GCxGC.tR1, pl_GCxGC.tR2, m=2, shape=:diamond)
	Plots.plot!(xlims=(fld(minimum(pl_GCxGC.tR1),sys.modules[3].PM)*sys.modules[3].PM - 10*sys.modules[3].PM, fld(maximum(pl_GCxGC.tR1),sys.modules[3].PM)*sys.modules[3].PM + 10*sys.modules[3].PM), ylims=(0,sys.modules[3].PM))
	#savefig("2D_chrom_3terpenes.svg")
end

# ╔═╡ 0cadf6e0-c0a7-4632-9d4f-e1754123b75a
itp_slice_mat = LinearInterpolation((t_D1[1:end-1], t_D2[1],), slice_mat, extrapolation_bc=Flat())

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
	offset_D1 = 0.0#15.0
	Plots.scatter!(meas.tR1[1:5].-offset_D1*PM, meas.tR2[1:5], markersize=4, c=:red, label="measured");# alcohols
	Plots.scatter!(meas.tR1[6:22].-offset_D1*PM, meas.tR2[6:22], markersize=4, c=:orange, label="measured") # terpenes
	Plots.scatter!(meas.tR1[23:28].-offset_D1*PM, meas.tR2[23:28], markersize=4, c=:yellow, label="measured"); # phenones
	Plots.scatter!(meas.tR1[29:35].-offset_D1*PM, meas.tR2[29:35], markersize=4, c=:lawngreen, label="measured"); # ketones

	# simulation
	offset_D2 = 0.0
	Plots.scatter!(pl_GCxGC.tR1.-offset_D1*PM, pl_GCxGC.tR2.+offset_D2, c=:lightblue, m=:diamond, markeralpha=1, msize=3, label="simulated", xlims=(0.0, 3500.0), ylims=(0.0, 4.0), xlabel="tR1 in s", ylabel="tR2 in s");
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
# ╠═9a6dea14-f8d4-4b79-9c02-268c0c799515
# ╠═e44bcb5f-353e-4ef1-bf4d-d226a864ad48
# ╠═7ccb0031-e7e0-43e9-952e-4f80501716ec
# ╠═95cc1f3d-eff3-4cc8-8ba1-7ee04e887b84
# ╠═ab143f29-0a96-4ae2-a178-6bf773e6cd26
# ╠═1a081250-c6a6-42a2-9839-f2fd6c6ea0a2
# ╠═0c1e27ab-3b5e-4f42-8352-73756042efa4
# ╠═196c80bb-55c9-4e71-8776-b99397930a4e
# ╠═54058244-e33a-42b6-badb-72cd3fc9696e
# ╠═9efd0ffb-99f7-4500-9ea2-762b91c82448
# ╠═5cf107e3-be39-44ae-a7a4-f9463828376a
# ╠═06fa229a-d4e5-4662-a8b9-33c3874ac3f0
# ╠═8d4391c3-3712-4aeb-9868-3d1d0fd09745
# ╠═f6ce117e-e382-439d-a943-862282714b90
# ╠═1f4ca07c-bac7-4d16-8570-e1b11222640d
# ╠═402ed0ed-2126-47fd-b044-9542f997a050
# ╟─d0adda7c-8980-4d74-8b79-c3c5edd8131f
# ╠═ab33d65f-f0ba-4e5d-b660-2ecd106d8d21
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
# ╠═1f62b350-332e-4acf-b907-75362157e4da
# ╠═f1bf21dd-07ee-4e4e-beb5-7d9225b3fc6b
# ╠═ceb6ba3f-aece-4ee3-b757-5904fa04c7ac
# ╠═3e8d2e43-1fa9-40e4-a132-7228664b7546
# ╟─bac5d026-0b84-412f-98bb-8ddedb2b92d9
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
# ╟─325a8e2e-ae2f-4dd9-a32b-00c7adf47461
# ╠═3388fdf6-dff9-4d44-9db3-5bf77a57980d
# ╠═4a20d35c-2156-4e5d-a94c-dc37df64545a
# ╠═e4172084-2dd7-4daf-b2b9-cc84891c06f2
# ╠═629d4ad1-9961-4d2c-bb41-3c3832bc20fd
# ╠═65c2d80b-152c-411d-b4e2-eff77e29acba
# ╠═f436afe9-0440-4149-be61-b426e243b34d
# ╠═91cf4b42-140d-490f-949a-b655c2254c83
# ╠═846c3cab-90d4-4e73-9ce5-c1b41f937cdf
# ╠═4418401d-c418-422c-a341-7f01e7a8a3dd
# ╠═3c22e8aa-376d-4eec-9650-3bb102d4a0b8
# ╠═ca0bbf4a-5cef-4827-b040-cdad9d1dacb8
# ╠═7a81f5de-a2bc-4828-b5df-92c9b7950991
# ╠═9c30bd27-c649-4ed9-9660-23de44b8ba06
# ╠═ccd40fbd-ee6b-4571-a4ad-d2bc6d07a987
# ╠═5ae945d5-831a-46b4-b0bb-c461b3301ace
# ╠═ed8d11a2-8583-4e42-8576-e3bdcb29b8d5
# ╠═4f9e4d5e-be9a-432d-8032-0176033975f6
# ╠═0fac473b-42c1-41ca-a673-4e3e48e825d8
# ╠═d8949e4a-77aa-485e-8911-e9faf7bbffc5
# ╠═40bd41a7-cf33-4da1-9869-332c16ba597c
# ╠═406b56cb-a215-4e89-a1ab-37233156ad9a
# ╠═c34d2c93-9a01-4674-91f4-37cfccb399e4
# ╠═bca11846-b79f-48b6-86e0-3fdb1b4c1fac
# ╠═4721202c-0b9c-43ac-9561-d2b040c5967e
# ╠═9e84e9e1-76f2-4598-ad62-950513d1b07b
# ╠═2b7e5fca-28d5-4801-beea-d60c0533e82f
# ╠═daf6224d-b893-43e0-9197-3b24232fccdb
# ╠═f19311a3-2719-41b7-8970-6f519213adcc
# ╠═a08277ae-7d69-4da0-8088-5f52c6ca6b73
# ╠═014a923d-95ce-4dd7-bcb2-97f1aee40e71
# ╠═2d4c7225-2e85-4a68-a9fa-5fd405d76ca4
# ╠═330def8f-12ac-4573-ad4b-373a6f8a6ad0
# ╠═4f7975da-85a8-4e42-a14d-086219952863
# ╠═4521feaf-80fc-41dd-8d92-ca96a3200278
# ╠═81957543-aede-40ee-aea0-3720dcfcf7f7
# ╠═30236719-29f5-48fa-8670-d753aebfc549
# ╠═c8f4396f-03bb-45e1-84bb-df4eea13c22d
# ╠═0bcfb830-7127-4eee-a696-e54596453ed1
# ╠═9e6a5920-6a48-4531-8fd6-1d2358c574bc
# ╠═8a41c6e8-49a0-4450-a977-27d5e658b20d
# ╠═4d4eb359-0396-4e16-8268-1f7ae6a02fc5
# ╠═296d9824-27e2-4d70-88f5-74f9b25a8348
# ╠═73e76fd5-a799-4340-b89f-1a452f96aa48
# ╠═d4b3911c-7ea8-4e02-8600-ae8b4c475e49
# ╠═c994a411-c4ad-4294-92d9-1806292f90d2
# ╠═db5c5f61-556a-4e7e-bdb1-751a74c7d370
# ╠═dbac14bf-2143-4de8-9bd1-915319690f3d
# ╠═fe255b07-c5cd-4778-bb95-080d29b80d2e
# ╠═47ddd70c-681c-44dc-aee2-4b1a38eb64f6
# ╠═f3a793ea-833e-4155-8d5a-ba73db5cb9ee
# ╠═3f384f8f-4b97-4cf0-b28f-e567402dcf93
# ╠═aa8dca86-a3d0-449d-88ef-a598ed469bc2
# ╠═6ce76586-68e3-4058-aebd-d721af6fe59d
# ╠═3c20b2b4-ca9a-4b90-abb4-4cc6fd25a733
# ╠═f087274e-5147-48f9-8200-4d6ed90397c9
# ╠═402e9e25-3d5b-46c2-8b5e-4a472fb357f7
# ╠═9c5ec167-212d-413a-a2bc-fb00d015af51
# ╠═740e977e-4cf8-496d-bb1a-000871bc4211
# ╠═06c3c74b-4ced-405e-89be-a05092fab027
# ╠═67a4dff7-2634-4ae2-adb0-515807988b17
# ╠═84111543-21e2-4891-9050-0c102cd46a52
# ╠═9c6573e8-bc0f-4c40-9eea-9733726fbd76
# ╠═095ee9ba-d974-4625-a47c-3634d91c6a5d
# ╠═d28d487a-10f0-45e3-8df3-e8c07692813c
# ╠═3fd0d9d9-740e-4c8f-8cc3-6c0a41242611
# ╠═0ea91a94-2e95-4e7b-814c-b6058894c608
# ╠═30feef76-872e-4f5c-b303-c661846e7fcc
# ╠═471878ef-b115-411c-9d9d-764d6db08069
# ╠═47bd51af-7f68-4967-aac2-c80e4ee73ff7
# ╠═4c2c3d26-951a-412f-833d-83b153b4354e
# ╠═f6093dce-30a7-4628-8362-0c5ed8c55ebc
# ╠═beb15214-b1bb-4b39-a2aa-86f0727bc2c4
# ╠═2c1fc71c-5f5d-4e18-84f5-42980d461753
# ╠═6a3f4043-63ba-44ef-b9ce-a391bc697fe6
# ╠═74121702-9f4a-410b-961b-7865ab3af941
# ╠═5d0cb969-8a98-4d28-aa5c-233723755f41
# ╠═96200092-8afd-4ce3-8287-79806a331656
# ╠═486273df-608e-4f21-ab97-d5466f38f6e8
# ╠═6b6e1ada-aa27-4156-967a-49e5f076727f
# ╠═093e0c79-5f78-4616-b765-48b46d64eec5
# ╠═66eaa639-416c-4167-9e00-14dcbf583a4c
# ╠═bb483a0e-72d7-405e-bf47-c5fa878e0994
# ╠═53aff8ef-4585-4d03-921b-7f19827dc04e
# ╠═b31fa530-922f-4862-b1cf-04d2ec34b1c1
# ╠═86274690-dfe8-4f71-8084-66fe5f03f932
# ╟─0d51d9c5-8a70-4d4d-b247-4c5fe07b7d1c
# ╠═e1ce8764-cec4-4c43-9354-ccceea375f3f
# ╠═269dde5f-9611-458e-a731-ce37cd5bdb97
# ╠═086c2c55-bd1f-4fb3-ad05-e7f6abbf87b8
# ╠═a05a4081-db0d-406c-bee6-f5399c49b123
# ╠═191670bd-952d-4521-b8eb-88202db5f98f
# ╠═6f727703-05c3-4c1b-afad-084bb59dfb9b
# ╠═cd1a2bed-0d5a-4530-8b47-51ca8e46709a
# ╠═723b771c-d727-4438-8092-06a9ad060acc
# ╠═becced0a-364f-4414-bcec-4b2b009d46b5
# ╠═0924d96b-c598-45a3-ba08-9d0ae95e4a3a
# ╠═2390b6a8-8927-466d-a392-f76c4cd2c7f2
# ╠═f071a5bc-922d-4a17-81c1-2e809ff380a5
# ╠═5ea43b31-b41b-4e0c-9f4e-dafa8560505e
# ╠═30c43181-3ccb-4ec1-bab9-31f35f2727ce
# ╠═854277d4-9a74-4abe-b827-b37ff532c6db
# ╠═bb24116c-5828-40ce-9e8e-673dceb19a66
# ╠═3699039a-62fe-43e5-be1d-b1b166e9afbe
# ╠═d5f9d884-0274-4f71-b932-0d39fc457b14
# ╠═71bc9925-13fb-4823-9f7c-4d9585114df5
# ╠═06771585-9592-4542-979e-e9a271028094
# ╠═3b8841e3-60dd-4cc2-8e93-9cc274a77225
# ╠═f7b718fb-488a-4eea-9e5d-b6d7b4cb29da
# ╠═ce48d436-4ebe-4887-b513-fb2562cd76d8
# ╠═914aa8c6-6ae1-4f56-a978-60496e3fbbc1
# ╠═0f213e62-ed46-4633-bdc6-3d3b729ef12a
# ╠═43d79322-59f5-4d07-9a42-2abb964cc0a9
# ╠═5355575b-d3a6-436d-9c68-93a969cad4aa
# ╠═3e95800f-69c3-445d-936c-1da337f77dbd
# ╠═9b3cd294-d6b3-42cd-81a1-07ad6119e012
# ╠═f845cfe5-bdde-4bc8-851f-3b05d14d75ad
# ╠═16630a1c-d7e5-405d-b3fc-c802b2f48c85
# ╠═44c7cc8e-a02b-437a-8da0-89f94cac878f
# ╠═96a0bd4f-cfe1-4dc3-8059-ee270e979ddb
# ╠═a205623b-bd86-4cec-923f-7a87ede6d408
# ╠═768be212-504b-4530-8456-8a719ff2707f
# ╠═e78ffe7d-3f9f-471c-925c-4feae33231b8
# ╠═f17a9a69-a9bf-4a7f-b6a6-39ed6350135c
# ╠═db04ef60-3490-49c3-a5c1-f198c4015b06
# ╠═ebc1b905-eef4-4c66-abb6-b070d19f17d9
# ╠═672e5876-0f0d-483a-9936-604b3fc7f97a
# ╠═231edaa6-8ef3-4a85-b48e-e0e1829d124d
# ╠═64e635f8-f5fd-4935-b44e-b02ebdd96b64
# ╠═78dbaf44-0468-417a-8c69-90c139e2e2ad
# ╠═0400ee95-e5f2-48f3-bc57-8f3eb24a2e21
# ╠═bb01e67d-fbca-461e-b592-508bc5480b94
# ╠═71875c85-58b3-46bc-ad5a-787bfb6efcc3
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
# ╠═3a938dd0-1799-450c-9341-6373af448192
# ╠═3361d83c-8beb-465e-8af6-25c0d78dc981
# ╠═0205e873-f3b5-4100-a49b-5824fc5eecf3
# ╠═aecb6572-cf17-428d-845f-ba3cf2f0a5d4
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
# ╠═85ecccaa-f38f-4d54-b63a-33abf5390dc5
# ╠═8361678e-60f9-4a0c-b027-71304fe52c4a
# ╠═4a7bf73a-7087-4f35-b4be-54ad38355a64
# ╠═af6e107a-6a36-4a23-8a1b-9c58dff9d22f
# ╠═0cadf6e0-c0a7-4632-9d4f-e1754123b75a
# ╠═622a22d0-ecac-402f-96a5-79c6dc1683d3
# ╠═e2eaa21a-fbac-4d00-8663-c69ee4e0a8a9
# ╠═e9f8df08-7538-4dd4-8bbb-4a6649ae5b1a
# ╠═f55d7b20-2966-4a52-95fd-2f11a10d387a
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
