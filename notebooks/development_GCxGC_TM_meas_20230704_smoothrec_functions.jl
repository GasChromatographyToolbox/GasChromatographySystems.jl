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

# ╔═╡ 95cc1f3d-eff3-4cc8-8ba1-7ee04e887b84
using SpecialFunctions

# ╔═╡ d7f045e5-47c7-44b6-8bce-ce9dc3154fff
using Integrals

# ╔═╡ 2165e2d3-bca9-41f1-acee-a72a0ad9d898
using OrdinaryDiffEq

# ╔═╡ 39ddb95f-40a1-4cbb-8317-3ef18fec4d0b
using Interpolations

# ╔═╡ 08d9257e-6d49-4cf6-97a9-f8b546d9e933
using LsqFit

# ╔═╡ d419921b-1472-4b66-bae4-469537259814
md"""
# Development of a thermal modulated GCxGC system in detail
"""

# ╔═╡ dff75c18-eee5-4f70-9509-0e8228932819
md"""
## Periodic smoothed rectangle temperature function
"""

# ╔═╡ e44bcb5f-353e-4ef1-bf4d-d226a864ad48
# copied to GasChromatographySystems
begin
	# definitions for the periodic smoothed rectangle function
	# attention to the used values

	# general smooth rectangle function
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

	# smooth rectangle with mitpoint of the rising flank at 'xstart', after 'width' the falling flank, minimum values 'min' and maximum values 'max'. 
	function smooth_rectangle(x, xstart, width, min, max; flank=20)
		a = xstart #- Δx₁/2
		b = xstart + width#/2
		m = width/flank
		val = (max-min)*smooth_rectangle(x, a, b, m) + min
		return val
	end

	# periodic repeated smoothed rectangle function with period 'PM', a shift by 'shift', 'ratio' of time of Tcold to time of Thot. A small shift is incorporated to move the falling flank from the beginning to the end
	function therm_mod(t, shift, PM, ratio, Tcold, Thot; flank=20) 
		width = (1-ratio)*PM
		tmod = mod(t+shift-width/2, PM)
		tstart = ratio*PM
		return smooth_rectangle.(tmod+2/3*width, tstart, width, Tcold, Thot; flank=flank)
	end

end

# ╔═╡ 1f5c1ca8-7b61-4dc4-91e6-9fb187073313
# del
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
	Plots.plot!(_t, therm_mod.(_t, 0.0, 5, 0.5, -30, 30))
	Plots.plot!(_t, therm_mod.(_t, 0.5, 5, 0.5, -30, 30))
	Plots.plot!(_t, therm_mod.(_t, -5/6*0.5, 5, 0.5, -30, 30))

	Plots.plot!(_t, therm_mod.(_t, 0.0, 6, 0.8, -60, -40))
	Plots.plot!(_t, therm_mod.(_t, 0.0, 6, 0.8, -60, -40; flank=10))
	Plots.plot!(_t, therm_mod.(_t, 0.0, 6, 0.8, -60, -40; flank=5))
	Plots.plot!(_t, therm_mod.(_t, 0.0, 6, 0.8, -60, -40; flank=7))
end

# ╔═╡ d0adda7c-8980-4d74-8b79-c3c5edd8131f
md"""
## `ModuleTM`
"""

# ╔═╡ ab33d65f-f0ba-4e5d-b660-2ecd106d8d21
# copied to GasChromatographySystems
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
# copied to GasChromatographySystems
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
# copied to GasChromatographySystems
function GCxGC_TM(L1, d1, df1, sp1, TP1, L2, d2, df2, sp2, TP2, LTL, dTL, dfTL, spTL, TPTL, LM::Array{Float64,1}, dM, dfM, spM, shiftM1, shiftM2, PM, ratioM, HotM, ColdM, TPM, F, pin, pout; opt=GasChromatographySystems.Options(ng=true))

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
	modules[3] = GasChromatographySystems.ModuleTM("TM1", LM[2], dM*1e-3, dfM*1e-6, spM, TPM, shiftM1, PM, ratioM, HotM, ColdM, NaN)
	modules[4] = GasChromatographySystems.ModuleColumn("mod loop", LM[3], dM*1e-3, dfM*1e-6, spM, TPM)
	modules[5] = GasChromatographySystems.ModuleTM("TM2", LM[4], dM*1e-3, dfM*1e-6, spM, TPM, shiftM2, PM, ratioM, HotM, ColdM, NaN)
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
	sys = GCxGC_TM(31.5, 0.25, 0.25, "ZB1ms", GCxGC_TP, 0.56, 0.1, 0.1, "Stabilwax", GCxGC_TP, 0.1, 0.1, 0.1, "Stabilwax", 280.0, [0.30, 0.01, 0.90, 0.01, 0.30], 0.1, 0.1, "Stabilwax", 0.0, 0.0, PM, coldjet/PM, 25.0, -120.0, GCxGC_TP, 0.8, NaN, 0.0)
#	sys = GCxGC_TM(30.6, 0.25, 0.25, "ZB1ms", GCxGC_TP, 0.72, 0.1, 0.1, "Stabilwax", GCxGC_TP, 0.235, 0.1, 0.1, "Stabilwax", 280.0, [0.30, 0.005, 0.90, 0.005, 0.30], 0.1, 0.1, "Stabilwax", 0.0, 0.0, PM, coldjet/PM, 25.0, -120.0, GCxGC_TP, 0.8, NaN, 0.0)
end

# ╔═╡ ca604ac6-fdb1-41d7-83d2-a0c9b2440dad
GasChromatographySystems.plot_graph(sys)

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

# ╔═╡ ac60e9eb-0dfb-4a8d-a6bb-5289110b40cd
# copied to GasChromatographySystems
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

	function module_temperature(module_::GasChromatographySystems.ModuleTM, sys; Tcold_abs=true, spatial=true)
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
		
		spot(x,t) = if spatial == true
			GasChromatographySystems.smooth_rectangle(x, 0.0, sys.modules[5].length, T_itp_(x, t), T_itp(x,t))
		else
			T_itp(x,t)
		end
		return time_steps, temp_steps, gf, a_gf, spot
	end
end

# ╔═╡ 54c85e34-941d-41a1-9780-cf156f82d765
# copied to GasChromatographySystems
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

# ╔═╡ 607265de-c784-4ed0-8e89-ebb206785b0b
begin
	plotly()
	TM_itp = module_temperature(sys.modules[3], sys)[5]
	Plots.plot(0.0:0.01:sys.modules[3].PM, TM_itp.(sys.modules[3].length/2, 0.0:0.01:sys.modules[3].PM).-273.115, title="temperature at modulator point", xlabel="time in s", ylabel="temperature in °C", label="")
end

# ╔═╡ dde2bb43-fb5e-4d5b-bf88-8dee81e69063
function module_temperature_(module_::GasChromatographySystems.ModuleTM, sys; Tcold_abs=true, spatial=true, flank=40)
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
		
		spot(x,t) = if spatial == true
			GasChromatographySystems.smooth_rectangle(x, 0.0, sys.modules[5].length, T_itp_(x, t), T_itp(x,t), flank=flank)
		else
			T_itp(x,t)
		end
		return time_steps, temp_steps, gf, a_gf, spot
	end

# ╔═╡ 804e57b0-8879-4888-b4b1-682ead7df322
TT = module_temperature_(sys.modules[5], sys; Tcold_abs=false, spatial=true, flank=60)[end]

# ╔═╡ ac13726f-5e4b-43eb-a45c-7fc3e628f005
begin
	Plots.plot(0.0:sys.modules[5].length/100.0:sys.modules[5].length, TT.(0.0:sys.modules[5].length/100.0:sys.modules[5].length, 0.0).-273.15, label="t=0.0s")
	Plots.plot!(0.0:sys.modules[5].length/100.0:sys.modules[5].length, TT.(0.0:sys.modules[5].length/100.0:sys.modules[5].length, 3.5).-273.15, label="t=3.5s")
	Plots.plot!(0.0:sys.modules[5].length/100.0:sys.modules[5].length, TT.(0.0:sys.modules[5].length/100.0:sys.modules[5].length, 3.6).-273.15, label="t=3.6s")
	Plots.plot!(0.0:sys.modules[5].length/100.0:sys.modules[5].length, TT.(0.0:sys.modules[5].length/100.0:sys.modules[5].length, 3.7).-273.15, label="t=3.7s")
	Plots.plot!(0.0:sys.modules[5].length/100.0:sys.modules[5].length, TT.(0.0:sys.modules[5].length/100.0:sys.modules[5].length, 3.8).-273.15, label="t=3.8s")
	Plots.plot!(0.0:sys.modules[5].length/100.0:sys.modules[5].length, TT.(0.0:sys.modules[5].length/100.0:sys.modules[5].length, 3.9).-273.15, label="t=3.9s")
	Plots.plot!(0.0:sys.modules[5].length/100.0:sys.modules[5].length, TT.(0.0:sys.modules[5].length/100.0:sys.modules[5].length, 4.0).-273.15, label="t=4.0s")
end

# ╔═╡ 7277a646-fd70-4032-9691-569eddeafec6
# copied to GasChromatographySystems
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
		time_steps, temp_steps, gf, a_gf, T_itp = module_temperature(sys.modules[i], sys; Tcold_abs=Tcold_abs)
		prog = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)

		# substance parameters
		sub = GasChromatographySimulator.load_solute_database(db_dataframe, sys.modules[i].stationary_phase, sys.options.mobile_phase, selected_solutes, NaN.*ones(length(selected_solutes)), NaN.*ones(length(selected_solutes)))

		# option parameters
		opt = if typeof(sys.modules[i]) == GasChromatographySystems.ModuleTM
			GasChromatographySimulator.Options(alg=sys.options.alg, abstol=sys.options.abstol, reltol=sys.options.reltol, Tcontrol=sys.options.Tcontrol, odesys=sys.options.odesys, ng=false, vis=sys.options.vis, control=sys.options.control, k_th=sys.options.k_th)
		else
			GasChromatographySimulator.Options(alg=sys.options.alg, abstol=sys.options.abstol, reltol=sys.options.reltol, Tcontrol=sys.options.Tcontrol, odesys=sys.options.odesys, ng=sys.options.ng, vis=sys.options.vis, control=sys.options.control, k_th=sys.options.k_th)
		end

		parameters[i] = GasChromatographySimulator.Parameters(col, prog, sub, opt)
	end
	return parameters
end

# ╔═╡ 793560f7-d364-4c68-81ec-994441a41059
par = graph_to_parameters(sys, db, selected_solutes; Tcold_abs=false)

# ╔═╡ 6993afed-4519-4c8a-9cfc-87c72723c444
begin
	gr()#plotly()
	p_TM_Prog = Plots.plot(0.0:0.01:800.0, par[3].prog.T_itp.(sys.modules[3].length/2, 0.0:0.01:800.0).-273.15)
	Plots.plot!(p_TM_Prog, 0.0:0.01:800.0, par[4].prog.T_itp.(0.0, 0.0:0.01:800.0).-273.15)
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
### Simulation using the new function

The new function will reproduce the step-by-step simulation.
"""

# ╔═╡ 8b79cb42-f2b9-42b2-afb7-4e7fd98d2807
paths_ = [paths[1][1:8]]

# ╔═╡ ef88c521-3802-4b8f-a1c1-4c85c249770f
par

# ╔═╡ 1f8538e8-fc1f-44c8-ab14-9f3a3f01db29
# after GC1

# ╔═╡ ad3a2208-69e3-4cb0-8b00-c155392b2a3d
sim[2][1][1]

# ╔═╡ ff46e4e4-63e0-444f-84fc-f5539e867a70
# after mod in

# ╔═╡ f38735a7-9672-4807-8bba-221ea8a7c723
sim[2][1][2]

# ╔═╡ f05077a4-e287-4be4-a217-28a8832fd7e3
# after TM1

# ╔═╡ 25dc6864-d493-4d2b-98df-c62ced2848a0
# after loop

# ╔═╡ ea7a4510-609e-4bff-9601-90f2e9652aba
# after TM2

# ╔═╡ 0bc059a1-2b66-44f0-ab25-cef2482e87af
# after mod out

# ╔═╡ bb77b0ce-65a1-47fc-9e57-058c242dd8e0
# after GC2

# ╔═╡ 60a48a4a-510f-418e-802d-7981c0869ef3
# after TL

# ╔═╡ 3e8d2e43-1fa9-40e4-a132-7228664b7546
sys.modules[3]

# ╔═╡ 049eb024-9002-4f74-9ede-cbf19ba987e7
md"""
#### Test simulation of system with other, simple systems
"""

# ╔═╡ e7866add-cfce-4663-9352-a8353e4474a2
ex_series = GasChromatographySystems.SeriesSystem([10.0, 5.0, 2.0, 1.0], [0.53, 0.32, 0.25, 0.1], [0.53, 0.32, 0.25, 0.1], ["ZB1ms", "ZB1ms", "Stabilwax", "ZB1ms"], [GasChromatographySystems.default_TP(), GasChromatographySystems.default_TP(), GasChromatographySystems.default_TP(), GasChromatographySystems.default_TP()], NaN, 300.0, 0.0)

# ╔═╡ af9b5e44-1afc-45fa-a574-5dbcf6abd0d8
par_series = graph_to_parameters(ex_series, db, selected_solutes) 

# ╔═╡ 8b19420d-00e9-4606-aef4-47df49bcacf5
ex_split = GasChromatographySystems.SplitSystem([10.0, 1.0, 5.0], [0.25, 0.1, 0.25], [0.25, 0.0, 0.0], ["ZB1ms", "Stabilwax", ""], [GasChromatographySystems.default_TP(), 300.0, 300.0], [1.0, NaN, NaN], NaN, 0.0, 101.3)

# ╔═╡ e08de85a-8e1e-436d-876c-38ee0003aa8b
par_split = graph_to_parameters(ex_split, db, selected_solutes) 

# ╔═╡ 95dcd097-2842-4b96-9f3b-4b1dae33c35c
GasChromatographySystems.all_paths(ex_split.g, 2)[2]

# ╔═╡ 6e24a0f4-8d56-4798-99cb-d9e8ceabad9c
path = [GasChromatographySystems.all_paths(ex_split.g, 2)[2][2][1:2]]

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

# ╔═╡ 70597ee8-760d-4c21-af6c-981922208d10
md"""
### 1st Modulator spot
"""

# ╔═╡ 89e09993-ed50-480d-ab1d-deb745e3c43f
par[2]

# ╔═╡ b83330dd-4df3-433a-ab29-273e39f8b32c
begin
	# increase the peak width befor the modulator by 50%
#	sim[2][1][2][!,:τR] = 1.5.*sim[2][1][2][!,:τR]
	sim[2][1][2]
end

# ╔═╡ fc7f8b92-f0ee-4a89-ab7c-0c8f5ba1bb0f
begin
	p_chrom_before_mark = Plots.plot(t_before, c_before, xlabel="time in s", label="")
	for i=1:length(c_befores)
		Plots.plot!(p_chrom_before_mark, t_befores[i], c_befores[i], label=sim[2][1][2].Name[i])
	end
	Plots.plot!(p_chrom_before_mark, xticks=0:sys.modules[3].PM:3000)
	

	#add marker for PM and tcold, what about the shift?
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
#=begin
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
end=#

# ╔═╡ 3388fdf6-dff9-4d44-9db3-5bf77a57980d
# copied to GasChromatographySystems

# this splicing function is with separating focussed and unfocussed peak segments

# area as an input is needed 'AR'
# A_focussed calculated in relation to it
# A_focussed is the complete area during a modulation period
# not focussed segment is already included in the focussed segment, assuming it will be focussed in the 2nd modulation in a multi stage modulation
function slicing(tR, τR, AR, PM, ratio, shift, par::GasChromatographySimulator.Parameters, prev_par::GasChromatographySimulator.Parameters; nτ=6, τ₀=zeros(length(τR)), abstol=1e-8, reltol=1e-6, alg=OwrenZen5())
	tcold = PM*ratio
	thot = PM*(1-ratio)
	init_t_start = (fld.(tR.-nτ.*τR, PM)).*PM .+ shift # start time of the peaks, rounded to multiple of PM
	init_t_end = (fld.(tR.+nτ.*τR, PM)).*PM .+ shift # end time of the peaks, rounded to multiple of PM
	n_slice = Int.((init_t_end .- init_t_start)./PM .+ 1) # number of slices for every substance
	sub_TM_focussed = Array{GasChromatographySimulator.Substance}(undef, sum(n_slice))
	A_focussed = Array{Float64}(undef, sum(n_slice))
	A_unfocussed = Array{Float64}(undef, sum(n_slice))
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
			
			sub_TM_focussed[ii] = GasChromatographySimulator.Substance(prev_par.sub[i].name, prev_par.sub[i].CAS, prev_par.sub[i].Tchar, prev_par.sub[i].θchar, prev_par.sub[i].ΔCp, prev_par.sub[i].φ₀, "s$(j)_"*prev_par.sub[i].ann, prev_par.sub[i].Cag, t₀, τ₀[i])

			# Integrals:
			p = [tR[i], τR[i]]
			# approximated integrals
			prob_focussed = IntegralProblem(g, init_t_start[i]+(j-1)*PM, init_t_start[i]+(j-1)*PM+tcold, p)
			prob_unfocussed = IntegralProblem(g, init_t_start[i]+(j-1)*PM+tcold, init_t_start[i]+(j-1)*PM+tcold+thot, p)
			A_focussed[ii] = solve(prob_focussed, QuadGKJL(); reltol = 1e-8, abstol = 1e-14).u * AR[i]
			A_unfocussed[ii] = solve(prob_unfocussed, QuadGKJL(); reltol = 1e-8, abstol = 1e-14).u * AR[i]
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

	df_A_foc = DataFrame(Name=Name, CAS=CAS, Annotations=Ann_focussed, A=A_focussed+A_unfocussed, t0=t0_foc)
	
	return newpar_focussed, df_A_foc
end

# ╔═╡ b40149ac-a07f-4ec4-913a-a4f83bfc8169
CAS_par = [par[3].sub[i].CAS for i in 1:length(par[3].sub)]

# ╔═╡ 0b7251a0-6dfb-4ca6-990f-163e606df5f2
# copied to GasChromatographySystems
# fixed version
function slicing(pl, PM, ratio, shift, par::GasChromatographySimulator.Parameters; nτ=6, τ₀=zeros(length(pl.τR)), abstol=1e-8, reltol=1e-6, alg=OwrenZen5())
	tR = pl.tR
	τR = pl.τR
	AR = pl.A
	tcold = PM*ratio
	thot = PM*(1-ratio)
	init_t_start = (fld.(tR.-nτ.*τR, PM)).*PM .+ shift # start time of the peaks, rounded to multiple of PM
	init_t_end = (fld.(tR.+nτ.*τR, PM)).*PM .+ shift # end time of the peaks, rounded to multiple of PM
	n_slice = Int.((init_t_end .- init_t_start)./PM .+ 1) # number of slices for every substance
	sub_TM_focussed = Array{GasChromatographySimulator.Substance}(undef, sum(n_slice))
	A_focussed = Array{Float64}(undef, sum(n_slice))
	A_unfocussed = Array{Float64}(undef, sum(n_slice))
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
			prob_focussed = IntegralProblem(g, init_t_start[i]+(j-1)*PM, init_t_start[i]+(j-1)*PM+tcold, p)
			prob_unfocussed = IntegralProblem(g, init_t_start[i]+(j-1)*PM+tcold, init_t_start[i]+(j-1)*PM+tcold+thot, p)
			A_focussed[ii] = solve(prob_focussed, QuadGKJL(); reltol = 1e-8, abstol = 1e-30).u * AR[i]
			A_unfocussed[ii] = solve(prob_unfocussed, QuadGKJL(); reltol = 1e-8, abstol = 1e-30).u * AR[i]
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

	df_A_foc = DataFrame(Name=Name, CAS=CAS, Annotations=Ann_focussed, A=A_focussed+A_unfocussed, t0=t0_foc)
	
	return newpar_focussed, df_A_foc
end

# ╔═╡ 76f2e4ad-638a-4232-aa09-ba7d42da79c3
slicing(sim[2][1][2], sys.modules[3].PM, sys.modules[3].ratio, sys.modules[3].shift, par[3]; nτ=6, abstol=1e-10, reltol=1e-8, alg=Vern9())

# ╔═╡ e4172084-2dd7-4daf-b2b9-cc84891c06f2
par_1st_TM = slicing(sim[2][1][2], sys.modules[3].PM, sys.modules[3].ratio, sys.modules[3].shift, par[3]; nτ=6, abstol=1e-10, reltol=1e-8, alg=Vern9()) # competly slicing up the peaks at the modulator

# ╔═╡ 3c22e8aa-376d-4eec-9650-3bb102d4a0b8
# copied to GasChromatographySystems
# identfy the common index of a substances CAS number and Annotations in a peaklist 
function common_index(pl, CAS, Annotation)
	ii_CAS = findall(CAS.==pl.CAS)
	ii_ann = findall(occursin.(Annotation, pl.Annotations))
	ii = intersect(ii_CAS, ii_ann)[1]
	return ii
end

# ╔═╡ 846c3cab-90d4-4e73-9ce5-c1b41f937cdf
# copied to GasChromatographySystems
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

# ╔═╡ 4418401d-c418-422c-a341-7f01e7a8a3dd
function add_t0_to_pl!(pl, df_t0)
	# add the initial time in the dataframe df_t0 to the peaklist, same CAS and Annotations are used to identify the correct solute
	sort_t0 = Array{Float64}(undef, length(df_t0.t0))
	for i=1:length(df_t0.t0)
		ii = common_index(pl, df_t0.CAS[i], join(split(df_t0.Annotations[i], "_")[1:end-1], "_"))
		#ii = common_index(pl, df_A.CAS[i], split(df_A.Annotations[i], ", ")[1])
		sort_t0[ii] = df_t0.t0[i]
	end
	pl[!, :t0] = sort_t0
	return pl
end

# ╔═╡ f436afe9-0440-4149-be61-b426e243b34d
begin
	sol_1st_TM = GasChromatographySimulator.simulate(par_1st_TM[1])
	
	add_A_to_pl!(sol_1st_TM[1], par_1st_TM[2])
	add_t0_to_pl!(sol_1st_TM[1], par_1st_TM[2])
end

# ╔═╡ ceb6ba3f-aece-4ee3-b757-5904fa04c7ac
sol_1st_TM[1]

# ╔═╡ 47ddd70c-681c-44dc-aee2-4b1a38eb64f6
md"""
#### Sliced, not focussed rectangular peaks
"""

# ╔═╡ 82c44295-495e-4033-a669-aa6868fc4bc0
#sol_1st_TM_unfoc = DataFrame(Name=sol_1st_TM[1].Name, CAS=sol_1st_TM[1].CAS, tR=sol_1st_TM[1].t0.+sys.modules[3].ratio*sys.modules[3].PM, τR=(1-sys.modules[3].ratio)*sys.modules[3].PM, A=sol_1st_TM[1].A_unfoc, Annotations=sol_1st_TM[1].Annotations)

# ╔═╡ 3c20b2b4-ca9a-4b90-abb4-4cc6fd25a733
md"""
#### Misc
"""

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

# ╔═╡ 9c5ec167-212d-413a-a2bc-fb00d015af51
begin
	plotly()
	p_Tt = Plots.plot(xlabel="time in s", ylabel="temperature in °C", legend=false)
	for i=1:length(sol_1st_TM[2])
		trace = traces(sol_1st_TM[2], par_1st_TM[1], i)
		Plots.plot!(p_Tt, trace.t, trace.T.-273.15)
	end
	p_Tt
end

# ╔═╡ 740e977e-4cf8-496d-bb1a-000871bc4211
begin
	p_τt = Plots.plot(xlabel="time in s", ylabel="τ in s", legend=false)
	for i=1:length(sol_1st_TM[2])
		trace = traces(sol_1st_TM[2], par_1st_TM[1], i)
		Plots.plot!(p_τt, trace.t, sqrt.(trace.τ²))
	end
	p_τt
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
function chrom_after_modulation(pl, par; nτ=6)
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
	t_sum = collect(t1:minimum(pl.τR)/10:t2)
	c_sum = fill(0.0, length(t_sum))
	for i=1:length(pl.Name)
		if isnan(tstart[i])
			c[i] = fill(NaN,length(t_sum))
		elseif "A_foc" in names(pl)
			c[i] = GasChromatographySimulator.chromatogram(t_sum, [pl.tR[i]], [pl.τR[i]])*pl.A_foc[i]
			c_sum = c_sum .+ c[i]
		else
			c[i] = GasChromatographySimulator.chromatogram(t_sum, [pl.tR[i]], [pl.τR[i]])*pl.A[i]
			c_sum = c_sum .+ c[i]
		end
	end
	return t_sum, c, c_sum
end

# ╔═╡ 3c97b7aa-64a0-43ea-938d-149dbbd9cc63
"A_foc" in names(sol_1st_TM[1])

# ╔═╡ 095ee9ba-d974-4625-a47c-3634d91c6a5d
begin
	gr()
	# focussed peaks
	t_TM1, c_TM1s, c_TM1 = chrom_after_modulation(sol_1st_TM[1], par_1st_TM[1])
	p_chrom_TM1 = Plots.plot(t_TM1, c_TM1, xlabel="time in s")
	# unfocussed peaks 
#	t_unfoc, c_unfoc, h_unfoc, t_unfoc_sum, c_unfoc_sum = rec_chrom_after_modulation(sol_1st_TM_unfoc)
	#for i=1:length(c_TM1s)
	#	Plots.plot!(p_chrom_TM1, t_TM1, c_TM1s[i], label=sol_1st_TM[1].Name[i])
	#end
	Plots.plot!(p_chrom_TM1, xticks=0:4:3000, legend=false)
#	Plots.plot!(p_chrom_TM1, t_unfoc_sum, c_unfoc_sum)
#	Plots.plot!(p_chrom_TM1)#, xlims=(1899.5, 1900.2), ylims=(-0.01, 0.5))
	p_chrom_TM1
end

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

# ╔═╡ f6093dce-30a7-4628-8362-0c5ed8c55ebc
md"""
### Loop between Modulator spots
"""

# ╔═╡ beb15214-b1bb-4b39-a2aa-86f0727bc2c4
md"""
#### Sliced focused peaks
"""

# ╔═╡ 6a3f4043-63ba-44ef-b9ce-a391bc697fe6
# copied to GasChromatographySystems
function change_initial(par_, prev_par, prev_pl)
# change initials, using the substance parameters from the previous parameter set (parameter set of the previous segment/module, should have the same stationary phase!)
# this version is needed for an increase of the number of substances due to slicing at a modulator relative the the 
	new_sub = Array{GasChromatographySimulator.Substance}(undef, length(prev_par.sub))
	for i=1:length(prev_par.sub)
		# filter for correct CAS and annotation (slice number)
		ii_ = common_index(prev_pl, prev_par.sub[i].CAS, join(split(prev_par.sub[i].ann, ", ")[1:end-1], ", "))
		new_sub[ii_] = GasChromatographySimulator.Substance(prev_par.sub[i].name, prev_par.sub[i].CAS, prev_par.sub[i].Tchar, prev_par.sub[i].θchar, prev_par.sub[i].ΔCp, prev_par.sub[i].φ₀, prev_par.sub[i].ann, prev_par.sub[i].Cag, prev_pl.tR[ii_], prev_pl.τR[ii_])
	end
	# here changes of options could be applied
	new_par = GasChromatographySimulator.Parameters(par_.col, par_.prog, new_sub, par_.opt)
	return new_par
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

# ╔═╡ 66eaa639-416c-4167-9e00-14dcbf583a4c
md"""
#### Chromatogram after loop
"""

# ╔═╡ 82d9ac61-136c-4503-8031-baf0cbf4e548
#sol_loop_unfoc = DataFrame(Name=sol_loop[1].Name, CAS=sol_loop[1].CAS, tR=sol_loop[1].tR, τR=(1-sys.modules[3].ratio)*sys.modules[3].PM, A=sol_loop[1].A_unfoc, Annotations=sol_loop[1].Annotations)

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

# ╔═╡ 2390b6a8-8927-466d-a392-f76c4cd2c7f2
function warning_2nd_TM(ΔtR_loop, τR_foc_loop, PM, ratio)
	if all(ΔtR_loop .- nτ .* τR_foc .> (1-sys.modules[5].ratio)*sys.modules[5].PM) == false
		@warn "Time in loop is shorter than time of hot jet. Breaktrough in the modulator possible."
	elseif all(ΔtR_loop .+ 1.05*(1-sys.modules[5].ratio)*sys.modules[5].PM .< sys.modules[5].ratio*sys.modules[5].PM) == false
		@warn "Time in loop (plus width of unfocussed part) is longer than time of cold jet. Breaktrough in the modulator possible."
	end
end

# ╔═╡ db04ef60-3490-49c3-a5c1-f198c4015b06
md"""
#### Chromatogram after 2nd TM
"""

# ╔═╡ ebc1b905-eef4-4c66-abb6-b070d19f17d9
p_chrom_TM1

# ╔═╡ 231edaa6-8ef3-4a85-b48e-e0e1829d124d
#savefig(p_chrom_TM2, "Chrom_after_2nd_TM_new.svg")

# ╔═╡ 64e635f8-f5fd-4935-b44e-b02ebdd96b64
md"""
### Connection to GC2
"""

# ╔═╡ b3b66a27-cea8-49a8-8e05-2d75178b954a
# fixed function
# copied to GasChromatographySystems
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

# ╔═╡ bac5d026-0b84-412f-98bb-8ddedb2b92d9
# fixed version
# copied to GasChromatographySystems
function simulate_along_paths(sys, paths, par_sys; t₀=zeros(length(par_sys[1].sub)), τ₀=zeros(length(par_sys[1].sub)), nτ=6, abstol=1e-10, reltol=1e-8, refocus=falses(ne(sys.g)), τ₀_focus=zeros(length(par_sys[1].sub)))
	
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
						
						new_par_sys[i_par[j]], df_A = slicing(peaklists_[j-1], sys.modules[i_par[j]].PM, sys.modules[i_par[j]].ratio, sys.modules[i_par[j]].shift, par_sys[i_par[j]]; nτ=nτ, τ₀=τ₀, abstol=abstol, reltol=reltol, alg=Vern9())
						
						peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(new_par_sys[i_par[j]])
						add_A_to_pl!(peaklists_[j], df_A)
					else # ModuleColumn
						#if length(par_sys[i_par[j]].sub) < length(new_par_sys[i_par[j-1]].sub)
							# number of peaks increased previously due to slicing in a modulator
							new_par_sys[i_par[j]] = change_initial(par_sys[i_par[j]], peaklists_[j-1])
						#else
						#	new_par_sys[i_par[j]] = GasChromatographySystems.change_initial(par_sys[i_par[j]], peaklists_[j-1])
						#end
						peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(new_par_sys[i_par[j]])
						add_A_to_pl!(peaklists_[j], peaklists_[j-1])
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
sim_ = simulate_along_paths(sys, paths_, par)

# ╔═╡ c667f7e8-d588-451e-b0c6-3422aa140324
sim_[4]

# ╔═╡ 1f62b350-332e-4acf-b907-75362157e4da
sim_[2][1][1]

# ╔═╡ 3e6687f5-a81e-403a-a877-179fbc6e923f
sim_[2][1][2]

# ╔═╡ f1bf21dd-07ee-4e4e-beb5-7d9225b3fc6b
sim_[2][1][3]

# ╔═╡ 33c2bf36-cf9d-4200-ab0b-71f80ce60a5b
sim_[2][1][4]

# ╔═╡ 843f5cba-a29b-473c-85aa-5f1ce3d3138f
sim_[2][1][5]

# ╔═╡ 1a113ae2-70fa-430e-9ff7-b58efadc9da9
sim_[2][1][6]

# ╔═╡ 3c6520f7-847a-4115-b17b-b2f980165446
sim_[2][1][7]

# ╔═╡ 9b5a69f5-79a5-4704-823c-90c916718b24
sim_[2][1][8]

# ╔═╡ 10a84b28-87f2-47dd-8d4e-feff196b5df5
sim_series = simulate_along_paths(ex_series, GasChromatographySystems.all_paths(ex_series.g, 1)[2], par_series)

# ╔═╡ 7a6772ac-1ce0-469b-966b-59003dee4bf1
sim_series[2]

# ╔═╡ 805dd021-f381-4be6-9236-ee86500c0af2
sim_split = simulate_along_paths(ex_split, GasChromatographySystems.all_paths(ex_split.g, 2)[2], par_split)

# ╔═╡ 2c1fc71c-5f5d-4e18-84f5-42980d461753
newpar_loop = change_initial(par[4], par_1st_TM[1], sol_1st_TM[1])

# ╔═╡ 74121702-9f4a-410b-961b-7865ab3af941
begin
	sol_loop = GasChromatographySimulator.simulate(newpar_loop)
	add_A_to_pl!(sol_loop[1], sol_1st_TM[1])
end

# ╔═╡ d5e17234-fc50-4863-b916-e9ad7df94cf9
sol_loop[1]

# ╔═╡ 96200092-8afd-4ce3-8287-79806a331656
# time in loop
ΔtR_loop = sol_loop[1].tR .- sol_1st_TM[1].tR

# ╔═╡ bb483a0e-72d7-405e-bf47-c5fa878e0994
begin
	plotly()
	
	# focussed peaks
	t_loop, c_loops, c_loop = chrom_after_loop(sol_loop[1])
	# unfocussed peaks 
#	t_loop_unfoc, c_loop_unfoc, h_loop_unfoc, t_loop_unfoc_sum, c_loop_unfoc_sum = rec_chrom_after_modulation(sol_loop_unfoc)
	
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
	
#	Plots.plot!(p_chrom_loop, t_loop_unfoc_sum, c_loop_unfoc_sum, linestyle=:dash, c=1)
#	for i=1:length(c_loop_unfoc)
#		if sol_loop_unfoc.Name[i] == selected_solutes[1]
#			color = :blue
#		elseif sol_loop_unfoc.Name[i] == selected_solutes[2]
#			color = :red
#		else
#			color = :orange
#		end
#		Plots.plot!(p_chrom_loop, t_loop_unfoc[i], c_loop_unfoc[i], label=sol_loop_unfoc.Name[i], c=color, linestyle=:dash)
#	end
	p_chrom_loop
end

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
	
#	Plots.plot!(p_chrom_loop_mark, t_loop_unfoc_sum, c_loop_unfoc_sum, linestyle=:dash, c=1)
#	for i=1:length(c_loop_unfoc)
#		
#		Plots.plot!(p_chrom_loop_mark, t_loop_unfoc[i], c_loop_unfoc[i], label=sol_loop_unfoc.Name[i], c=dot_color(sol_loop_unfoc.Name[i], selected_solutes), linestyle=:dash)
#	end

	#add marker for PM and tcold
	n = unique(fld.(t_loop, sys.modules[5].PM))
	for i=1:length(n)
		Plots.plot!(p_chrom_loop_mark, [n[i]*sys.modules[5].PM, n[i]*sys.modules[5].PM], [0.0, maximum(c_loop)], c=:black)
		Plots.plot!(p_chrom_loop_mark, [n[i]*sys.modules[5].PM+sys.modules[5].PM*sys.modules[5].ratio, n[i]*sys.modules[5].PM+sys.modules[5].PM*sys.modules[5].ratio], [0.0, maximum(c_loop)], c=:black, linestyle=:dash)
	end
	p_chrom_loop_mark
end

# ╔═╡ 1cb5a03e-b438-4afa-bc4d-bb66e4d5511e
par_2nd_TM_ = slicing(sol_loop[1], sys.modules[5].PM, sys.modules[5].ratio, sys.modules[5].shift, par[5]; nτ=6, abstol=1e-10, reltol=1e-8, alg=Vern9()) # competly slicing up the peaks at the modulator

# ╔═╡ 1e08231c-c890-4bcf-8b79-638918602634
par_2nd_TM_[2]

# ╔═╡ e78ffe7d-3f9f-471c-925c-4feae33231b8
begin
	sol_2nd_TM = GasChromatographySimulator.simulate(par_2nd_TM_[1])
	add_A_to_pl!(sol_2nd_TM[1], par_2nd_TM_[2])
	add_t0_to_pl!(sol_2nd_TM[1], par_2nd_TM_[2])
	sol_2nd_TM[1]
end

# ╔═╡ a821704e-af96-4431-a231-479a9daf2609
sol_2nd_TM[1]

# ╔═╡ 6f5421fb-97e4-43fc-9fc0-e5486b5958a9
maximum(diff(sol_2nd_TM[2][1].t))

# ╔═╡ 189ae09d-1b47-44ac-b657-b12bfa694912
sol_2nd_TM[1]

# ╔═╡ 081c75dc-4207-4ef7-aada-472f8908df47
begin
	plotly()
	p_τt_ = Plots.plot(xlabel="time in s", ylabel="τ in s", legend=false)
	for i=1:length(sol_2nd_TM[2])
		trace = traces(sol_2nd_TM[2], par_2nd_TM_[1], i)
		Plots.plot!(p_τt_, trace.t, sqrt.(trace.τ²))
	end
	p_τt_
end

# ╔═╡ 672e5876-0f0d-483a-9936-604b3fc7f97a
begin
	plotly()
	# focussed peaks
	t_TM2, c_TM2s, c_TM2 = chrom_after_modulation(sol_2nd_TM[1], par_2nd_TM_[1])
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

# ╔═╡ 78dbaf44-0468-417a-8c69-90c139e2e2ad
newpar_mod_out = change_initial(par[6], sol_2nd_TM[1])

# ╔═╡ 0400ee95-e5f2-48f3-bc57-8f3eb24a2e21
sol_2nd_TM[1]

# ╔═╡ bb01e67d-fbca-461e-b592-508bc5480b94
begin
	sol_mod_out = GasChromatographySimulator.simulate(newpar_mod_out)
	add_A_to_pl!(sol_mod_out[1], sol_2nd_TM[1])
end

# ╔═╡ b8960aa6-44d7-492f-b04f-1b664cdcc3a6
sol_mod_out[1]

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
newpar_gc2 = change_initial(par[7], newpar_mod_out, sol_mod_out[1])

# ╔═╡ caef921f-f06c-471b-a56a-8415ac2d4fa6
begin
	sol_gc2 = GasChromatographySimulator.simulate(newpar_gc2)
	add_A_to_pl!(sol_gc2[1], sol_mod_out[1])
end 

# ╔═╡ a7c7f2a5-1f70-4d94-91e0-7df6e532981e
sol_gc2[1]

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
newpar_tl = change_initial(par[8], newpar_gc2, sol_gc2[1])

# ╔═╡ 32e40a19-ebab-47a6-b4ca-9300c2b04dab
begin
	sol_tl = GasChromatographySimulator.simulate(newpar_tl)
	add_A_to_pl!(sol_tl[1], sol_gc2[1])
end 

# ╔═╡ f7d37712-b63f-4e6d-951e-5fcc555a7dc7
sol_tl[1]

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
	pl_GCxGC = DataFrame(Name=name, tR1=tR1s.-tR2s, tR2=tR2s)
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
# correct?
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

# ╔═╡ 17c4fbd1-4a9d-4f7a-922a-a72eb97db33a
md"""
# End
"""

# ╔═╡ Cell order:
# ╠═23c71f14-efdf-11ed-33f2-95304516114d
# ╠═95cc1f3d-eff3-4cc8-8ba1-7ee04e887b84
# ╠═d7f045e5-47c7-44b6-8bce-ce9dc3154fff
# ╠═2165e2d3-bca9-41f1-acee-a72a0ad9d898
# ╠═39ddb95f-40a1-4cbb-8317-3ef18fec4d0b
# ╠═d419921b-1472-4b66-bae4-469537259814
# ╠═dff75c18-eee5-4f70-9509-0e8228932819
# ╠═1f5c1ca8-7b61-4dc4-91e6-9fb187073313
# ╠═e44bcb5f-353e-4ef1-bf4d-d226a864ad48
# ╟─d0adda7c-8980-4d74-8b79-c3c5edd8131f
# ╠═ab33d65f-f0ba-4e5d-b660-2ecd106d8d21
# ╠═60d86bdd-5527-49da-8c35-ecd508171e6c
# ╠═b63ea500-c2bb-49c3-9915-032fd3e33e14
# ╠═7bb430cd-6bc5-4c48-85f8-3ffce3220d50
# ╠═80150c58-f9b9-43d9-bd7b-e38388f29fd8
# ╠═522a7582-0ff1-42d6-960a-adb952a225d2
# ╠═ca604ac6-fdb1-41d7-83d2-a0c9b2440dad
# ╠═54c85e34-941d-41a1-9780-cf156f82d765
# ╠═607265de-c784-4ed0-8e89-ebb206785b0b
# ╠═bdbd8cc9-29fb-43ce-850b-8d18980808d3
# ╠═4adbf7e9-232d-4889-92c6-75e8bde0d88d
# ╠═fc36ffe3-b1df-42c8-ad14-e3b90f27a272
# ╠═43f27a81-dca0-4188-baf4-3f8140d9e57c
# ╠═ac60e9eb-0dfb-4a8d-a6bb-5289110b40cd
# ╠═dde2bb43-fb5e-4d5b-bf88-8dee81e69063
# ╠═804e57b0-8879-4888-b4b1-682ead7df322
# ╠═ac13726f-5e4b-43eb-a45c-7fc3e628f005
# ╠═7277a646-fd70-4032-9691-569eddeafec6
# ╠═793560f7-d364-4c68-81ec-994441a41059
# ╠═6993afed-4519-4c8a-9cfc-87c72723c444
# ╠═a14d1ca2-46ef-4296-aa1d-06590eca97cf
# ╟─69598506-cf92-4626-b6e1-d8f57bf8388f
# ╠═1a0695c5-0a67-4158-a3bd-c135be31e408
# ╠═6556a858-65ef-49a6-bdaa-562a2b568a70
# ╠═8b79cb42-f2b9-42b2-afb7-4e7fd98d2807
# ╠═ef88c521-3802-4b8f-a1c1-4c85c249770f
# ╠═fcdc5015-31d0-4808-8218-f33901437c0b
# ╠═c667f7e8-d588-451e-b0c6-3422aa140324
# ╠═1f8538e8-fc1f-44c8-ab14-9f3a3f01db29
# ╠═1f62b350-332e-4acf-b907-75362157e4da
# ╠═ad3a2208-69e3-4cb0-8b00-c155392b2a3d
# ╠═ff46e4e4-63e0-444f-84fc-f5539e867a70
# ╠═3e6687f5-a81e-403a-a877-179fbc6e923f
# ╠═f38735a7-9672-4807-8bba-221ea8a7c723
# ╠═f05077a4-e287-4be4-a217-28a8832fd7e3
# ╠═f1bf21dd-07ee-4e4e-beb5-7d9225b3fc6b
# ╠═ceb6ba3f-aece-4ee3-b757-5904fa04c7ac
# ╠═25dc6864-d493-4d2b-98df-c62ced2848a0
# ╠═33c2bf36-cf9d-4200-ab0b-71f80ce60a5b
# ╠═d5e17234-fc50-4863-b916-e9ad7df94cf9
# ╠═ea7a4510-609e-4bff-9601-90f2e9652aba
# ╠═843f5cba-a29b-473c-85aa-5f1ce3d3138f
# ╠═a821704e-af96-4431-a231-479a9daf2609
# ╠═0bc059a1-2b66-44f0-ab25-cef2482e87af
# ╠═1a113ae2-70fa-430e-9ff7-b58efadc9da9
# ╠═b8960aa6-44d7-492f-b04f-1b664cdcc3a6
# ╠═bb77b0ce-65a1-47fc-9e57-058c242dd8e0
# ╠═3c6520f7-847a-4115-b17b-b2f980165446
# ╠═a7c7f2a5-1f70-4d94-91e0-7df6e532981e
# ╠═60a48a4a-510f-418e-802d-7981c0869ef3
# ╠═9b5a69f5-79a5-4704-823c-90c916718b24
# ╠═f7d37712-b63f-4e6d-951e-5fcc555a7dc7
# ╠═3e8d2e43-1fa9-40e4-a132-7228664b7546
# ╠═bac5d026-0b84-412f-98bb-8ddedb2b92d9
# ╠═049eb024-9002-4f74-9ede-cbf19ba987e7
# ╠═e7866add-cfce-4663-9352-a8353e4474a2
# ╠═af9b5e44-1afc-45fa-a574-5dbcf6abd0d8
# ╠═10a84b28-87f2-47dd-8d4e-feff196b5df5
# ╠═7a6772ac-1ce0-469b-966b-59003dee4bf1
# ╠═8b19420d-00e9-4606-aef4-47df49bcacf5
# ╠═e08de85a-8e1e-436d-876c-38ee0003aa8b
# ╠═95dcd097-2842-4b96-9f3b-4b1dae33c35c
# ╠═6e24a0f4-8d56-4798-99cb-d9e8ceabad9c
# ╠═805dd021-f381-4be6-9236-ee86500c0af2
# ╠═e3e34deb-9249-4811-b74a-44a4c2d06ac2
# ╠═54e20666-b217-415d-b2fd-305980588668
# ╠═ac3db64b-82c8-455a-97aa-6bbddcc03830
# ╠═ee2e3313-9d05-4815-97a2-5b292dc47b68
# ╠═82c69944-c4e0-4699-87c0-29500b33ba78
# ╠═70597ee8-760d-4c21-af6c-981922208d10
# ╠═89e09993-ed50-480d-ab1d-deb745e3c43f
# ╠═b83330dd-4df3-433a-ab29-273e39f8b32c
# ╠═fc7f8b92-f0ee-4a89-ab7c-0c8f5ba1bb0f
# ╠═00886bb1-39ff-48f2-9329-2a195f705a30
# ╠═325a8e2e-ae2f-4dd9-a32b-00c7adf47461
# ╠═3388fdf6-dff9-4d44-9db3-5bf77a57980d
# ╠═b40149ac-a07f-4ec4-913a-a4f83bfc8169
# ╠═0b7251a0-6dfb-4ca6-990f-163e606df5f2
# ╠═76f2e4ad-638a-4232-aa09-ba7d42da79c3
# ╠═e4172084-2dd7-4daf-b2b9-cc84891c06f2
# ╠═f436afe9-0440-4149-be61-b426e243b34d
# ╠═846c3cab-90d4-4e73-9ce5-c1b41f937cdf
# ╠═4418401d-c418-422c-a341-7f01e7a8a3dd
# ╠═3c22e8aa-376d-4eec-9650-3bb102d4a0b8
# ╠═47ddd70c-681c-44dc-aee2-4b1a38eb64f6
# ╠═82c44295-495e-4033-a669-aa6868fc4bc0
# ╠═3c20b2b4-ca9a-4b90-abb4-4cc6fd25a733
# ╠═9c5ec167-212d-413a-a2bc-fb00d015af51
# ╠═740e977e-4cf8-496d-bb1a-000871bc4211
# ╠═06c3c74b-4ced-405e-89be-a05092fab027
# ╠═67a4dff7-2634-4ae2-adb0-515807988b17
# ╠═84111543-21e2-4891-9050-0c102cd46a52
# ╠═9c6573e8-bc0f-4c40-9eea-9733726fbd76
# ╠═3c97b7aa-64a0-43ea-938d-149dbbd9cc63
# ╠═095ee9ba-d974-4625-a47c-3634d91c6a5d
# ╠═47bd51af-7f68-4967-aac2-c80e4ee73ff7
# ╠═4c2c3d26-951a-412f-833d-83b153b4354e
# ╠═f6093dce-30a7-4628-8362-0c5ed8c55ebc
# ╠═beb15214-b1bb-4b39-a2aa-86f0727bc2c4
# ╠═2c1fc71c-5f5d-4e18-84f5-42980d461753
# ╠═6a3f4043-63ba-44ef-b9ce-a391bc697fe6
# ╠═74121702-9f4a-410b-961b-7865ab3af941
# ╠═5d0cb969-8a98-4d28-aa5c-233723755f41
# ╠═96200092-8afd-4ce3-8287-79806a331656
# ╠═66eaa639-416c-4167-9e00-14dcbf583a4c
# ╠═82d9ac61-136c-4503-8031-baf0cbf4e548
# ╠═bb483a0e-72d7-405e-bf47-c5fa878e0994
# ╠═53aff8ef-4585-4d03-921b-7f19827dc04e
# ╠═b31fa530-922f-4862-b1cf-04d2ec34b1c1
# ╠═86274690-dfe8-4f71-8084-66fe5f03f932
# ╠═0d51d9c5-8a70-4d4d-b247-4c5fe07b7d1c
# ╠═2390b6a8-8927-466d-a392-f76c4cd2c7f2
# ╠═1cb5a03e-b438-4afa-bc4d-bb66e4d5511e
# ╠═1e08231c-c890-4bcf-8b79-638918602634
# ╠═e78ffe7d-3f9f-471c-925c-4feae33231b8
# ╠═6f5421fb-97e4-43fc-9fc0-e5486b5958a9
# ╠═081c75dc-4207-4ef7-aada-472f8908df47
# ╠═db04ef60-3490-49c3-a5c1-f198c4015b06
# ╠═ebc1b905-eef4-4c66-abb6-b070d19f17d9
# ╠═672e5876-0f0d-483a-9936-604b3fc7f97a
# ╠═231edaa6-8ef3-4a85-b48e-e0e1829d124d
# ╠═64e635f8-f5fd-4935-b44e-b02ebdd96b64
# ╠═189ae09d-1b47-44ac-b657-b12bfa694912
# ╠═b3b66a27-cea8-49a8-8e05-2d75178b954a
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
# ╠═17c4fbd1-4a9d-4f7a-922a-a72eb97db33a
