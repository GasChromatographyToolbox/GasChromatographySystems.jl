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
	TableOfContents()
end

# ╔═╡ d7f045e5-47c7-44b6-8bce-ce9dc3154fff
using Integrals

# ╔═╡ 927e9354-4130-47cd-bd94-8eb37da15b47
using Waveforms

# ╔═╡ 8e38d02b-1161-46f7-92ad-5d6ef1e8ede4
using ForwardDiff

# ╔═╡ 8ca0c042-7e2e-4248-a183-58235e909f59
using DifferentialEquations

# ╔═╡ d419921b-1472-4b66-bae4-469537259814
md"""
# Development of a thermal modulated GCxGC system in detail
"""

# ╔═╡ dff75c18-eee5-4f70-9509-0e8228932819
md"""
## Periodic rectangle temperature function
"""

# ╔═╡ b4a1fe1b-f44a-4844-b041-1c9eabcd7ce2
function rectangular(x, a)
	y = if abs(x) > a/2
		0
	elseif abs(x) == a/2
		1//2
	elseif abs(x) < a/2
		1
	end
	return y
end

# ╔═╡ cbe42d0d-e765-4a3f-b7dc-ba5687ac77e2
x = -10.0:0.1:10.0

# ╔═╡ 2777203c-c074-469c-9194-efa526676002
a = 1.0

# ╔═╡ cbc53a06-1efc-4543-9151-52f9e95e6ae3
begin
	plotly()
	Plots.plot(x, rectangular.(x, a))
end

# ╔═╡ b283ebdd-38f5-4b8c-a1b7-34476d21dff2
begin
	plotly()
	Plots.plot(x, cos.(x*π/a))
end

# ╔═╡ 5cc7905e-477f-46ee-816a-f1327f914fd8
begin
	plotly()
	Plots.plot(x, rectangular.(x, a).*cos.(x*π/a))
end

# ╔═╡ 58441de1-11cd-4d44-8664-fbc0b9039ed5
begin
	plotly()
	Plots.plot(x, cos.(x*π/a))
	Plots.plot!(x, rectangular.(x, a).*cos.(x*π/a))
end

# ╔═╡ b25f6d02-173a-4d17-adcf-c54ad854a871
begin
	plotly()
	Plots.plot(x, sin.(x*π/a))
	Plots.plot!(x, squarewave.(x*π/a))
end

# ╔═╡ 6e5eba86-f046-4cf9-b6a5-c2b0622bf1fc
def_TP = GasChromatographySystems.default_TP()

# ╔═╡ 5aa38a4e-d9fa-4a16-b6e1-c9ee2eb7f2d9
TP_itp = GasChromatographySimulator.temperature_interpolation(def_TP.timesteps, def_TP.temperaturesteps, def_TP.gradient_function, 1.0)

# ╔═╡ b1a8b34e-49f4-4678-91dc-6d60bcb743e3
ΔT = 80.0

# ╔═╡ 14865482-42d6-4d13-b738-54a1a67a65bc
begin # symetrical hot and cold jet
	Plots.plot(0.0:0.1:1800.0, (80.0.*squarewave.((0.0:0.1:1800.0).*π/10)).+TP_itp.(0.0,0.0:0.1:1800.0).-273.15, label="TM1") # 20s period
	Plots.plot!(0.0:0.1:1800.0, (80.0.*squarewave.(((0.0:0.1:1800.0).-10).*π/10)).+TP_itp.(0.0,0.0:0.1:1800.0).-273.15, label="TM2") # 20s period, shifted by half period
end

# ╔═╡ 735300f0-b41a-41ce-8a96-707078dcd808
Plots.plot(0.0:0.1:1800.0, (80.0.*squarewave1.((0.0:0.1:1800.0)./20)).+TP_itp.(0.0,0.0:0.1:1800.0).-273.15)# 20s period

# ╔═╡ 0080eaa8-89ef-4c50-b9e8-a6630d088532
# variable length of hot and cold jet?

# ╔═╡ 14f47b21-6d1f-4e46-9bcb-7240f61bc492
squarewave1_new(x, hw) = ifelse(mod(x, 1) < hw, 1.0, -1.0)

# ╔═╡ 9031c353-5dbf-4296-bacc-d3f6338e9d8e
Plots.plot(0.0:0.1:1800.0, (80.0.*squarewave1_new.((0.0:0.1:1800.0)./20, 0.2)))

# ╔═╡ 02299b66-0054-4b77-ad36-d7f0d4f6fd21
f(x) = squarewave1_new(x, 0.2)

# ╔═╡ 1231d383-532b-4f01-8975-f65f574def0c
∂f(x) = ForwardDiff.derivative(f, x)

# ╔═╡ 615d984f-2452-4234-b9be-7eb97f3de94a
x0 = 0.3

# ╔═╡ 5c056b8c-304c-4348-a49b-24ef86561418
f(x0), ∂f(x0)

# ╔═╡ f635daa6-5bc5-4831-8320-3abf4a9dcc14
Plots.plot(0.0:0.001:2.0, f.(0.0:0.001:2.0))

# ╔═╡ c0a92e37-2592-4cf5-a572-228def68b8a1
Plots.plot(0.0:0.001:2.0, ∂f.(0.0:0.001:2.0))

# ╔═╡ 264e80c1-d8f9-4513-b599-ac25a021f321
Plots.plot(0.0:0.1:1800.0, (80.0.*squarewave1_new.((0.0:0.1:1800.0)./20, 0.8)).+TP_itp.(0.0,0.0:0.1:1800.0).-273.15)

# ╔═╡ 3c545fbb-2732-4ac5-b628-7cfda7bc6770
ff(t) = 80.0.*squarewave1_new.(t./20, 0.8).+TP_itp.(0.0,t).-273.15

# ╔═╡ 2842aa10-55a8-4e53-b582-d0693a91d0c2
∂ff(t) = ForwardDiff.derivative(ff, t)

# ╔═╡ 7f19f5df-6626-48c7-b4d9-7eb897c12bde
fff(t) = TP_itp.(0.0,t)

# ╔═╡ b38300b5-9b16-45e9-9d5f-19e0e825349f
∂fff(t) = ForwardDiff.derivative(fff, t)

# ╔═╡ 3722291a-8f23-4bdd-a86d-8ae2056ee809
begin
	Plots.plot(0.0:0.1:1800.0, ∂ff.(0.0:0.1:1800.0))
	Plots.plot!(0.0:0.1:1800.0, ∂fff.(0.0:0.1:1800.0))
end

# ╔═╡ 324b1803-dab5-4b68-a9a9-385b562c1c3f
squarewave1_ratio(x, ratio) = ifelse(mod(x, 1) < ratio, 1.0, -1.0)

# ╔═╡ 326d2249-dfcc-40b6-b8c0-d9ab5ff88627
function therm_mod(t, shift, PM, ratio, Thot, Tcold) 
	return ifelse(mod(t+shift, PM) < ratio*PM, Thot, Tcold)
end

# ╔═╡ 673fc2e1-43c3-429c-a7df-4930c407bd8f
begin
	Plots.plot(0.0:0.01:100.0, therm_mod.(0.0:0.01:100.0, 0.0, 6.0, 0.5, 100.0, -80.0).+TP_itp.(0.0,0.0:0.01:100.0).-273.15)
	Plots.plot!(0.0:0.01:100.0, therm_mod.(0.0:0.01:100.0, 3.0, 6.0, 0.5, 100.0, -80.0).+TP_itp.(0.0,0.0:0.01:100.0).-273.15)
end

# ╔═╡ 33a88b81-a65b-4ba3-b148-7339a8ff97a9
begin
	Plots.plot(0.0:0.01:100.0, therm_mod.(0.0:0.01:100.0, 0.0, 6.0, 0.5, TP_itp.(0.0,0.0:0.01:100.0).+ 100.0.-273.15, -80.0))
	Plots.plot!(0.0:0.01:100.0, therm_mod.(0.0:0.01:100.0, 3.0, 6.0, 0.5, 100.0, -80.0).+TP_itp.(0.0,0.0:0.01:100.0).-273.15)
end

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

# ╔═╡ 10427844-f557-4178-807b-f7cd194548c3
ModuleTM("TM1", 0.02, 0.25e-3, 0.25e-6, "Wax", GasChromatographySystems.default_TP(), 0.0, 6.0, 0.5, 80.0, -80.0, NaN)

# ╔═╡ b63ea500-c2bb-49c3-9915-032fd3e33e14
GCxGC_TP = GasChromatographySystems.default_TP()

# ╔═╡ f5a4bfba-1379-4b4f-a694-79426695d254
tsteps = [0.0, 1800.0]

# ╔═╡ 60d86bdd-5527-49da-8c35-ecd508171e6c
md"""
## Definition System
"""

# ╔═╡ e44cf22a-4fbe-43b7-a3ca-2d17a41a3f43
# system with two cold jets and two hot jets  

# ╔═╡ 27d0f80e-52b9-4942-a519-fcd787fdfaae
GasChromatographySystems.Options(ng=true, abstol=1e-10, reltol=1e-7)

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

# ╔═╡ 3fde3b56-e5a4-4a28-95a6-0d8b36fe37d1
begin
	PM = 6.0
	g = SimpleDiGraph(8)
	add_edge!(g, 1, 2) # 1st-D GC
	add_edge!(g, 2, 3) # modulator
	add_edge!(g, 3, 4) # hot/cold 1
	add_edge!(g, 4, 5) # modulator
	add_edge!(g, 5, 6) # hot/cold 2 
	add_edge!(g, 6, 7) # modulator
	add_edge!(g, 7, 8) # 2nd-D GC
	#add_edge!(g2, 3, 4) # 2nd TL column inlet -> TL column -> Det
	# pressure points
	pp = Array{GasChromatographySystems.PressurePoint}(undef, nv(g))
	pp[1] = GasChromatographySystems.PressurePoint("p₁", tsteps, [NaN, NaN]) # inlet 
	for i=2:(nv(g)-1)
		pp[i] = GasChromatographySystems.PressurePoint("p$(i)", tsteps, [NaN, NaN])
	end
	pp[end] = GasChromatographySystems.PressurePoint("p$(nv(g))", tsteps, [eps(Float64), eps(Float64)]) # outlet
	# modules
	modules = Array{GasChromatographySystems.AbstractModule}(undef, ne(g))
	modules[1] = GasChromatographySystems.ModuleColumn("GC column 1", 30.0, 0.25e-3, 0.25e-6, "Rxi5ms", GCxGC_TP, 1.0/60e6)
	modules[2] = GasChromatographySystems.ModuleColumn("mod in", 0.1, 0.25e-3, 0.25e-6, "Rxi17SilMS", GCxGC_TP)
	modules[3] = ModuleTM("TM1", 0.005, 0.25e-3, 0.25e-6, "Rxi17SilMS", GCxGC_TP, 0.1*PM, PM, 0.1, 80.0, -80.0, NaN)
	modules[4] = GasChromatographySystems.ModuleColumn("mod loop", 0.1, 0.25e-3, 0.25e-6, "Rxi17SilMS", GCxGC_TP)
	modules[5] = ModuleTM("TM2", 0.005, 0.25e-3, 0.25e-6, "Rxi17SilMS", GCxGC_TP, 0.1*PM, PM, 0.1, 80.0, -80.0, NaN)
	modules[6] = GasChromatographySystems.ModuleColumn("mod out", 0.1, 0.25e-3, 0.25e-6, "Rxi17SilMS", GCxGC_TP)
	modules[7] = GasChromatographySystems.ModuleColumn("GC column 2", 2.0, 0.25e-3, 0.25e-6, "Rxi17SilMS", GCxGC_TP)
	# system
	sys = update_system(GasChromatographySystems.System(g, pp, modules, GasChromatographySystems.Options(ng=true, abstol=1e-12, reltol=1e-10)))
	#sys = update_system(GasChromatographySystems.System(g, pp, modules, GasChromatographySystems.Options(ng=true, abstol=1e-12, reltol=1e-10, alg=AutoVern9(Rodas5P()))))
	#sys = update_system(GasChromatographySystems.System(g, pp, modules, GasChromatographySystems.Options(ng=true, abstol=1e-12, reltol=1e-10, alg=Vern6())))
end

# ╔═╡ c7e771b6-7ec7-41bf-940f-2e2617af3c02
ModuleTM("TM2", 0.02, 0.25e-3, 0.25e-6, "Rxi17SilMS", GCxGC_TP, PM/2, PM, 0.5, 80.0, -80.0, NaN)

# ╔═╡ 80150c58-f9b9-43d9-bd7b-e38388f29fd8
function GCxGC_TM(L1, d1, df1, sp1, TP1, L2, d2, df2, sp2, TP2, LM::Array{Float64,1}, dM, dfM, spM, shiftM, PM, ratioM, HotM, ColdM, TPM, F, pin, pout; opt=GasChromatographySystems.Options(ng=true))

	TPs = [TP1, TP2, TPM]
	
	g = SimpleDiGraph(8)
	add_edge!(g, 1, 2) # 1st-D GC
	add_edge!(g, 2, 3) # modulator
	add_edge!(g, 3, 4) # hot/cold 1
	add_edge!(g, 4, 5) # modulator
	add_edge!(g, 5, 6) # hot/cold 2 
	add_edge!(g, 6, 7) # modulator
	add_edge!(g, 7, 8) # 2nd-D GC
	
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
	modules[3] = ModuleTM("TM1", LM[2], dM*1e-3, dfM*1e-6, spM, TPM, 0.0, PM, ratioM, HotM, ColdM, NaN)
	modules[4] = GasChromatographySystems.ModuleColumn("mod loop", LM[3], dM*1e-3, dfM*1e-6, spM, TPM)
	modules[5] = ModuleTM("TM2", LM[4], dM*1e-3, dfM*1e-6, spM, TPM, shiftM, PM, ratioM, HotM, ColdM, NaN)
	modules[6] = GasChromatographySystems.ModuleColumn("mod out", LM[5], dM*1e-3, dfM*1e-6, spM, TPM, NaN)
	modules[7] = GasChromatographySystems.ModuleColumn("GC column 2", L2, d2*1e-3, df2*1e-6, sp2, TP2, NaN)
	# system
	sys = update_system(GasChromatographySystems.System(g, pp, modules, opt))
	return sys
end

# ╔═╡ 7bb430cd-6bc5-4c48-85f8-3ffce3220d50
test_sys = GCxGC_TM(30.0, 0.25, 0.25, "Rxi5ms", GCxGC_TP, 2.0, 0.25, 0.25, "Rxi17SilMS", GCxGC_TP, [0.1, 0.02, 0.1, 0.02, 0.1], 0.25, 0.25, "Rxi17SilMS", 0.0, 6.0, 0.4, 80.0, -80.0, GCxGC_TP, 1.0, NaN, 0.0)

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
	T_itp(x, t) = therm_mod(t, module_.shift, module_.PM, module_.ratio, T_itp_(x, t) .+ module_.Thot .- 273.15, module_.Tcold) .+ 273.15
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

# ╔═╡ 7fc154cd-d1a4-4f14-8750-e92a58da35f0
TM_itp(x, t) = modulator_temperature(x, t, sys.modules[3])

# ╔═╡ daf5e89f-343c-4076-80fd-2dd52ba14f29
TM_itp(0.0,0.0)

# ╔═╡ 051e9ad6-51de-40eb-a3cf-2b929d594400
TM_itp(0.0,3.0)

# ╔═╡ cf03f9b5-62e6-46af-924e-3a16638c768d
TM_itp(0.0,3.1)

# ╔═╡ 607265de-c784-4ed0-8e89-ebb206785b0b
Plots.plot(0.0:0.01:100.0, TM_itp.(0.0, 0.0:0.01:100.0).-273.115)

# ╔═╡ 2e4612e7-ce0f-4ae7-97b0-7076d9bf5651
κ = flow_restrictions(sys)

# ╔═╡ c366df65-59a4-427a-beee-bf0bef5fb409
Plots.plot(0.0:0.01:100.0, κ[3].(0.0:0.01:100.0))

# ╔═╡ c1cb420c-7869-4da8-a240-1583545bde7c
GasChromatographySystems.plot_pressure_over_time(sys)

# ╔═╡ 93f19f34-5766-4863-bfba-fde045afcdb9
pf = GasChromatographySystems.pressure_functions(sys)

# ╔═╡ 5428a416-48e4-4ebd-a493-54fd30cd2556
Plots.plot(10.0:0.01:20.0, pf[3].(10.0:0.01:20.0))

# ╔═╡ 6378cb76-57fd-4e59-a183-fd2c745d0f69
Ff = GasChromatographySystems.flow_functions(sys)

# ╔═╡ bdbd8cc9-29fb-43ce-850b-8d18980808d3
md"""
## Simulation
"""

# ╔═╡ 4adbf7e9-232d-4889-92c6-75e8bde0d88d
db = filter([:CAS] => x -> !ismissing(x), DataFrame(CSV.File(joinpath(dirname(pwd()), "data", "database_Rxi17SilMS_Rxi5SilMS_.csv"), header=1, silencewarnings=true)))

# ╔═╡ fc36ffe3-b1df-42c8-ad14-e3b90f27a272
selected_solutes_ = GasChromatographySystems.common_solutes(db, sys).Name

# ╔═╡ 43f27a81-dca0-4188-baf4-3f8140d9e57c
selected_solutes = selected_solutes_[1:13]

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
			time_steps = common_timesteps(sys)
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
	plotly()
	Plots.plot(0.0:0.01:100.0, par[3].prog.T_itp.(0.0, 0.0:0.01:100.0))
	Plots.plot!(0.0:0.01:100.0, par[4].prog.T_itp.(0.0, 0.0:0.01:100.0))
	Plots.plot!(0.0:0.01:100.0, par[5].prog.T_itp.(0.0, 0.0:0.01:100.0))
end

# ╔═╡ f21bca1b-6d3a-4ca8-a6c6-3107f9f5dbb9
test_par = graph_to_parameters(test_sys, db, selected_solutes)

# ╔═╡ 175b7e36-71f3-4f6b-9467-b3c9515da114
begin
	plotly()
	Plots.plot(0.0:0.01:100.0, test_par[3].prog.T_itp.(0.0, 0.0:0.01:100.0))
	Plots.plot!(0.0:0.01:100.0, test_par[4].prog.T_itp.(0.0, 0.0:0.01:100.0))
	Plots.plot!(0.0:0.01:100.0, test_par[5].prog.T_itp.(0.0, 0.0:0.01:100.0))
end

# ╔═╡ a14d1ca2-46ef-4296-aa1d-06590eca97cf
paths = GasChromatographySystems.all_paths(sys.g, 1)[2]

# ╔═╡ 69598506-cf92-4626-b6e1-d8f57bf8388f
# simulation with the old function results in ???

# ╔═╡ 1a0695c5-0a67-4158-a3bd-c135be31e408
sim = GasChromatographySystems.simulate_along_paths(sys, paths, par)

# ╔═╡ 6556a858-65ef-49a6-bdaa-562a2b568a70
# simulation using the new function

# ╔═╡ fcdc5015-31d0-4808-8218-f33901437c0b
#sim_ = simulate_along_paths(sys, paths, par, abstol=1e-6, reltol=1e-3)

# ╔═╡ 1f62b350-332e-4acf-b907-75362157e4da
sim_[2][1][3]

# ╔═╡ f1bf21dd-07ee-4e4e-beb5-7d9225b3fc6b
sim_[2][1][5]

# ╔═╡ e3e34deb-9249-4811-b74a-44a4c2d06ac2
md"""
## Modulation
"""

# ╔═╡ ce317588-2f42-477e-a0a6-cc40047e1506
# reconstruction of function simulate_along_paths(sys, paths, par_sys; t₀=zeros(length(par_sys[1].sub)), τ₀=zeros(length(par_sys[1].sub)), refocus=falses(ne(sys.g)), τ₀_focus=zeros(length(par_sys[1].sub)), nτ=6):

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
sim[2][1][2].tR

# ╔═╡ 325a8e2e-ae2f-4dd9-a32b-00c7adf47461
begin
	# if next module == moduleTM
	par_prev = par[2]
	# initial peak width?
	τ0 = 0.0 # or PM?
#	init_t_start = (fld.(sim[2][1][2].tR.-6*sim[2][1][2].τR, PM).-1).*PM # is the -1 here correct?????
#	init_t_end = (fld.(sim[2][1][2].tR.+6*sim[2][1][2].τR, PM).-1).*PM 
	nτ = 20
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
	newopt = GasChromatographySimulator.Options(par[3].opt.alg, 1e-14, 1e-9, par[3].opt.Tcontrol, par[3].opt.odesys, par[3].opt.ng, par[3].opt.vis, par[3].opt.control, par[3].opt.k_th)
	#newopt = GasChromatographySimulator.Options(BS3(), 1e-6, 1e-3, par[3].opt.Tcontrol, par[3].opt.odesys, par[3].opt.ng, par[3].opt.vis, par[3].opt.control, par[3].opt.k_th)
	newpar = GasChromatographySimulator.Parameters(par[3].col, par[3].prog, sub_TM, newopt)
	
end

# ╔═╡ c626e4de-3c65-4230-bc23-2fa6ea9bbbd9
par[3]

# ╔═╡ 03aeeaac-0a83-4d34-944f-d3457904d69d
newpar

# ╔═╡ e4172084-2dd7-4daf-b2b9-cc84891c06f2
slicing(sim[2][1][2].tR, sim[2][1][2].τR, PM, par[3], par[2]; nτ=6)

# ╔═╡ 6faa498d-08de-49a6-94c5-9252da61b5a8
zeros(length(sim[2][1][2].τR))

# ╔═╡ d3c212b5-4a80-4d0f-ab9f-94893dec552e
init_t_start

# ╔═╡ 3fab9b05-a578-4198-b0c3-850760b89bc3
init_t_end

# ╔═╡ c072add2-9159-456c-beb2-54a081132bf6
newpar.sub

# ╔═╡ 532b0a0e-9595-47aa-add8-53dfdf54f568


# ╔═╡ 3057e8b2-e23c-4bb6-94ad-403a9fb94e61
# from GasChromatographySimulator
function solving_odesystem_r(col::GasChromatographySimulator.Column, prog::GasChromatographySimulator.Program, sub::GasChromatographySimulator.Substance, opt::GasChromatographySimulator.Options; dt=1e-10)#mod
    t₀ = [sub.t₀; sub.τ₀^2]
    zspan = (0.0,col.L)
	p = [col, prog, sub, opt]
    prob = ODEProblem(GasChromatographySimulator.odesystem_r!, t₀, zspan, p)

    solution = solve(prob, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol, dt=dt)#mod

    if solution.t[end]<col.L
        solution = solve(prob, alg=opt.alg, abstol=opt.abstol,reltol=opt.reltol, dt=col.L/1000000)
    end
    return solution
end

# ╔═╡ e8a0ca14-f566-4d7e-b49b-34934d54aa3c
# from GasChromatographySimulator
function solve_system_multithreads(par; dt=1e-10)#mod
	n = length(par.sub)
	sol = Array{Any}(undef, n)
	Threads.@threads for i=1:n
		sol[i] = solving_odesystem_r(par.col, par.prog, par.sub[i], par.opt, dt=dt)#mod
	end
	return sol
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
    df = sort!(DataFrame(Name = Name, CAS = CAS, tR = tR, τR = τR, TR=TR, σR = σR, uR = uR, kR = kR, ), [:tR])
    Threads.@threads for i=1:n-1
        Res[i] = (df.tR[i+1] - df.tR[i])/(2*(df.τR[i+1] + df.τR[i]))
        Δs[i] = (df.tR[i+1] - df.tR[i])/(df.τR[i+1] - df.τR[i]) * log(df.τR[i+1]/df.τR[i])
    end
    df[!, :Res] = Res
    df[!, :Δs] = Δs 
	df[!, :ann] = ann
    return df
end

# ╔═╡ f48d790f-6ef0-4e31-a80d-8acaae8e454f
# from GasChromatographySimulator
function simulate(par; dt=1e-10)#mod
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
#						peaklists_[j][!,:A] = ones(length(peaklists_[j].Name))
						# increase peakwidth of 1st dimension by 100%
#						peaklists_[j][!,:τR] = 2.0.*peaklists_[j][!,:τR]
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
soltest = simulate(newpar)

# ╔═╡ f087274e-5147-48f9-8200-4d6ed90397c9
par[3].col.L ./ soltest[1].uR

# ╔═╡ 90d32247-3e00-41ee-bd2e-b8e698c2613c
soltest[1].ann

# ╔═╡ 2b89c457-0723-47cc-bb69-7cd83c84e9ad
Plots.plot(soltest[2][8])

# ╔═╡ 8cd33b6b-68c5-4237-81fa-7a26f37e76ae
Plots.plot(738.0:0.1:750.0, newpar.prog.T_itp.(0.0, 738.0:0.1:750.0).-273.15)

# ╔═╡ d3d40b69-d659-4604-b18a-7722f3e7446c
init_t_start

# ╔═╡ f6093dce-30a7-4628-8362-0c5ed8c55ebc
md"""
### between Modulator spots
"""

# ╔═╡ beb15214-b1bb-4b39-a2aa-86f0727bc2c4
begin
	# loop between the modulator spots
	par_prev_ = newpar
	# initial peak width?
	# spacial width of the previous modulator spot * solute residency
	τ0_ = par_prev_.col.L ./ soltest[1].uR
	sub_loop = Array{GasChromatographySimulator.Substance}(undef, length(sub_TM))
	for i=1:length(sub_TM)
		# should check with CAS/name and ann for the same solute
		sub_loop[i] = GasChromatographySimulator.Substance(par_prev_.sub[i].name, par_prev_.sub[i].CAS, par_prev_.sub[i].Tchar, par_prev_.sub[i].θchar, par_prev_.sub[i].ΔCp, par_prev_.sub[i].φ₀, par_prev_.sub[i].ann, par_prev_.sub[i].Cag, soltest[1].tR[i], τ0_[i])
	end
	newpar_loop = GasChromatographySimulator.Parameters(par[4].col, par[4].prog, sub_loop, par[4].opt)
end

# ╔═╡ 74121702-9f4a-410b-961b-7865ab3af941
soltest_loop = simulate(newpar_loop)

# ╔═╡ 86274690-dfe8-4f71-8084-66fe5f03f932
md"""
### 2nd Modulator spot
"""

# ╔═╡ 12c5ce26-4a91-4b1b-b86e-0465c1c155a3
begin
	# if next module == moduleTM
	par_prev__ = newpar_loop
	# initial peak width?
	τ0__ = soltest_loop[1].τR
	init_t_start__ = (fld.(soltest_loop[1].tR.-6*soltest_loop[1].τR, PM)).*PM
	init_t_end__ = (fld.(soltest_loop[1].tR.+6*soltest_loop[1].τR, PM)).*PM 
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
	newopt__ = GasChromatographySimulator.Options(par[5].opt.alg, 1e-14, 1e-10, par[5].opt.Tcontrol, par[5].opt.odesys, par[5].opt.ng, par[5].opt.vis, par[5].opt.control, par[5].opt.k_th)
	newpar__ = GasChromatographySimulator.Parameters(par[5].col, par[5].prog, sub_TM__, newopt__)
	
end

# ╔═╡ f9385263-08dc-49d1-991d-637e053f4f99
begin
	# loop between the modulator spots
	par_prev___ = newpar_loop
	# initial peak width?
	# spacial width of the previous modulator spot * solute residency
	τ0___ = soltest_loop[1].τR
	sub_TM___ = Array{GasChromatographySimulator.Substance}(undef, length(sub_TM))
	for i=1:length(sub_TM)
		# should check with CAS/name and ann for the same solute
		sub_TM___[i] = GasChromatographySimulator.Substance(par_prev_.sub[i].name, par_prev_.sub[i].CAS, par_prev_.sub[i].Tchar, par_prev_.sub[i].θchar, par_prev_.sub[i].ΔCp, par_prev_.sub[i].φ₀, par_prev_.sub[i].ann, par_prev_.sub[i].Cag, (1+fld(soltest_loop[1].tR[i], PM))*PM, τ0___[i])
	end #                      +1 PM ????
	newopt___ = GasChromatographySimulator.Options(par[5].opt.alg, 1e-13, 1e-10, par[5].opt.Tcontrol, par[5].opt.odesys, par[5].opt.ng, par[5].opt.vis, par[5].opt.control, par[5].opt.k_th)
	newpar___ = GasChromatographySimulator.Parameters(par[5].col, par[5].prog, sub_TM___, newopt___)
end

# ╔═╡ 4c85cfd1-cb1b-449d-a101-e537abe25e86
soltest_2 = simulate(newpar___)

# ╔═╡ 6c53f8a8-f6f3-4af1-84f4-e03dc161aca6
soltest_loop[1].tR.- soltest[1].tR

# ╔═╡ 17c9a0d3-7904-49d9-8ef6-88bf74c1ba59
soltest[1]

# ╔═╡ dc4f9bec-695b-46b6-90ed-54ac2d3648bb
soltest_loop[1]

# ╔═╡ 06dceade-48b4-459d-b71f-1149891d1a6c
soltest_2[1]

# ╔═╡ 64e635f8-f5fd-4935-b44e-b02ebdd96b64
md"""
### connection to GC2
"""

# ╔═╡ f388ca37-0167-4459-858b-619f1d2ce2d7
begin
	# loop between the modulator spots
	par_prev_c = newpar___
	# initial peak width?
	# spacial width of the previous modulator spot * solute residency
	τ0_c = par_prev_c.col.L ./ soltest_2[1].uR
	sub_c = Array{GasChromatographySimulator.Substance}(undef, length(sub_TM))
	for i=1:length(sub_TM)
		# should check with CAS/name and ann for the same solute
		sub_c[i] = GasChromatographySimulator.Substance(par_prev_c.sub[i].name, par_prev_c.sub[i].CAS, par_prev_c.sub[i].Tchar, par_prev_c.sub[i].θchar, par_prev_c.sub[i].ΔCp, par_prev_c.sub[i].φ₀, par_prev_c.sub[i].ann, par_prev_c.sub[i].Cag, soltest_2[1].tR[i], τ0_c[i])
	end
	newpar_c = GasChromatographySimulator.Parameters(par[6].col, par[6].prog, sub_c, par[6].opt)
end

# ╔═╡ bb01e67d-fbca-461e-b592-508bc5480b94
soltest_c = simulate(newpar_c)

# ╔═╡ 14bf7364-ce16-467e-894b-43d5ad1d91dc
md"""
### GC2
"""

# ╔═╡ 98713f9d-123a-4eda-8a8c-72c6ada1a3a4
begin
	# loop between the modulator spots
	par_prev_gc2 = newpar_c
	# initial peak width?
	# spacial width of the previous modulator spot * solute residency
	τ0_gc2 = soltest_c[1].τR
	sub_gc2 = Array{GasChromatographySimulator.Substance}(undef, length(sub_TM))
	for i=1:length(sub_TM)
		# should check with CAS/name and ann for the same solute
		sub_gc2[i] = GasChromatographySimulator.Substance(par_prev_gc2.sub[i].name, par_prev_gc2.sub[i].CAS, par_prev_gc2.sub[i].Tchar, par_prev_gc2.sub[i].θchar, par_prev_gc2.sub[i].ΔCp, par_prev_gc2.sub[i].φ₀, par_prev_gc2.sub[i].ann, par_prev_gc2.sub[i].Cag, soltest_c[1].tR[i], τ0_gc2[i])
	end
	newpar_gc2 = GasChromatographySimulator.Parameters(par[7].col, par[7].prog, sub_gc2, par[7].opt)
end

# ╔═╡ caef921f-f06c-471b-a56a-8415ac2d4fa6
soltest_gc2 = simulate(newpar_gc2)

# ╔═╡ 915d8972-5423-4dde-8039-9a912bd0ae61
slice_tR1 = fld.(soltest_gc2[1].tR, PM).*PM

# ╔═╡ 556c157c-321d-4eb8-b2ef-ee8c4fe0f038
slice_tR2 = soltest_gc2[1].tR .- fld.(soltest_gc2[1].tR, PM).*PM

# ╔═╡ a76b6b3e-a41c-413e-9b4d-920823e9edda
md"""
## Chromatograms
"""

# ╔═╡ ac2e34ed-a519-48fc-a87e-fb8dde6f3ff5
md"""
### Result from step-by-step
"""

# ╔═╡ ae9bb1ec-fc67-46de-91fe-6e3f47844ce0
begin
	t_ = 0.0:0.01:sum(GCxGC_TP.timesteps)
	c_ = GasChromatographySimulator.chromatogram(collect(t_), soltest_gc2[1].tR, soltest_gc2[1].τR)
	p_chrom_1 = Plots.plot(t_, c_)
end

# ╔═╡ 5af9ebbe-b126-4ad4-9a44-c934d0d18cb3
md"""
### Result from new function
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

# ╔═╡ 947a0624-2484-470e-a6c1-afbcc4e92d9c
md"""
### **Problem**
At the thermal modulator spot the solution of the migration can become wrong. The periodic heating up with the hot jet seems to be "over looked" by the solving algorithm and only to a later hot jet moment the substance migrates through the modulator spot.
"""

# ╔═╡ 569114c5-9cae-44b0-9983-c5e36b7723ad
sim_[3][1][3][3]

# ╔═╡ 7d383c07-6a67-42cb-8659-a5d21e24854a
begin
	Δt_TM1 = Array{Float64}(undef, length(sim_[3][1][3]))
	for i=1:length(sim_[3][1][3])
		Δt_TM1[i] = sim_[3][1][3][i].u[end][1] - sim_[3][1][3][i].u[1][1]
	end
	Δt_TM1
end

# ╔═╡ c513358a-a18e-4678-9adb-59352592f77b
begin
	Δt_TM2 = Array{Float64}(undef, length(sim_[3][1][5]))
	for i=1:length(sim_[3][1][5])
		Δt_TM2[i] = sim_[3][1][5][i].u[end][1] - sim_[3][1][5][i].u[1][1]
	end
	Δt_TM2
end

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
	c_slices, t_D1, t_D2 = chrom_slicing(t_, c_, 6.0)
	p_chrom_2 = Plots.plot(t_D2[1], c_slices[1])
	for i=2:length(c_slices)
		Plots.plot!(p_chrom_2, t_D2[i], c_slices[i].+i*10)
	end
	p_chrom_2
end

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
sys_simp = GCxGC_TM_simp(30.0, 0.25, 0.25, "Rxi5ms", GCxGC_TP, 2.34, 0.25, 0.25, "Rxi17SilMS", GCxGC_TP, 1.0, NaN, 0.0)

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
	offset_D2 = 0.0
	p_sliced = Plots.scatter(slice_tR1.-2*PM, slice_tR2.+offset_D2, ylims=(0.0, 6.0))
	Plots.scatter!(p_sliced, sim[2][1][1].tR, 2.0.*ones(length(sim[2][1][1].tR)))
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
mod.(pl.tR2.-pl.tR1, 6.0)

# ╔═╡ 2ca1ac9f-a4f7-48d3-ad4a-43dc01d97b18
(pl.tR2.-pl.tR1) .- mod.(pl.tR2.-pl.tR1, 6.0)

# ╔═╡ ffaca8b2-5d88-4aeb-a8a2-b7fdbd661cad
sim

# ╔═╡ 709868c9-4f3d-4445-8831-bcdda9cacaf4
sim[3][1][1][1]

# ╔═╡ 82d56178-1c4f-43b1-992a-bbdc90ddbf78
sim[3][1][2][1]

# ╔═╡ 1edc472e-23be-460a-94b4-2c6d776b3890
sim[3][1][3][1]

# ╔═╡ 39d38ca1-0cfa-423f-8493-e2fdd8912350
sim[3][1][3][1].u[end][1]-sim[3][1][3][1].u[1][1]

# ╔═╡ 50c61ee0-37d2-474a-b971-9469d11027a4
sim[3][1][4][1]

# ╔═╡ f3e9d229-5cde-4062-9af6-a478f94013cb
sim[3][1][5][1]

# ╔═╡ fdb6a9a3-720d-4712-90c3-5037fb6ba5b1
sim[3][1][5][1].u[end][1]-sim[3][1][5][1].u[1][1]

# ╔═╡ 657f0ab1-8d7b-4398-9e2f-ef636681b7c1
sim[3][1][6][1]

# ╔═╡ 65f5b190-9c52-4d7c-82b6-0758a4c9db14
sim[3][1][7][1]

# ╔═╡ c4366067-3f97-41dd-8288-38fb03b9c16d
begin
	Δt_TM = Array{Float64}(undef, length(selected_solutes), 2)
	for i=1:length(selected_solutes)
		Δt_TM[i,1] = sim[3][1][3][i].u[end][1] - sim[3][1][3][i].u[1][1]
		Δt_TM[i,2] = sim[3][1][5][i].u[end][1] - sim[3][1][5][i].u[1][1] 
	end
	pl[!,:Δt_TM1] = Δt_TM[:,1]
	pl[!,:Δt_TM2] = Δt_TM[:,2]
	pl
end

# ╔═╡ e2eaa21a-fbac-4d00-8663-c69ee4e0a8a9
md"""
## Traces
"""

# ╔═╡ e2459cb1-3aa2-4192-9de7-c57b997041bc
function traces(sim, i_select)
	#i_select = 10
	z = Float64[]
	t = Float64[]
	τ² = Float64[]
	T = Float64[]
	k = Float64[]
	for i=1:ne(sys.g)
		if i == 1
			append!(z, sim[3][1][i][i_select].t)
		else
			append!(z, sim[3][1][i][i_select].t .+ z[end])
		end
		#z0 = z[end]
		tt = Array{Float64}(undef, length(sim[3][1][i][i_select].t))
		ττ = Array{Float64}(undef, length(sim[3][1][i][i_select].t))
		TT = Array{Float64}(undef, length(sim[3][1][i][i_select].t))
		kk = Array{Float64}(undef, length(sim[3][1][i][i_select].t))
		for j=1:length(sim[3][1][i][i_select].t)
			tt[j] = sim[3][1][i][i_select].u[j][1]
			ττ[j] = sim[3][1][i][i_select].u[j][2]
			TT[j] = par[i].prog.T_itp(sim[3][1][i][i_select].t[j], sim[3][1][i][i_select].u[j][1])
			kk[j] = GasChromatographySimulator.retention_factor(sim[3][1][i][i_select].t[j], sim[3][1][i][i_select].u[j][1], par[i].col, par[i].prog, par[i].sub[i_select], par[i].opt)
		end
		append!(t, tt)
		append!(τ², ττ)
		append!(T, TT)
		append!(k, kk)
	end
	trace = DataFrame(z=z, t=t, τ²=τ², T=T, k=k)
end

# ╔═╡ 4dbe3105-5ff0-4d16-a163-81f4ae6c9889
trace = traces(sim, 10)

# ╔═╡ 5e5a5472-25e0-4914-9086-a288f81ba8b6
Plots.scatter(trace.z, trace.t)

# ╔═╡ 9067aa11-0b8f-4453-bb48-fe75eaf9800b
Plots.scatter(trace.z, trace.τ²)

# ╔═╡ 9e7a5d52-68c7-4edb-915e-b5091acebeba
Plots.scatter(trace.z, trace.T.-273.15)

# ╔═╡ ed599e7c-f221-43ba-b56d-8fcec83c4bd9
Plots.scatter(trace.z, trace.k)

# ╔═╡ 17c4fbd1-4a9d-4f7a-922a-a72eb97db33a
md"""
# End
"""

# ╔═╡ Cell order:
# ╠═23c71f14-efdf-11ed-33f2-95304516114d
# ╠═d7f045e5-47c7-44b6-8bce-ce9dc3154fff
# ╠═d419921b-1472-4b66-bae4-469537259814
# ╠═dff75c18-eee5-4f70-9509-0e8228932819
# ╠═b4a1fe1b-f44a-4844-b041-1c9eabcd7ce2
# ╠═cbe42d0d-e765-4a3f-b7dc-ba5687ac77e2
# ╠═2777203c-c074-469c-9194-efa526676002
# ╠═cbc53a06-1efc-4543-9151-52f9e95e6ae3
# ╠═b283ebdd-38f5-4b8c-a1b7-34476d21dff2
# ╠═5cc7905e-477f-46ee-816a-f1327f914fd8
# ╠═58441de1-11cd-4d44-8664-fbc0b9039ed5
# ╠═927e9354-4130-47cd-bd94-8eb37da15b47
# ╠═b25f6d02-173a-4d17-adcf-c54ad854a871
# ╠═6e5eba86-f046-4cf9-b6a5-c2b0622bf1fc
# ╠═5aa38a4e-d9fa-4a16-b6e1-c9ee2eb7f2d9
# ╠═b1a8b34e-49f4-4678-91dc-6d60bcb743e3
# ╠═14865482-42d6-4d13-b738-54a1a67a65bc
# ╠═735300f0-b41a-41ce-8a96-707078dcd808
# ╠═0080eaa8-89ef-4c50-b9e8-a6630d088532
# ╠═14f47b21-6d1f-4e46-9bcb-7240f61bc492
# ╠═9031c353-5dbf-4296-bacc-d3f6338e9d8e
# ╠═02299b66-0054-4b77-ad36-d7f0d4f6fd21
# ╠═8e38d02b-1161-46f7-92ad-5d6ef1e8ede4
# ╠═1231d383-532b-4f01-8975-f65f574def0c
# ╠═615d984f-2452-4234-b9be-7eb97f3de94a
# ╠═5c056b8c-304c-4348-a49b-24ef86561418
# ╠═f635daa6-5bc5-4831-8320-3abf4a9dcc14
# ╠═c0a92e37-2592-4cf5-a572-228def68b8a1
# ╠═264e80c1-d8f9-4513-b599-ac25a021f321
# ╠═3c545fbb-2732-4ac5-b628-7cfda7bc6770
# ╠═2842aa10-55a8-4e53-b582-d0693a91d0c2
# ╠═3722291a-8f23-4bdd-a86d-8ae2056ee809
# ╠═7f19f5df-6626-48c7-b4d9-7eb897c12bde
# ╠═b38300b5-9b16-45e9-9d5f-19e0e825349f
# ╠═324b1803-dab5-4b68-a9a9-385b562c1c3f
# ╠═326d2249-dfcc-40b6-b8c0-d9ab5ff88627
# ╠═673fc2e1-43c3-429c-a7df-4930c407bd8f
# ╠═33a88b81-a65b-4ba3-b148-7339a8ff97a9
# ╟─d0adda7c-8980-4d74-8b79-c3c5edd8131f
# ╠═ab33d65f-f0ba-4e5d-b660-2ecd106d8d21
# ╠═10427844-f557-4178-807b-f7cd194548c3
# ╠═b63ea500-c2bb-49c3-9915-032fd3e33e14
# ╠═f5a4bfba-1379-4b4f-a694-79426695d254
# ╠═60d86bdd-5527-49da-8c35-ecd508171e6c
# ╠═e44cf22a-4fbe-43b7-a3ca-2d17a41a3f43
# ╠═8ca0c042-7e2e-4248-a183-58235e909f59
# ╠═3fde3b56-e5a4-4a28-95a6-0d8b36fe37d1
# ╠═7bb430cd-6bc5-4c48-85f8-3ffce3220d50
# ╠═27d0f80e-52b9-4942-a519-fcd787fdfaae
# ╠═80150c58-f9b9-43d9-bd7b-e38388f29fd8
# ╠═c7e771b6-7ec7-41bf-940f-2e2617af3c02
# ╠═522a7582-0ff1-42d6-960a-adb952a225d2
# ╠═ca604ac6-fdb1-41d7-83d2-a0c9b2440dad
# ╠═b9e3acb8-a7ec-40ef-ba58-932f79b08fb0
# ╠═95f57b5b-652c-40c2-b66d-137bb7b2d878
# ╠═9605b6ad-7cdc-47ee-9edf-c3ce531ae9a2
# ╠═06487ba0-8752-40a8-b5d8-07fb08514ffe
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
# ╠═7277a646-fd70-4032-9691-569eddeafec6
# ╠═793560f7-d364-4c68-81ec-994441a41059
# ╠═d40361fe-9e0c-4df8-a09c-0c5bf143dbf6
# ╠═5d63fd50-208e-4db1-9609-4f08c493d90c
# ╠═6993afed-4519-4c8a-9cfc-87c72723c444
# ╠═f21bca1b-6d3a-4ca8-a6c6-3107f9f5dbb9
# ╠═175b7e36-71f3-4f6b-9467-b3c9515da114
# ╠═a14d1ca2-46ef-4296-aa1d-06590eca97cf
# ╠═69598506-cf92-4626-b6e1-d8f57bf8388f
# ╠═1a0695c5-0a67-4158-a3bd-c135be31e408
# ╠═6556a858-65ef-49a6-bdaa-562a2b568a70
# ╠═fcdc5015-31d0-4808-8218-f33901437c0b
# ╠═1f62b350-332e-4acf-b907-75362157e4da
# ╠═f1bf21dd-07ee-4e4e-beb5-7d9225b3fc6b
# ╠═e3e34deb-9249-4811-b74a-44a4c2d06ac2
# ╠═ce317588-2f42-477e-a0a6-cc40047e1506
# ╠═6e8f2101-8245-483c-9c05-e4b16c260607
# ╠═bac5d026-0b84-412f-98bb-8ddedb2b92d9
# ╟─ce958d12-3dad-4e3c-ab57-f84136610870
# ╠═70597ee8-760d-4c21-af6c-981922208d10
# ╠═b83330dd-4df3-433a-ab29-273e39f8b32c
# ╠═325a8e2e-ae2f-4dd9-a32b-00c7adf47461
# ╠═c626e4de-3c65-4230-bc23-2fa6ea9bbbd9
# ╠═03aeeaac-0a83-4d34-944f-d3457904d69d
# ╠═e4172084-2dd7-4daf-b2b9-cc84891c06f2
# ╠═6faa498d-08de-49a6-94c5-9252da61b5a8
# ╠═d3c212b5-4a80-4d0f-ab9f-94893dec552e
# ╠═3fab9b05-a578-4198-b0c3-850760b89bc3
# ╠═c072add2-9159-456c-beb2-54a081132bf6
# ╠═f436afe9-0440-4149-be61-b426e243b34d
# ╠═f087274e-5147-48f9-8200-4d6ed90397c9
# ╠═532b0a0e-9595-47aa-add8-53dfdf54f568
# ╠═90d32247-3e00-41ee-bd2e-b8e698c2613c
# ╠═f48d790f-6ef0-4e31-a80d-8acaae8e454f
# ╠═e8a0ca14-f566-4d7e-b49b-34934d54aa3c
# ╠═3057e8b2-e23c-4bb6-94ad-403a9fb94e61
# ╠═7708f828-0ea3-4568-8614-ff0a724b9eea
# ╠═2b89c457-0723-47cc-bb69-7cd83c84e9ad
# ╠═8cd33b6b-68c5-4237-81fa-7a26f37e76ae
# ╠═d3d40b69-d659-4604-b18a-7722f3e7446c
# ╠═f6093dce-30a7-4628-8362-0c5ed8c55ebc
# ╠═beb15214-b1bb-4b39-a2aa-86f0727bc2c4
# ╠═74121702-9f4a-410b-961b-7865ab3af941
# ╠═86274690-dfe8-4f71-8084-66fe5f03f932
# ╠═12c5ce26-4a91-4b1b-b86e-0465c1c155a3
# ╠═f9385263-08dc-49d1-991d-637e053f4f99
# ╠═4c85cfd1-cb1b-449d-a101-e537abe25e86
# ╠═6c53f8a8-f6f3-4af1-84f4-e03dc161aca6
# ╠═17c9a0d3-7904-49d9-8ef6-88bf74c1ba59
# ╠═dc4f9bec-695b-46b6-90ed-54ac2d3648bb
# ╠═06dceade-48b4-459d-b71f-1149891d1a6c
# ╠═64e635f8-f5fd-4935-b44e-b02ebdd96b64
# ╠═f388ca37-0167-4459-858b-619f1d2ce2d7
# ╠═bb01e67d-fbca-461e-b592-508bc5480b94
# ╠═14bf7364-ce16-467e-894b-43d5ad1d91dc
# ╠═98713f9d-123a-4eda-8a8c-72c6ada1a3a4
# ╠═caef921f-f06c-471b-a56a-8415ac2d4fa6
# ╠═915d8972-5423-4dde-8039-9a912bd0ae61
# ╠═556c157c-321d-4eb8-b2ef-ee8c4fe0f038
# ╠═a76b6b3e-a41c-413e-9b4d-920823e9edda
# ╠═ac2e34ed-a519-48fc-a87e-fb8dde6f3ff5
# ╠═9c3223ab-1dcc-40b6-8996-c8dbc22a2c6a
# ╠═ae9bb1ec-fc67-46de-91fe-6e3f47844ce0
# ╠═5a510b13-491d-48d3-84b1-1335aaf6cd06
# ╠═5af9ebbe-b126-4ad4-9a44-c934d0d18cb3
# ╠═18cc22ce-0f7c-4683-ab2a-cae5a25ad9c4
# ╠═f2522fe2-e5e4-44bc-8e5a-528b86dfc290
# ╠═20d112f4-f5cc-4e6b-89c5-a0b3ad4ff4c5
# ╠═fd4026bf-7085-4a63-a074-c1bbf49bde2a
# ╠═947a0624-2484-470e-a6c1-afbcc4e92d9c
# ╠═569114c5-9cae-44b0-9983-c5e36b7723ad
# ╠═7d383c07-6a67-42cb-8659-a5d21e24854a
# ╠═c513358a-a18e-4678-9adb-59352592f77b
# ╠═622a22d0-ecac-402f-96a5-79c6dc1683d3
# ╠═5558e8e6-bb35-497b-9d61-ba0f3518e75e
# ╠═21d1659a-e50f-4908-8b7b-bcee33f7acae
# ╠═0ccf11db-c934-4592-95ae-d28bd856f97f
# ╠═9bb3302a-16c7-442c-b6a1-9ee1b9d20587
# ╠═04ce95ea-affb-45e1-b09c-dfd412a7291b
# ╠═b360f18a-fb8a-4e66-b13f-67e14d008238
# ╠═3e0f5a36-4db1-4ded-863b-2f0bfe01108a
# ╠═690988a3-e9f0-4c01-adb1-399ffbe2f093
# ╠═1f1fe7e1-d27d-485c-ad65-022b0862cc16
# ╠═b5598366-a6d2-49b3-b71f-1aefe5b9ef26
# ╠═6ad69424-5cd5-4ac2-ae15-c3220b6ce449
# ╠═5c0103b4-5d18-4100-8fa6-a7e8fdbc7ef1
# ╠═f74033df-d0d0-4641-9109-f0ff7b07abd9
# ╠═6e970b15-ff40-427c-9dd2-517e7ede00f8
# ╠═e5d6f0da-f249-4c09-b259-427fe240af34
# ╠═2ca1ac9f-a4f7-48d3-ad4a-43dc01d97b18
# ╠═ffaca8b2-5d88-4aeb-a8a2-b7fdbd661cad
# ╠═709868c9-4f3d-4445-8831-bcdda9cacaf4
# ╠═82d56178-1c4f-43b1-992a-bbdc90ddbf78
# ╠═1edc472e-23be-460a-94b4-2c6d776b3890
# ╠═39d38ca1-0cfa-423f-8493-e2fdd8912350
# ╠═50c61ee0-37d2-474a-b971-9469d11027a4
# ╠═f3e9d229-5cde-4062-9af6-a478f94013cb
# ╠═fdb6a9a3-720d-4712-90c3-5037fb6ba5b1
# ╠═657f0ab1-8d7b-4398-9e2f-ef636681b7c1
# ╠═65f5b190-9c52-4d7c-82b6-0758a4c9db14
# ╠═c4366067-3f97-41dd-8288-38fb03b9c16d
# ╠═e2eaa21a-fbac-4d00-8663-c69ee4e0a8a9
# ╠═e2459cb1-3aa2-4192-9de7-c57b997041bc
# ╠═4dbe3105-5ff0-4d16-a163-81f4ae6c9889
# ╠═5e5a5472-25e0-4914-9086-a288f81ba8b6
# ╠═9067aa11-0b8f-4453-bb48-fe75eaf9800b
# ╠═9e7a5d52-68c7-4edb-915e-b5091acebeba
# ╠═ed599e7c-f221-43ba-b56d-8fcec83c4bd9
# ╠═17c4fbd1-4a9d-4f7a-922a-a72eb97db33a
