### A Pluto.jl notebook ###
# v0.19.17

using Markdown
using InteractiveUtils

# ╔═╡ 80af05a0-8c11-11ed-1008-cb6a5443493d
begin
	using PlutoUI, Graphs, GraphMakie, Symbolics, NetworkLayout, CairoMakie, GasChromatographySimulator, Interpolations, CSV, DataFrames, Plots
	TableOfContents()
end

# ╔═╡ 8ace0fa0-52ef-4db8-92d0-7f99587bb065
md"""
# Development of Packackage
"""

# ╔═╡ a6d8c5bc-32ea-4f7b-89f5-56372243d27d
md"""
## [x] Constants
"""

# ╔═╡ 524d119e-9de7-43f0-b6dc-d719c71b7131
Tst = 273.15            # K

# ╔═╡ c8d84ca7-b160-4a47-aa6c-720ccd46090b
R = 8.31446261815324    # J mol⁻¹ K⁻¹

# ╔═╡ f0541435-7af0-408a-9d1a-378e57336989
Tn = 25.0 + Tst         # K

# ╔═╡ c2d60366-f54c-4b10-9803-d80e19a0855f
pn = 101300             # Pa

# ╔═╡ a693515b-e24d-4609-b97d-2dbcb6248414
md"""
## [x] Structures
"""

# ╔═╡ 46d410e9-b45a-40a0-932d-a025bba99a89
begin
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
end

# ╔═╡ d3107cc8-aa83-424c-be15-40de49084623
md"""
### [x] Module structures
"""

# ╔═╡ 932ffdbb-0dd6-48ca-98b2-af70f229d4f5
abstract type AbstractModule end

# ╔═╡ 73d0e4d5-76c6-4860-92f5-5888ae605861
struct ModuleTransferline<:AbstractModule
    # Module
	# transferline, constant temperature
	name::String
	length::Float64
	diameter::Float64
	film_thickness::Float64
	stationary_phase::String
	temperature::Float64
end

# ╔═╡ 471193bb-c2e1-4d27-9f0f-a9d65820ef7e
md"""
### [x] Temperature program structure
"""

# ╔═╡ a177a0a9-7e67-4644-896e-8863d4bd15c6
begin
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
end

# ╔═╡ 8a9c2704-5854-4aa4-9782-943e82cd38f6
begin
	struct ModuleColumn<:AbstractModule
	    # Module
		# GC column, gradients are possible
		name::String
		length::Float64
		diameter#::Fd # Function
	    a_diameter::Array{Float64,1} # Parameters of diameter function, just for information
		film_thickness#::Fdf # Function
	    a_film_thickness::Array{Float64,1} # Parameters of film_thickness function, just for information
		stationary_phase::String
		temperature_program::TemperatureProgram
	end

	function ModuleColumn(name, L, d, df, sp, tp)
	    # function to construct the Column structure
	    # for the case of constant diameter and constant film thickness
	    #d_func(x) = gradient(x, [d])
	    #df_func(x) = gradient(x, [df])
	    col = ModuleColumn(name, L, d, [d], df, [df], sp, tp)
	    return col
	end
end

# ╔═╡ 291ed3ff-e500-48e1-b1ec-2d28b5cbfa54
default_TP = TemperatureProgram([0.0, 1800.0], [40.0, 300.0])

# ╔═╡ 40663239-9c5b-466d-9f43-76ab913d5750
md"""
### [x] Pressure point structure
"""

# ╔═╡ f11b32e2-c447-4f54-8f42-c1e81c9fb489
struct PressurePoint
	# Pressure program, same structure for inlet and outlet
	name::String
	timesteps::Array{Float64,1}
	pressure_steps::Array{Float64,1}
	PressurePoint(n,ts,ps) = length(ts)!=length(ps) ? error("Mismatch between length(timesteps) = $(length(ts)), length(pressure_steps) = $(length(ps))") : new(n,ts,ps)
end

# ╔═╡ e7181f3b-e148-4327-a420-0091826b327c
md"""
### [x] System structure
"""

# ╔═╡ c5c010e6-80c6-445b-b977-0bf8559af0c4
struct System
	g::SimpleDiGraph{Int}
	pressurepoints::Array{PressurePoint}
	modules::Array{AbstractModule}
	options::Options
	System(g_,pressurepoints_,modules_,options_) = nv(g_)!=length(pressurepoints_) || ne(g_)!=length(modules_) ? error("Mismatch between number of nodes ($(nv(g_))) and number of pressure points ($(length(pressurepoints_))) and/or mismatch between number of edges ($(ne(g_))) and number of modules ($(length(modules_))).") : new(g_,pressurepoints_,modules_,options_)
end

# ╔═╡ 6f74b758-46a7-45cb-a177-8d766a8c42e8
md"""
## [x] Some Graphs
"""

# ╔═╡ 411fc6ec-cb7e-4ec9-b2f5-83b3b6347b52
function plot_graph(g, edge_labels, node_labels; lay = Stress(), color=:lightblue, node_size=40)
	fig, ax, p = GraphMakie.graphplot(g, 
						layout=lay,
						nlabels=node_labels, 
						nlabels_align=(:center,:center),
						node_size = [node_size for i in 1:nv(g)],
						node_color = [color for i in 1:nv(g)],
						elabels = edge_labels
					)
	hidedecorations!(ax)
	hidespines!(ax)
	ax.aspect = DataAspect()
	return fig
end

# ╔═╡ 4d8c7f3a-582c-4fe9-a692-5966492cc31f
function plot_graph(sys; lay = Stress(), color=:lightblue, node_size=40)
	fig, ax, p = GraphMakie.graphplot(sys.g, 
						layout=lay,
						nlabels = [sys.pressurepoints[i].name for i in 1:nv(sys.g)],
						nlabels_align=(:center,:center),
						node_size = [node_size for i in 1:nv(sys.g)],
						node_color = [color for i in 1:nv(sys.g)],
						elabels = [sys.modules[i].name for i in 1:ne(sys.g)]
					)
	hidedecorations!(ax)
	hidespines!(ax)
	ax.aspect = DataAspect()
	return fig
end

# ╔═╡ 8915304d-d902-4988-98c7-943f4f29478b
md"""
### [x] One Column
"""

# ╔═╡ 7cea41a4-a94d-4b44-8f32-04b324206d36
begin
	g0 = SimpleDiGraph(2)
	add_edge!(g0, 1, 2) # Inj -> GC column -> Det
	# pressure points:
	pp0 = Array{PressurePoint}(undef, nv(g0))
	pp0[1] = PressurePoint("p₁", [0.0, 1800.0], [200000.0, 200000.0]) # inlet 
	pp0[2] = PressurePoint("p₂", [0.0, 1800.0], [101300.0, 101300.0]) # outlet
	# modules
	modules0 = Array{AbstractModule}(undef, ne(g0))
	modules0[1] = ModuleColumn("GC column", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	# system
	sys0 = System(g0, pp0, modules0, Options())
end

# ╔═╡ 5a1c12c4-d1d3-4a6e-8832-2fb8cdf598b4
nv(g0)

# ╔═╡ 3d52db0f-92e6-425d-8618-7b9422d89b43
ne(g0)

# ╔═╡ e7f1e626-b002-4848-aac4-891e5efc2d63
md"""
### [x] Series of columns
"""

# ╔═╡ a1e4f8b0-f704-484e-87b0-b8c9a17a608f
UndefArray = Array{Real}(undef, 2)

# ╔═╡ d5991b35-16c2-4969-8b41-9255451ea843
# ╠═╡ disabled = true
#=╠═╡
a = missing
  ╠═╡ =#

# ╔═╡ 775eaf12-eb3d-4dc1-891a-a1a2be066be4
begin
	g_series = SimpleDiGraph(4)
	add_edge!(g_series, 1, 2) # Inj -> TL column -> GC column inlet
	add_edge!(g_series, 2, 3) # GC column inlet -> GC column -> 2nd TL column inlet
	add_edge!(g_series, 3, 4) # 2nd TL column inlet -> TL column -> Det
	# pressure points
	pp_series = Array{PressurePoint}(undef, nv(g_series))
	pp_series[1] = PressurePoint("p₁", [0.0, 600.0, 1200.0], [200000.0, 200000.0, 250000.0]) # inlet 
	pp_series[2] = PressurePoint("p₂", [0.0, 1800.0], [NaN, NaN]) 
	pp_series[3] = PressurePoint("p₃", [0.0, 1800.0], [NaN, NaN]) 
	pp_series[4] = PressurePoint("p₄", [0.0, 1800.0], [101300.0, 101300.0]) # outlet
	# modules
	modules_series = Array{AbstractModule}(undef, ne(g_series))
	modules_series[1] = ModuleTransferline("TL column 1", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", 300.0)
	modules_series[2] = ModuleColumn("GC column", 0.5, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	modules_series[3] = ModuleTransferline("TL column 2", 1.0, 0.25e-3, 0.25e-6, "SLB5ms", 300.0)
	# system
	sys_series = System(g_series, pp_series, modules_series, Options())
end

# ╔═╡ 70a6aff1-5547-46be-8884-17984c41fc64
md"""
### [x] Split
"""

# ╔═╡ 4c79d09f-5cd0-4788-ae06-a2be74bbc580
begin
	g_split = SimpleDiGraph(4)
	add_edge!(g_split, 1, 2) # Inj -> GC column -> Split point
	add_edge!(g_split, 2, 3) # Split point -> TL column -> Det 1
	add_edge!(g_split, 2, 4) # Split point -> TL column -> Det 2
	# pressure points
	pp_split = Array{PressurePoint}(undef, nv(g_split))
	pp_split[1] = PressurePoint("p₁", [0.0, 1800.0], [200000.0, 250000.0]) # inlet 
	pp_split[2] = PressurePoint("p₂", [0.0, 1800.0], [NaN, NaN]) # 
	pp_split[3] = PressurePoint("p₃", [0.0, 1800.0], [0.0, 0.0]) # outlet 1 
	pp_split[4] = PressurePoint("p₄", [0.0, 1800.0], [101300.0, 101300.0]) # outlet 2
	# modules
	modules_split = Array{AbstractModule}(undef, ne(g_split))
	modules_split[1] = ModuleColumn("GC column", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	modules_split[2] = ModuleTransferline("TL column 1", 0.5, 0.15e-3, 0.0e-6, "", 300.0)
	modules_split[3] = ModuleTransferline("TL column 2", 3.0, 0.25e-3, 0.25e-6, "SLB5ms", 300.0)
	# system
	sys_split = System(g_split, pp_split, modules_split, Options())
end

# ╔═╡ 3c476391-74bd-45e5-847a-9c239c3e2f0b
md"""
### [x] Complex
"""

# ╔═╡ fb417960-388d-44f8-8a53-39b549ae3ec2
begin
	# setup of the Graph
	g = SimpleDiGraph(9)
	add_edge!(g, 1, 2) # Inj->GC1->pressure regulator 
	add_edge!(g, 2, 3) # pressure regulator -> GC2 -> split
	add_edge!(g, 3, 4) # Split -> TL1 -> Det1
	add_edge!(g, 3, 5) # Split -> TL2 -> Modulator inlet
	add_edge!(g, 5, 6) # Modulator inlet -> Modulator -> GC3 inlet
	add_edge!(g, 6, 7) # GC3 inlet -> GC3 -> Split
	add_edge!(g, 7, 8) # Split -> TL3 -> Det2
	add_edge!(g, 7, 9) # Split -> TL4 -> Det3
	#edge_labels = ["GC₁", "GC₂", "TL₁", "TL₂", "Modulator", "GC₃", "TL₃", "TL₄"]
	#node_labels = ["p₁", "p₂", "p₃", "p₄", "p₅", "p₆", "p₇", "p₈", "p₉"]
	# pressure points
	pp = Array{PressurePoint}(undef, nv(g))
	pp[1] = PressurePoint("p₁", [0.0, 1800.0], [200000.0, 200000.0]) # inlet 
	pp[2] = PressurePoint("p₂", [0.0, 1800.0], [NaN, NaN]) # 
	pp[3] = PressurePoint("p₃", [0.0, 1800.0], [NaN, NaN]) #  
	pp[4] = PressurePoint("p₄", [0.0, 1800.0], [101300.0, 101300.0]) # outlet 
	pp[5] = PressurePoint("p₅", [0.0, 1800.0], [NaN, NaN]) #
	pp[6] = PressurePoint("p₆", [0.0, 1800.0], [NaN, NaN]) # 
	pp[7] = PressurePoint("p₇", [0.0, 1800.0], [NaN, NaN]) #
	pp[8] = PressurePoint("p₈", [0.0, 1800.0], [0.0, 0.0]) # outlet 2
	pp[9] = PressurePoint("p₉", [0.0, 1800.0], [101300.0, 101300.0]) # outlet 3
	# modules
	modules = Array{AbstractModule}(undef, ne(g))
	modules[1] = ModuleColumn("GC column 1", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	modules[2] = ModuleColumn("GC column 2", 30.0, 0.25e-3, 0.25e-6, "SPB50", default_TP)
	modules[3] = ModuleTransferline("TL column 1", 3.0, 0.25e-3, 0.25e-6, "SPB50", 300.0)
	modules[4] = ModuleTransferline("TL column 2", 1.0, 0.25e-3, 0.25e-6, "SPB50", 300.0)
	modules[5] = ModuleTransferline("Modulator", 0.3, 0.1e-3, 0.1e-6, "Wax", 300.0) # replace later with ModuleModulator
	modules[6] = ModuleColumn("GC column 3", 4.0, 0.1e-3, 0.1e-6, "Wax", default_TP)
	modules[7] = ModuleTransferline("TL column 3", 0.5, 0.1e-3, 0.1e-6, "Wax", 300.0)
	modules[8] = ModuleTransferline("TL column 4", 1.5, 0.1e-3, 0.1e-6, "Wax", 300.0)
	# system
	sys = System(g, pp, modules, Options())
end

# ╔═╡ 24de1d54-5696-491d-943e-107f7099c7ef
typeof(g)

# ╔═╡ 0415d6bc-f28f-456c-a30f-3c2d48c90a9e
md"""
### [x] Loop
"""

# ╔═╡ e6a9dd2f-48b6-4420-8f30-1921b9013a0a
begin
	# setup of the Graph -> ATTENTION to the right order of the modules
	# look at `collect(edges(sys_loop.g))`
	g_loop = SimpleDiGraph(6)
	add_edge!(g_loop, 1, 2) # Inj->TL1->Split 
	add_edge!(g_loop, 2, 3) # Split -> TL2
	add_edge!(g_loop, 2, 5) # GC2
	add_edge!(g_loop, 3, 4) # GC1
	add_edge!(g_loop, 4, 5) # TL3
	add_edge!(g_loop, 5, 6) # TL4 -> Det
	# pressure points
	pp_loop = Array{PressurePoint}(undef, nv(g_loop))
	pp_loop[1] = PressurePoint("p₁", [0.0, 1800.0], [200000.0, 200000.0]) # inlet 
	pp_loop[2] = PressurePoint("p₂", [0.0, 1800.0], [NaN, NaN]) # 
	pp_loop[3] = PressurePoint("p₃", [0.0, 1800.0], [NaN, NaN]) #  
	pp_loop[4] = PressurePoint("p₄", [0.0, 1800.0], [NaN, NaN]) #
	pp_loop[5] = PressurePoint("p₅", [0.0, 1800.0], [NaN, NaN]) #
	pp_loop[6] = PressurePoint("p₆", [0.0, 1800.0], [101300.0, 101300.0]) # outlet
	# modules
	modules_loop = Array{AbstractModule}(undef, ne(g_loop))
	modules_loop[1] = ModuleTransferline("1 TL1", 0.5, 0.25e-3, 0.25e-6, "Rxi5ms", 300.0)
	modules_loop[2] = ModuleTransferline("2 TL2", 1.0, 0.25e-3, 0.25e-6, "Rxi15ms", 300.0)
	modules_loop[3] = ModuleColumn("5 GC2", 2.0, 0.1e-3, 0.1e-6, "Rxi17SilMS", default_TP) 
	modules_loop[4] = ModuleColumn("3 GC1", 10.0, 0.25e-3, 0.25e-6, "Rxi5ms", default_TP)
	modules_loop[5] = ModuleColumn("4 TL3", 1.0, 0.1e-3, 0.1e-6, "Rxi5ms", default_TP)
	
	modules_loop[6] = ModuleTransferline("6 TL4", 1.5, 0.1e-3, 0.1e-6, "", 300.0)
	# system
	sys_loop = System(g_loop, pp_loop, modules_loop, Options())
end

# ╔═╡ 117717f5-efa7-4e8a-aaa3-5cb245dec9c9
collect(edges(sys_loop.g))

# ╔═╡ 755a4eb3-50f4-4188-8804-be4dcf7f547b
md"""
## [x] Common programs
"""

# ╔═╡ 4ee06f0e-4c8d-44e8-b8ce-14c0e75e0ed2
md"""
### [x] Common time steps
"""

# ╔═╡ 6f0ab9f0-27a4-4a84-94b6-ac4f155ef182
function common_timesteps(sys)
	com_timesteps = []
	for i=1:nv(sys.g)
		com_timesteps = GasChromatographySimulator.common_time_steps(com_timesteps, sys.pressurepoints[i].timesteps)
	end
	for i=1:ne(sys.g)
		if typeof(sys.modules[i]) == ModuleColumn
			com_timesteps = GasChromatographySimulator.common_time_steps(com_timesteps, sys.modules[i].temperature_program.timesteps)
		end
	end
	return com_timesteps
end

# ╔═╡ 380e891d-ce28-4b26-8d54-99ff1ad2621f
md"""
### [x] Index of modules with temperature program
"""

# ╔═╡ b8fe811d-350c-4806-99ef-481c40cd7445
function index_modules_with_temperature_program(sys)
	i_tempprog = Int[]
	for i=1:ne(sys.g)
		if typeof(sys.modules[i])==ModuleColumn
			push!(i_tempprog, i)
		end
	end
	return i_tempprog
end

# ╔═╡ 031d0b0e-6780-46ec-8ee7-398ebd4a62e9
md"""
### [x] Matching programs
"""

# ╔═╡ e727f789-e08e-435d-aede-59867c373e67
function match_programs(sys)
	com_times = common_timesteps(sys)
	new_press_steps = Array{Array{Float64,1}}(undef, nv(sys.g))
	for i=1:nv(sys.g)
		new_press_steps[i] = GasChromatographySimulator.new_value_steps(sys.pressurepoints[i].pressure_steps, sys.pressurepoints[i].timesteps, com_times)
	end
	i_tempprog = index_modules_with_temperature_program(sys)
	new_temp_steps = Array{Array{Float64,1}}(undef, length(i_tempprog))
	for i=1:length(i_tempprog)
		new_temp_steps[i] = GasChromatographySimulator.new_value_steps(sys.modules[i_tempprog[i]].temperature_program.temperaturesteps, sys.modules[i_tempprog[i]].temperature_program.timesteps, com_times)
	end
	# add for gradient new_a_gf
	return com_times, new_press_steps, new_temp_steps, i_tempprog
end

# ╔═╡ 870948f5-e827-4315-9786-e9b3d092642f
ts, ps, Ts, iT = match_programs(sys_series)

# ╔═╡ 20a1fe9a-4aca-4a16-869a-cbda38df6db4
md"""
### [x] Updating system variable
"""

# ╔═╡ d0f3d270-71b6-40a5-b9bc-49a48657cf22
function update_system(sys)
	new_timesteps, new_pressuresteps, new_temperaturesteps, index_module_tempprog = match_programs(sys)
	new_pp = Array{PressurePoint}(undef, nv(sys.g))
	for i=1:nv(sys.g)
		new_pp[i] = PressurePoint(sys.pressurepoints[i].name, new_timesteps, new_pressuresteps[i])
	end
	new_modules = Array{AbstractModule}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		if typeof(sys.modules[i]) == ModuleTransferline
			new_modules[i] = sys.modules[i]
		elseif typeof(sys.modules[i]) == ModuleColumn
			# add/modify for gradient
			ii = findfirst(index_module_tempprog.==i)
			new_tp = TemperatureProgram(new_timesteps, new_temperaturesteps[ii])
			new_modules[i] = ModuleColumn(sys.modules[i].name, sys.modules[i].length, sys.modules[i].diameter, sys.modules[i].film_thickness, sys.modules[i].stationary_phase, new_tp)
		end
	end
	new_sys = System(sys.g, new_pp, new_modules, sys.options)
end

# ╔═╡ ed4fa21c-fce3-4e12-89fd-73a780af0250
nsys0 = update_system(sys0)

# ╔═╡ 2b87a94a-d9ac-4a56-83a2-e70abcbd0a98
nsys_series = update_system(sys_series)

# ╔═╡ 5159c3d6-9bc7-4377-9839-9f0c26ce0276
nsys_split = update_system(sys_split)

# ╔═╡ 1fae0996-4183-445a-880c-330e52405265
nsys = update_system(sys)

# ╔═╡ 23e7b7f7-5f24-428f-8faa-b3778b5ba7e1
nsys_loop = update_system(sys_loop)

# ╔═╡ eebff747-334b-402d-a0ce-107f98207f30
md"""
## [x] Calculate unkown pressures
"""

# ╔═╡ ed23953d-214c-49a1-ac91-5b31abd935d9
md"""
### [x] Setup the balance equations
"""

# ╔═╡ 437fb458-4d5c-49a1-b029-3e48e6e92e5c
function flow_balance(g, i_n, P², κ)
	#@variables P²[1:nv(g)], κ[1:ne(g)]
	E = collect(edges(g))
	srcE = src.(E)
	dstE = dst.(E)
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
end

# ╔═╡ 7e4e16be-0972-43a7-bf03-f342bdc975af
function flow_balance(sys)
	@variables P²[1:nv(sys.g)], κ[1:ne(sys.g)]
	inner_V = findall(Graphs.degree(sys.g).>1)
	bal_eq = Array{Symbolics.Equation}(undef, length(inner_V))
	for i=1:length(inner_V)
		bal_eq[i] = flow_balance(sys.g, inner_V[i], P², κ) ~ 0
	end
	return bal_eq
end

# ╔═╡ 2cfd8772-1812-4af5-8ac7-fd360ee798f1
flow_balance(nsys0)

# ╔═╡ 04970711-8673-41bf-9d55-14aefabb6ed4
flow_balance(nsys_series)

# ╔═╡ f50d8baa-071c-4e3c-a130-45b504f1fa42
flow_balance(nsys_split)

# ╔═╡ ea7bc95b-5207-4096-9145-e9c339dec50f
flow_balance(nsys)

# ╔═╡ 790a86b7-3e2d-45b1-8a06-64fb100080f2
md"""
### [x] Determine the unknown pressures
"""

# ╔═╡ 042fa3b0-36c2-48b0-80f7-e0e01edd13a5
function unknown_p(sys)
	i_unknown_pressures = Int[]
	for i=1:nv(sys.g)
		if isnan(sys.pressurepoints[i].pressure_steps[1])
			push!(i_unknown_pressures, i)
		end
	end
	return i_unknown_pressures
end

# ╔═╡ f08c5dfe-e554-4c2d-a0e0-0826c166fa21
unknown_p(nsys0)

# ╔═╡ 2da5e914-1e07-411c-b942-f4bbba14f835
unknown_p(nsys_series)

# ╔═╡ ba57147b-a813-478c-93a8-dc429fc24a7c
unknown_p(nsys_split)

# ╔═╡ 509829b9-8edb-4495-9ea0-c3de242f0436
i_unknown_p = unknown_p(nsys)

# ╔═╡ 39b179a2-849a-4466-916c-9abffe084f00
md"""
### [x] Solve the balance equations for the unknown pressures
"""

# ╔═╡ 9e1706f7-919e-4c26-b356-0f2834c3f2a5
function solve_balance(sys)
	@variables P²[1:nv(sys.g)], κ[1:ne(sys.g)]
	i_unknown_p = unknown_p(sys)
	sol = Symbolics.solve_for(flow_balance(sys), [P²[i_unknown_p[i]] for i=1:length(i_unknown_p)])
	return sol
end

# ╔═╡ e38a0de6-34ca-412b-9367-989cdae07788
solve_balance(nsys0)

# ╔═╡ 19247928-adb9-45b4-a8e7-4b0edce080d9
solve_balance(nsys_series)

# ╔═╡ 7b7e7cff-491a-46f1-943a-ee5d2b195bdb
solve_balance(nsys_split)

# ╔═╡ 299615c8-5c97-4bda-a249-9c55197c3371
solve_balance(nsys)

# ╔═╡ f04645f1-b1ba-4235-ab91-0fda740a931e
md"""
### [x] Calculate flow restrictions κ
"""

# ╔═╡ 55c016a4-d09c-4dfa-b886-c650baf0b645
function flow_restriction(t, module_, mobile_phase)
	if typeof(module_) == ModuleColumn
		T_itp = GasChromatographySimulator.temperature_interpolation(module_.temperature_program.timesteps, module_.temperature_program.temperaturesteps, module_.temperature_program.gradient_function, module_.length)
	elseif typeof(module_) == ModuleTransferline
		gf(x) = [zero(x), zero(x)]
		T_itp = GasChromatographySimulator.temperature_interpolation([0.0, 3600.0], [module_.temperature, module_.temperature], gf, module_.length)
	end
	κ = GasChromatographySimulator.flow_restriction(module_.length, t, T_itp, module_.diameter, mobile_phase)
	return κ
end

# ╔═╡ c3f170cf-2776-4827-a447-5d4ace14fe2e
function flow_restrictions(sys)
	kappas = Array{Function}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		f(t) = flow_restriction(t, sys.modules[i], sys.options.mobile_phase)
		kappas[i] = f
	end
	return kappas
end

# ╔═╡ 7359bb10-a2f0-4010-a5c7-25822327b96d
κ_split = flow_restrictions(nsys_split)

# ╔═╡ 18ca6d06-b94a-453a-9c61-255033998525
md"""
### [x] Calculate squared pressures $p^2$
"""

# ╔═╡ 9b7209c7-384c-414e-b198-13daa3d47d68
function pressures_squared(sys)
	#p² = Array{Interpolations.Extrapolation}(undef, nv(sys.g))
	p² = Array{Any}(undef, nv(sys.g))
	for i=1:nv(sys.g)
		f = if isnan(sys.pressurepoints[i].pressure_steps[1])
			GasChromatographySimulator.steps_interpolation(sys.pressurepoints[i].timesteps, NaN.*ones(length(sys.pressurepoints[i].timesteps)))
		else
			GasChromatographySimulator.steps_interpolation(sys.pressurepoints[i].timesteps, identity.(sys.pressurepoints[i].pressure_steps).^2)
		end
		p²[i] = f
	end
	return p²
end

# ╔═╡ 701bb671-39a6-480c-b715-6304a09058f0
p²_split = pressures_squared(nsys_split)

# ╔═╡ b4a3ad9e-bff2-44a2-9275-48dff2be585f
p²_split[2](0)

# ╔═╡ f1cdee20-c959-4f7a-b6d0-7a0ff7b15d86
md"""
### [x] Substitut values
"""

# ╔═╡ 1a3c53b3-840a-4130-b40e-856f4a64006f
function solve_pressure(sys)
	@variables P²[1:nv(sys.g)], κ[1:ne(sys.g)]
	balance = flow_balance(sys)
	unknowns = unknown_p(sys)
	solutions = solve_balance(sys)
	κs = flow_restrictions(sys)
	p²s = pressures_squared(sys)
	p_solution = Array{Function}(undef, length(unknowns))
	for i=1:length(unknowns)
		sub_dict(t) = merge(Dict((P²[j] => p²s[j](t) for j=setdiff(1:nv(sys.g), unknowns))), Dict((κ[j] => κs[j](t) for j=1:ne(sys.g))))
		f(t) = sqrt(substitute(solutions[i], sub_dict(t)))
		p_solution[i] = f
	end
	return p_solution, unknowns
end

# ╔═╡ 30853cdb-06c7-4599-9ac2-76e78992922b
solve_pressure(nsys0)

# ╔═╡ 6bbdde5a-4522-45a6-aaec-5657371d6ca3
solve_pressure(nsys_series)

# ╔═╡ f529973b-5338-4987-aa59-0d342918a4aa
solve_pressure(nsys_split)

# ╔═╡ ac054b76-0e44-41cd-bb91-afcbb897bf8b
solve_pressure(nsys)

# ╔═╡ 1b9949f1-48e9-4232-858a-4712f362b4d3
md"""
### [x] Plot the pressures
"""

# ╔═╡ 3611cfbc-b987-46f5-9909-4bb345e28c54
pres, unk = solve_pressure(nsys)

# ╔═╡ ae40c445-e1e6-4b5c-b23a-190fab0e0158
p²s = pressures_squared(nsys)

# ╔═╡ e0a28475-9195-44cd-ae90-38e7092b87c4
com_timesteps = common_timesteps(nsys)

# ╔═╡ b8f09841-fe72-4f89-befd-625758d9a354
t = 0:sum(com_timesteps)/100:sum(com_timesteps)

# ╔═╡ 30fcbde8-70e9-453c-b621-f277209bd609
function pressure_functions(sys)
	pres, unk = solve_pressure(sys)
	p²s = pressures_squared(sys)
	p_func = Array{Any}(undef, nv(sys.g))
	for i=1:nv(sys.g)
		f(t) = if i in unk
			pres[findfirst(unk.==i)](t)
		else
			sqrt(p²s[i](t))
		end
	p_func[i] = f
	end
	return p_func
end

# ╔═╡ 7d032b2c-114d-4caa-84bd-de7e2ac4ac0d
p_func = pressure_functions(nsys)

# ╔═╡ 61500d0e-ca1f-49ab-9782-296186808108
begin
	plotly()
	p_pres = Plots.plot(xlabel="time in s", ylabel="pressure in Pa(a)")
	for i=1:nv(sys.g)
		#if i in unk
		#	ii = findfirst(unk.==i)
		#	Plots.plot!(p_pres, t, pres[ii].(t), label="p$(i)")
		#else
		#	Plots.plot!(p_pres, t, sqrt.(p²s[i].(t)), label="p$(i)")
		#end
		Plots.plot!(p_pres, t, p_func[i].(t), label="p$(i)")
	end
	p_pres
end

# ╔═╡ 80c260d7-f7c3-45a2-8880-960e3ece5c54
md"""
### [x] Plot flows
"""

# ╔═╡ e838d374-cdfe-4338-b850-6f0fa11a897d
function flow_functions(sys)
	p_func = pressure_functions(sys)
	F_func = Array{Function}(undef, ne(sys.g))
	E = collect(edges(sys.g))
	srcE = src.(E)
	dstE = dst.(E)
	for i=1:ne(sys.g)
		pin(t) = p_func[srcE[i]](t)
		pout(t) = p_func[dstE[i]](t)
		if typeof(sys.modules[i]) == ModuleColumn
			T_itp = GasChromatographySimulator.temperature_interpolation(sys.modules[i].temperature_program.timesteps, sys.modules[i].temperature_program.temperaturesteps, sys.modules[i].temperature_program.gradient_function, sys.modules[i].length)
		elseif typeof(sys.modules[i]) == ModuleTransferline
			gf(x) = [zero(x), zero(x)]
			T_itp = GasChromatographySimulator.temperature_interpolation([0.0, 3600.0], [sys.modules[i].temperature, sys.modules[i].temperature], gf, sys.modules[i].length)
		end
		f(t) = GasChromatographySimulator.flow(t, T_itp, pin, pout, sys.modules[i].length, sys.modules[i].diameter, sys.options.mobile_phase; ng=sys.options.ng, vis=sys.options.vis, control=sys.options.control)
		F_func[i] = f
	end
	return F_func
end

# ╔═╡ 44b8bd75-5893-4e8b-9eca-297fd0f7c137
F_func = flow_functions(nsys)

# ╔═╡ 70930065-e396-4418-b6ce-80a30cf85c07
begin
	plotly()
	p_flow = Plots.plot(xlabel="time in s", ylabel="flow in mL/min")
	for i=1:ne(nsys.g)
		Plots.plot!(p_flow, t, F_func[i].(t)*60*1e6, label="F_$(nsys.modules[i].name)")
	end
	p_flow
end

# ╔═╡ fc2ea452-c71c-4408-adc9-51e0a633f14e
function plot_graph(sys, t; lay = Stress(), color=:lightblue, node_size=80)
	p_func = pressure_functions(sys)
	F_func = flow_functions(sys)
	fig, ax, p = GraphMakie.graphplot(sys.g, 
						layout=lay,
						nlabels = [sys.pressurepoints[i].name*"\n$(trunc(Int, p_func[i](t)/1000))kPa" for i in 1:nv(sys.g)],
						nlabels_align=(:center,:center),
						node_size = [node_size for i in 1:nv(sys.g)],
						node_color = [color for i in 1:nv(sys.g)],
						elabels = [sys.modules[i].name*"\n $(round(F_func[i](t)*60*1e6; sigdigits=2)) mL/min" for i in 1:ne(sys.g)],
						elabels_distance = 20
					)
	hidedecorations!(ax)
	hidespines!(ax)
	ax.aspect = DataAspect()
	return fig
end

# ╔═╡ b85e57b4-645b-452b-bae7-e55ee6cedc0e
plot_graph(sys0; color=:pink, node_size=20)

# ╔═╡ e923dc1f-c289-46f7-9d97-b2d4d8e04488
plot_graph(sys_series; color=:red)

# ╔═╡ 6d3fbb35-3f06-45c1-aa78-1c80997c1e21
plot_graph(sys_split; color=:green)

# ╔═╡ ab0f33cd-55bd-48e4-a1f2-ecd5e4720805
plot_graph(sys)

# ╔═╡ fc64bb5e-08e7-4a65-bc03-79e6ebaa616f
plot_graph(nsys)

# ╔═╡ 85ed2383-4038-44ac-9e49-c71d005a9b27
plot_graph(nsys, 0)

# ╔═╡ f7c369a1-32ce-41a1-9170-fa39c18f82a8
plot_graph(nsys, 1200)

# ╔═╡ 5be9c1cb-5e40-4e59-8cf7-b013a896f2ae
plot_graph(nsys0, 0)

# ╔═╡ 8108835b-1763-4df5-860b-86eebb86260f
plot_graph(nsys_series, 0)

# ╔═╡ 1e7d2b50-dc43-40fb-9e2e-00e3f3f00763
plot_graph(nsys_split, 0)

# ╔═╡ 6ba5df8a-caeb-452e-bff0-4604f922fdf3
plot_graph(nsys_loop, 0)

# ╔═╡ 71c150ce-91a7-40ca-9c15-bad5593f1f4a
begin
	plotly()
	F_func_ = flow_functions(nsys_split)
	p_flow_ = Plots.plot(xlabel="time in s", ylabel="flow in mL/min")
	for i=1:ne(nsys_split.g)
		Plots.plot!(p_flow_, t, F_func_[i].(t)*60*1e6, label="F_$(nsys_split.modules[i].name)")
	end
	p_flow_
end

# ╔═╡ ed6a1cea-c9c7-44a7-8c3d-471fe1605721
md"""
## [] Transform system to GasChromatographySimulator structures
"""

# ╔═╡ 4942e000-421e-47bd-9829-92914de7d1bd
md"""
- every module becomes one set of parameter
- have to define some solutes
- initial time t and peak variance τ² is NaN (except for the first module, there they are 0)
"""

# ╔═╡ 59cabe32-fe06-4589-a403-cdc86501aa73
md"""
### [x] Load database
"""

# ╔═╡ 29486f2f-3c72-4226-83ad-890ed9d83f66
db_path = joinpath(dirname(pwd()), "data", "database_Kcentric.csv")

# ╔═╡ 8d69770e-c243-428b-a634-4da77c544235
db_dataframe = DataFrame(CSV.File(db_path, header=1, silencewarnings=true))

# ╔═╡ 34ca8b05-3b6b-4636-9386-0f0ebe836e52
md"""
### [x] Default solutes
"""

# ╔═╡ fdde9769-0822-47bb-890e-55612adda3b3
function all_stationary_phases(sys)
	stat_phases = Array{String}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		stat_phases[i] = sys.modules[i].stationary_phase
	end
	return stat_phases
end

# ╔═╡ af8f790b-afbd-41e2-a703-40a6aa89e951
function common_solutes(db, sys)
	# gives the soultes with common stationary phases between the database 'db' and the
	# GC-System 'GCsys'
	usp = setdiff(unique(all_stationary_phases(sys)), [""])
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
		common_solutes = DataFrame(Name=common_db.Name, CAS=common_db.CAS)
	end	
	return common_solutes
end

# ╔═╡ 6991ce74-f615-462d-8313-b5ee52e015a0
unique(all_stationary_phases(nsys0))

# ╔═╡ 4bf9f619-17b8-421e-ae59-02d04c70a8c2
common_solutes(db_dataframe, nsys0)

# ╔═╡ adb2715d-5cb4-4642-9236-600388ad4e5f
unique(all_stationary_phases(nsys_series))

# ╔═╡ fde79b23-ec4b-45c0-b874-9b59d96682b9
common_solutes(db_dataframe, nsys_series)

# ╔═╡ cb7b96d1-5f7b-4099-a020-a067c580900d
unique(all_stationary_phases(nsys_split))

# ╔═╡ 535ead23-da4c-4ab1-b981-d59c79685773
common_solutes(db_dataframe, nsys_split)

# ╔═╡ e57591c8-c66e-4544-9674-cc0d0337c593
unique(all_stationary_phases(nsys))

# ╔═╡ 4771a723-bea4-4f1f-be27-c83fee15b284
common_solutes(db_dataframe, nsys)

# ╔═╡ 53493ee2-1719-4a8f-a1ac-66bb0ad0cc49
selected_solutes = ["Octane", "Decane", "Dodecane", "Tetradecane", "2-Octanone", "2-Octanol", "5-Nonanol", "2-Nonanol"]

# ╔═╡ da7f7047-cb19-4a2d-aff2-7c03cea49c31
GasChromatographySimulator.load_solute_database(db_dataframe, "Rxi5ms", nsys.options.mobile_phase, selected_solutes, NaN.*ones(length(selected_solutes)), NaN.*ones(length(selected_solutes)))

# ╔═╡ b29455bc-05d7-4d1d-8b52-0ad2836da6a0
md"""
### [] Modules to parameters
"""

# ╔═╡ 54183ac7-d928-4250-aeaf-75adf5c0bc42
parameters = Array{GasChromatographySimulator.Parameters}(undef, ne(nsys.g))

# ╔═╡ 1ebe2956-3706-45b4-91f5-804e541bcc91
nsys

# ╔═╡ 036f3a55-2fd5-48e1-adfd-1dd2ab3d03eb
# delete the definition of the type as `Interpolations.Extrapolation` for T_itp, Fpin_itp and pout_itp in GasChromatographySimulator.Program

# ╔═╡ 887ec730-dde8-4fa1-8671-f0815af0c466
begin
	E = collect(edges(nsys.g))
	srcE = src.(E)
	dstE = dst.(E)
	for i=1:ne(nsys.g)
		col = GasChromatographySimulator.Column(nsys.modules[i].length, nsys.modules[i].diameter, [nsys.modules[i].diameter], nsys.modules[i].film_thickness, [nsys.modules[i].film_thickness], nsys.modules[i].stationary_phase, nsys.options.mobile_phase)

		pin_steps = nsys.pressurepoints[srcE[i]].pressure_steps
		pout_steps = nsys.pressurepoints[dstE[i]].pressure_steps
		pin_itp = p_func[srcE[i]]
		pout_itp = p_func[dstE[i]]	
		if typeof(nsys.modules[i]) == ModuleColumn
			time_steps = nsys.modules[i].temperature_program.timesteps
			temp_steps = nsys.modules[i].temperature_program.temperaturesteps	
			gf = nsys.modules[i].temperature_program.gradient_function
			a_gf = nsys.modules[i].temperature_program.a_gradient_function
			T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, nsys.modules[i].length)		
		elseif typeof(nsys.modules[i]) == ModuleTransferline
			time_steps = common_timesteps(nsys)
			temp_steps = nsys.modules[i].temperature.*ones(length(time_steps))
			gf(x) = zero(x).*ones(length(time_steps))
			a_gf = [zeros(length(time_steps)) zeros(length(time_steps)) ones(length(time_steps)) zeros(length(time_steps))]
			T_itp = GasChromatographySimulator.temperature_interpolation([0.0, 3600.0], [sys.modules[i].temperature, sys.modules[i].temperature], gf, sys.modules[i].length)
		end
		prog = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)

		sub = GasChromatographySimulator.load_solute_database(db_dataframe, nsys.modules[i].stationary_phase, nsys.options.mobile_phase, selected_solutes, NaN.*ones(length(selected_solutes)), NaN.*ones(length(selected_solutes)))

		opt = GasChromatographySimulator.Options(nsys.options.alg, nsys.options.abstol, nsys.options..reltol, nsys.options.Tcontrol, nsys.options.odesys, nsys.options.ng, nsys.options.vis, nsys.options.control, nsys.options.k_th)

		parameters[i] = GasChromatographySimulator.Parameters(col, prog, sub, opt)
	end
end

# ╔═╡ 44ebd883-232c-4a3c-a2a3-f7c8c85ad8e3
nsys.options

# ╔═╡ 8d7befb0-6225-404c-8dd2-12834978bde1
md"""
## [] Simulate the system
"""

# ╔═╡ 6dd088fa-5907-400d-bf13-4890a224401e
md"""
# End
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
GasChromatographySimulator = "dd82b6e2-56ef-419d-b271-0be268cb65f5"
GraphMakie = "1ecd5474-83a3-4783-bb4f-06765db800d2"
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
NetworkLayout = "46757867-2c16-5918-afeb-47bfcb05e46a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
CSV = "~0.10.8"
CairoMakie = "~0.10.1"
DataFrames = "~1.4.4"
GasChromatographySimulator = "~0.3.14"
GraphMakie = "~0.5.0"
Graphs = "~1.7.4"
Interpolations = "~0.14.7"
NetworkLayout = "~0.4.4"
Plots = "~1.38.0"
PlutoUI = "~0.7.49"
Symbolics = "~4.14.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "14b710f2694e0ef9409042ed6640dee874dc4401"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Markdown", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "7772df04fda9bc25a44c9ef61e9dc7c92bb35d86"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.27.7"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "69f7020bd72f069c219b5e8c236c1fa90d2cb409"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.2.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "52b3b436f8f73133d7bc3a6c71ee7ed6ab2ab754"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.3"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterface]]
deps = ["ArrayInterfaceCore", "Compat", "IfElse", "LinearAlgebra", "Static"]
git-tree-sha1 = "6d0918cb9c0d3db7fe56bea2bc8638fc4014ac35"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "6.0.24"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "14c3f84a763848906ac681f94cf469a851601d92"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.28"

[[deps.ArrayInterfaceGPUArrays]]
deps = ["Adapt", "ArrayInterfaceCore", "GPUArraysCore", "LinearAlgebra"]
git-tree-sha1 = "fc114f550b93d4c79632c2ada2924635aabfa5ed"
uuid = "6ba088a2-8465-4c0a-af30-387133b534db"
version = "0.2.2"

[[deps.ArrayInterfaceOffsetArrays]]
deps = ["ArrayInterface", "OffsetArrays", "Static"]
git-tree-sha1 = "3d1a9a01976971063b3930d1aed1d9c4af0817f8"
uuid = "015c0d05-e682-4f19-8f0a-679ce4c54826"
version = "0.1.7"

[[deps.ArrayInterfaceStaticArrays]]
deps = ["Adapt", "ArrayInterface", "ArrayInterfaceCore", "ArrayInterfaceStaticArraysCore", "LinearAlgebra", "Static", "StaticArrays"]
git-tree-sha1 = "f12dc65aef03d0a49650b20b2fdaf184928fd886"
uuid = "b0d46f97-bff5-4637-a19a-dd75974142cd"
version = "0.1.5"

[[deps.ArrayInterfaceStaticArraysCore]]
deps = ["Adapt", "ArrayInterfaceCore", "LinearAlgebra", "StaticArraysCore"]
git-tree-sha1 = "93c8ba53d8d26e124a5a8d4ec914c3a16e6a0970"
uuid = "dd5226c6-a4d4-4bc7-8575-46859f9c95b9"
version = "0.1.3"

[[deps.Arrow]]
deps = ["ArrowTypes", "BitIntegers", "CodecLz4", "CodecZstd", "DataAPI", "Dates", "LoggingExtras", "Mmap", "PooledArrays", "SentinelArrays", "Tables", "TimeZones", "UUIDs", "WorkerUtilities"]
git-tree-sha1 = "3d04ab3584ece56c39397e01b55e1bd4fb8f0b30"
uuid = "69666777-d1a9-59fb-9406-91d4454c9d45"
version = "2.4.1"

[[deps.ArrowTypes]]
deps = ["UUIDs"]
git-tree-sha1 = "a0633b6d6efabf3f76dacd6eb1b3ec6c42ab0552"
uuid = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
version = "1.2.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AutoHashEquals]]
git-tree-sha1 = "45bb6705d93be619b81451bb2006b7ee5d4e4453"
uuid = "15f4f7f2-30c1-5605-9d31-71845cf9641f"
version = "0.2.0"

[[deps.Automa]]
deps = ["Printf", "ScanByte", "TranscodingStreams"]
git-tree-sha1 = "d50976f217489ce799e366d9561d56a98a30d7fe"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "0.8.2"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "1dd4d9f5beebac0c03446918741b1a03dc5e5788"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.6"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "Future", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables", "ZygoteRules"]
git-tree-sha1 = "7fe6d92c4f281cf4ca6f2fba0ce7b299742da7ca"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.37"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.Bijections]]
git-tree-sha1 = "fe4f8c5ee7f76f2198d5c2a06d3961c249cce7bd"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.4"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.BitIntegers]]
deps = ["Random"]
git-tree-sha1 = "fc54d5837033a170f3bad307f993e156eefc345f"
uuid = "c3b6d118-76ef-56ca-8cc7-ebb389d030a1"
version = "0.2.7"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "0c5f81f47bbbcf4aea7b2959135713459170798b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.5"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "Static"]
git-tree-sha1 = "a7157ab6bcda173f533db4c93fc8a27a48843757"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.1.30"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "SnoopPrecompile", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "8c73e96bd6817c2597cfd5615b91fca5deccf1af"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.8"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["Base64", "Cairo", "Colors", "FFTW", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "SHA", "SnoopPrecompile"]
git-tree-sha1 = "439517f69683932a078b2976ca040e21dd18598c"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.10.1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRules]]
deps = ["Adapt", "ChainRulesCore", "Compat", "Distributed", "GPUArraysCore", "IrrationalConstants", "LinearAlgebra", "Random", "RealDot", "SparseArrays", "Statistics", "StructArrays"]
git-tree-sha1 = "99a39b0f807499510e2ea14b0eef8422082aa372"
uuid = "082447d4-558c-5d27-93f4-14fc19e9eca2"
version = "1.46.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e7ff6cadf743c098e08fca25c91103ee4303c9bb"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.6"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.ChemicalIdentifiers]]
deps = ["Arrow", "Downloads", "Preferences", "Scratch", "UUIDs", "Unicode"]
git-tree-sha1 = "f55715e75fb14cbe72808d43b75c4d66e5200daf"
uuid = "fa4ea961-1416-484e-bda2-883ee1634ba5"
version = "0.1.7"

[[deps.CloseOpenIntervals]]
deps = ["ArrayInterface", "Static"]
git-tree-sha1 = "d61300b9895f129f4bd684b2aff97cf319b6c493"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.11"

[[deps.CodecLz4]]
deps = ["Lz4_jll", "TranscodingStreams"]
git-tree-sha1 = "59fe0cb37784288d6b9f1baebddbf75457395d40"
uuid = "5ba52731-8f18-5e0d-9241-30f10d1ec561"
version = "0.4.0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.CodecZstd]]
deps = ["CEnum", "TranscodingStreams", "Zstd_jll"]
git-tree-sha1 = "849470b337d0fa8449c21061de922386f32949d9"
uuid = "6b39b394-51ab-5f42-8807-6242bab2b4c2"
version = "0.7.2"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random", "SnoopPrecompile"]
git-tree-sha1 = "aa3edc8f8dea6cbfa176ee12f7c2fc82f0608ed3"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.20.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "9441451ee712d1aec22edad62db1a9af3dc8d852"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.3"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "00a2cccc7f098ff3b66806862d275ca3db9e6e5a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.5.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "02d2316b7ffceff992f3096ae48c7829a8aa0638"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.3"

[[deps.CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d4f69885afa5e6149d0cab3818491565cf41446d"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.4.4"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffEqBase]]
deps = ["ArrayInterfaceCore", "ChainRulesCore", "DataStructures", "Distributions", "DocStringExtensions", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "MuladdMacro", "Parameters", "PreallocationTools", "Printf", "RecursiveArrayTools", "Reexport", "Requires", "SciMLBase", "Setfield", "SimpleNonlinearSolve", "SparseArrays", "Static", "StaticArrays", "Statistics", "Tricks", "ZygoteRules"]
git-tree-sha1 = "29777943a9e73c7d6b47d93830038ebdaacc18db"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.113.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "c5b6685d53f933c11404a3ae9822afe30d522494"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.12.2"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "a7756d098cbabec6b3ac44f369f74915e8cfd70a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.79"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "988e2db482abeb69efc76ae8b6eba2e93805ee70"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.15"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "d0fa82f39c2a5cdb3ee385ad52bc05c42cb4b9f0"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.4.5"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceGPUArrays", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "9837d3f3a904c7a7ab9337759c0093d3abea1d81"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.22.0"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "90630efff0894f8142308e334473eba54c433549"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.5.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "LinearAlgebra", "Polyester", "Static", "StrideArraysCore"]
git-tree-sha1 = "4bef892787c972913d4d84e7255400759bb650e5"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.2.4"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastLapackInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "7fbaf9f73cd4c8561702ea9b16acf3f99d913fe4"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "1.2.8"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "7be5f99f7d15578798f338f5433b6c432ea8037b"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "e27c4ebe80e8699540f2d6c805cc12203b614f12"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.20"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "9a0472ec2f5409db243160a8b030f94c380167a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.6"

[[deps.FiniteDiff]]
deps = ["ArrayInterfaceCore", "LinearAlgebra", "Requires", "Setfield", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "04ed1f0029b6b3af88343e439b995141cb0d0b8d"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.17.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "a69dd6db8a809f78846ff259298678f0d6212180"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.34"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "cabd77ab6a6fdff49bfd24af2ebe76e6e018a2b4"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.0.0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "38a92e40157100e796690421e34a11c107205c86"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "a5e6e7f12607e90d71b09e6ce2c965e41b337968"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.1"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GPUArrays]]
deps = ["Adapt", "GPUArraysCore", "LLVM", "LinearAlgebra", "Printf", "Random", "Reexport", "Serialization", "Statistics"]
git-tree-sha1 = "45d7deaf05cbb44116ba785d147c518ab46352d7"
uuid = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7"
version = "8.5.0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "6872f5ec8fd1a38880f027a26739d42dcda6691f"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.2"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "bcc737c4c3afc86f3bbc55eb1b9fabcee4ff2d81"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.71.2"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "64ef06fa8f814ff0d09ac31454f784c488e22b29"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.71.2+0"

[[deps.GasChromatographySimulator]]
deps = ["CSV", "ChemicalIdentifiers", "DataFrames", "ForwardDiff", "HypertextLiteral", "Integrals", "Interpolations", "OrdinaryDiffEq", "Plots", "PlutoUI", "Reexport", "UrlDownload"]
git-tree-sha1 = "9cc55b54adefd6746b3befca489f24f4332de84b"
uuid = "dd82b6e2-56ef-419d-b271-0be268cb65f5"
version = "0.3.14"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "fb69b2a645fa69ba5f474af09221b9308b160ce6"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.3"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "fb28b5dc239d0174d7297310ef7b84a11804dfab"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.0.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "fe9aea4ed3ec6afdfbeb5a4f39a2208909b162a6"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.5"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.GraphMakie]]
deps = ["GeometryBasics", "Graphs", "LinearAlgebra", "Makie", "NetworkLayout", "StaticArrays"]
git-tree-sha1 = "da596204780670d848c5bf35aff1f8580b885e09"
uuid = "1ecd5474-83a3-4783-bb4f-06765db800d2"
version = "0.5.0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "ba2d094a88b6b287bd25cfa86f301e7693ffae2f"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.7.4"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "678d136003ed5bceaab05cf64519e3f956ffa4ba"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.9.1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "Logging", "MultivariatePolynomials", "Primes", "Random"]
git-tree-sha1 = "47f0f03eddecd7ad59c42b1dd46d5f42916aff63"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.2.11"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9e1a5e9f3b81ad6a5c613d181664a0efc6fe6dd7"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.0"

[[deps.HCubature]]
deps = ["Combinatorics", "DataStructures", "LinearAlgebra", "QuadGK", "StaticArrays"]
git-tree-sha1 = "e95b36755023def6ebc3d269e6483efa8b2f7f65"
uuid = "19dc6840-f33b-545b-b366-655c7e3ffd49"
version = "1.5.1"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "fd9861adba6b9ae4b42582032d0936d456c8602d"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.6.3"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "f64b890b2efa4de81520d2b0fbdc9aadb65bdf53"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.13"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IRTools]]
deps = ["InteractiveUtils", "MacroTools", "Test"]
git-tree-sha1 = "2e99184fca5eb6f075944b04c22edec29beb4778"
uuid = "7869d1d1-7146-5819-86e3-90919afe41df"
version = "0.4.7"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "c54b581a83008dc7f292e205f4c409ab5caa0f04"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.10"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "342f789fd041a55166764c351da1710db97ce0e0"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.6"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "36cbaebed194b292590cba2593da27b34763804a"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.8"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "87f7662e03a649cffa2e05bf19c303e168732d3e"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "0cf92ec945125946352f3d46c96976ab972bde6f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.3.2"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "f366daebdfb079fd1fe4e3d560f99a0c892e15bc"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.0"

[[deps.Integrals]]
deps = ["ChainRulesCore", "CommonSolve", "Distributions", "ForwardDiff", "HCubature", "LinearAlgebra", "MonteCarloIntegration", "QuadGK", "Reexport", "ReverseDiff", "SciMLBase", "Zygote", "ZygoteRules"]
git-tree-sha1 = "6a1426c218f1dd1042f84892331f0baea6c2cc00"
uuid = "de52edbc-65ea-441a-8357-d3a637375a31"
version = "3.4.0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "3f91cd3f56ea48d4d2a75c2a65455c5fc74fa347"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.3"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.InvertedIndices]]
git-tree-sha1 = "82aec7a3dd64f4d9584659dc0b62ef7db2ef3e19"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.2.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "a77b273f1ddec645d1b7c4fd5fb98c8f90ad10a5"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "764164ed65c30738750965d55652db9c94c59bfe"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.4.0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "9816b296736292a80b9a3200eb7fbb57aaa3917a"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.5"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "dd90aacbfb622f898a97c2a4411ac49101ebab8a"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.0"

[[deps.KrylovKit]]
deps = ["ChainRulesCore", "GPUArraysCore", "LinearAlgebra", "Printf"]
git-tree-sha1 = "1a5e1d9941c783b0119897d29f2eb665d876ecf3"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.6.0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVM]]
deps = ["CEnum", "LLVMExtra_jll", "Libdl", "Printf", "Unicode"]
git-tree-sha1 = "088dd02b2797f0233d92583562ab669de8517fd1"
uuid = "929cbde3-209d-540e-8aea-75f648917ca0"
version = "4.14.1"

[[deps.LLVMExtra_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg", "TOML"]
git-tree-sha1 = "771bfe376249626d3ca12bcd58ba243d3f961576"
uuid = "dad2f222-ce93-54a1-a47d-0025e8a3acab"
version = "0.0.16+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterfaceCore", "ArrayInterfaceStaticArrays", "ArrayInterfaceStaticArraysCore", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "dae002226b59701dbafd7e2dd757df1bd83442fd"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.12.5"

[[deps.LambertW]]
git-tree-sha1 = "2d9f4009c486ef676646bca06419ac02061c088e"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.5"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "ab9aa169d2160129beb241cb2750ca499b4e90e9"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.17"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "ArrayInterfaceOffsetArrays", "ArrayInterfaceStaticArrays", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static"]
git-tree-sha1 = "7e34177793212f6d64d045ee47d2883f09fffacc"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.12"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearSolve]]
deps = ["ArrayInterfaceCore", "DocStringExtensions", "FastLapackInterface", "GPUArraysCore", "IterativeSolvers", "KLU", "Krylov", "KrylovKit", "LinearAlgebra", "Preferences", "RecursiveFactorization", "Reexport", "SciMLBase", "Setfield", "SnoopPrecompile", "SparseArrays", "Sparspak", "SuiteSparse", "UnPack"]
git-tree-sha1 = "cf1227e369513687658476e466a5b73a7c3dfa1f"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "1.33.0"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "946607f84feb96220f480e0422d3484c49c00239"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.19"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "5d4d2d9904227b8bd66386c1138cf4d5ffa826bf"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.9"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "ArrayInterfaceOffsetArrays", "ArrayInterfaceStaticArrays", "CPUSummary", "ChainRulesCore", "CloseOpenIntervals", "DocStringExtensions", "ForwardDiff", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDDualNumbers", "SIMDTypes", "SLEEFPirates", "SnoopPrecompile", "SpecialFunctions", "Static", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "0897c2d2138f090064559a3d90b12b8725d118c2"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.143"

[[deps.Lz4_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "5d494bc6e85c4c9b626ee0cab05daa4085486ab1"
uuid = "5ced341a-0733-55b8-9ab6-a4889d929147"
version = "1.9.3+0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "2ce8695e1e699b68702c03402672a69f54b8aca9"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Makie]]
deps = ["Animations", "Base64", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG", "FileIO", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "InteractiveUtils", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MakieCore", "Markdown", "Match", "MathTeXEngine", "MiniQhull", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "Printf", "Random", "RelocatableFolders", "Setfield", "Showoff", "SignedDistanceFields", "SnoopPrecompile", "SparseArrays", "StableHashTraits", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun"]
git-tree-sha1 = "20f42c8f4d70a795cb7927d7312b98a255209155"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.19.1"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "c5b3ce048ee73a08bbca1b9f4a776e64257611d5"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.6.1"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Match]]
git-tree-sha1 = "1d9bc5c1a6e7ee24effb93f175c9342f9154d97f"
uuid = "7eb4fadd-790c-5f42-8a69-bfa0b872bfbf"
version = "1.2.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "Test", "UnicodeFun"]
git-tree-sha1 = "f04120d9adf4f49be242db0b905bea0be32198d1"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.5.4"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Metatheory]]
deps = ["AutoHashEquals", "DataStructures", "Dates", "DocStringExtensions", "Parameters", "Reexport", "TermInterface", "ThreadsX", "TimerOutputs"]
git-tree-sha1 = "0f39bc7f71abdff12ead4fc4a7d998fb2f3c171f"
uuid = "e9d8d322-4543-424a-9be4-0cc815abe26c"
version = "1.3.5"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "4d5917a26ca33c66c8e5ca3247bd163624d35493"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.3"

[[deps.MiniQhull]]
deps = ["QhullMiniWrapper_jll"]
git-tree-sha1 = "9dc837d180ee49eeb7c8b77bb1c860452634b0d1"
uuid = "978d7f02-9e05-4691-894f-ae31a51d76ca"
version = "0.4.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "c272302b22479a24d1cf48c114ad702933414f80"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.5"

[[deps.MonteCarloIntegration]]
deps = ["Distributions", "Random"]
git-tree-sha1 = "3f78ebce296c927d5c854e83cccdb5dcb1845629"
uuid = "4886b29c-78c9-11e9-0a6e-41e1f4161f7b"
version = "0.0.3"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "393fc4d82a73c6fe0e2963dd7c882b09257be537"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.4.6"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "aa532179d4a643d4bd9f328589ca01fa20a0d197"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.1.0"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "5ae7ca23e13855b3aba94550f26146c01d259267"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.0"

[[deps.NetworkLayout]]
deps = ["GeometryBasics", "LinearAlgebra", "Random", "Requires", "SparseArrays"]
git-tree-sha1 = "cac8fc7ba64b699c678094fa630f49b80618f625"
uuid = "46757867-2c16-5918-afeb-47bfcb05e46a"
version = "0.4.4"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ArrayInterfaceCore", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SnoopPrecompile", "SparseArrays", "SparseDiffTools", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "7142ca5ab9bd7452cafb29f7d51f574a09d69052"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "1.1.1"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "f71d8950b724e9ff6110fc948dff5a329f901d64"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "923319661e9a22712f24596ce81c54fc0366f304"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.1+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "df6830e37943c7aaa10023471ca47fb3065cc3c4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.2"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6e9dba33f9f2c44e08a020b0caf6903be540004"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.19+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.OrdinaryDiffEq]]
deps = ["Adapt", "ArrayInterface", "ArrayInterfaceCore", "ArrayInterfaceGPUArrays", "ArrayInterfaceStaticArrays", "ArrayInterfaceStaticArraysCore", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "LinearAlgebra", "LinearSolve", "Logging", "LoopVectorization", "MacroTools", "MuladdMacro", "NLsolve", "NonlinearSolve", "Polyester", "PreallocationTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLNLSolve", "SimpleNonlinearSolve", "SnoopPrecompile", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "e1563399318752a2df41d08ab1033a772bd0fa4b"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.36.2"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "f809158b27eba0c18c269cf2a2be6ed751d3e81d"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.17"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "ec3edfe723df33528e085e632414499f26650501"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.0"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "84a314e3926ba9ec66ac097e3635e270986b0f10"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.9+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "6466e524967496866901a78fca3f2e9ea445a559"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.2"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f6cf8e7944e50901594838951729a1861e668cb8"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.2"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "5b7690dd212e026bbab1860016a6601cb077ab66"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.2"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SnoopPrecompile", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "513084afca53c9af3491c94224997768b9af37e8"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eadad7b14cf046de6eb41f13c9275e5aa2711ab6"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.49"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "a5071cd52fc3fc0a960b825ddeb64e352fdf41e1"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.6.20"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "43883d15c7cf16f340b9367c645cf88372f55641"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.1.13"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ForwardDiff"]
git-tree-sha1 = "758f3283aba57c53960c8e1900b4c724bf24ba74"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.8"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "LaTeXStrings", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "96f6db03ab535bdb901300f88335257b0018689d"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.2"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "311a2aa90a64076ea0fac2ad7492e914e6feeb81"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QhullMiniWrapper_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Qhull_jll"]
git-tree-sha1 = "607cf73c03f8a9f83b36db0b86a3a9c14179621f"
uuid = "460c41e3-6112-5d7f-b78c-b6823adb3f2d"
version = "1.0.0+1"

[[deps.Qhull_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "238dd7e2cc577281976b9681702174850f8d4cbc"
uuid = "784f63db-0788-585a-bace-daefebcd302b"
version = "8.0.1001+0"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "97aa253e65b784fd13e83774cadc95b38011d734"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.6.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "062986376ce6d394b23d5d90f01d81426113a3c9"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.3"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "18c35ed630d7229c5584b945641a73ca83fb5213"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.2"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase", "SnoopPrecompile"]
git-tree-sha1 = "e974477be88cb5e3040009f3767611bc6357846f"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.11"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceStaticArraysCore", "ChainRulesCore", "DocStringExtensions", "FillArrays", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "ZygoteRules"]
git-tree-sha1 = "66e6a85fd5469429a3ac30de1bd491e48a6bac00"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.34.1"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "SnoopPrecompile", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "2979cbb21580760431d2afb9b8f0f522899542f7"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.13"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Referenceables]]
deps = ["Adapt"]
git-tree-sha1 = "e681d3bfa49cd46c3c161505caddf20f0e62aaa9"
uuid = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"
version = "0.1.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.ReverseDiff]]
deps = ["ChainRulesCore", "DiffResults", "DiffRules", "ForwardDiff", "FunctionWrappers", "LinearAlgebra", "LogExpFunctions", "MacroTools", "NaNMath", "Random", "SpecialFunctions", "StaticArrays", "Statistics"]
git-tree-sha1 = "afc870db2b2c2df1ba3f7b199278bb071e4f6f90"
uuid = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
version = "1.14.4"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "50314d2ef65fce648975a8e80ae6d8409ebbf835"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.5"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "34edfe91375e5883875987e740c554b92c48fc41"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.3"

[[deps.SIMDDualNumbers]]
deps = ["ForwardDiff", "IfElse", "SLEEFPirates", "VectorizationBase"]
git-tree-sha1 = "dd4195d308df24f33fb10dde7c22103ba88887fa"
uuid = "3cdde19b-5bb0-4aaf-8931-af3e248e098b"
version = "0.1.1"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "c8679919df2d3c71f74451321f1efea6433536cc"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.37"

[[deps.ScanByte]]
deps = ["Libdl", "SIMD"]
git-tree-sha1 = "2436b15f376005e8790e318329560dcc67188e84"
uuid = "7b38b023-a4d7-4c5e-8d43-3f3097f304eb"
version = "0.3.3"

[[deps.SciMLBase]]
deps = ["ArrayInterfaceCore", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Preferences", "RecipesBase", "RecursiveArrayTools", "RuntimeGeneratedFunctions", "StaticArraysCore", "Statistics", "Tables"]
git-tree-sha1 = "fe89a8113ea445bcff9ee570077830674babb534"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.81.0"

[[deps.SciMLNLSolve]]
deps = ["LineSearches", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "b35d1f5d8afeee44e24915bb767e34fae867502f"
uuid = "e9a6253c-8580-4d32-9898-8661bb511710"
version = "0.1.1"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "efd23b378ea5f2db53a55ae53d3133de4e080aa9"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.16"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ArrayInterfaceCore", "FiniteDiff", "ForwardDiff", "Reexport", "SciMLBase", "SnoopPrecompile", "StaticArraysCore"]
git-tree-sha1 = "fc4b9f81a033cf6879c91bb7f5b3ff59008c7dd2"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "0.1.4"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

[[deps.SnoopPrecompile]]
git-tree-sha1 = "f604441450a3c0569830946e5b33b78c928e1a85"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SparseDiffTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceStaticArrays", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "4245283bee733122a9cb4545748d64e0c63337c0"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.30.0"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "dbe7d12e530416a26135d9566f00946aea1b2474"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.2"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

[[deps.StableHashTraits]]
deps = ["CRC32c", "Compat", "Dates", "SHA", "Tables", "TupleTools", "UUIDs"]
git-tree-sha1 = "0b8b801b8f03a329a4e86b44c5e8a7d7f4fe10a3"
uuid = "c5dd0088-6c3f-4803-b00e-f31a60c170fa"
version = "0.3.1"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "c35b107b61e7f34fa3f124026f2a9be97dea9e1c"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.3"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "6954a456979f23d05085727adb17c4551c19ecd1"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.12"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "ab6083f09b3e617e34a956b43e9d51b824206932"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.1.1"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "SIMDTypes", "Static", "ThreadingUtilities"]
git-tree-sha1 = "8e91e5618bbca975312313c39ff827ea8f802fe3"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.4.4"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "b03a3b745aa49b566f128977a7dd1be8711c5e71"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.14"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+0"

[[deps.SymbolicIndexingInterface]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "6b764c160547240d868be4e961a5037f47ad7379"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.2.1"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "Metatheory", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TermInterface", "TimerOutputs"]
git-tree-sha1 = "027b43d312f6d52187bb16c2d4f0588ddb8c4bb2"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "0.19.11"

[[deps.Symbolics]]
deps = ["ArrayInterfaceCore", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Markdown", "Metatheory", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TermInterface", "TreeViews"]
git-tree-sha1 = "111fbf43883d95989577133aeeb889f2040d0aea"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "4.14.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TermInterface]]
git-tree-sha1 = "7aa601f12708243987b88d1b453541a75e3d8c7a"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "0.2.3"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "f8629df51cab659d70d2e5618a430b4d3f37f2c3"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.0"

[[deps.ThreadsX]]
deps = ["ArgCheck", "BangBang", "ConstructionBase", "InitialValues", "MicroCollections", "Referenceables", "Setfield", "SplittablesBase", "Transducers"]
git-tree-sha1 = "34e6bcf36b9ed5d56489600cf9f3c16843fa2aa2"
uuid = "ac1d9e8a-700a-412c-b207-f0111f4b6c0d"
version = "0.1.11"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "7e6b0e3e571be0b4dd4d2a9a3a83b65c04351ccc"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.3"

[[deps.TimeZones]]
deps = ["Dates", "Downloads", "InlineStrings", "LazyArtifacts", "Mocking", "Printf", "RecipesBase", "Scratch", "Unicode"]
git-tree-sha1 = "a92ec4466fc6e3dd704e2668b5e7f24add36d242"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.9.1"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f2fd3f288dfc6f507b0c3a2eb3bac009251e548b"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.22"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "e4bdc63f5c6d62e80eb1c0043fcc0360d5950ff7"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.10"

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "c42fa452a60f022e9e087823b47e5a5f8adc53d5"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.75"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "SnoopPrecompile", "Static", "VectorizationBase"]
git-tree-sha1 = "766f5b1fc80f667ad6eab787af21e19cd65cb1b7"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.16"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.TupleTools]]
git-tree-sha1 = "3c712976c47707ff893cf6ba4354aa14db1d8938"
uuid = "9d95972d-f1c8-5527-a6e0-b4b365fa01f6"
version = "1.3.0"

[[deps.URIs]]
git-tree-sha1 = "ac00576f90d8a259f2c9d823e91d1de3fd44d348"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.UrlDownload]]
deps = ["HTTP", "ProgressMeter"]
git-tree-sha1 = "758aeefcf0bfdd04cd88e6e3bc9feb9e8c7d6d70"
uuid = "856ac37a-3032-4c1c-9122-f86d88358c8b"
version = "1.0.1"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static"]
git-tree-sha1 = "fc79d0f926592ecaeaee164f6a4ca81b51115c3b"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.56"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "ed8d92d9774b077c53e1da50fd81a36af3744c1c"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.Zygote]]
deps = ["AbstractFFTs", "ChainRules", "ChainRulesCore", "DiffRules", "Distributed", "FillArrays", "ForwardDiff", "GPUArrays", "GPUArraysCore", "IRTools", "InteractiveUtils", "LinearAlgebra", "LogExpFunctions", "MacroTools", "NaNMath", "Random", "Requires", "SparseArrays", "SpecialFunctions", "Statistics", "ZygoteRules"]
git-tree-sha1 = "a6f1287943ac05fae56fa06049d1a7846dfbc65f"
uuid = "e88e6eb3-aa80-5325-afca-941959d7151f"
version = "0.6.51"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╠═80af05a0-8c11-11ed-1008-cb6a5443493d
# ╟─8ace0fa0-52ef-4db8-92d0-7f99587bb065
# ╠═a6d8c5bc-32ea-4f7b-89f5-56372243d27d
# ╠═524d119e-9de7-43f0-b6dc-d719c71b7131
# ╠═c8d84ca7-b160-4a47-aa6c-720ccd46090b
# ╠═f0541435-7af0-408a-9d1a-378e57336989
# ╠═c2d60366-f54c-4b10-9803-d80e19a0855f
# ╠═a693515b-e24d-4609-b97d-2dbcb6248414
# ╠═46d410e9-b45a-40a0-932d-a025bba99a89
# ╠═d3107cc8-aa83-424c-be15-40de49084623
# ╠═932ffdbb-0dd6-48ca-98b2-af70f229d4f5
# ╠═73d0e4d5-76c6-4860-92f5-5888ae605861
# ╠═8a9c2704-5854-4aa4-9782-943e82cd38f6
# ╠═471193bb-c2e1-4d27-9f0f-a9d65820ef7e
# ╠═a177a0a9-7e67-4644-896e-8863d4bd15c6
# ╠═291ed3ff-e500-48e1-b1ec-2d28b5cbfa54
# ╠═40663239-9c5b-466d-9f43-76ab913d5750
# ╠═f11b32e2-c447-4f54-8f42-c1e81c9fb489
# ╠═e7181f3b-e148-4327-a420-0091826b327c
# ╠═24de1d54-5696-491d-943e-107f7099c7ef
# ╠═c5c010e6-80c6-445b-b977-0bf8559af0c4
# ╠═5a1c12c4-d1d3-4a6e-8832-2fb8cdf598b4
# ╠═3d52db0f-92e6-425d-8618-7b9422d89b43
# ╠═6f74b758-46a7-45cb-a177-8d766a8c42e8
# ╠═411fc6ec-cb7e-4ec9-b2f5-83b3b6347b52
# ╠═4d8c7f3a-582c-4fe9-a692-5966492cc31f
# ╠═8915304d-d902-4988-98c7-943f4f29478b
# ╠═7cea41a4-a94d-4b44-8f32-04b324206d36
# ╠═b85e57b4-645b-452b-bae7-e55ee6cedc0e
# ╠═e7f1e626-b002-4848-aac4-891e5efc2d63
# ╠═a1e4f8b0-f704-484e-87b0-b8c9a17a608f
# ╠═d5991b35-16c2-4969-8b41-9255451ea843
# ╠═775eaf12-eb3d-4dc1-891a-a1a2be066be4
# ╠═e923dc1f-c289-46f7-9d97-b2d4d8e04488
# ╠═70a6aff1-5547-46be-8884-17984c41fc64
# ╠═4c79d09f-5cd0-4788-ae06-a2be74bbc580
# ╠═6d3fbb35-3f06-45c1-aa78-1c80997c1e21
# ╠═3c476391-74bd-45e5-847a-9c239c3e2f0b
# ╠═fb417960-388d-44f8-8a53-39b549ae3ec2
# ╠═ab0f33cd-55bd-48e4-a1f2-ecd5e4720805
# ╠═0415d6bc-f28f-456c-a30f-3c2d48c90a9e
# ╠═e6a9dd2f-48b6-4420-8f30-1921b9013a0a
# ╠═117717f5-efa7-4e8a-aaa3-5cb245dec9c9
# ╠═755a4eb3-50f4-4188-8804-be4dcf7f547b
# ╠═4ee06f0e-4c8d-44e8-b8ce-14c0e75e0ed2
# ╠═6f0ab9f0-27a4-4a84-94b6-ac4f155ef182
# ╠═380e891d-ce28-4b26-8d54-99ff1ad2621f
# ╠═b8fe811d-350c-4806-99ef-481c40cd7445
# ╠═031d0b0e-6780-46ec-8ee7-398ebd4a62e9
# ╠═e727f789-e08e-435d-aede-59867c373e67
# ╠═870948f5-e827-4315-9786-e9b3d092642f
# ╠═20a1fe9a-4aca-4a16-869a-cbda38df6db4
# ╠═d0f3d270-71b6-40a5-b9bc-49a48657cf22
# ╠═ed4fa21c-fce3-4e12-89fd-73a780af0250
# ╠═2b87a94a-d9ac-4a56-83a2-e70abcbd0a98
# ╠═5159c3d6-9bc7-4377-9839-9f0c26ce0276
# ╠═1fae0996-4183-445a-880c-330e52405265
# ╠═23e7b7f7-5f24-428f-8faa-b3778b5ba7e1
# ╠═eebff747-334b-402d-a0ce-107f98207f30
# ╠═ed23953d-214c-49a1-ac91-5b31abd935d9
# ╠═437fb458-4d5c-49a1-b029-3e48e6e92e5c
# ╠═7e4e16be-0972-43a7-bf03-f342bdc975af
# ╠═2cfd8772-1812-4af5-8ac7-fd360ee798f1
# ╠═04970711-8673-41bf-9d55-14aefabb6ed4
# ╠═f50d8baa-071c-4e3c-a130-45b504f1fa42
# ╠═ea7bc95b-5207-4096-9145-e9c339dec50f
# ╠═790a86b7-3e2d-45b1-8a06-64fb100080f2
# ╠═042fa3b0-36c2-48b0-80f7-e0e01edd13a5
# ╠═f08c5dfe-e554-4c2d-a0e0-0826c166fa21
# ╠═2da5e914-1e07-411c-b942-f4bbba14f835
# ╠═ba57147b-a813-478c-93a8-dc429fc24a7c
# ╠═509829b9-8edb-4495-9ea0-c3de242f0436
# ╠═39b179a2-849a-4466-916c-9abffe084f00
# ╠═9e1706f7-919e-4c26-b356-0f2834c3f2a5
# ╠═e38a0de6-34ca-412b-9367-989cdae07788
# ╠═19247928-adb9-45b4-a8e7-4b0edce080d9
# ╠═7b7e7cff-491a-46f1-943a-ee5d2b195bdb
# ╠═299615c8-5c97-4bda-a249-9c55197c3371
# ╠═f04645f1-b1ba-4235-ab91-0fda740a931e
# ╠═55c016a4-d09c-4dfa-b886-c650baf0b645
# ╠═c3f170cf-2776-4827-a447-5d4ace14fe2e
# ╠═7359bb10-a2f0-4010-a5c7-25822327b96d
# ╠═18ca6d06-b94a-453a-9c61-255033998525
# ╠═9b7209c7-384c-414e-b198-13daa3d47d68
# ╠═701bb671-39a6-480c-b715-6304a09058f0
# ╠═b4a3ad9e-bff2-44a2-9275-48dff2be585f
# ╠═f1cdee20-c959-4f7a-b6d0-7a0ff7b15d86
# ╠═1a3c53b3-840a-4130-b40e-856f4a64006f
# ╠═30853cdb-06c7-4599-9ac2-76e78992922b
# ╠═6bbdde5a-4522-45a6-aaec-5657371d6ca3
# ╠═f529973b-5338-4987-aa59-0d342918a4aa
# ╠═ac054b76-0e44-41cd-bb91-afcbb897bf8b
# ╠═1b9949f1-48e9-4232-858a-4712f362b4d3
# ╠═3611cfbc-b987-46f5-9909-4bb345e28c54
# ╠═ae40c445-e1e6-4b5c-b23a-190fab0e0158
# ╠═e0a28475-9195-44cd-ae90-38e7092b87c4
# ╠═b8f09841-fe72-4f89-befd-625758d9a354
# ╠═7d032b2c-114d-4caa-84bd-de7e2ac4ac0d
# ╠═61500d0e-ca1f-49ab-9782-296186808108
# ╠═30fcbde8-70e9-453c-b621-f277209bd609
# ╠═80c260d7-f7c3-45a2-8880-960e3ece5c54
# ╠═e838d374-cdfe-4338-b850-6f0fa11a897d
# ╠═44b8bd75-5893-4e8b-9eca-297fd0f7c137
# ╠═70930065-e396-4418-b6ce-80a30cf85c07
# ╠═fc64bb5e-08e7-4a65-bc03-79e6ebaa616f
# ╠═fc2ea452-c71c-4408-adc9-51e0a633f14e
# ╠═85ed2383-4038-44ac-9e49-c71d005a9b27
# ╠═f7c369a1-32ce-41a1-9170-fa39c18f82a8
# ╠═5be9c1cb-5e40-4e59-8cf7-b013a896f2ae
# ╠═8108835b-1763-4df5-860b-86eebb86260f
# ╠═1e7d2b50-dc43-40fb-9e2e-00e3f3f00763
# ╠═6ba5df8a-caeb-452e-bff0-4604f922fdf3
# ╠═71c150ce-91a7-40ca-9c15-bad5593f1f4a
# ╠═ed6a1cea-c9c7-44a7-8c3d-471fe1605721
# ╠═4942e000-421e-47bd-9829-92914de7d1bd
# ╠═59cabe32-fe06-4589-a403-cdc86501aa73
# ╠═29486f2f-3c72-4226-83ad-890ed9d83f66
# ╠═8d69770e-c243-428b-a634-4da77c544235
# ╠═34ca8b05-3b6b-4636-9386-0f0ebe836e52
# ╠═fdde9769-0822-47bb-890e-55612adda3b3
# ╠═af8f790b-afbd-41e2-a703-40a6aa89e951
# ╠═6991ce74-f615-462d-8313-b5ee52e015a0
# ╠═4bf9f619-17b8-421e-ae59-02d04c70a8c2
# ╠═adb2715d-5cb4-4642-9236-600388ad4e5f
# ╠═fde79b23-ec4b-45c0-b874-9b59d96682b9
# ╠═cb7b96d1-5f7b-4099-a020-a067c580900d
# ╠═535ead23-da4c-4ab1-b981-d59c79685773
# ╠═e57591c8-c66e-4544-9674-cc0d0337c593
# ╠═4771a723-bea4-4f1f-be27-c83fee15b284
# ╠═53493ee2-1719-4a8f-a1ac-66bb0ad0cc49
# ╠═da7f7047-cb19-4a2d-aff2-7c03cea49c31
# ╠═b29455bc-05d7-4d1d-8b52-0ad2836da6a0
# ╠═54183ac7-d928-4250-aeaf-75adf5c0bc42
# ╠═1ebe2956-3706-45b4-91f5-804e541bcc91
# ╠═036f3a55-2fd5-48e1-adfd-1dd2ab3d03eb
# ╠═887ec730-dde8-4fa1-8671-f0815af0c466
# ╠═44ebd883-232c-4a3c-a2a3-f7c8c85ad8e3
# ╠═8d7befb0-6225-404c-8dd2-12834978bde1
# ╠═6dd088fa-5907-400d-bf13-4890a224401e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
