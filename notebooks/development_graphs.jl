### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ 80af05a0-8c11-11ed-1008-cb6a5443493d
begin
	using PlutoUI, Graphs, GraphMakie, Symbolics, NetworkLayout, CairoMakie, GasChromatographySimulator, Interpolations, CSV, DataFrames, Plots
	TableOfContents()
end

# ╔═╡ d8b77441-b471-4efe-9e0f-93f369be63b8
using Roots

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
### [] Module structures
"""

# ╔═╡ 932ffdbb-0dd6-48ca-98b2-af70f229d4f5
abstract type AbstractModule end

# ╔═╡ 68c043ba-ce1e-4693-a9fd-f761d1edec58
# combine Transferline and Column in one module. `temperature` can be a number or a `TemperatureProgram`

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
default_TP = TemperatureProgram([0.0, 1800.0], [40.0, 240.0])

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
	g::Graphs.SimpleDiGraph{Int}
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
						node_size = [node_size for i in 1:Graphs.nv(g)],
						node_color = [color for i in 1:Graphs.nv(g)],
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
						nlabels = [sys.pressurepoints[i].name for i in 1:Graphs.nv(sys.g)],
						nlabels_align=(:center,:center),
						node_size = [node_size for i in 1:Graphs.nv(sys.g)],
						node_color = [color for i in 1:Graphs.nv(sys.g)],
						elabels = [sys.modules[i].name for i in 1:Graphs.ne(sys.g)]
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

# ╔═╡ a206d81b-86b6-41cd-9a5b-a291d361642d
101300+137300

# ╔═╡ 7cea41a4-a94d-4b44-8f32-04b324206d36
begin
	g0 = Graphs.SimpleDiGraph(2)
	Graphs.add_edge!(g0, 1, 2) # Inj -> GC column -> Det
	# pressure points:
	pp0 = Array{PressurePoint}(undef, Graphs.nv(g0))
	pp0[1] = PressurePoint("p₁", [0.0, 1800.0], [174700.0, 238600.0]) # inlet 
	pp0[2] = PressurePoint("p₂", [0.0, 1800.0], [101300.0, 101300.0]) # outlet
	# modules
	modules0 = Array{AbstractModule}(undef, Graphs.ne(g0))
	modules0[1] = ModuleColumn("GC column", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	# system
	sys0 = System(g0, pp0, modules0, Options(mobile_phase="H2", ng=true))
end

# ╔═╡ e7f1e626-b002-4848-aac4-891e5efc2d63
md"""
### [x] Series of columns
"""

# ╔═╡ 775eaf12-eb3d-4dc1-891a-a1a2be066be4
begin
	g_series = SimpleDiGraph(4)
	add_edge!(g_series, 1, 2) # Inj -> TL column -> GC column inlet
	add_edge!(g_series, 2, 3) # GC column inlet -> GC column -> 2nd TL column inlet
	add_edge!(g_series, 3, 4) # 2nd TL column inlet -> TL column -> Det
	# pressure points
	pp_series = Array{PressurePoint}(undef, nv(g_series))
	pp_series[1] = PressurePoint("p₁", [0.0, 600.0, 1200.0], [200000.0, 200000.0, 200000.0]) # inlet 
	pp_series[2] = PressurePoint("p₂", [0.0, 1800.0], [NaN, NaN]) 
	pp_series[3] = PressurePoint("p₃", [0.0, 1800.0], [NaN, NaN]) 
	pp_series[4] = PressurePoint("p₄", [0.0, 1800.0], [101300.0, 101300.0]) # outlet
	# modules
	modules_series = Array{AbstractModule}(undef, ne(g_series))
	modules_series[1] = ModuleColumn("GC column", 10.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	modules_series[2] = ModuleColumn("GC column", 10.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	modules_series[3] = ModuleColumn("GC column", 10.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	# system
	sys_series = System(g_series, pp_series, modules_series, Options(ng=true))
end

# ╔═╡ bc7a259b-9d43-4b7e-9fca-755f391ed890
adjacency_matrix(sys_series.g) 

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
	pp_split[3] = PressurePoint("p₃", [0.0, 1800.0], [eps(Float64), eps(Float64)]) # outlet 1 
	pp_split[4] = PressurePoint("p₄", [0.0, 1800.0], [101300.0, 101300.0]) # outlet 2
	# modules
	modules_split = Array{AbstractModule}(undef, ne(g_split))
	modules_split[1] = ModuleColumn("GC column", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	modules_split[2] = ModuleTransferline("TL column 1", 0.5, 0.15e-3, 0.0e-6, "", 300.0)
	modules_split[3] = ModuleTransferline("TL column 2", 1.0, 0.25e-3, 0.25e-6, "SLB5ms", 300.0)
	# system
	sys_split = System(g_split, pp_split, modules_split, Options(ng=true))
end

# ╔═╡ f89e34e1-b10a-42c5-b526-e69517095159
adjacency_matrix(sys_split.g) 

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
	pp[8] = PressurePoint("p₈", [0.0, 1800.0], [eps(Float64), eps(Float64)])#[0.0, 0.0]) # outlet 2
	pp[9] = PressurePoint("p₉", [0.0, 1800.0], [101300.0, 101300.0]) # outlet 3
	# modules
	modules = Array{AbstractModule}(undef, ne(g))
	modules[1] = ModuleColumn("GC column 1", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	modules[2] = ModuleColumn("GC column 2", 30.0, 0.25e-3, 0.25e-6, "SPB50", default_TP)
	modules[3] = ModuleTransferline("TL column 1", 3.0, 0.25e-3, 0.25e-6, "SPB50", 300.0)
	modules[4] = ModuleTransferline("TL column 2", 1.0, 0.25e-3, 0.25e-6, "SPB50", 300.0)
	modules[5] = ModuleTransferline("Modulator", 0.3, 0.1e-3, 0.1e-6, "Wax", 300.0) # replace later with ModuleModulator
	modules[6] = ModuleColumn("GC column 3", 4.0, 0.1e-3, 0.1e-6, "Wax", default_TP)
	modules[7] = ModuleTransferline("TL column 3", 2.5, 0.1e-3, 0.1e-6, "", 300.0)
	modules[8] = ModuleTransferline("TL column 4", 1.5, 0.1e-3, 0.1e-6, "", 300.0)
	# system
	sys = System(g, pp, modules, Options(ng=true))
end

# ╔═╡ 627e0139-70da-490b-a85d-13024d76730f
eps(Float64)

# ╔═╡ 69dbe34c-2075-4d0e-bbd0-f7b13cce4ed8
adjacency_matrix(sys.g) 

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
	modules_loop[1] = ModuleTransferline("1 TL1", 0.5, 0.25e-3, 0.25e-6, "Wax", 300.0)
	modules_loop[2] = ModuleTransferline("2 TL2", 1.0, 0.25e-3, 0.25e-6, "Wax", 300.0)
	modules_loop[3] = ModuleColumn("5 GC2", 2.0, 0.1e-3, 0.1e-6, "SLB5ms", default_TP) 
	modules_loop[4] = ModuleColumn("3 GC1", 10.0, 0.25e-3, 0.25e-6, "Wax", default_TP)
	modules_loop[5] = ModuleColumn("4 TL3", 1.0, 0.1e-3, 0.1e-6, "Wax", default_TP)
	
	modules_loop[6] = ModuleTransferline("6 TL4", 1.5, 0.1e-3, 0.1e-6, "", 300.0)
	# system
	sys_loop = System(g_loop, pp_loop, modules_loop, Options(ng=true))
end

# ╔═╡ 117717f5-efa7-4e8a-aaa3-5cb245dec9c9
collect(edges(sys_loop.g))

# ╔═╡ 3d792c30-4770-4eef-be13-dfd1ccb04e5c
adjacency_matrix(sys_loop.g) 

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

# ╔═╡ 8517a63b-198d-4658-a015-aab8ac7969ce
begin
	g_test = SimpleDiGraph(6)
	add_edge!(g_test, 1, 2) # Inj->TL1->Split 
	add_edge!(g_test, 1, 3) # Split -> TL2
	add_edge!(g_test, 2, 5) # GC2
	add_edge!(g_test, 3, 4) # GC1
	add_edge!(g_test, 4, 5) # TL3
	add_edge!(g_test, 5, 6)
end

# ╔═╡ 0645c482-7a21-4e2e-b574-6f733b7e1545
Graphs.degree(g_test)

# ╔═╡ 3e1a13da-b7e3-4d40-82a2-2b91f58ed191
neighbors(g_test, 1)

# ╔═╡ ab162833-59f4-4b52-b3e2-0500aefd365a
inneighbors(g_test, 1)

# ╔═╡ 0e42d2f7-7b5a-4acb-bc8c-2ec7a7e890ff
outneighbors(g_test, 1)

# ╔═╡ e6cd895b-3d69-41ff-8733-6cfaa50452f6
inneighbors(g_test, 6)

# ╔═╡ 42af3ef1-8969-4428-94e1-3195e03913fc
outneighbors(g_test, 6)

# ╔═╡ ef85d984-b494-4b3f-8243-512d9eca4864
vertices(g_test)

# ╔═╡ cd9eb028-97be-4e5a-983a-7c0f09c3669b
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

# ╔═╡ d83b9980-7d28-414c-910d-8314d613baf3
inlet_vertices(g_test)

# ╔═╡ 540bb955-da9d-4628-8eaf-f86f918da1ed
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

# ╔═╡ b141c05e-cccc-43d1-8bd5-495b2b14215e
outlet_vertices(g_test)

# ╔═╡ 86572e5d-179f-4257-a96e-43fc075ccbab
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

# ╔═╡ b0a7d30b-4e2c-49c8-acf6-428836c5d7c1
inner_vertices(g_test)

# ╔═╡ ea162e25-2803-435a-aeb6-137f4365b3f7
findall(Graphs.degree(g_test).>1)

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
function flow_restriction(t, module_, options)
	if typeof(module_) == ModuleColumn
		T_itp = GasChromatographySimulator.temperature_interpolation(module_.temperature_program.timesteps, module_.temperature_program.temperaturesteps, module_.temperature_program.gradient_function, module_.length)
	elseif typeof(module_) == ModuleTransferline
		gf(x) = [zero(x), zero(x)]
		T_itp = GasChromatographySimulator.temperature_interpolation([0.0, 3600.0], [module_.temperature, module_.temperature], gf, module_.length)
	end
	κ = GasChromatographySimulator.flow_restriction(module_.length, t, T_itp, module_.diameter, options.mobile_phase; ng=options.ng, vis=options.vis)
	return κ
end

# ╔═╡ c3f170cf-2776-4827-a447-5d4ace14fe2e
function flow_restrictions(sys)
	kappas = Array{Function}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		f(t) = flow_restriction(t, sys.modules[i], sys.options)
		kappas[i] = f
	end
	return kappas
end

# ╔═╡ 7359bb10-a2f0-4010-a5c7-25822327b96d
κ_split = flow_restrictions(nsys_split)

# ╔═╡ 2d4366c7-bd2c-40fb-92c8-a03f32bf5fd0
κ_0 = flow_restrictions(nsys0)

# ╔═╡ ed5d3089-bb8a-427a-a3ee-75986e6cc583
κ_series = flow_restrictions(nsys_series)

# ╔═╡ f9a69795-bee2-46f5-80e9-9a262fa5b6b1
begin
	plotly()
	p_κ = Plots.plot(xlabel="t in s", ylabel="κ", legend=:topleft)
	Plots.scatter!(0.0:100:1800.0, κ_0[1].(0.0:100:1800.0))
	for i=1:3
		Plots.scatter!(0.0:100:1800, κ_series[i].(0.0:100:1800))
	end
	Plots.scatter!(0.0:100:1800, κ_series[1].(0.0:100:1800).+κ_series[2].(0.0:100:1800).+κ_series[3].(0.0:100:1800))
	p_κ
end

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

# ╔═╡ de33be5f-458b-4fdf-8d37-4a733ff41320
plot_graph(sys_loop)

# ╔═╡ 47413580-d2c0-4ec1-9b85-7710ef24d4e2
plot_graph(g_test, ["V₁ => V₂", "V1 => V₃", "V₂ => V5", "V3 => V4", "V4 => V5", "V5 => V6"], ["V₁", "V₂", "V₃", "V₄", "V5", "V6"])

# ╔═╡ fc64bb5e-08e7-4a65-bc03-79e6ebaa616f
plot_graph(nsys)

# ╔═╡ 85ed2383-4038-44ac-9e49-c71d005a9b27
plot_graph(nsys, 0)

# ╔═╡ f7c369a1-32ce-41a1-9170-fa39c18f82a8
plot_graph(nsys, 1200)

# ╔═╡ 5be9c1cb-5e40-4e59-8cf7-b013a896f2ae
plot_graph(nsys0, 0)

# ╔═╡ 8108835b-1763-4df5-860b-86eebb86260f
plot_graph(nsys_series, 1200)

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
## [x] Transform system to GasChromatographySimulator structures
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
selected_solutes = ["Octane", "Decane"]#, "Dodecane", "Tetradecane", "2-Octanone", "2-Octanol", "5-Nonanol", "2-Nonanol"]

# ╔═╡ da7f7047-cb19-4a2d-aff2-7c03cea49c31
GasChromatographySimulator.load_solute_database(db_dataframe, "Rxi5ms", nsys.options.mobile_phase, selected_solutes, NaN.*ones(length(selected_solutes)), NaN.*ones(length(selected_solutes)))

# ╔═╡ b29455bc-05d7-4d1d-8b52-0ad2836da6a0
md"""
### [x] Modules to parameters
"""

# ╔═╡ 887ec730-dde8-4fa1-8671-f0815af0c466
function graph_to_parameters(sys, db_dataframe, selected_solutes)
	E = collect(edges(sys.g))
	srcE = src.(E)
	dstE = dst.(E)
	p_func = pressure_functions(sys)
	parameters = Array{GasChromatographySimulator.Parameters}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		col = GasChromatographySimulator.Column(sys.modules[i].length, sys.modules[i].diameter, [sys.modules[i].diameter], sys.modules[i].film_thickness, [sys.modules[i].film_thickness], sys.modules[i].stationary_phase, sys.options.mobile_phase)

		pin_steps = sys.pressurepoints[srcE[i]].pressure_steps
		pout_steps = sys.pressurepoints[dstE[i]].pressure_steps
		pin_itp = p_func[srcE[i]]
		pout_itp = p_func[dstE[i]]	
		if typeof(sys.modules[i]) == ModuleColumn
			time_steps = sys.modules[i].temperature_program.timesteps
			temp_steps = sys.modules[i].temperature_program.temperaturesteps	
			gf = sys.modules[i].temperature_program.gradient_function
			a_gf = sys.modules[i].temperature_program.a_gradient_function
			T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, sys.modules[i].length)		
		elseif typeof(sys.modules[i]) == ModuleTransferline
			time_steps = common_timesteps(sys)
			temp_steps = sys.modules[i].temperature.*ones(length(time_steps))
			gf(x) = zero(x).*ones(length(time_steps))
			a_gf = [zeros(length(time_steps)) zeros(length(time_steps)) ones(length(time_steps)) zeros(length(time_steps))]
			T_itp = GasChromatographySimulator.temperature_interpolation([0.0, 3600.0], [sys.modules[i].temperature, sys.modules[i].temperature], gf, sys.modules[i].length)
		end
		prog = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)

		sub = GasChromatographySimulator.load_solute_database(db_dataframe, sys.modules[i].stationary_phase, sys.options.mobile_phase, selected_solutes, NaN.*ones(length(selected_solutes)), NaN.*ones(length(selected_solutes)))

		opt = GasChromatographySimulator.Options(sys.options.alg, sys.options.abstol, sys.options.reltol, sys.options.Tcontrol, sys.options.odesys, sys.options.ng, sys.options.vis, sys.options.control, sys.options.k_th)

		parameters[i] = GasChromatographySimulator.Parameters(col, prog, sub, opt)
	end
	return parameters
end

# ╔═╡ 44ebd883-232c-4a3c-a2a3-f7c8c85ad8e3
graph_to_parameters(nsys, db_dataframe, selected_solutes)

# ╔═╡ 8d7befb0-6225-404c-8dd2-12834978bde1
md"""
## [] Simulate the system
"""

# ╔═╡ 9eedabe6-811c-41a0-a90e-f2f62fddb170
md"""
1. determine number of paths from injector to detector/detectors
- in graph without loop/cycle the number of paths from injector to detectors is the same as the number of detectors
2. determine the order of the simulations according to flow directions
3. starting from inlet
"""

# ╔═╡ 01d9e464-3d6c-4ffe-a92d-dee1d88258bb
md"""
### [x] Paths for graph
"""

# ╔═╡ e3fcbb11-cce6-4d5f-aa7c-13d1d11174b2
function inlets_outlets(g)
# outlets = nodes, which are only destinations of edges but no source of an edge
# inlets = nodes, which are only source of edges but no destination of an edge
	E = collect(edges(g))
	outer_V = findall(degree(g).==1) # outer nodes
	inlets_V = Int[]
	outlets_V = Int[]
	for i=1:nv(g)
		if i in src.(E) && !(i in dst.(E))
			push!(inlets_V, i)
		elseif !(i in src.(E)) && i in dst.(E)
			push!(outlets_V, i)
		end
	end
	inlets_V, outlets_V
end

# ╔═╡ 93cc1041-8883-4c8e-88e5-a26389a28ee2
function split_merge(g)
	# split has multiple edges comming in
	# merge has multiple edges going out
	E = collect(edges(g))
	split_V = []
	merge_V = []
	sum_out = 0
	sum_in = 0
	split_merge_V = findall(degree(g).>2)
	for i=1:length(split_merge_V)
		src_ = findall(src.(E).==split_merge_V[i])
		if length(src_)>1
			push!(split_V, split_merge_V[i])
		end
		dst_ = findall(dst.(E).==split_merge_V[i])
		if length(dst_)>1
			push!(merge_V, split_merge_V[i])
		end
		
		sum_out = sum_out + length(src_)
		sum_in = sum_in + length(dst_)
	end
	return split_V, merge_V, sum_in, sum_out
end

# ╔═╡ 70c82cf4-012e-4137-a0fd-7060cda73fb5
function number_of_paths(g)# wrong
	split_V = findall(degree(g).>2) # split nodes
	inlets_V, outlets_V = inlets_outlets(g)
	num_paths = length(outlets_V) + trunc(Int, length(split_V)/2)
	return num_paths
end

# ╔═╡ 30cdb150-fad1-43fa-8ef0-a776fc34b4d8
function paths(g) # works only for graphs without `cycles` (split of path and later merging of the paths)
	inlets_V, outlets_V = inlets_outlets(g)
	paths = Array{Array{Graphs.SimpleGraphs.SimpleEdge{Int64},1}}(undef, length(outlets_V))
	for i=1:length(outlets_V) # only one inlet at node `1` assumed
		paths[i] = a_star(g, 1, outlets_V[i])
	end
	return paths
end

# ╔═╡ 7290bbc2-76d4-4b4f-9f46-e369bff7afae
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

# ╔═╡ bce24317-737e-4b43-8ad0-4b82ffe9b973
function all_paths(g, num_paths)
	rand_paths = Any[]
	while length(unique(rand_paths))<num_paths
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

# ╔═╡ d8ee7f45-67a0-4b11-a401-4d2b1cf31e10
paths(nsys.g)

# ╔═╡ dbf731b1-c897-4be0-a6fe-3358f0782cfd
all_paths(nsys.g)

# ╔═╡ 0c3eda88-5d13-43b2-b083-7577af04cf84
all_paths(nsys0.g)

# ╔═╡ 203bf2a5-64eb-4f6d-a6f6-1bf12b19a968
all_paths(nsys_series.g)

# ╔═╡ 318075b8-76f3-4e9e-8655-5424cb42fe47
all_paths(nsys_split.g)

# ╔═╡ afc385df-0184-4c97-8201-12a6ad714fd2
all_paths(nsys_loop.g)

# ╔═╡ 047bb299-8391-4da8-9c0e-2c4bca80371d
md"""
### [] Simulate along the path
"""

# ╔═╡ 7b13e4c3-be31-446e-8edd-b3e1c75cde4d
paths_V, paths_E = all_paths(nsys.g)

# ╔═╡ 97b19b85-841d-41a5-9a47-eaf497adfee2
function index_parameter(g, path)
	i_par = Array{Int}(undef, length(path))
	for i=1:length(path)
		i_par[i] = findall(src.(path)[i].==src.(edges(g)))[findfirst(findall(src.(path)[i].==src.(edges(g))).==findall(dst.(path)[i].==dst.(edges(g))))]
	end
	return i_par
end

# ╔═╡ 9b16093f-2db4-485e-9f6e-353dc8db886d
index_parameter(nsys.g, paths_E[1])

# ╔═╡ aeb33e75-f64d-4b36-8c49-ec73e5dd4011
index_parameter(nsys.g, paths_E[2])

# ╔═╡ 47d4cdd7-b009-4f3d-834b-6b1809a5c5d2
index_parameter(nsys.g, paths_E[3])

# ╔═╡ bdc69707-46d1-47ea-ba8d-6855c3b68251
paths_E

# ╔═╡ 25696111-910c-42fb-942d-b0fb382d1334
# common edges between two paths

# ╔═╡ 02ccb0b0-2494-4733-b1e5-55bf863064a5
function common_edges(path1, path2)
	return path2[findall(x->x in path2, path1)]
end

# ╔═╡ 6e173c0a-fe02-455b-ab66-74e5cb6f0164
common_edges(paths_E[1], paths_E[2])

# ╔═╡ dc1a71ad-3bbe-48b5-926f-514bb8bdf7e3
common_edges(paths_E[2], paths_E[1])

# ╔═╡ 8095da90-dc1d-46db-b383-9980e0485733
common_edges(paths_E[1], paths_E[3])

# ╔═╡ 4d3c608c-2d60-454f-8c4f-2c01fb780cf0
common_edges(paths_E[3], paths_E[2])

# ╔═╡ 033c19b4-f696-44ac-a23e-9a6f1e8368af
common_edges(paths_E[2], paths_E[3])

# ╔═╡ a91f5cd2-adb4-488f-a936-e2f8c48d8f83
# simulations on common edges, connected to the inlet, can be reused for the other paths

# ╔═╡ e652a548-4c03-471a-b6a9-717bead6b355
# check if the flow is always positive over the whole runtime in the used edges

# ╔═╡ 92f7549c-ee98-4074-934c-b6b69bde8bd9
function positive_flow(sys, path)
	F_func = flow_functions(sys)
	tend = sum(common_timesteps(sys))
	t = 0:tend/100:tend
	ipar = index_parameter(sys.g, path)
	pos_Flow = Array{Bool}(undef, length(ipar))
	for i=1:length(ipar)
		if isempty(findall(F_func[ipar[i]].(t).<=0))
			pos_Flow[i] = true
		else
			pos_Flow[i] = false
		end
	end
	return pos_Flow
end

# ╔═╡ 284553a9-5e63-4f47-ad9b-7a66c75c12b7
positive_flow(nsys, paths_E[1])

# ╔═╡ 8b9f4f16-0bbf-4319-bada-15f0275ad78d
positive_flow(nsys, paths_E[2])

# ╔═╡ f4bac904-83ee-42a1-8169-cd335745e1f7
positive_flow(nsys, paths_E[3])

# ╔═╡ 07eb454e-f035-4419-8689-626d120bf012
find_zero(F_func[1], 1000.0)

# ╔═╡ a067f417-99e0-4844-bec6-dbbf9166205e
function zero_flow(sys, path)
	F_func = flow_functions(sys)
	tend = sum(common_timesteps(sys))

	ipar = index_parameter(sys.g, path)
	zero_flow = Array{Float64}(undef, length(ipar))
	for i=1:length(ipar)
		zero_flow[i] = find_zero(F_func[ipar[i]], tend/2)
	end
	return zero_flow
end

# ╔═╡ f63a4d32-1877-466c-9df9-24da42365c08
zero_flow(nsys, paths_E[1])

# ╔═╡ 84fe5a51-3011-4c91-b070-c71f2a64656e
zero_flow(nsys, paths_E[2])

# ╔═╡ cb46b45f-9a50-480e-99ea-0bbed4a4ea81
zero_flow(nsys, paths_E[3])

# ╔═╡ 4605cc1e-6885-4c3c-ad57-ae5643324c5f
function zero_flow_(sys, F_func)
	#F_func = flow_functions(sys)
	#tend = sum(common_timesteps(sys))

	zero_flow = Array{Float64}(undef, length(F_func))
	for i=1:length(F_func)
		zero_flow[i] = find_zero(F_func[i], 0.0)
	end
	return zero_flow
end

# ╔═╡ 6d66210f-e504-48a3-8832-190fd7fb4cbb
zero_flow_(nsys, F_func)

# ╔═╡ 6a84d685-0e8d-41fc-b116-57b95de99814
function positive_flow_(sys)
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

# ╔═╡ 7bf076d8-3ae2-4c0a-ae82-db1164d7af4e
positive_flow_(nsys)

# ╔═╡ be21c243-2aa9-4c9f-8fdf-3656717082cc
findall(positive_flow_(nsys))

# ╔═╡ 60f69e16-60e3-4f17-ad79-74af60d74491
function path_possible(sys, paths)
	F_func = flow_functions(sys)
	i_E = index_parameter(sys.g, paths)
	if length(i_E) == length(findall(positive_flow_(nsys)[i_E]))
		possible = true
	else
		possible = false
	end
	return possible
end

# ╔═╡ 5c738d5d-12ea-4ca9-94ee-da4c56d73318
path_possible(nsys, paths_E[1])

# ╔═╡ 2ddab7b9-2d0c-4d07-9a5b-0314bedbe9de
path_possible(nsys, paths_E[2])

# ╔═╡ 5e1f9b33-6bd3-4acb-9132-09551e13f8f8
path_possible(nsys, paths_E[3])

# ╔═╡ 1e055bcd-0713-4992-83c5-f9621376f891
function change_initial(par::GasChromatographySimulator.Parameters, init_t, init_τ)
	# copys the parameters `par` and changes the values of par.sub[i].τ₀ and par.sub[i].t₀ to init_τ[i] resp. init_t[i]
	newsub = Array{GasChromatographySimulator.Substance}(undef, length(par.sub))
	for i=1:length(par.sub)
		newsub[i] = GasChromatographySimulator.Substance(par.sub[i].name, par.sub[i].CAS, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀, par.sub[i].ann, par.sub[i].Cag, init_t[i], init_τ[i])
	end
	newpar = GasChromatographySimulator.Parameters(par.col, par.prog, newsub, par.opt)
	return newpar
end

# ╔═╡ 5e9c922b-5eca-4c5b-872d-9fb73dc20389
function simulate_along_paths(sys, paths, db_dataframe, selected_solutes; t₀=zeros(length(selected_solutes)), τ₀=zeros(length(selected_solutes)))
	par_sys = graph_to_parameters(sys, db_dataframe, selected_solutes)
	E = collect(edges(sys.g))
	peaklists = Array{Array{DataFrame,1}}(undef, length(paths))
	solutions = Array{Array{Any,1}}(undef, length(paths))
	path_pos = Array{String}(undef, length(paths))
	new_par_sys = Array{GasChromatographySimulator.Parameters}(undef, length(par_sys))
	for i=1:length(paths)
		i_par = index_parameter(sys.g, paths[i])
		if path_possible(sys, paths[i]) == true
			path_pos[i] = "path is possible"
			peaklists_ = Array{DataFrame}(undef, length(i_par))
			solutions_ = Array{Any}(undef, length(i_par))
			for j=1:length(i_par)
				if (isassigned(new_par_sys, i_par[j]) == true) && (i > 1) # new parameters allready assigned, than the simulation for this edge was already done, assing the results ... 
					# find the edge E[i_par[j]] in the previous paths and use the results from there

					#!!!! BUT ONLY RESULTS DIRECTLY CONNECTED WITH ALLREADY SIMULATED MODULES CAN BE USED, see sys_loop 
					i_prev_path = Int[]
					for k=1:i-1
						if in(i_par[j], index_parameter(sys.g, paths[k]))
							push!(i_prev_path, k)
						end
					end
					println("i=$(i), j=$(j), i_prev_path=$(i_prev_path)")
					i_par_prev = index_parameter(sys.g, paths[i_prev_path[end]])
					peaklists_[i_par[j]] = peaklists[i_prev_path[end]][findfirst(i_par[j].==i_par_prev)]
					solutions_[i_par[j]] = solutions[i_prev_path[end]][findfirst(i_par[j].==i_par_prev)]
				elseif j == 1 #src(E[i_par[j]]) == 1 # edge starts at node 1, which is assumed to be the injector
				##if j == 1
					new_par_sys[i_par[j]] = change_initial(par_sys[i_par[j]], t₀, τ₀)
					peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(new_par_sys[i_par[j]])
				else
					new_par_sys[i_par[j]] = change_initial(par_sys[i_par[j]], peaklists_[j-1].tR, peaklists_[j-1].τR)
					peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(new_par_sys[i_par[j]])
				end
			end
			peaklists[i] = peaklists_
			solutions[i] = solutions_
		else
			neg_flow_modules = sys.modules[findall(paths[i][findall(positive_flow_(sys)[i_par].==false)].==edges(sys.g))]
			str_neg_flow = "path not possible: "
			for j=1:length(neg_flow_modules)
				str_neg_flow = str_neg_flow*"Flow in module >$(neg_flow_modules[j].name)< becomes negative during the program. "
			end
			path_pos[i] = str_neg_flow
		end
	end
	return path_pos, peaklists, solutions, new_par_sys
end

# ╔═╡ a1a2ace8-fbd8-4674-b26a-c29bfe2678f1
function simulate_along_all_paths(sys, db_dataframe, selected_solutes; t₀=zeros(length(selected_solutes)), τ₀=zeros(length(selected_solutes)))
	nsys = update_system(sys)
	allpaths = all_paths(nsys.g)[2]
	path_pos, peaklists, solutions, new_par_sys = simulate_along_paths(nsys, allpaths, db_dataframe, selected_solutes; t₀=t₀, τ₀=τ₀)
	return allpaths, path_pos, peaklists, solutions, new_par_sys
end

# ╔═╡ 61a93ebd-0680-4b5e-9cee-102915f7ed72
function simulate_along_paths_brute_force(sys, paths, db_dataframe, selected_solutes; t₀=zeros(length(selected_solutes)), τ₀=zeros(length(selected_solutes)))
	par_sys = graph_to_parameters(sys, db_dataframe, selected_solutes)
	E = collect(edges(sys.g))
	peaklists = Array{Array{DataFrame,1}}(undef, length(paths))
	solutions = Array{Array{Any,1}}(undef, length(paths))
	path_pos = Array{String}(undef, length(paths))
	new_par_sys = Array{GasChromatographySimulator.Parameters}(undef, length(par_sys))
	for i=1:length(paths)
		i_par = index_parameter(sys.g, paths[i])
		if path_possible(nsys, paths[i]) == true
			path_pos[i] = "path is possible"
			peaklists_ = Array{DataFrame}(undef, length(i_par))
			solutions_ = Array{Any}(undef, length(i_par))
			for j=1:length(i_par)
				if j == 1
					new_par_sys[i_par[j]] = change_initial(par_sys[i_par[j]], t₀, τ₀)
					peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(new_par_sys[i_par[j]])
				else
					new_par_sys[i_par[j]] = change_initial(par_sys[i_par[j]], peaklists_[j-1].tR, peaklists_[j-1].τR)
					peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(new_par_sys[i_par[j]])
				end
			end
			peaklists[i] = peaklists_
			solutions[i] = solutions_
		else
			neg_flow_modules = sys.modules[findall(paths[i][findall(positive_flow_(sys)[i_par].==false)].==edges(sys.g))]
			str_neg_flow = "path not possible: "
			for j=1:length(neg_flow_modules)
				str_neg_flow = str_neg_flow*"Flow in module >$(neg_flow_modules[j].name)< becomes negative during the program. "
			end
			path_pos[i] = str_neg_flow
		end
	end
	return path_pos, peaklists, solutions, new_par_sys
end

# ╔═╡ f212b6d9-909a-480d-aacb-113a00a7fa99
vertices(sys.g)

# ╔═╡ 1cedb103-11af-42a3-a933-e0ee0652fe3e
collect(edges(sys.g))

# ╔═╡ b02713df-0bf3-4ddf-bdfc-15d2890a912c
function simulate_along_all_paths_brute_force(sys, db_dataframe, selected_solutes; t₀=zeros(length(selected_solutes)), τ₀=zeros(length(selected_solutes)))
	nsys = update_system(sys)
	allpaths = all_paths(nsys.g)[2]
	path_pos, peaklists, solutions, new_par_sys = simulate_along_paths_brute_force(nsys, allpaths, db_dataframe, selected_solutes; t₀=t₀, τ₀=τ₀)
	return allpaths, path_pos, peaklists, solutions, new_par_sys
end

# ╔═╡ e11e520b-c226-4077-8526-7793aa9a76f6
sim_nsys = simulate_along_paths(nsys, paths_E, db_dataframe, selected_solutes)

# ╔═╡ e6b8a2dd-41e8-4d66-97c3-ccddae0e109e
sim_nsys0 = simulate_along_paths(nsys0, all_paths(nsys0.g)[2], db_dataframe, selected_solutes)

# ╔═╡ 4f1080d6-8a11-4280-802b-0f5dc80fba08
sim_sys_series = simulate_along_all_paths(sys_series, db_dataframe, selected_solutes)

# ╔═╡ bf5b09d0-384c-4222-9f5c-ab2c869a564c
path_possible(nsys_split, all_paths(nsys_split.g)[2][1])

# ╔═╡ 05d97331-ef89-4c4b-a26d-ae61c50bcf03
path_possible(nsys_split, all_paths(nsys_split.g)[2][2])

# ╔═╡ dde7e5a0-9073-48a7-a28a-1575072067c2
index_parameter(nsys_split.g, all_paths(nsys_split.g)[2][1])

# ╔═╡ 4c39462a-1495-494c-b2c4-8a5d35c78f0c
index_parameter(nsys_split.g, all_paths(nsys_split.g)[2][2])

# ╔═╡ d659b94f-27d0-4b2d-a07a-7f77e0dde77e
sim_sys_split = simulate_along_paths_brute_force(nsys_split, all_paths(nsys_split.g)[2], db_dataframe, selected_solutes)

# ╔═╡ 36b6623f-5c1e-4d33-b818-5f52f6e43eb4
#=begin
	t₀=zeros(length(selected_solutes))
	τ₀=zeros(length(selected_solutes))
	paths_split = all_paths(nsys_split.g)[2]
	par_sys_split = graph_to_parameters(nsys_split, db_dataframe, selected_solutes)
	E_split = collect(edges(nsys_split.g))
	peaklists_split = Array{Array{DataFrame,1}}(undef, length(paths_split))
	solutions_split = Array{Array{Any,1}}(undef, length(paths_split))
	path_pos_split = Array{String}(undef, length(paths_split))
	new_par_sys_split = Array{GasChromatographySimulator.Parameters}(undef, length(par_sys_split))
	for i=1:length(paths_split)
		i_par = index_parameter(nsys_split.g, paths_split[i])
		if path_possible(nsys_split, paths_split[i]) == true
			path_pos_split[i] = "path is possible"
			peaklists_split_ = Array{DataFrame}(undef, length(i_par))
			solutions_split_ = Array{Any}(undef, length(i_par))
			for j=1:length(i_par)
				if j == 1
					new_par_sys_split[i_par[j]] = change_initial(par_sys_split[i_par[j]], t₀, τ₀)
					peaklists_split_[j], solutions_split_[j] = GasChromatographySimulator.simulate(new_par_sys_split[i_par[j]])
				else
					new_par_sys_split[i_par[j]] = change_initial(par_sys_split[i_par[j]], peaklists_split_[j-1].tR, peaklists_split_[j-1].τR)
					peaklists_split_[j], solutions_split_[j] = GasChromatographySimulator.simulate(new_par_sys_split[i_par[j]])
				end
			end
			peaklists_split[i] = peaklists_split_
			solutions_split[i] = solutions_split_
		else
			neg_flow_modules = nsys_split.modules[findall(paths_split[i][findall(positive_flow_(nsys_split)[i_par].==false)].==edges(nsys_split.g))]
			str_neg_flow = "path not possible: "
			for j=1:length(neg_flow_modules)
				str_neg_flow = str_neg_flow*"Flow in module >$(neg_flow_modules[j].name)< becomes negative during the program. "
			end
			path_pos_split[i] = str_neg_flow
		end
	end
	path_pos_split, peaklists_split, solutions_split, new_par_sys_split
end=#

# ╔═╡ bf9779bd-b365-48c5-8dcc-44139397898a
t₀=zeros(length(selected_solutes))

# ╔═╡ 0f65afa8-c9df-4c93-8882-337cfc4fab5b
τ₀=zeros(length(selected_solutes))

# ╔═╡ 4e495ff0-1ee1-45e6-abb3-8e6bd09424fa
paths_split = all_paths(nsys_split.g)[2]

# ╔═╡ 3565aac6-cdaa-4788-a73c-59a15b54711d
par_sys_split = graph_to_parameters(nsys_split, db_dataframe, selected_solutes)

# ╔═╡ 72928529-a1aa-4b6a-b040-1c67f37f435a
E_split = collect(edges(nsys_split.g))

# ╔═╡ 2d402085-68b5-485a-b5ac-9d82369e1311
peaklists_split = Array{Array{DataFrame,1}}(undef, length(paths_split))

# ╔═╡ 5dfa6a11-c56a-4fc4-9c16-75cf85b1fb8a
solutions_split = Array{Array{Any,1}}(undef, length(paths_split))

# ╔═╡ bbab45a5-3b66-4b5d-9c52-7e0a4f1e0266
path_pos_split = Array{String}(undef, length(paths_split))

# ╔═╡ da476db9-4603-4a51-b7ea-8f9bafd62d50
new_par_sys_split = Array{GasChromatographySimulator.Parameters}(undef, length(par_sys_split))

# ╔═╡ 50a2fc50-34ed-4ba4-a014-9a619bf7c071
length(paths_split)

# ╔═╡ 4a07f579-b9fc-47fc-ad18-11a9b12d9d0d
i = 1

# ╔═╡ 7c768889-9473-486c-a6cc-f9b76c62a53b
i_par = index_parameter(nsys_split.g, paths_split[i])

# ╔═╡ 836ff6dc-4166-4096-8e3f-cc97856908ee
path_possible(nsys_split, paths_split[i]) == true

# ╔═╡ 00dde33e-fccb-411a-ba75-4dda4be9b035
path_pos_split[i] = "path is possible"

# ╔═╡ 0575d3fb-101b-46d7-819d-c79a2125520a
peaklists_split_ = Array{DataFrame}(undef, length(i_par))

# ╔═╡ b23ed445-e112-415a-8c3d-952cd9aa69ae
solutions_split_ = Array{Any}(undef, length(i_par))

# ╔═╡ c67a3d24-164d-4205-8306-682f98bc6cbe
j = 1

# ╔═╡ 0ca097e1-f821-48b2-bc07-a87944f16bef
j == 1

# ╔═╡ 80337eea-b187-4a07-bc0a-c7e7b43e1b47
i_par[j]

# ╔═╡ 76f4e987-163e-47af-8689-436b98c69348
new_par_sys_split[i_par[j]] = change_initial(par_sys_split[i_par[j]], t₀, τ₀)

# ╔═╡ 31d3fe45-bde0-4038-8157-c677997f2ae3
peaklists_split_[j], solutions_split_[j] = GasChromatographySimulator.simulate(new_par_sys_split[i_par[j]])

# ╔═╡ ee039325-c34c-439c-8d1e-8b8938a6a670
#j=2

# ╔═╡ 7660bc48-5e8c-40b4-a6dd-1a41c93923f0
j2 = 2

# ╔═╡ 86b145f7-712e-4147-a49a-23d59a11246d
j2 == 1

# ╔═╡ 4f87c3bc-fdcc-4211-9818-99ed006d5ca0
new_par_sys_split[i_par[j2]] = change_initial(par_sys_split[i_par[j2]], peaklists_split_[j2-1].tR, peaklists_split_[j2-1].τR)

# ╔═╡ f473d963-30c6-4786-aad8-7df9b7d48511
peaklists_split_[j2], solutions_split_[j2] = GasChromatographySimulator.simulate(new_par_sys_split[i_par[j2]])

# ╔═╡ b275a4e9-fd84-4ba3-96dd-8c298ea8a711
peaklists_split[i] = peaklists_split_

# ╔═╡ 233844fc-4beb-44db-a4a8-5eec12c28dd2
solutions_split[i] = solutions_split_

# ╔═╡ 3255a840-05ac-4490-809f-736a5e27353b
#i=2

# ╔═╡ 8743715a-d966-40a6-ba92-a794b4b57c1c
i2 = 2

# ╔═╡ 56d864dc-fc52-430b-8478-10aa32f95561
i2_par = index_parameter(nsys_split.g, paths_split[i2])

# ╔═╡ 0e9c7b54-34d0-4784-8f67-b9f907ab0e93
path_possible(nsys_split, paths_split[i2]) == true

# ╔═╡ 70829b28-8efc-47dd-b12f-49f9cb29bf5a
path_pos_split[i2] = "path is possible"

# ╔═╡ 3213765f-4df8-4c1d-8189-ed29c8a5e667
peaklists_split_2 = Array{DataFrame}(undef, length(i2_par))

# ╔═╡ e1ac1411-6bf2-45d4-9e33-76a62e1d6ec3
solutions_split_2 = Array{Any}(undef, length(i2_par))

# ╔═╡ dd284c83-7825-44cb-9e01-928aa49897cb
length(i2_par)

# ╔═╡ e5bbe44b-0f0d-4a66-a15f-0303a3f57978
_j = 1

# ╔═╡ 6babaa9e-103a-4525-afe9-ea9fd4f7b347
_j == 1

# ╔═╡ 79d90de0-79d4-4b88-bdc2-3ae2f172e695
new_par_sys_split[i2_par[_j]] = change_initial(par_sys_split[i2_par[_j]], t₀, τ₀)

# ╔═╡ 6ee58b3c-ebfd-4b9a-ae84-8b88916913d3
peaklists_split_2[_j], solutions_split_2[_j] = GasChromatographySimulator.simulate(new_par_sys_split[i2_par[_j]])

# ╔═╡ 96447c64-a8e6-4b59-af54-785b5bfe1803
#_j=2

# ╔═╡ d8c2a71a-9854-4cee-a8d0-7654db698525
_j2 = 2

# ╔═╡ 7d238885-276c-4d3b-90e4-261b09534cde
_j2 == 1

# ╔═╡ a0978e46-704e-43b7-8519-2e8ff698aae0
new_par_sys_split[i2_par[_j2]] = change_initial(par_sys_split[i2_par[_j2]], peaklists_split_2[_j2-1].tR, peaklists_split_2[_j2-1].τR)

# ╔═╡ 45656162-bcdd-4395-af01-7b6f69eed043
peaklists_split_2[_j2], solutions_split_2[_j2] = GasChromatographySimulator.simulate(new_par_sys_split[i2_par[_j2]])

# ╔═╡ 1e1911d1-f783-4f14-9202-a65762cfd8bd
####

# ╔═╡ 6dbccda5-8459-46b6-9966-431447d34b93
graph_to_parameters(nsys_split, db_dataframe, selected_solutes)

# ╔═╡ 669edc63-042b-4c55-99cc-0f19263de4e1
all_paths(nsys_loop.g)[2]

# ╔═╡ 6c482076-52b9-46a7-9b5f-fc0e80a59b2e
sim_sys_loop = simulate_along_all_paths_brute_force(sys_loop, db_dataframe, selected_solutes)

# ╔═╡ 8d63f820-8966-47f3-a60a-c784f51add52
md"""
### Manually simulate path1 of nsys

Simulate along the several paths, reuse of allready simulated modules.
"""

# ╔═╡ 22b1221d-050d-4ff4-9a55-24c5ae6eaec8
i_E = index_parameter(nsys.g, paths_E[3])

# ╔═╡ f765b518-ca59-42a5-b500-e3bf3fbece76
nsys.modules[findall(paths_E[3][findall(positive_flow_(nsys)[i_E].==false)].==edges(nsys.g))]

# ╔═╡ 1cb665cf-694c-49b7-8115-bb24bc98af27
nsys

# ╔═╡ 26da194d-0c9c-40f5-a5aa-d5054f1458dc
new_par_sys = Array{GasChromatographySimulator.Parameters}(undef, 3)

# ╔═╡ 479101b8-83b0-421f-925c-dd95a5df1286
isassigned(new_par_sys, 4)

# ╔═╡ 0d6e7d0c-cee2-4fdb-8782-e77bd6cd8445
path1 = paths_E[1]

# ╔═╡ fa311bb6-3edb-434d-bca3-6efd9138231c
path_possible(nsys, path1)

# ╔═╡ 5bc7a5d1-7106-484d-9078-f1d6acdd8c83
par_nsys = graph_to_parameters(nsys, db_dataframe, selected_solutes)

# ╔═╡ 7729714d-f0be-4be4-b9ed-39c17a6c5ffb
new_par_sys[2] = par_nsys[1]

# ╔═╡ 4dc19ebe-e5c8-4301-8f93-f5cdb155927b
ipar1 = index_parameter(nsys.g, path1)

# ╔═╡ b54430da-c90c-40dc-9eca-d6cd1d10d0a6
par_nsys[1].sub

# ╔═╡ 5cc2f15a-d26d-4c79-889a-d398ead57f7e
par_path1 = Array{GasChromatographySimulator.Parameters}(undef, length(ipar1))

# ╔═╡ 8e172e01-4f9f-4f23-8893-e9c8179558a9
par_path1[1] = change_initial(par_nsys[ipar1[1]], zeros(length(par_nsys[ipar1[1]].sub)), zeros(length(par_nsys[ipar1[1]].sub))) 

# ╔═╡ ecb880c6-c912-4de3-93b9-fe67dcb9eaec
pl1 = Array{DataFrame}(undef, length(ipar1))

# ╔═╡ ea08662a-177f-422b-ab47-ec9b4963a762
sol1 = Array{Any}(undef, length(ipar1))

# ╔═╡ 13da4e88-1977-4322-923f-9c9669f85172
pl1[1], sol1[1] = GasChromatographySimulator.simulate(par_path1[1])

# ╔═╡ e40184c6-fbe9-4bce-83a3-4dc21c6fbd0a
par_path1[2] = change_initial(par_nsys[ipar1[2]], pl1[2-1].tR, pl1[2-1].τR) 

# ╔═╡ a48bbb8b-a42c-4d2e-9e04-45bca6688589
pl1[2], sol1[2] = GasChromatographySimulator.simulate(par_path1[2])

# ╔═╡ 740edb41-6bde-4fd3-850a-1ff2f4026f79
par_path1[3] = change_initial(par_nsys[ipar1[3]], pl1[3-1].tR, pl1[3-1].τR) 

# ╔═╡ 85d0ac9b-02ed-4cf8-9a48-e8701a88e79d
pl1[3], sol1[3] = GasChromatographySimulator.simulate(par_path1[3])

# ╔═╡ 5ec0fae0-f3f4-4379-bf72-b7be1aba1250
md"""
### Manually simulate path2 of nsys

Simulate along the several paths, reuse of allready simulated modules.
"""

# ╔═╡ e4377d5c-d9a8-4d8f-bb38-3091bbb094c2
path2 = paths_E[2]

# ╔═╡ 3e359060-f8de-41ee-b390-64bf7c6ce95c
path_possible(nsys, path2)

# ╔═╡ e17378a3-a43c-4d48-a785-19b9dad838a3
ipar2 = index_parameter(nsys.g, path2)

# ╔═╡ d366fe51-d833-43db-902b-18d3e9270b24
par_path2 = Array{GasChromatographySimulator.Parameters}(undef, length(ipar2))

# ╔═╡ 2d1b40a4-926c-4943-977c-62bc1a59c55d
pl2 = Array{DataFrame}(undef, length(ipar2))

# ╔═╡ dd73e349-51b0-4698-9ce4-b40a4210dcdd
sol2 = Array{Any}(undef, length(ipar2))

# ╔═╡ 271bef3b-0fc3-4199-8a87-4a04ff7424f9
common_edges(path2, path1)

# ╔═╡ 6ef9edcc-562f-4b6d-9849-53c2ab040ac6
common_ipar12 = index_parameter(nsys.g, common_edges(path2, path1))

# ╔═╡ 1c52ca3a-2d70-47f9-9377-193d122d51f1
par_path2[common_ipar12[1]] = par_path1[common_ipar12[1]]

# ╔═╡ 085b42eb-9c5c-4131-9614-b4fa18aa48ce
par_path2[common_ipar12[2]] = par_path1[common_ipar12[2]]

# ╔═╡ e3aae89f-a268-4b89-8934-0849c162b887
pl2[common_ipar12[1]] = pl1[common_ipar12[1]]

# ╔═╡ c6532915-7440-4502-9079-02902eab4ffb
pl2[common_ipar12[2]] = pl1[common_ipar12[2]]

# ╔═╡ 0b8c055c-e51e-49e3-bd63-e1845cf7d2ad
sol2[common_ipar12[1]] = sol1[common_ipar12[1]]

# ╔═╡ 06cd4272-567b-4967-ae50-6452b4d0fe98
sol2[common_ipar12[2]] = sol1[common_ipar12[2]]

# ╔═╡ 5e1c84bc-ea00-4aeb-abc6-bc7b286f6874
par_path2[3] = change_initial(par_nsys[ipar2[3]], pl2[3-1].tR, pl2[3-1].τR) 

# ╔═╡ 88896515-691b-4578-bd97-919c9fb45123
pl2[3], sol2[3] = GasChromatographySimulator.simulate(par_path2[3])

# ╔═╡ fc5f1f29-6b3d-4359-abd5-4f4b219018c5
par_path2[4] = change_initial(par_nsys[ipar2[4]], pl2[4-1].tR, pl2[4-1].τR) 

# ╔═╡ 8326c50e-9541-40ac-be81-2993d77ef8c4
pl2[4], sol2[4] = GasChromatographySimulator.simulate(par_path2[4])

# ╔═╡ b9ac802e-f3a1-4664-a0e4-9549a0b36a47
par_path2[5] = change_initial(par_nsys[ipar2[5]], pl2[5-1].tR, pl2[5-1].τR) 

# ╔═╡ 6f891b87-93a1-4f83-99ba-aecc68e67807
pl2[5], sol2[5] = GasChromatographySimulator.simulate(par_path2[5])

# ╔═╡ 49eae5d8-b740-4dbd-ad15-1c24edef59e3
par_path2[6] = change_initial(par_nsys[ipar2[6]], pl2[6-1].tR, pl2[6-1].τR)

# ╔═╡ 1ba35120-f42d-4a0e-a610-bc6860abc010
pl2[6], sol2[6] = GasChromatographySimulator.simulate(par_path2[6])

# ╔═╡ be1847b1-08bb-490f-a8ea-021dceb6c2a0
md"""
### Simulate from edge to edge
"""

# ╔═╡ 48188dab-d497-4384-bf13-7ae320ac9b35
par_nsys4 = change_initial(par_nsys[4], sol2[1].tR, sol2[1].τR) # first split point

# ╔═╡ 604084b7-4c3d-4c69-a30f-0ba20d892233
sol4 = GasChromatographySimulator.simulate(par_nsys4)

# ╔═╡ 94362d7b-8996-4f66-bd2a-07f8bdf32941
par_nsys5 = change_initial(par_nsys[5], sol4[1].tR, sol4[1].τR) # initial values from the result of the module before (src of this module = dst of the previous module) 

# ╔═╡ 1a8fdaba-6e72-4b7b-ad62-760b521e0571
sol5 = GasChromatographySimulator.simulate(par_nsys5)

# ╔═╡ 3837d27e-5d2e-4008-8cec-1fc47ecc9cc0
par_nsys6 = change_initial(par_nsys[6], sol5[1].tR, sol5[1].τR)

# ╔═╡ 09e0e77a-da3c-40e4-b30f-7ce6653d899c
sol6 = GasChromatographySimulator.simulate(par_nsys6)

# ╔═╡ dd3a144c-3e4e-431e-8f5e-22188bfd0f8f
par_nsys7 = change_initial(par_nsys[7], sol1[1].tR, sol1[1].τR)#change_initial(par_nsys[7], sol6[1].tR, sol6[1].τR)

# ╔═╡ 86959362-d31d-4d7e-aef9-7228a38c24ea
par_nsys7.prog.Fpin_itp(0)

# ╔═╡ 29a8d9f0-c5c0-4973-ae75-e90781afb1f1
par_nsys7.prog.pout_itp(0)

# ╔═╡ ca1caa5a-7396-4f24-b7f2-d9c1def63ec5
sol7 = GasChromatographySimulator.simulate(par_nsys7) # what happens here? -> DomainError

# ╔═╡ b8e845ba-138b-4b0b-a7e0-b06fb0041f61
#GasChromatographySimulator.simulate(GasChromatographySimulator.Parameters(par_nsys7.col, par_nsys7.prog, [par_nsys7.sub[end]], par_nsys7.opt))

# ╔═╡ 8cffd15a-a23c-4af6-b7ef-a25d86324adf
plot_graph(nsys, 3000)

# ╔═╡ c94bfe43-b1c8-4206-a438-8aadf9f7fb2b
par_nsys8 = change_initial(par_nsys[8], sol6[1].tR, sol6[1].τR) # 2nd split point

# ╔═╡ 18a7afb2-51ed-4032-ad99-b427b18a3800
sol8 = GasChromatographySimulator.simulate(par_nsys8)

# ╔═╡ 453ea126-7af2-400a-a8a9-b1d1651a8222
# reduction of retention time and negative velocity???? -> here we have the negative flow -> this result is not realistic

# ╔═╡ 07ae180a-17e0-4caf-8251-a6eb17bcdcc8
md"""
### Identity of `One Column` and `Series of columns`
#### One Column
"""

# ╔═╡ e1f42a3e-54a1-40e9-8d97-3e28c9798a84
par_nsys0 = graph_to_parameters(nsys0, db_dataframe, selected_solutes)

# ╔═╡ 62f8329f-bcee-4601-8d55-e86038d96588
par_nsys0_ = change_initial(par_nsys0[1], zeros(length(par_nsys0[1].sub)), zeros(length(par_nsys0[1].sub))) 

# ╔═╡ 56a6a1fb-ec84-4720-9c91-2292830648cd
sol_nsys0 = GasChromatographySimulator.simulate(par_nsys0_)

# ╔═╡ 2cdbcc48-10d0-4c1c-85f8-9fea1d539e85
md"""
#### Series of Column
"""

# ╔═╡ 3f8318f7-fbe9-43a3-ac74-3dfbecc885d2
par_nsys_series = graph_to_parameters(nsys_series, db_dataframe, selected_solutes)

# ╔═╡ 51bad243-d5d0-412e-bb67-e639ecd15f05
par_nsys_series1 = change_initial(par_nsys_series[1], zeros(length(par_nsys_series[1].sub)), zeros(length(par_nsys_series[1].sub))) 

# ╔═╡ 7636abfe-14f4-4abd-802a-23c26a3b3b78
sol_nsys_series1 = GasChromatographySimulator.simulate(par_nsys_series1)

# ╔═╡ eed83b74-6f19-493f-9e71-49d3879fa496
par_nsys_series2 = change_initial(par_nsys_series[2], sol_nsys_series1[1].tR, sol_nsys_series1[1].τR)

# ╔═╡ 73e17898-54c9-4ff2-b44e-839f89e3a185
#sol_nsys_series2 = GasChromatographySimulator.simulate(par_nsys_series2)

# ╔═╡ 80cd1c90-6668-4a85-8f00-3d12a192683b
#par_nsys_series3 = change_initial(par_nsys_series[3], sol_nsys_series2[1].tR, sol_nsys_series2[1].τR)

# ╔═╡ 894855c4-3285-42a4-a79d-fa845f41d330
#sol_nsys_series3 = GasChromatographySimulator.simulate(par_nsys_series3)

# ╔═╡ c0e0228e-649b-44b7-8c50-a8c8d61172db
sol_nsys0[2][2]

# ╔═╡ 7ad6f427-ab41-4633-9c4d-da2e9a7df0d2
begin
	plotly()
	p_trace = Plots.plot(xlabel="x in m", ylabel="t in s", legend=:bottomright)
	t0 = Array{Float64}(undef, length(sol_nsys0[2][2].t))
	for i=1:length(sol_nsys0[2][2].t)
		t0[i] = sol_nsys0[2][2].u[i][1]
	end
	Plots.scatter!(p_trace, sol_nsys0[2][2].t, t0)
	
	t1 = Array{Float64}(undef, length(sol_nsys_series1[2][2].t))
	for i=1:length(sol_nsys_series1[2][2].t)
		t1[i] = sol_nsys_series1[2][2].u[i][1]
	end
	Plots.scatter!(p_trace, sol_nsys_series1[2][2].t, t1)

	t2 = Array{Float64}(undef, length(sol_nsys_series2[2][2].t))
	for i=1:length(sol_nsys_series2[2][2].t)
		t2[i] = sol_nsys_series2[2][2].u[i][1]
	end
	Plots.scatter!(p_trace, sol_nsys_series2[2][2].t.+10.0, t2)

	t3 = Array{Float64}(undef, length(sol_nsys_series3[2][2].t))
	for i=1:length(sol_nsys_series3[2][2].t)
		t3[i] = sol_nsys_series3[2][2].u[i][1]
	end
	Plots.scatter!(p_trace, sol_nsys_series3[2][2].t.+20.0, t3)
	
	p_trace
end

# ╔═╡ 6d3d88d5-30d2-448d-81ff-fc1738bda199
begin
	plotly()
	F_func0 = flow_functions(nsys0)
	F_func_series = flow_functions(nsys_series)
	p_flow0_series = Plots.plot(xlabel="time in s", ylabel="flow in mL/min")
	for i=1:ne(nsys0.g)
		Plots.plot!(p_flow0_series, t, F_func0[i].(t)*60*1e6, label="F_$(nsys0.modules[i].name)")
	end
	for i=1:ne(nsys_series.g)
		Plots.plot!(p_flow0_series, t, F_func_series[i].(t)*60*1e6, label="F_$(nsys_series.modules[i].name)")
	end
	p_flow0_series
end

# ╔═╡ 2178a42a-0845-4d30-abe6-c056f6cc14e8
GasChromatographySimulator.local_plots("z", "u", sol_nsys0[2], par_nsys0_)

# ╔═╡ 54409689-a7c0-4a16-8035-ee0a25d227de
GasChromatographySimulator.local_plots("z", "u", sol_nsys_series1[2], par_nsys_series1)

# ╔═╡ 92710af2-3516-4d59-9621-19ba7e645295
GasChromatographySimulator.local_plots("z", "u", sol_nsys_series2[2], par_nsys_series2)

# ╔═╡ 43134a36-64c8-46c3-a5b3-3f5a632a1ea3
GasChromatographySimulator.local_plots("z", "u", sol_nsys_series3[2], par_nsys_series3)

# ╔═╡ debf351e-60a7-4f70-8bd9-5f4f83cb4239
pM0(x, t) = GasChromatographySimulator.pressure(x, t, par_nsys0_)

# ╔═╡ d673ee10-a534-48f8-9601-f02a7f35334f
pM_series1(x, t) = GasChromatographySimulator.pressure(x, t, par_nsys_series1)

# ╔═╡ 3399ed30-d00e-4a9c-ada0-33f66f0899a6
pM_series2(x, t) = GasChromatographySimulator.pressure(x, t, par_nsys_series2)

# ╔═╡ c87dad2a-5356-4eb6-b837-e9f199d2b4e7
pM_series3(x, t) = GasChromatographySimulator.pressure(x, t, par_nsys_series3)

# ╔═╡ c566dbb9-4525-470b-b499-e13a34a7c53b
par_nsys0_.prog.Fpin_itp(0)

# ╔═╡ c5f207cf-8374-4ba5-93ed-d0d99545fe5a
par_nsys0_.prog.pout_itp(1800)

# ╔═╡ 9161612d-5b82-44b1-9ca9-7ab0ad61e9d3
begin
	plotly()
	p_series, unknown_series = solve_pressure(nsys_series)
	pres_plot = Plots.plot(xlabel="x in m", ylabel="p in Pa", legend=:topright)
	Plots.plot!(pres_plot, 0:30, pM0.(0:30,0))
	Plots.plot!(pres_plot, 0:10, pM_series1.(0:10,0))
	Plots.plot!(pres_plot, 10:20, pM_series2.(0:10,0))
	Plots.plot!(pres_plot, 20:30, pM_series3.(0:10,0))
	Plots.scatter!((10, p_series[1](0)))
	Plots.scatter!((20, p_series[2](0)))
	pres_plot
end

# ╔═╡ f6d8ba16-b151-42ff-90b9-4259afd15a73
md"""
### [] Other Examples
"""

# ╔═╡ 11a2927c-69ef-4a46-b4a6-6818855fd867
begin
	g2 = SimpleDiGraph(14)
	add_edge!(g2, 1, 2)
	add_edge!(g2, 2, 3)
	add_edge!(g2, 2, 5)
	add_edge!(g2, 2, 6)
	add_edge!(g2, 3, 4)
	add_edge!(g2, 4, 10)
	add_edge!(g2, 5, 8)
	add_edge!(g2, 6, 9)
	add_edge!(g2, 6, 10)
	add_edge!(g2, 7, 11)
	add_edge!(g2, 9, 12)
	add_edge!(g2, 10, 13)
	add_edge!(g2, 11, 13)
	add_edge!(g2, 13, 14)

	fig2, ax2, p2 = GraphMakie.graphplot(g2, 
						layout=Stress(),
						nlabels=["$(i)" for i in 1:nv(g2)], 
						nlabels_align=(:center,:center),
						node_size = [40 for i in 1:nv(g2)],
						node_color = [:lightblue for i in 1:nv(g2)]
					)
	hidedecorations!(ax2)
	hidespines!(ax2)
	ax2.aspect = DataAspect()
	fig2
end

# ╔═╡ 70ebec42-dbb6-416f-8262-71a21eaac15e
all_paths(g2)

# ╔═╡ d7ea0d14-0353-421b-abaf-bf3718276cdd
md"""
### [] Path estimation 2
"""

# ╔═╡ a3f4b9fe-b129-4a84-98d7-02ab363c57fc
begin
	g3 = SimpleDiGraph(7)
	add_edge!(g3, 1, 2)
	add_edge!(g3, 2, 3)
	add_edge!(g3, 2, 4)
	add_edge!(g3, 4, 5)
	add_edge!(g3, 4, 6)
	add_edge!(g3, 5, 6)
	add_edge!(g3, 6, 7)

	fig3, ax3, p3 = GraphMakie.graphplot(g3, 
						layout=Stress(),
						nlabels=["$(i)" for i in 1:nv(g3)], 
						nlabels_align=(:center,:center),
						node_size = [40 for i in 1:nv(g3)],
						node_color = [:lightblue for i in 1:nv(g3)]
					)
	hidedecorations!(ax3)
	hidespines!(ax3)
	ax3.aspect = DataAspect()
	fig3
end

# ╔═╡ 37e8f35a-acf3-43a1-b4af-3dc0e46f7883
all_paths(g3)

# ╔═╡ 64a6b2a5-2fec-4582-9c7d-25a3cfa5adc6
md"""
- Path 1: 1->2, 2->3
- Path 2: 1->2, 2->4, 4->5, 5->6, 6->7
- Path 3: 1->2, 2->4, 4->6, 6->7
"""

# ╔═╡ 8ab53481-e145-41d8-952d-f621f4348afd
E3 = collect(edges(g3))

# ╔═╡ bda07d80-47cf-4c51-acba-bcc8ab4b7d32
split_merge(g3)

# ╔═╡ b6eea36f-7b68-4420-bbad-a129f5ce53a3
inlets_V3, outlets_V3 = inlets_outlets(g3)

# ╔═╡ 04265cc6-19e9-4f77-8873-63795eea6801


# ╔═╡ ab069891-a6b3-4171-98a4-9e915e9400dd


# ╔═╡ b218b461-a772-4a06-a223-638c6a64c46e


# ╔═╡ 80de0103-1133-48b0-9fa2-479b9c7b1cc6
all_paths(g3)

# ╔═╡ a027bbfd-acc6-422b-a844-aedad9a19b5d
all_paths(nsys.g, 3)

# ╔═╡ 81d681ae-c32b-49a7-9e90-15b7a2dd3f51
all_paths(nsys.g)

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
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
CSV = "~0.10.9"
CairoMakie = "~0.10.1"
DataFrames = "~1.4.4"
GasChromatographySimulator = "~0.3.15"
GraphMakie = "~0.5.1"
Graphs = "~1.7.4"
Interpolations = "~0.14.7"
NetworkLayout = "~0.4.4"
Plots = "~1.38.2"
PlutoUI = "~0.7.49"
Roots = "~2.0.8"
Symbolics = "~5.0.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "3e21177599949928fbde7429c956f26f8ec9220b"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Markdown", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "df23d15b1090a3332a09a7a51da45bd9f0a07f92"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.27.8"

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
git-tree-sha1 = "e97bdb5e241bb57f628968fd56efd9590078ada4"
uuid = "69666777-d1a9-59fb-9406-91d4454c9d45"
version = "2.4.2"

[[deps.ArrowTypes]]
deps = ["UUIDs"]
git-tree-sha1 = "563d60f89fcb730668bd568ba3e752ee71dde023"
uuid = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
version = "2.0.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

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

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

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
git-tree-sha1 = "5b735f654bdfd7b6c18c49f1d3ebff34b4b8af43"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.1"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "SnoopPrecompile", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "c700cce799b51c9045473de751e9319bdd1c6e94"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.9"

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
git-tree-sha1 = "c46adabdd0348f0ee8de91142cfc4a72a613ac0a"
uuid = "082447d4-558c-5d27-93f4-14fc19e9eca2"
version = "1.46.1"

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
git-tree-sha1 = "b410f0b8a52752e1c1723b4316382203f914672c"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.113.1"

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
git-tree-sha1 = "74911ad88921455c6afcad1eefa12bd7b1724631"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.80"

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
git-tree-sha1 = "387d2b8b3ca57b791633f0993b31d8cb43ea3292"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.71.3"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "5982b5e20f97bff955e9a2343a14da96a746cd8c"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.71.3+0"

[[deps.GasChromatographySimulator]]
deps = ["CSV", "ChemicalIdentifiers", "DataFrames", "ForwardDiff", "HypertextLiteral", "Integrals", "Interpolations", "OrdinaryDiffEq", "Plots", "PlutoUI", "Reexport", "UrlDownload"]
git-tree-sha1 = "f5e10e50f23366488b89a1584b5ad16ca5eae960"
uuid = "dd82b6e2-56ef-419d-b271-0be268cb65f5"
version = "0.3.15"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "fb69b2a645fa69ba5f474af09221b9308b160ce6"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.3"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "e315c4f9d43575cf6b4e511259433803c15ebaa2"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.1.0"

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
git-tree-sha1 = "3e2a15c851ea53cc28501c600f3df30647e3885b"
uuid = "1ecd5474-83a3-4783-bb4f-06765db800d2"
version = "0.5.1"

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
git-tree-sha1 = "eb5aa5e3b500e191763d35198f859e4b40fff4a6"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.7.3"

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
git-tree-sha1 = "2422f47b34d4b127720a18f86fa7b1aa2e141f29"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.18"

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
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "ArrayInterfaceOffsetArrays", "ArrayInterfaceStaticArrays", "CPUSummary", "ChainRulesCore", "CloseOpenIntervals", "DocStringExtensions", "ForwardDiff", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDDualNumbers", "SIMDTypes", "SLEEFPirates", "SnoopPrecompile", "SpecialFunctions", "Static", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "155132d68bc33c826dbdeb452c5d0a79e2d0e586"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.146"

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
git-tree-sha1 = "6503b77492fd7fcb9379bf73cd31035670e3c509"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.3"

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
git-tree-sha1 = "ca0c8939dbd3617ae3fdca13374d0b7501a2dd28"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.37.0"

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
git-tree-sha1 = "8175fc2b118a3755113c8e68084dc1a9e63c61ee"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.3"

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
git-tree-sha1 = "a99bbd3664bb12a775cda2eba7f3b2facf3dad94"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eadad7b14cf046de6eb41f13c9275e5aa2711ab6"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.49"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "7f8dd47630b265df9e1d117137ee1894b195e032"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.1"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "5d0a598c95f67ee0787723e38745cb954d143684"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.0"

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
git-tree-sha1 = "261dddd3b862bd2c940cf6ca4d1c8fe593e457c8"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.3"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase", "SnoopPrecompile"]
git-tree-sha1 = "e974477be88cb5e3040009f3767611bc6357846f"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.11"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceStaticArraysCore", "ChainRulesCore", "DocStringExtensions", "FillArrays", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "ZygoteRules"]
git-tree-sha1 = "fcf0962b399f3bc0fa13ae7274db7879c3ef9f1e"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.35.0"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "SnoopPrecompile", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "664aba1c5259821356f2ef771eabc502d67a8f0d"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.16"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

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

[[deps.Roots]]
deps = ["ChainRulesCore", "CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "a3db467ce768343235032a1ca0830fc64158dadf"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.0.8"

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
git-tree-sha1 = "8b20084a97b004588125caebf418d8cab9e393d1"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.4"

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
git-tree-sha1 = "c02bd3c9c3fc8463d3591a62a378f90d2d8ab0f3"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.17"

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
deps = ["ArrayInterfaceCore", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "Reexport", "SciMLBase", "SnoopPrecompile", "StaticArraysCore"]
git-tree-sha1 = "61b8ffdb22453132e02a10c5638dfb42834c776b"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "0.1.5"

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
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

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
git-tree-sha1 = "4149bb50ccf65bc9c9d55e7001ca7aa4a7649603"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.4"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

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
git-tree-sha1 = "70b6ee0e5cc1745a28dd9ba040b8e5ee28fffc69"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.4.5"

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
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TimerOutputs", "Unityper"]
git-tree-sha1 = "348ad5af9c916b6e1641c74378fac8bb49236688"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.0.1"

[[deps.Symbolics]]
deps = ["ArrayInterfaceCore", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Markdown", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TreeViews"]
git-tree-sha1 = "eb34c1e20b225c1de5adeeff3a085c9e985df532"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.0.0"

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

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "f8629df51cab659d70d2e5618a430b4d3f37f2c3"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.0"

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
git-tree-sha1 = "94f38103c984f89cf77c402f2a68dbd870f8165f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.11"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "SnoopPrecompile", "Static", "VectorizationBase"]
git-tree-sha1 = "6cca884e0fe17916da63c62dc1bf5896ce5d723e"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.17"

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

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "d5f4ec8c22db63bd3ccb239f640e895cfde145aa"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.2"

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
git-tree-sha1 = "6b1dc4fc039d273abc247eba675ac1299380e5d9"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.57"

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
git-tree-sha1 = "e76246c67099856cc9a85a3a809a57980eef2eae"
uuid = "e88e6eb3-aa80-5325-afca-941959d7151f"
version = "0.6.54"

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
# ╠═68c043ba-ce1e-4693-a9fd-f761d1edec58
# ╠═73d0e4d5-76c6-4860-92f5-5888ae605861
# ╠═8a9c2704-5854-4aa4-9782-943e82cd38f6
# ╠═471193bb-c2e1-4d27-9f0f-a9d65820ef7e
# ╠═a177a0a9-7e67-4644-896e-8863d4bd15c6
# ╠═291ed3ff-e500-48e1-b1ec-2d28b5cbfa54
# ╠═40663239-9c5b-466d-9f43-76ab913d5750
# ╠═f11b32e2-c447-4f54-8f42-c1e81c9fb489
# ╠═e7181f3b-e148-4327-a420-0091826b327c
# ╠═c5c010e6-80c6-445b-b977-0bf8559af0c4
# ╠═6f74b758-46a7-45cb-a177-8d766a8c42e8
# ╠═411fc6ec-cb7e-4ec9-b2f5-83b3b6347b52
# ╠═4d8c7f3a-582c-4fe9-a692-5966492cc31f
# ╠═8915304d-d902-4988-98c7-943f4f29478b
# ╠═a206d81b-86b6-41cd-9a5b-a291d361642d
# ╠═7cea41a4-a94d-4b44-8f32-04b324206d36
# ╠═b85e57b4-645b-452b-bae7-e55ee6cedc0e
# ╠═e7f1e626-b002-4848-aac4-891e5efc2d63
# ╠═775eaf12-eb3d-4dc1-891a-a1a2be066be4
# ╠═bc7a259b-9d43-4b7e-9fca-755f391ed890
# ╠═e923dc1f-c289-46f7-9d97-b2d4d8e04488
# ╠═70a6aff1-5547-46be-8884-17984c41fc64
# ╠═4c79d09f-5cd0-4788-ae06-a2be74bbc580
# ╠═f89e34e1-b10a-42c5-b526-e69517095159
# ╠═6d3fbb35-3f06-45c1-aa78-1c80997c1e21
# ╠═3c476391-74bd-45e5-847a-9c239c3e2f0b
# ╠═fb417960-388d-44f8-8a53-39b549ae3ec2
# ╠═627e0139-70da-490b-a85d-13024d76730f
# ╠═69dbe34c-2075-4d0e-bbd0-f7b13cce4ed8
# ╠═ab0f33cd-55bd-48e4-a1f2-ecd5e4720805
# ╠═0415d6bc-f28f-456c-a30f-3c2d48c90a9e
# ╠═e6a9dd2f-48b6-4420-8f30-1921b9013a0a
# ╠═117717f5-efa7-4e8a-aaa3-5cb245dec9c9
# ╠═3d792c30-4770-4eef-be13-dfd1ccb04e5c
# ╠═de33be5f-458b-4fdf-8d37-4a733ff41320
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
# ╠═8517a63b-198d-4658-a015-aab8ac7969ce
# ╠═47413580-d2c0-4ec1-9b85-7710ef24d4e2
# ╠═0645c482-7a21-4e2e-b574-6f733b7e1545
# ╠═3e1a13da-b7e3-4d40-82a2-2b91f58ed191
# ╠═ab162833-59f4-4b52-b3e2-0500aefd365a
# ╠═0e42d2f7-7b5a-4acb-bc8c-2ec7a7e890ff
# ╠═e6cd895b-3d69-41ff-8733-6cfaa50452f6
# ╠═42af3ef1-8969-4428-94e1-3195e03913fc
# ╠═ef85d984-b494-4b3f-8243-512d9eca4864
# ╠═cd9eb028-97be-4e5a-983a-7c0f09c3669b
# ╠═d83b9980-7d28-414c-910d-8314d613baf3
# ╠═540bb955-da9d-4628-8eaf-f86f918da1ed
# ╠═b141c05e-cccc-43d1-8bd5-495b2b14215e
# ╠═86572e5d-179f-4257-a96e-43fc075ccbab
# ╠═b0a7d30b-4e2c-49c8-acf6-428836c5d7c1
# ╠═ea162e25-2803-435a-aeb6-137f4365b3f7
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
# ╠═2d4366c7-bd2c-40fb-92c8-a03f32bf5fd0
# ╠═ed5d3089-bb8a-427a-a3ee-75986e6cc583
# ╠═f9a69795-bee2-46f5-80e9-9a262fa5b6b1
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
# ╠═887ec730-dde8-4fa1-8671-f0815af0c466
# ╠═44ebd883-232c-4a3c-a2a3-f7c8c85ad8e3
# ╠═8d7befb0-6225-404c-8dd2-12834978bde1
# ╠═9eedabe6-811c-41a0-a90e-f2f62fddb170
# ╠═01d9e464-3d6c-4ffe-a92d-dee1d88258bb
# ╠═e3fcbb11-cce6-4d5f-aa7c-13d1d11174b2
# ╠═93cc1041-8883-4c8e-88e5-a26389a28ee2
# ╠═70c82cf4-012e-4137-a0fd-7060cda73fb5
# ╠═30cdb150-fad1-43fa-8ef0-a776fc34b4d8
# ╠═7290bbc2-76d4-4b4f-9f46-e369bff7afae
# ╠═bce24317-737e-4b43-8ad0-4b82ffe9b973
# ╠═d8ee7f45-67a0-4b11-a401-4d2b1cf31e10
# ╠═dbf731b1-c897-4be0-a6fe-3358f0782cfd
# ╠═0c3eda88-5d13-43b2-b083-7577af04cf84
# ╠═203bf2a5-64eb-4f6d-a6f6-1bf12b19a968
# ╠═318075b8-76f3-4e9e-8655-5424cb42fe47
# ╠═afc385df-0184-4c97-8201-12a6ad714fd2
# ╠═37e8f35a-acf3-43a1-b4af-3dc0e46f7883
# ╠═70ebec42-dbb6-416f-8262-71a21eaac15e
# ╠═047bb299-8391-4da8-9c0e-2c4bca80371d
# ╠═7b13e4c3-be31-446e-8edd-b3e1c75cde4d
# ╠═97b19b85-841d-41a5-9a47-eaf497adfee2
# ╠═9b16093f-2db4-485e-9f6e-353dc8db886d
# ╠═aeb33e75-f64d-4b36-8c49-ec73e5dd4011
# ╠═47d4cdd7-b009-4f3d-834b-6b1809a5c5d2
# ╠═bdc69707-46d1-47ea-ba8d-6855c3b68251
# ╠═25696111-910c-42fb-942d-b0fb382d1334
# ╠═02ccb0b0-2494-4733-b1e5-55bf863064a5
# ╠═6e173c0a-fe02-455b-ab66-74e5cb6f0164
# ╠═dc1a71ad-3bbe-48b5-926f-514bb8bdf7e3
# ╠═8095da90-dc1d-46db-b383-9980e0485733
# ╠═4d3c608c-2d60-454f-8c4f-2c01fb780cf0
# ╠═033c19b4-f696-44ac-a23e-9a6f1e8368af
# ╠═a91f5cd2-adb4-488f-a936-e2f8c48d8f83
# ╠═e652a548-4c03-471a-b6a9-717bead6b355
# ╠═92f7549c-ee98-4074-934c-b6b69bde8bd9
# ╠═284553a9-5e63-4f47-ad9b-7a66c75c12b7
# ╠═8b9f4f16-0bbf-4319-bada-15f0275ad78d
# ╠═f4bac904-83ee-42a1-8169-cd335745e1f7
# ╠═d8b77441-b471-4efe-9e0f-93f369be63b8
# ╠═07eb454e-f035-4419-8689-626d120bf012
# ╠═a067f417-99e0-4844-bec6-dbbf9166205e
# ╠═f63a4d32-1877-466c-9df9-24da42365c08
# ╠═84fe5a51-3011-4c91-b070-c71f2a64656e
# ╠═cb46b45f-9a50-480e-99ea-0bbed4a4ea81
# ╠═4605cc1e-6885-4c3c-ad57-ae5643324c5f
# ╠═6d66210f-e504-48a3-8832-190fd7fb4cbb
# ╠═6a84d685-0e8d-41fc-b116-57b95de99814
# ╠═7bf076d8-3ae2-4c0a-ae82-db1164d7af4e
# ╠═be21c243-2aa9-4c9f-8fdf-3656717082cc
# ╠═60f69e16-60e3-4f17-ad79-74af60d74491
# ╠═5c738d5d-12ea-4ca9-94ee-da4c56d73318
# ╠═2ddab7b9-2d0c-4d07-9a5b-0314bedbe9de
# ╠═5e1f9b33-6bd3-4acb-9132-09551e13f8f8
# ╠═1e055bcd-0713-4992-83c5-f9621376f891
# ╠═5e9c922b-5eca-4c5b-872d-9fb73dc20389
# ╠═a1a2ace8-fbd8-4674-b26a-c29bfe2678f1
# ╠═61a93ebd-0680-4b5e-9cee-102915f7ed72
# ╠═f212b6d9-909a-480d-aacb-113a00a7fa99
# ╠═1cedb103-11af-42a3-a933-e0ee0652fe3e
# ╠═b02713df-0bf3-4ddf-bdfc-15d2890a912c
# ╠═e11e520b-c226-4077-8526-7793aa9a76f6
# ╠═e6b8a2dd-41e8-4d66-97c3-ccddae0e109e
# ╠═4f1080d6-8a11-4280-802b-0f5dc80fba08
# ╠═bf5b09d0-384c-4222-9f5c-ab2c869a564c
# ╠═05d97331-ef89-4c4b-a26d-ae61c50bcf03
# ╠═dde7e5a0-9073-48a7-a28a-1575072067c2
# ╠═4c39462a-1495-494c-b2c4-8a5d35c78f0c
# ╠═d659b94f-27d0-4b2d-a07a-7f77e0dde77e
# ╠═36b6623f-5c1e-4d33-b818-5f52f6e43eb4
# ╠═bf9779bd-b365-48c5-8dcc-44139397898a
# ╠═0f65afa8-c9df-4c93-8882-337cfc4fab5b
# ╠═4e495ff0-1ee1-45e6-abb3-8e6bd09424fa
# ╠═3565aac6-cdaa-4788-a73c-59a15b54711d
# ╠═72928529-a1aa-4b6a-b040-1c67f37f435a
# ╠═2d402085-68b5-485a-b5ac-9d82369e1311
# ╠═5dfa6a11-c56a-4fc4-9c16-75cf85b1fb8a
# ╠═bbab45a5-3b66-4b5d-9c52-7e0a4f1e0266
# ╠═da476db9-4603-4a51-b7ea-8f9bafd62d50
# ╠═50a2fc50-34ed-4ba4-a014-9a619bf7c071
# ╠═4a07f579-b9fc-47fc-ad18-11a9b12d9d0d
# ╠═7c768889-9473-486c-a6cc-f9b76c62a53b
# ╠═836ff6dc-4166-4096-8e3f-cc97856908ee
# ╠═00dde33e-fccb-411a-ba75-4dda4be9b035
# ╠═0575d3fb-101b-46d7-819d-c79a2125520a
# ╠═b23ed445-e112-415a-8c3d-952cd9aa69ae
# ╠═c67a3d24-164d-4205-8306-682f98bc6cbe
# ╠═0ca097e1-f821-48b2-bc07-a87944f16bef
# ╠═80337eea-b187-4a07-bc0a-c7e7b43e1b47
# ╠═76f4e987-163e-47af-8689-436b98c69348
# ╠═31d3fe45-bde0-4038-8157-c677997f2ae3
# ╠═ee039325-c34c-439c-8d1e-8b8938a6a670
# ╠═7660bc48-5e8c-40b4-a6dd-1a41c93923f0
# ╠═86b145f7-712e-4147-a49a-23d59a11246d
# ╠═4f87c3bc-fdcc-4211-9818-99ed006d5ca0
# ╠═f473d963-30c6-4786-aad8-7df9b7d48511
# ╠═b275a4e9-fd84-4ba3-96dd-8c298ea8a711
# ╠═233844fc-4beb-44db-a4a8-5eec12c28dd2
# ╠═3255a840-05ac-4490-809f-736a5e27353b
# ╠═8743715a-d966-40a6-ba92-a794b4b57c1c
# ╠═56d864dc-fc52-430b-8478-10aa32f95561
# ╠═0e9c7b54-34d0-4784-8f67-b9f907ab0e93
# ╠═70829b28-8efc-47dd-b12f-49f9cb29bf5a
# ╠═3213765f-4df8-4c1d-8189-ed29c8a5e667
# ╠═e1ac1411-6bf2-45d4-9e33-76a62e1d6ec3
# ╠═dd284c83-7825-44cb-9e01-928aa49897cb
# ╠═e5bbe44b-0f0d-4a66-a15f-0303a3f57978
# ╠═6babaa9e-103a-4525-afe9-ea9fd4f7b347
# ╠═79d90de0-79d4-4b88-bdc2-3ae2f172e695
# ╠═6ee58b3c-ebfd-4b9a-ae84-8b88916913d3
# ╠═96447c64-a8e6-4b59-af54-785b5bfe1803
# ╠═d8c2a71a-9854-4cee-a8d0-7654db698525
# ╠═7d238885-276c-4d3b-90e4-261b09534cde
# ╠═a0978e46-704e-43b7-8519-2e8ff698aae0
# ╠═45656162-bcdd-4395-af01-7b6f69eed043
# ╠═1e1911d1-f783-4f14-9202-a65762cfd8bd
# ╠═6dbccda5-8459-46b6-9966-431447d34b93
# ╠═669edc63-042b-4c55-99cc-0f19263de4e1
# ╠═6c482076-52b9-46a7-9b5f-fc0e80a59b2e
# ╠═8d63f820-8966-47f3-a60a-c784f51add52
# ╠═22b1221d-050d-4ff4-9a55-24c5ae6eaec8
# ╠═f765b518-ca59-42a5-b500-e3bf3fbece76
# ╠═1cb665cf-694c-49b7-8115-bb24bc98af27
# ╠═26da194d-0c9c-40f5-a5aa-d5054f1458dc
# ╠═7729714d-f0be-4be4-b9ed-39c17a6c5ffb
# ╠═479101b8-83b0-421f-925c-dd95a5df1286
# ╠═0d6e7d0c-cee2-4fdb-8782-e77bd6cd8445
# ╠═fa311bb6-3edb-434d-bca3-6efd9138231c
# ╠═5bc7a5d1-7106-484d-9078-f1d6acdd8c83
# ╠═4dc19ebe-e5c8-4301-8f93-f5cdb155927b
# ╠═b54430da-c90c-40dc-9eca-d6cd1d10d0a6
# ╠═5cc2f15a-d26d-4c79-889a-d398ead57f7e
# ╠═8e172e01-4f9f-4f23-8893-e9c8179558a9
# ╠═ecb880c6-c912-4de3-93b9-fe67dcb9eaec
# ╠═ea08662a-177f-422b-ab47-ec9b4963a762
# ╠═13da4e88-1977-4322-923f-9c9669f85172
# ╠═e40184c6-fbe9-4bce-83a3-4dc21c6fbd0a
# ╠═a48bbb8b-a42c-4d2e-9e04-45bca6688589
# ╠═740edb41-6bde-4fd3-850a-1ff2f4026f79
# ╠═85d0ac9b-02ed-4cf8-9a48-e8701a88e79d
# ╠═5ec0fae0-f3f4-4379-bf72-b7be1aba1250
# ╠═e4377d5c-d9a8-4d8f-bb38-3091bbb094c2
# ╠═3e359060-f8de-41ee-b390-64bf7c6ce95c
# ╠═e17378a3-a43c-4d48-a785-19b9dad838a3
# ╠═d366fe51-d833-43db-902b-18d3e9270b24
# ╠═2d1b40a4-926c-4943-977c-62bc1a59c55d
# ╠═dd73e349-51b0-4698-9ce4-b40a4210dcdd
# ╠═271bef3b-0fc3-4199-8a87-4a04ff7424f9
# ╠═6ef9edcc-562f-4b6d-9849-53c2ab040ac6
# ╠═1c52ca3a-2d70-47f9-9377-193d122d51f1
# ╠═085b42eb-9c5c-4131-9614-b4fa18aa48ce
# ╠═e3aae89f-a268-4b89-8934-0849c162b887
# ╠═c6532915-7440-4502-9079-02902eab4ffb
# ╠═0b8c055c-e51e-49e3-bd63-e1845cf7d2ad
# ╠═06cd4272-567b-4967-ae50-6452b4d0fe98
# ╠═5e1c84bc-ea00-4aeb-abc6-bc7b286f6874
# ╠═88896515-691b-4578-bd97-919c9fb45123
# ╠═fc5f1f29-6b3d-4359-abd5-4f4b219018c5
# ╠═8326c50e-9541-40ac-be81-2993d77ef8c4
# ╠═b9ac802e-f3a1-4664-a0e4-9549a0b36a47
# ╠═6f891b87-93a1-4f83-99ba-aecc68e67807
# ╠═49eae5d8-b740-4dbd-ad15-1c24edef59e3
# ╠═1ba35120-f42d-4a0e-a610-bc6860abc010
# ╠═be1847b1-08bb-490f-a8ea-021dceb6c2a0
# ╠═48188dab-d497-4384-bf13-7ae320ac9b35
# ╠═604084b7-4c3d-4c69-a30f-0ba20d892233
# ╠═94362d7b-8996-4f66-bd2a-07f8bdf32941
# ╠═1a8fdaba-6e72-4b7b-ad62-760b521e0571
# ╠═3837d27e-5d2e-4008-8cec-1fc47ecc9cc0
# ╠═09e0e77a-da3c-40e4-b30f-7ce6653d899c
# ╠═dd3a144c-3e4e-431e-8f5e-22188bfd0f8f
# ╠═86959362-d31d-4d7e-aef9-7228a38c24ea
# ╠═29a8d9f0-c5c0-4973-ae75-e90781afb1f1
# ╠═ca1caa5a-7396-4f24-b7f2-d9c1def63ec5
# ╠═b8e845ba-138b-4b0b-a7e0-b06fb0041f61
# ╠═8cffd15a-a23c-4af6-b7ef-a25d86324adf
# ╠═c94bfe43-b1c8-4206-a438-8aadf9f7fb2b
# ╠═18a7afb2-51ed-4032-ad99-b427b18a3800
# ╠═453ea126-7af2-400a-a8a9-b1d1651a8222
# ╠═07ae180a-17e0-4caf-8251-a6eb17bcdcc8
# ╠═e1f42a3e-54a1-40e9-8d97-3e28c9798a84
# ╠═62f8329f-bcee-4601-8d55-e86038d96588
# ╠═56a6a1fb-ec84-4720-9c91-2292830648cd
# ╠═2cdbcc48-10d0-4c1c-85f8-9fea1d539e85
# ╠═3f8318f7-fbe9-43a3-ac74-3dfbecc885d2
# ╠═51bad243-d5d0-412e-bb67-e639ecd15f05
# ╠═7636abfe-14f4-4abd-802a-23c26a3b3b78
# ╠═eed83b74-6f19-493f-9e71-49d3879fa496
# ╠═73e17898-54c9-4ff2-b44e-839f89e3a185
# ╠═80cd1c90-6668-4a85-8f00-3d12a192683b
# ╠═894855c4-3285-42a4-a79d-fa845f41d330
# ╠═c0e0228e-649b-44b7-8c50-a8c8d61172db
# ╠═7ad6f427-ab41-4633-9c4d-da2e9a7df0d2
# ╠═6d3d88d5-30d2-448d-81ff-fc1738bda199
# ╠═2178a42a-0845-4d30-abe6-c056f6cc14e8
# ╠═54409689-a7c0-4a16-8035-ee0a25d227de
# ╠═92710af2-3516-4d59-9621-19ba7e645295
# ╠═43134a36-64c8-46c3-a5b3-3f5a632a1ea3
# ╠═debf351e-60a7-4f70-8bd9-5f4f83cb4239
# ╠═d673ee10-a534-48f8-9601-f02a7f35334f
# ╠═3399ed30-d00e-4a9c-ada0-33f66f0899a6
# ╠═c87dad2a-5356-4eb6-b837-e9f199d2b4e7
# ╠═c566dbb9-4525-470b-b499-e13a34a7c53b
# ╠═c5f207cf-8374-4ba5-93ed-d0d99545fe5a
# ╠═9161612d-5b82-44b1-9ca9-7ab0ad61e9d3
# ╠═f6d8ba16-b151-42ff-90b9-4259afd15a73
# ╠═11a2927c-69ef-4a46-b4a6-6818855fd867
# ╠═d7ea0d14-0353-421b-abaf-bf3718276cdd
# ╠═a3f4b9fe-b129-4a84-98d7-02ab363c57fc
# ╟─64a6b2a5-2fec-4582-9c7d-25a3cfa5adc6
# ╠═8ab53481-e145-41d8-952d-f621f4348afd
# ╠═bda07d80-47cf-4c51-acba-bcc8ab4b7d32
# ╠═b6eea36f-7b68-4420-bbad-a129f5ce53a3
# ╠═04265cc6-19e9-4f77-8873-63795eea6801
# ╠═ab069891-a6b3-4171-98a4-9e915e9400dd
# ╠═b218b461-a772-4a06-a223-638c6a64c46e
# ╠═80de0103-1133-48b0-9fa2-479b9c7b1cc6
# ╠═a027bbfd-acc6-422b-a844-aedad9a19b5d
# ╠═81d681ae-c32b-49a7-9e90-15b7a2dd3f51
# ╠═6dd088fa-5907-400d-bf13-4890a224401e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
