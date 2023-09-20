### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 7f5b6a3c-ed83-11ed-2af0-4b696e429dc4
begin
	import Pkg
    # activate the shared project environment
	Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()

	using CSV, DataFrames
	using Plots, CairoMakie, GraphMakie
	using Graphs, NetworkLayout, Symbolics
	using GasChromatographySimulator
	using PlutoUI
	include(joinpath(dirname(pwd()), "src", "GasChromatographySystems.jl"))
	TableOfContents()
end

# ╔═╡ 60ed17af-1494-4dd4-ab5a-9ed76ed82c11
md"""
# SNAT modulation
## Splitter-based non-cryogenic artificial trapping modulation
according to Tungkijanansin2022
"""

# ╔═╡ af37fc60-57a4-4af9-b78e-fe467f3e9164
default_TP = GasChromatographySystems.TemperatureProgram([0.0, 1500.0], [40.0, 240.0])

# ╔═╡ 71e7d4bb-faa8-43da-9bab-41046ca38ade
begin
	Ldratio² = [10e7, 8e7, 6e7, 4e7]
	d = [0.1e-3, 0.15e-3, 0.18e-3, 0.25e-3]
	L = sqrt.(Ldratio².*d.^2)
end

# ╔═╡ 5548429b-da44-4abd-bb4d-00d9ab66649d
L./d.^4

# ╔═╡ 07fb0cea-e5aa-43be-900e-c7ff91f6e639
begin
	Ld4 = 1e16
	L_ = d.^4.0*Ld4
end

# ╔═╡ 31190533-2581-42f9-b352-23ae1cd54f10
#=begin
	g = SimpleDiGraph(10)
	add_edge!(g, 1, 2) # GC1
	add_edge!(g, 2, 3) # Split 1a
	add_edge!(g, 2, 4) # Split 1b
	add_edge!(g, 3, 5) # Split 2b1
	add_edge!(g, 3, 7) # Split 2a
	add_edge!(g, 4, 6) # Split 3b1
	add_edge!(g, 4, 8) # Split 3a
	add_edge!(g, 5, 7) # Split 2b2
	add_edge!(g, 6, 8) # Split 3b2
	add_edge!(g, 7, 9) # Split 4a
	add_edge!(g, 8, 9) # Split 4b
	add_edge!(g, 9, 10) # GC2
	# pressure points
	pp = Array{GasChromatographySystems.PressurePoint}(undef, nv(g))
	pp[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 1800.0], [NaN, NaN]) # inlet 
	pp[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 1800.0], [NaN, NaN]) #
	pp[3] = GasChromatographySystems.PressurePoint("p₃", [0.0, 1800.0], [NaN, NaN])
	pp[4] = GasChromatographySystems.PressurePoint("p₄", [0.0, 1800.0], [NaN, NaN])
	pp[5] = GasChromatographySystems.PressurePoint("p₅", [0.0, 1800.0], [NaN, NaN])
	pp[6] = GasChromatographySystems.PressurePoint("p₆", [0.0, 1800.0], [NaN, NaN])
	pp[7] = GasChromatographySystems.PressurePoint("p₇", [0.0, 1800.0], [NaN, NaN])
	pp[8] = GasChromatographySystems.PressurePoint("p₈", [0.0, 1800.0], [NaN, NaN])
	pp[9] = GasChromatographySystems.PressurePoint("p₉", [0.0, 1800.0], [NaN, NaN])
	pp[10] = GasChromatographySystems.PressurePoint("p₁₀", [0.0, 1800.0], [eps(Float64), eps(Float64)]) # outlet 1
	# modules
	modules = Array{GasChromatographySystems.AbstractModule}(undef, ne(g))
	modules[1] = GasChromatographySystems.ModuleColumn("GC1", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP, 1.0/60e6)
	modules[2] = GasChromatographySystems.ModuleColumn("Sp1a", 0.5, 0.25e-3, 0.0e-6, "", default_TP, NaN)
	modules[3] = GasChromatographySystems.ModuleColumn("Sp1b", 0.5, 0.25e-3, 0.0e-6, "", default_TP, NaN)
	modules[4] = GasChromatographySystems.ModuleColumn("Sp2b1", L[1]/2, d[1], 0.0e-6, "", default_TP, NaN)
	modules[5] = GasChromatographySystems.ModuleColumn("Sp2a", L[3], d[3], 0.0e-6, "", default_TP, NaN)#3
	modules[6] = GasChromatographySystems.ModuleColumn("Sp3b1", L[2]/2, d[2], 0.0e-6, "", default_TP, NaN)#2
	modules[7] = GasChromatographySystems.ModuleColumn("Sp3a", L[4], d[4], 0.0e-6, "", default_TP, NaN)#4
	modules[8] = GasChromatographySystems.ModuleColumn("Sp2b2", L[1]/2, d[1], 0.0e-6, "", default_TP, NaN)
	modules[9] = GasChromatographySystems.ModuleColumn("Sp3b2", L[2]/2, d[2], 0.0e-6, "", default_TP, NaN)#2
	modules[10] = GasChromatographySystems.ModuleColumn("Sp4a", 0.5, 0.25e-3, 0.0e-6, "", default_TP, NaN)
	modules[11] = GasChromatographySystems.ModuleColumn("Sp4b", 0.5, 0.25e-3, 0.0e-6, "", default_TP, NaN)
	modules[12] = GasChromatographySystems.ModuleColumn("GC2", 30.0, 0.25e-3, 0.25e-6, "Wax", default_TP, NaN)
	# system
	sys = GasChromatographySystems.update_system(GasChromatographySystems.System(g, pp, modules, GasChromatographySystems.Options(ng=true, vis="HP")))
end=#

# ╔═╡ deb6198d-6e0a-471b-9ab3-30bd21ca37a1
begin # column config from Tungkijanansin2022
	g = SimpleDiGraph(10)
	add_edge!(g, 1, 2) # GC1
	add_edge!(g, 2, 3) # Split 1a
	add_edge!(g, 2, 4) # Split 1b
	add_edge!(g, 3, 5) # Split 2b1
	add_edge!(g, 3, 7) # Split 2a
	add_edge!(g, 4, 6) # Split 3b1
	add_edge!(g, 4, 8) # Split 3a
	add_edge!(g, 5, 7) # Split 2b2
	add_edge!(g, 6, 8) # Split 3b2
	add_edge!(g, 7, 9) # Split 4a
	add_edge!(g, 8, 9) # Split 4b
	add_edge!(g, 9, 10) # GC2
	# pressure points
	pp = Array{GasChromatographySystems.PressurePoint}(undef, nv(g))
	pp[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 1800.0], [475000, 475000]) # inlet 
	pp[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 1800.0], [NaN, NaN]) #
	pp[3] = GasChromatographySystems.PressurePoint("p₃", [0.0, 1800.0], [NaN, NaN])
	pp[4] = GasChromatographySystems.PressurePoint("p₄", [0.0, 1800.0], [NaN, NaN])
	pp[5] = GasChromatographySystems.PressurePoint("p₅", [0.0, 1800.0], [NaN, NaN])
	pp[6] = GasChromatographySystems.PressurePoint("p₆", [0.0, 1800.0], [NaN, NaN])
	pp[7] = GasChromatographySystems.PressurePoint("p₇", [0.0, 1800.0], [NaN, NaN])
	pp[8] = GasChromatographySystems.PressurePoint("p₈", [0.0, 1800.0], [NaN, NaN])
	pp[9] = GasChromatographySystems.PressurePoint("p₉", [0.0, 1800.0], [NaN, NaN])
	pp[10] = GasChromatographySystems.PressurePoint("p₁₀", [0.0, 1800.0], [eps(Float64), eps(Float64)]) # outlet 1
	# modules
	modules = Array{GasChromatographySystems.AbstractModule}(undef, ne(g))
	modules[1] = GasChromatographySystems.ModuleColumn("GC1", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP, NaN)
	modules[2] = GasChromatographySystems.ModuleColumn("Sp1a", 0.5, 0.25e-3, 0.0e-6, "", default_TP, NaN)
	modules[3] = GasChromatographySystems.ModuleColumn("Sp1b", 0.5, 0.25e-3, 0.0e-6, "", default_TP, NaN)
	modules[4] = GasChromatographySystems.ModuleColumn("Sp2b1", 1.42, 0.25e-3, 0.0e-6, "", default_TP, NaN)
	modules[5] = GasChromatographySystems.ModuleColumn("Sp2a", 0.14, 0.1e-3, 0.0e-6, "", default_TP, NaN)#3
	modules[6] = GasChromatographySystems.ModuleColumn("Sp3b1", 2.93, 0.25e-3, 0.0e-6, "", default_TP, NaN)#2
	modules[7] = GasChromatographySystems.ModuleColumn("Sp3a", 4.46, 0.25e-3, 0.0e-6, "", default_TP, NaN)#4
	modules[8] = GasChromatographySystems.ModuleColumn("Sp2b2", 0.095, 0.1e-3, 0.0e-6, "", default_TP, NaN)
	modules[9] = GasChromatographySystems.ModuleColumn("Sp3b2", 0.05, 0.1e-3, 0.0e-6, "", default_TP, NaN)#2
	modules[10] = GasChromatographySystems.ModuleColumn("Sp4a", 0.5, 0.25e-3, 0.0e-6, "", default_TP, NaN)
	modules[11] = GasChromatographySystems.ModuleColumn("Sp4b", 0.5, 0.25e-3, 0.0e-6, "", default_TP, NaN)
	modules[12] = GasChromatographySystems.ModuleColumn("GC2", 30.0, 0.25e-3, 0.25e-6, "Wax", default_TP, NaN)
	# system
	sys = GasChromatographySystems.update_system(GasChromatographySystems.System(g, pp, modules, GasChromatographySystems.Options(ng=true, vis="HP")))
end

# ╔═╡ 0e3164b6-d583-477f-843f-f6d3d8a390d4
#=begin # column config from Tungkijanansin2022 but const flow
	g = SimpleDiGraph(10)
	add_edge!(g, 1, 2) # GC1
	add_edge!(g, 2, 3) # Split 1a
	add_edge!(g, 2, 4) # Split 1b
	add_edge!(g, 3, 5) # Split 2b1
	add_edge!(g, 3, 7) # Split 2a
	add_edge!(g, 4, 6) # Split 3b1
	add_edge!(g, 4, 8) # Split 3a
	add_edge!(g, 5, 7) # Split 2b2
	add_edge!(g, 6, 8) # Split 3b2
	add_edge!(g, 7, 9) # Split 4a
	add_edge!(g, 8, 9) # Split 4b
	add_edge!(g, 9, 10) # GC2
	# pressure points
	pp = Array{GasChromatographySystems.PressurePoint}(undef, nv(g))
	pp[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 1800.0], [NaN, NaN]) # inlet 
	pp[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 1800.0], [NaN, NaN]) #
	pp[3] = GasChromatographySystems.PressurePoint("p₃", [0.0, 1800.0], [NaN, NaN])
	pp[4] = GasChromatographySystems.PressurePoint("p₄", [0.0, 1800.0], [NaN, NaN])
	pp[5] = GasChromatographySystems.PressurePoint("p₅", [0.0, 1800.0], [NaN, NaN])
	pp[6] = GasChromatographySystems.PressurePoint("p₆", [0.0, 1800.0], [NaN, NaN])
	pp[7] = GasChromatographySystems.PressurePoint("p₇", [0.0, 1800.0], [NaN, NaN])
	pp[8] = GasChromatographySystems.PressurePoint("p₈", [0.0, 1800.0], [NaN, NaN])
	pp[9] = GasChromatographySystems.PressurePoint("p₉", [0.0, 1800.0], [NaN, NaN])
	pp[10] = GasChromatographySystems.PressurePoint("p₁₀", [0.0, 1800.0], [eps(Float64), eps(Float64)]) # outlet 1
	# modules
	modules = Array{GasChromatographySystems.AbstractModule}(undef, ne(g))
	modules[1] = GasChromatographySystems.ModuleColumn("GC1", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP, 1/60e6)
	modules[2] = GasChromatographySystems.ModuleColumn("Sp1a", 0.5, 0.25e-3, 0.0e-6, "", default_TP, NaN)
	modules[3] = GasChromatographySystems.ModuleColumn("Sp1b", 0.5, 0.25e-3, 0.0e-6, "", default_TP, NaN)
	modules[4] = GasChromatographySystems.ModuleColumn("Sp2b1", 1.42, 0.25e-3, 0.0e-6, "", default_TP, NaN)
	modules[5] = GasChromatographySystems.ModuleColumn("Sp2a", 0.14, 0.1e-3, 0.0e-6, "", default_TP, NaN)#3
	modules[6] = GasChromatographySystems.ModuleColumn("Sp3b1", 2.93, 0.25e-3, 0.0e-6, "", default_TP, NaN)#2
	modules[7] = GasChromatographySystems.ModuleColumn("Sp3a", 4.46, 0.25e-3, 0.0e-6, "", default_TP, NaN)#4
	modules[8] = GasChromatographySystems.ModuleColumn("Sp2b2", 0.095, 0.1e-3, 0.0e-6, "", default_TP, NaN)
	modules[9] = GasChromatographySystems.ModuleColumn("Sp3b2", 0.05, 0.1e-3, 0.0e-6, "", default_TP, NaN)#2
	modules[10] = GasChromatographySystems.ModuleColumn("Sp4a", 0.5, 0.25e-3, 0.0e-6, "", default_TP, NaN)
	modules[11] = GasChromatographySystems.ModuleColumn("Sp4b", 0.5, 0.25e-3, 0.0e-6, "", default_TP, NaN)
	modules[12] = GasChromatographySystems.ModuleColumn("GC2", 30.0, 0.25e-3, 0.25e-6, "Wax", default_TP, NaN)
	# system
	sys = GasChromatographySystems.update_system(GasChromatographySystems.System(g, pp, modules, GasChromatographySystems.Options(ng=true, vis="HP")))
end=#

# ╔═╡ 87fd20bd-2e93-401a-8978-6e4714bd334d
GasChromatographySystems.plot_graph(sys)

# ╔═╡ 15faeed7-066a-4f66-b64f-9cd915437ca0
GasChromatographySystems.plot_graph_with_flow(sys,0)

# ╔═╡ 6e5cc5cd-7a17-4505-b2cb-caec7f4e39d1
GasChromatographySystems.plot_graph_with_flow(sys,1800)

# ╔═╡ efb6cb49-811b-4a29-8eac-93c79a10b765
GasChromatographySystems.flow_balance(sys)

# ╔═╡ 2bb93f4d-3d7b-47ad-86e6-3a958a915a7e
GasChromatographySystems.substitute_unknown_flows(sys)

# ╔═╡ 02d8a178-410c-42f5-8784-596df16218b1
# ATTENETION Notebook chrashes with this GasChromatographySystems.solve_balance(sys)

# ╔═╡ 92e7ea3c-64a4-4282-bdb2-596986a87278
sol_bal = GasChromatographySystems.solve_balance(sys);

# ╔═╡ 16857ed2-4794-4346-94b7-114ad2a4a294
typeof(sol_bal)

# ╔═╡ 7fb06668-5d56-4304-baea-fa8d32a36d82
#sol_bal[1]

# ╔═╡ 556e42d6-d0b0-4c0a-9044-f5462fd7de23
function plot_graph_with_flow(sys, t; lay = Spring(), color=:lightblue, node_size=80, arrow_size=20, arrow_shift=0.8, dataaspect=false, nlabels_fontsize=14, elabels_fontsize=14, elabels_distance = 20)
	p_func = GasChromatographySystems.pressure_functions(sys)
	F_func = GasChromatographySystems.flow_functions(sys)
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

# ╔═╡ d6e02463-f9df-4658-ba6c-865f2640c3b9
function plot_graph_with_flow(sys, F_func, p_func, t; lay = Spring(), color=:lightblue, node_size=80, arrow_size=20, arrow_shift=0.8, dataaspect=false, nlabels_fontsize=14, elabels_fontsize=14, elabels_distance = 20)
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

# ╔═╡ 803485e0-fd65-4794-98ae-b79fb13d9518
p_func = GasChromatographySystems.pressure_functions(sys)

# ╔═╡ 72f68999-82e6-4660-b7eb-2114a16898af
F_func = GasChromatographySystems.flow_functions(sys)

# ╔═╡ 8fcbf371-3ca5-41d8-8149-7af890a6b90b
plot_graph_with_flow(sys,F_func, p_func, 0)

# ╔═╡ b04bd252-fea8-41ed-a249-cc564241a4ec
function flow_functions(sys, p_func)
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

# ╔═╡ 2b332378-2933-48b0-bd69-9c8687f13965
flow_functions(sys, p_func)

# ╔═╡ 2aa1f5f1-2f60-4eaf-ae82-949df1532fd9
function plot_flow_over_time(sys, F_func; dt=60.0)
	com_timesteps = GasChromatographySystems.common_timesteps(sys)
	trange = 0:dt:sum(com_timesteps)
	p_flow = Plots.plot(xlabel="time in s", ylabel="flow in mL/min")
	for i=1:ne(sys.g)
		Plots.plot!(p_flow, trange, F_func[i].(trange).*60e6, label="F_$(sys.modules[i].name)")
	end
	return p_flow
end

# ╔═╡ 523c7c64-84a7-4612-99c0-c0fdd23591df
begin
	gr()
	#p_flow = plot_flow_over_time(sys, F_func; dt=60.0)
end

# ╔═╡ 55f06025-db1c-4b42-9da8-5cb2d95e0234
#savefig(p_flow, "MultipSplit_flow.svg")

# ╔═╡ 714b3bdd-2b2b-44af-8672-f21245c57634
plotly()

# ╔═╡ 9df9bdcf-7efb-4cf3-a379-b89a4e81f25d
begin
	gr()
	#p_pres = GasChromatographySystems.plot_pressure_over_time(sys; dt=60.0)
	#Plots.plot!(p_pres, ylims=(1.5e5, 3.5e5))
	#p_pres
end

# ╔═╡ 60680135-ada1-48ae-bbef-b587cb81c50c
#savefig(p_pres, "MultipSplit_pres.svg")

# ╔═╡ 22e80b28-7670-4d80-a65b-942135d63a8a
GasChromatographySystems.flow_restrictions(sys)

# ╔═╡ d828afa5-8335-434a-96be-041cbd53a5a8
begin
	Ld_ratios = Array{Float64}(undef, length(sys.modules))
	for i=1:length(sys.modules)
		Ld_ratios[i] = sys.modules[i].length/sys.modules[i].diameter
	end
	Ld_ratios
end

# ╔═╡ 9eecf202-d3e7-405a-99d7-6a17a12d481a
md"""
## ----
"""

# ╔═╡ a1dcef3e-3d8f-418a-875a-4d05258d39a3
db = filter([:CAS] => x -> !ismissing(x), DataFrame(CSV.File("/Users/janleppert/Documents/GitHub/RetentionData/Databases/GCSim_database_nonflag.csv", header=1, silencewarnings=true)))

# ╔═╡ bd81cafc-208f-4afb-94af-990c6697335b
com_solutes = GasChromatographySystems.common_solutes(db, sys).Name

# ╔═╡ 70bc6f4d-a4aa-47e7-92f0-671d8f230e6a
md"""
## Graph to Parameters
"""

# ╔═╡ 54151540-4150-49c0-a260-d755d663462b
par = GasChromatographySystems.graph_to_parameters(sys, db, com_solutes[1:4]; interp=true, dt=60)

# ╔═╡ b531f3cf-1d0a-453a-ad7b-16c113991ba2
paths = GasChromatographySystems.all_paths(sys.g, 4)[2]

# ╔═╡ bd9955af-9d6e-4092-8a8e-8cb62b1e97a9
md"""
## Holdup times
"""

# ╔═╡ 94b05220-e3c2-4ccb-8193-64aa3b224915
GasChromatographySimulator.holdup_time(0.0, par[1])

# ╔═╡ efda5295-ca67-4746-a84b-6de830e87e09
GasChromatographySystems.index_parameter(sys.g, paths[1])

# ╔═╡ 399d3c27-4db3-4165-af03-4512df3262d7
function holdup_times_of_paths(t, paths, sys)
	tMs = Array{Float64,1}(undef, length(paths))
	tM = Array{Array{Float64,1},1}(undef, length(paths))
	for i=1:length(paths)
		i_par = GasChromatographySystems.index_parameter(sys.g, paths[i])
		tM_ = Array{Float64}(undef, length(i_par))
		for j=1:length(i_par)
			tM_[j] = GasChromatographySimulator.holdup_time(t, par[i_par[j]])
		end
		tMs[i] = sum(tM_)
		tM[i] = tM_
	end
	return tMs, tM
end

# ╔═╡ 2d856b07-7739-4b9b-a7f2-a2dedf5e882e
tMs = holdup_times_of_paths(0.0, paths, sys)[1]

# ╔═╡ 4bd3f0ea-03dd-4dd7-a8d8-99f20bcdb27a
tMs[1]-tMs[3]

# ╔═╡ a0314b36-7840-4fbe-8084-12d7c4fd5f9c
tMs[3]-tMs[2]

# ╔═╡ 9596cd2b-2e05-4247-a889-e64186fcca9a
tMs[2]-tMs[4]

# ╔═╡ 20b6a285-6449-4ec5-b0d3-a7c443978c0f
tMs[4]-tMs[3]

# ╔═╡ 82bb824c-a0c9-46a2-8aec-603fce1bf66d
tMs[3]-tMs[1]

# ╔═╡ 4864997a-b053-40bf-81e9-ff0a8b28867f
tMs[1]-tMs[2]

# ╔═╡ 80da63de-c461-4727-b57f-34856be6a71e
begin
	t = 0.0:1.0:1800.0
	tM = Array{Float64}(undef, length(t), 4)
	for i=1:length(t)
		tM_ = holdup_times_of_paths(t[i], paths, sys)[1]
		for j=1:4
			tM[i,j] = tM_[j]
		end
	end
	Plots.plot(t, tM[:,4].-tM[:,3], label="tM4-tM3")
	Plots.plot!(t, tM[:,3].-tM[:,1], label="tM3-tM1")
	Plots.plot!(t, tM[:,1].-tM[:,2], label="tM1-tM2")
	# these holdup time differences between the different paths should be the same, but could vary during the program
end

# ╔═╡ 5101c8d7-7e16-407a-b9cb-e726994034da
begin
	Plots.plot(t, tM[:,1], label="tM1")
	Plots.plot!(t, tM[:,2], label="tM2")
	Plots.plot!(t, tM[:,3], label="tM3")
	Plots.plot!(t, tM[:,4], label="tM4")
	# the holdup times between the different paths should be different
end

# ╔═╡ e13060a8-9ab7-46b5-807d-93f577bf7ce3
f(x)=0.0

# ╔═╡ a07e2feb-5ba9-4d35-a776-a01d46ab84bb
f_(x) = x

# ╔═╡ 83978df2-41f5-48f1-a3a8-88a346730490
f__(x) = f_(x) + f(x)

# ╔═╡ 685cd6f8-c3c6-4f53-b6c5-0814e633c17c
f__(10)

# ╔═╡ 06d4d265-0353-46ee-ad7f-6fd4df88d11c
sim_res = GasChromatographySystems.simulate_along_paths(sys, paths, par)

# ╔═╡ 33e4d5c0-42f8-4c3e-8465-fb41cc38336e
sim_res[3][1]

# ╔═╡ 4d1b2744-29da-439b-aa1a-e45138823e60
sim_res[3][2]

# ╔═╡ 414b08d5-8fa5-40d0-831c-0b63864f9cde
sim_res[2][1][end]

# ╔═╡ 360d72d3-4941-456d-9e75-e0a6d821acfa
sim_res[2][2][end]

# ╔═╡ 4122c986-ed3a-4303-a231-15094e4ebbc2
sim_res[2][3][end]

# ╔═╡ 2b2fb924-28ae-4e4b-98f9-5afb8144c206
sim_res[2][4][end]

# ╔═╡ f2e722fb-2edf-4159-9f9e-a0c841b4adf1
sim_res[2][1][end].tR[1] - sim_res[2][2][end].tR[1]

# ╔═╡ b9f39f74-85df-4b69-90fa-e033fa48114d
sim_res[2][3][end].tR[1] - sim_res[2][1][end].tR[1]

# ╔═╡ 5cd75ea5-8519-45da-bdf6-8eca2da64065
sim_res[2][4][end].tR[1] - sim_res[2][3][end].tR[1]

# ╔═╡ 85513069-ebee-4b1a-a2d2-b73d9dbdf848
sim_res[2][2][end].tR[1] - sim_res[2][2][1].tR[1] # 'time' after 1st D

# ╔═╡ ae818c15-ea60-45f8-a65f-4e2f6ad1b119
sim_res[2][1][end].tR[1] - sim_res[2][1][1].tR[1] # 'time' after 1st D

# ╔═╡ 410a3252-5541-4cbb-8254-f6a5d89d46cd
sim_res[2][3][end].tR[1] - sim_res[2][3][1].tR[1] # 'time' after 1st D

# ╔═╡ d6cb1be1-6df6-40f0-8374-398936cb17d3
sim_res[2][4][end].tR[1] - sim_res[2][4][1].tR[1] # 'time' after 1st D

# ╔═╡ 7636143f-79cd-4f1c-b939-0b031e1fbbe6
md"""
## Chromatograms
"""

# ╔═╡ efc2bd49-554b-4148-9516-d2668ef7d462
md"""
### Split of the flows
"""

# ╔═╡ 4a4e28a0-aeb1-4fc8-9113-711b7b4eb26d
# time of 1st split
tsplit1 = sim_res[2][1][1].tR[1]

# ╔═╡ 4657b632-4f49-4cd9-bfec-262805cb3c0b
F_func[2](tsplit1)

# ╔═╡ 145758e9-240a-471f-b113-aa33f9c3960c
F_func[3](tsplit1)

# ╔═╡ c7d1051e-c65f-4a95-8cc4-a29f86039fa7
split_ratio1a = F_func[2](tsplit1)/F_func[1](tsplit1)

# ╔═╡ 2794959c-9b5a-4e34-ae3d-dc3c687c7747
split_ratio1b = F_func[3](tsplit1)/F_func[1](tsplit1)

# ╔═╡ fa85bc90-401d-41c8-a19e-d68e12707b92
# time of second split on splitter 2
tsplit2 = sim_res[2][1][2].tR[1]

# ╔═╡ 63b4f098-c82d-4a0f-b83e-da92b5780d5c
split_ratio2b = F_func[4](tsplit2)/F_func[2](tsplit2)

# ╔═╡ abcc4a79-1ec9-4301-b4a0-3d1c2c41d2fc
split_ratio2a = F_func[5](tsplit2)/F_func[2](tsplit2)

# ╔═╡ 36cb9e74-e510-4305-87a2-86fac95fe041
# time of second split on splitter 3
tsplit3 = sim_res[2][3][2].tR[1]

# ╔═╡ a8d4b540-a9ce-44c0-94db-24931b68c9f3
split_ratio3b = F_func[6](tsplit3)/F_func[3](tsplit3)

# ╔═╡ 6dca011c-ae85-49da-9347-437a76635787
split_ratio3a = F_func[7](tsplit3)/F_func[3](tsplit3)

# ╔═╡ d858d61a-48f5-4104-8c56-bbc0295746bd
A1 = split_ratio1a*split_ratio2a

# ╔═╡ 2ef5bc9f-5a88-4fbe-8fa7-c92e1c44947a
A2 = split_ratio1a*split_ratio2b

# ╔═╡ 19b830a0-faff-457c-bd35-c2991184f6c0
A3 = split_ratio1b*split_ratio3a

# ╔═╡ 31719c59-ab86-4b3a-81bf-1e2d2841a64f
A4 = split_ratio1b*split_ratio3b

# ╔═╡ 27c4aeb2-7fe4-4003-b8ee-b9ee30775cd5
A1+A2+A3+A4

# ╔═╡ 7bb1a8d3-ebab-49ac-b1c3-de2206f441f9
plot1, t1, c1 = GasChromatographySimulator.plot_chromatogram(sim_res[2][2][end], (0, 1800.0))

# ╔═╡ fe919b98-bbe1-47ea-8ed8-d9d1acba6e40
plot2, t2, c2 = GasChromatographySimulator.plot_chromatogram(sim_res[2][1][end], (0, 1800.0))

# ╔═╡ 5623cb51-ceae-49ab-974b-c2a530d55487
plot3, t3, c3 = GasChromatographySimulator.plot_chromatogram(sim_res[2][4][end], (0, 1800.0))

# ╔═╡ 01bae1cb-a3f4-4a2c-a5d3-8201a776c562
plot4, t4, c4 = GasChromatographySimulator.plot_chromatogram(sim_res[2][3][end], (0, 1800.0))

# ╔═╡ afb77882-db2f-48a9-85c6-889920e6ab78
begin
	p_chrom = Plots.plot(xlabel="time in s", label="detector signal")
	Plots.plot!(p_chrom, t1, A1.*c1.+ 0.0, label="1", c=2, lw=2)
	Plots.plot!(p_chrom, t1, A2.*c2.+ 0.2, label="2", c=3, lw=2)
	Plots.plot!(p_chrom, t1, A3.*c3.+ 0.4, label="3", c=4, lw=2)
	Plots.plot!(p_chrom, t1, A4.*c4.+ 0.6, label="4", c=5, lw=2)
	Plots.plot!(p_chrom, t1, (A1.*c1 .+ A2.*c2 .+ A3.*c3 .+ A4.*c4).+ 0.8, label="detector signal", c=1, lw=2)
	Plots.plot!(p_chrom, xlims=(0.0, 520.0))
	p_chrom#_
end

# ╔═╡ a995c54c-6386-4b1f-a3e3-c6ab84f4544d
paths

# ╔═╡ f2f0188d-bc26-49da-9d40-03238ed0c3e6
#savefig(p_chrom, "chrom_SNAT_pconst.svg")

# ╔═╡ edb3cbb3-040f-45fe-85be-5ec90ab93376
md"""
## Slicing the chromatogram
"""

# ╔═╡ 17cbcc86-9360-4702-b019-7b8071d320fc
c = A1.*c1 .+ A2.*c2 .+ A3.*c3 .+ A4.*c4

# ╔═╡ f4ea54b0-b95f-44f3-a76e-067eae2727ca
Int.(fld.(collect(t1), 30.0))

# ╔═╡ 6a33ee2c-8ed9-4e10-a878-75c742a9332e
begin # slicing the 1D chrom to 2D
	PM = 25.0
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
end

# ╔═╡ 8231226c-7ab7-49aa-89ff-00394975049f
n

# ╔═╡ db7beb42-7f7d-49e7-9944-c81f606c31a9
t_D2

# ╔═╡ 20d5756a-9663-4ddc-a243-a6cae54481b6
begin
	p_slices = Plots.plot(t_D2[1], slices[1])
	for i=2:length(unique(n))
		Plots.plot!(p_slices, t_D2[i], slices[i].+i*0.1)
	end
	p_slices
end

# ╔═╡ 11f8842a-1c32-4a52-8762-209f80e3cbe5
function slicing(t1, c, PM)
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

# ╔═╡ 57d2199c-a1f7-45f0-aaf2-c1a4492967d5
@bind select_PM PlutoUI.Slider(0.0:0.1:60.0, default=30.0)

# ╔═╡ 2dc1893c-2c92-4453-a0a7-c85138b49482
select_PM

# ╔═╡ c00417e0-a227-4420-b9c8-a7770a2bdc3b
begin
	slices_, t_D1_, t_D2_ = slicing(t1, c, select_PM)
	p_slices_ = Plots.plot(t_D2_[1], slices_[1])
	for i=2:length(slices_)
		Plots.plot!(p_slices_, t_D2_[i], slices_[i].+i*0.1, legend=false)
	end
	p_slices_
end

# ╔═╡ f4f3a3a1-7ec6-4b87-959a-ae6a7fbc3164
t_D2_

# ╔═╡ Cell order:
# ╠═7f5b6a3c-ed83-11ed-2af0-4b696e429dc4
# ╠═60ed17af-1494-4dd4-ab5a-9ed76ed82c11
# ╠═af37fc60-57a4-4af9-b78e-fe467f3e9164
# ╠═71e7d4bb-faa8-43da-9bab-41046ca38ade
# ╠═5548429b-da44-4abd-bb4d-00d9ab66649d
# ╠═07fb0cea-e5aa-43be-900e-c7ff91f6e639
# ╟─31190533-2581-42f9-b352-23ae1cd54f10
# ╠═deb6198d-6e0a-471b-9ab3-30bd21ca37a1
# ╟─0e3164b6-d583-477f-843f-f6d3d8a390d4
# ╠═4bd3f0ea-03dd-4dd7-a8d8-99f20bcdb27a
# ╠═a0314b36-7840-4fbe-8084-12d7c4fd5f9c
# ╠═9596cd2b-2e05-4247-a889-e64186fcca9a
# ╠═87fd20bd-2e93-401a-8978-6e4714bd334d
# ╠═15faeed7-066a-4f66-b64f-9cd915437ca0
# ╠═6e5cc5cd-7a17-4505-b2cb-caec7f4e39d1
# ╠═efb6cb49-811b-4a29-8eac-93c79a10b765
# ╠═2bb93f4d-3d7b-47ad-86e6-3a958a915a7e
# ╠═02d8a178-410c-42f5-8784-596df16218b1
# ╠═92e7ea3c-64a4-4282-bdb2-596986a87278
# ╠═16857ed2-4794-4346-94b7-114ad2a4a294
# ╠═7fb06668-5d56-4304-baea-fa8d32a36d82
# ╠═556e42d6-d0b0-4c0a-9044-f5462fd7de23
# ╠═d6e02463-f9df-4658-ba6c-865f2640c3b9
# ╠═803485e0-fd65-4794-98ae-b79fb13d9518
# ╠═72f68999-82e6-4660-b7eb-2114a16898af
# ╠═2b332378-2933-48b0-bd69-9c8687f13965
# ╠═8fcbf371-3ca5-41d8-8149-7af890a6b90b
# ╠═b04bd252-fea8-41ed-a249-cc564241a4ec
# ╠═2aa1f5f1-2f60-4eaf-ae82-949df1532fd9
# ╠═523c7c64-84a7-4612-99c0-c0fdd23591df
# ╠═55f06025-db1c-4b42-9da8-5cb2d95e0234
# ╠═714b3bdd-2b2b-44af-8672-f21245c57634
# ╠═9df9bdcf-7efb-4cf3-a379-b89a4e81f25d
# ╠═60680135-ada1-48ae-bbef-b587cb81c50c
# ╠═22e80b28-7670-4d80-a65b-942135d63a8a
# ╠═d828afa5-8335-434a-96be-041cbd53a5a8
# ╠═9eecf202-d3e7-405a-99d7-6a17a12d481a
# ╠═a1dcef3e-3d8f-418a-875a-4d05258d39a3
# ╠═bd81cafc-208f-4afb-94af-990c6697335b
# ╠═70bc6f4d-a4aa-47e7-92f0-671d8f230e6a
# ╠═54151540-4150-49c0-a260-d755d663462b
# ╠═b531f3cf-1d0a-453a-ad7b-16c113991ba2
# ╠═bd9955af-9d6e-4092-8a8e-8cb62b1e97a9
# ╠═94b05220-e3c2-4ccb-8193-64aa3b224915
# ╠═efda5295-ca67-4746-a84b-6de830e87e09
# ╠═399d3c27-4db3-4165-af03-4512df3262d7
# ╠═2d856b07-7739-4b9b-a7f2-a2dedf5e882e
# ╠═20b6a285-6449-4ec5-b0d3-a7c443978c0f
# ╠═82bb824c-a0c9-46a2-8aec-603fce1bf66d
# ╠═4864997a-b053-40bf-81e9-ff0a8b28867f
# ╠═80da63de-c461-4727-b57f-34856be6a71e
# ╠═5101c8d7-7e16-407a-b9cb-e726994034da
# ╠═e13060a8-9ab7-46b5-807d-93f577bf7ce3
# ╠═a07e2feb-5ba9-4d35-a776-a01d46ab84bb
# ╠═83978df2-41f5-48f1-a3a8-88a346730490
# ╠═685cd6f8-c3c6-4f53-b6c5-0814e633c17c
# ╠═33e4d5c0-42f8-4c3e-8465-fb41cc38336e
# ╠═4d1b2744-29da-439b-aa1a-e45138823e60
# ╠═414b08d5-8fa5-40d0-831c-0b63864f9cde
# ╠═360d72d3-4941-456d-9e75-e0a6d821acfa
# ╠═4122c986-ed3a-4303-a231-15094e4ebbc2
# ╠═2b2fb924-28ae-4e4b-98f9-5afb8144c206
# ╠═f2e722fb-2edf-4159-9f9e-a0c841b4adf1
# ╠═b9f39f74-85df-4b69-90fa-e033fa48114d
# ╠═5cd75ea5-8519-45da-bdf6-8eca2da64065
# ╠═85513069-ebee-4b1a-a2d2-b73d9dbdf848
# ╠═ae818c15-ea60-45f8-a65f-4e2f6ad1b119
# ╠═410a3252-5541-4cbb-8254-f6a5d89d46cd
# ╠═d6cb1be1-6df6-40f0-8374-398936cb17d3
# ╠═06d4d265-0353-46ee-ad7f-6fd4df88d11c
# ╠═7636143f-79cd-4f1c-b939-0b031e1fbbe6
# ╠═efc2bd49-554b-4148-9516-d2668ef7d462
# ╠═4a4e28a0-aeb1-4fc8-9113-711b7b4eb26d
# ╠═4657b632-4f49-4cd9-bfec-262805cb3c0b
# ╠═145758e9-240a-471f-b113-aa33f9c3960c
# ╠═c7d1051e-c65f-4a95-8cc4-a29f86039fa7
# ╠═2794959c-9b5a-4e34-ae3d-dc3c687c7747
# ╠═fa85bc90-401d-41c8-a19e-d68e12707b92
# ╠═63b4f098-c82d-4a0f-b83e-da92b5780d5c
# ╠═abcc4a79-1ec9-4301-b4a0-3d1c2c41d2fc
# ╠═36cb9e74-e510-4305-87a2-86fac95fe041
# ╠═a8d4b540-a9ce-44c0-94db-24931b68c9f3
# ╠═6dca011c-ae85-49da-9347-437a76635787
# ╠═d858d61a-48f5-4104-8c56-bbc0295746bd
# ╠═2ef5bc9f-5a88-4fbe-8fa7-c92e1c44947a
# ╠═19b830a0-faff-457c-bd35-c2991184f6c0
# ╠═31719c59-ab86-4b3a-81bf-1e2d2841a64f
# ╠═27c4aeb2-7fe4-4003-b8ee-b9ee30775cd5
# ╠═7bb1a8d3-ebab-49ac-b1c3-de2206f441f9
# ╠═fe919b98-bbe1-47ea-8ed8-d9d1acba6e40
# ╠═5623cb51-ceae-49ab-974b-c2a530d55487
# ╠═01bae1cb-a3f4-4a2c-a5d3-8201a776c562
# ╠═afb77882-db2f-48a9-85c6-889920e6ab78
# ╠═a995c54c-6386-4b1f-a3e3-c6ab84f4544d
# ╠═f2f0188d-bc26-49da-9d40-03238ed0c3e6
# ╠═edb3cbb3-040f-45fe-85be-5ec90ab93376
# ╠═17cbcc86-9360-4702-b019-7b8071d320fc
# ╠═f4ea54b0-b95f-44f3-a76e-067eae2727ca
# ╠═6a33ee2c-8ed9-4e10-a878-75c742a9332e
# ╠═8231226c-7ab7-49aa-89ff-00394975049f
# ╠═db7beb42-7f7d-49e7-9944-c81f606c31a9
# ╠═20d5756a-9663-4ddc-a243-a6cae54481b6
# ╠═11f8842a-1c32-4a52-8762-209f80e3cbe5
# ╠═57d2199c-a1f7-45f0-aaf2-c1a4492967d5
# ╠═2dc1893c-2c92-4453-a0a7-c85138b49482
# ╠═f4f3a3a1-7ec6-4b87-959a-ae6a7fbc3164
# ╠═c00417e0-a227-4420-b9c8-a7770a2bdc3b
