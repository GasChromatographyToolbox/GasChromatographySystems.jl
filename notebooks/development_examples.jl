### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ c774e732-ee3f-11ed-3a6e-35b02ecb1068
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

# ╔═╡ e80f4b40-79af-4a8b-a720-0f6cd0262cf7
md"""
# Example GC-Systems as defined functions/modules
"""

# ╔═╡ 0dd386aa-1673-4eb9-965a-adeed1c7f4f8
default_TP = GasChromatographySystems.TemperatureProgram([0.0, 1800.0], [40.0, 340.0])

# ╔═╡ 1eacdd92-1ca2-419e-94a1-f0619198c6c1
db = filter([:CAS] => x -> !ismissing(x), DataFrame(CSV.File("/Users/janleppert/Documents/GitHub/RetentionData/Databases/GCSim_database_nonflag.csv", header=1, silencewarnings=true)))

# ╔═╡ d02056cf-1116-4396-9432-caa5011e67da
md"""
## Series of n columns
"""

# ╔═╡ 29921324-acf5-471d-96d0-9887646b9387
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

# ╔═╡ 19fefff3-29d4-4946-b443-c20487b9c742
example_SeriesSystem() = SeriesSystem([10.0, 5.0, 2.0, 1.0], [0.53, 0.32, 0.25, 0.1], [0.53, 0.32, 0.25, 0.1], ["Rxi17SilMS", "Rxi17SilMS", "Rxi17SilMS", "Rxi17SilMS"], [default_TP, default_TP, default_TP, default_TP], NaN, 300.0, 0.0)

# ╔═╡ 0d3776fe-0fc8-4262-9d11-e2b9fed1ad1e
SeriesSystem([30.0, 20.0, 10.0, 1.0], [0.53, 0.32, 0.25, 0.1], [0.53, 0.32, 0.25, 0.1], ["Rxi17SilMS", "Rxi17SilMS", "Rxi17SilMS", "Rxi17SilMS"], [default_TP, default_TP, default_TP, default_TP], NaN, 300.0, 0.0)

# ╔═╡ 28f5c420-0f0d-471f-aa5a-f7c4b42a9322
md"""
### Graph plots
"""

# ╔═╡ cb23a0dd-35f6-4a11-8ba0-67fd8257ef89
GasChromatographySystems.plot_graph(SeriesSystem([10.0, 5.0, 2.0, 1.0], [0.53, 0.32, 0.25, 0.1], [0.53, 0.32, 0.25, 0.1], ["Rxi17SilMS", "Rxi17SilMS", "Rxi17SilMS", "Rxi17SilMS"], [default_TP, default_TP, default_TP, default_TP], NaN, 300.0, 0.0))

# ╔═╡ 9a26363e-9fd7-4cae-8b8b-8104e861dd3d
GasChromatographySystems.plot_graph_with_flow(example_SeriesSystem(),1800)

# ╔═╡ 71fcaf11-b474-424d-8076-d7882183588a
md"""
### Simulation
"""

# ╔═╡ 06c64864-0d52-4317-a1b3-414cba3b3e9b
sys_series = example_SeriesSystem()

# ╔═╡ c0b83ac2-27dd-4119-84e2-8bd3b930dbf8
solutes_series = GasChromatographySystems.common_solutes(db, sys_series).Name[1:10]

# ╔═╡ 4d922d8f-8262-49de-ab5b-11d2fc0e14af
par_series = GasChromatographySystems.graph_to_parameters(sys_series, db, solutes_series)

# ╔═╡ a4de2a51-8fd7-4d9c-aef0-c300ed8952b0
paths_series = GasChromatographySystems.all_paths(sys_series.g, 1)[2]

# ╔═╡ 0fb1e8c2-2c84-430b-ad00-b8b8a2134310
sim_series = GasChromatographySystems.simulate_along_paths(sys_series, paths_series, par_series)

# ╔═╡ d54f6427-d941-4dab-bca2-61551868c396
p_series = GasChromatographySimulator.plot_chromatogram(sim_series[2][1][end], (0, 2000.0))[1]

# ╔═╡ ff0014cb-cea0-4a41-9028-d28150f45ffd
md"""
## Split
"""

# ╔═╡ e2021079-b4ba-46b1-99d7-8dfebe9eeb34
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

# ╔═╡ 0120987d-757a-404f-a00a-6c6bd26ac115
example_SplitSystem() = SplitSystem([10.0, 1.0, 5.0], [0.25, 0.1, 0.25], [0.25, 0.0, 0.0], ["Rxi17SilMS", "", ""], [default_TP, 300.0, 300.0], [1.0, NaN, NaN], NaN, 0.0, 101.3)

# ╔═╡ 9069f693-8a0d-4593-a87b-6be4c532560f
SplitSystem([10.0, 30.0, 5.0], [0.25, 0.25, 0.25], [0.25, 0.0, 0.0], ["Rxi17SilMS", "", ""], [300.0, 300.0, 300.0], [NaN, NaN, NaN], 300.0, 0.0, 101.3)

# ╔═╡ 66ed5c2a-d5f6-452b-bb69-5dd138ad388b
SplitSystem([10.0, 3.0, 5.0], [0.25, 0.1, 0.25], [0.25, 0.0, 0.0], ["Rxi17SilMS", "", ""], [default_TP, 300.0, 300.0], [1.0, NaN, NaN], NaN, 0.0, 101.3)

# ╔═╡ 099aa77a-39fd-484c-a8a9-654b6033a8b7
md"""
### Graph plots
"""

# ╔═╡ 0ca9155f-16e3-4a56-a128-c31c1ef43118
GasChromatographySystems.plot_graph(SplitSystem([10.0, 30.0, 5.0], [0.25, 0.25, 0.25], [0.25, 0.0, 0.0], ["Rxi17SilMS", "", ""], [300.0, 300.0, 300.0], [NaN, NaN, NaN], 300.0, 0.0, 101.3))

# ╔═╡ c9f120e0-94b9-43a5-b60d-3d84873b5c3d
GasChromatographySystems.plot_graph(SplitSystem([10.0, 3.0, 5.0], [0.25, 0.1, 0.25], [0.25, 0.0, 0.0], ["Rxi17SilMS", "", ""], [default_TP, 300.0, 300.0], [1.0, NaN, NaN], NaN, 0.0, 101.3))

# ╔═╡ 2ba6124a-0e10-4340-8023-5d5e7dec1899
GasChromatographySystems.plot_graph_with_flow(example_SplitSystem(),0)

# ╔═╡ 182a76d2-d2fe-47d0-9824-704b15103cfe
md"""
### Simulation
"""

# ╔═╡ f4af960a-6c7f-4348-8faf-e5b205484057
sys_split = example_SplitSystem()

# ╔═╡ cca1fcbc-2988-45cf-b6d3-4076a0aaeda1
ne(sys_split.g)

# ╔═╡ 0820d89c-70da-4e2b-bfc9-1099e277bc1d
solutes_split = GasChromatographySystems.common_solutes(db, sys_split).Name[1:10]

# ╔═╡ ee8ba306-eac0-4ff2-b509-6411e5d138d8
par_split = GasChromatographySystems.graph_to_parameters(sys_split, db, solutes_split)

# ╔═╡ d2936bc3-08b1-4fb9-8eb5-38d63d0c02a9
paths_split = GasChromatographySystems.all_paths(sys_series.g, 2)[2]

# ╔═╡ 3a48baff-0cec-495d-8888-41161c4b2d7e
sim_split = GasChromatographySystems.simulate_along_paths(sys_split, paths_split, par_split)

# ╔═╡ e63232a0-1ed0-4a48-8069-069d2b737526
p_split = GasChromatographySimulator.plot_chromatogram(sim_split[2][1][end], (0, 2000.0))[1]

# ╔═╡ 7ca5d83f-7f9d-40ee-b5dd-f2a227f8fa27
md"""
## Simplified GCxGC thermal modulator
"""

# ╔═╡ 5a9b2a5c-0ca1-4d95-989f-7ad43390c84e
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

# ╔═╡ 82122a6e-8504-44da-80c2-6002106c0fa3
example_GCxGC_TM_simp() = GCxGC_TM_simp(30.0, 0.25, 0.25, "Rxi17SilMS", default_TP, 2.0, 0.25, 0.25, "Wax", default_TP, NaN, 200.0, 0.0)

# ╔═╡ e2d02c1c-4659-44ad-b7ce-47acdba2454f
test = example_GCxGC_TM_simp()

# ╔═╡ 4a8e2f78-80f1-4fc1-863e-6d66d1d19b30
md"""
### Graph plots
"""

# ╔═╡ 286da99e-7610-483b-8fd1-434f59a16de4
GasChromatographySystems.plot_graph(GCxGC_TM_simp(30.0, 0.25, 0.25, "Rxi17SilMS", default_TP, 2.0, 0.25, 0.25, "Wax", default_TP, NaN, 300.0, 0.0))

# ╔═╡ 00fae38e-572a-4ebe-992a-a43c0a7ab1b4
GasChromatographySystems.plot_graph_with_flow(example_GCxGC_TM_simp(),0)

# ╔═╡ 3e599175-5caf-44ba-be79-89957568bcd6
md"""
### Simulation
"""

# ╔═╡ 076f8fd5-3c15-479b-90a0-bf5e8693845a
sys_GCxGC_TM = example_GCxGC_TM_simp()

# ╔═╡ 4bc17bb0-ba03-48ee-81dd-13fc52bc28bb
solutes_GCxGC_TM = GasChromatographySystems.common_solutes(db, sys_GCxGC_TM).Name#[1:15]

# ╔═╡ fd59e29d-4648-4310-b4d9-a0a66ac28fe3
par_GCxGC_TM = GasChromatographySystems.graph_to_parameters(sys_GCxGC_TM, db, solutes_GCxGC_TM)

# ╔═╡ 7002a7ec-447d-4937-9e17-d6adbfa91b28
paths_GCxGC_TM = GasChromatographySystems.all_paths(sys_GCxGC_TM.g, 2)[2]

# ╔═╡ 92af9dc0-c997-4819-8617-532d81f6606f
sim_GCxGC_TM = GasChromatographySystems.simulate_along_paths(sys_GCxGC_TM, paths_GCxGC_TM, par_GCxGC_TM; refocus=trues(ne(sys_GCxGC_TM.g)), τ₀_focus=0.1.*ones(length(par_GCxGC_TM[1].sub)))

# ╔═╡ 6a3adb74-fcb6-4ec2-9a91-4d595cc29050
pl_GCxGC_TM = GasChromatographySystems.peaklist_GCxGC(sim_GCxGC_TM[2][1][1], sim_GCxGC_TM[2][1][2])

# ╔═╡ 124ab4c6-7a1d-463b-98e9-8ebfe45dd45e
p_GCxGC_TM = GasChromatographySystems.plot_GCxGC(pl_GCxGC_TM, sys_GCxGC_TM)

# ╔═╡ 5ab9c8d2-3ff1-47c7-82ff-a8d6ebe82437
md"""
## Simplified GCxGC flow modulator
"""

# ╔═╡ f06ec257-ec20-4f90-a0a2-488ac648e491
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

# ╔═╡ 2cd02dce-9e24-4342-ba4f-f1b583e7d407
example_GCxGC_FM_simp() = GCxGC_FM_simp(30.0, 0.25, 0.25, "Rxi17SilMS", default_TP, 1.0, 2.0, 0.25, 0.25, "Wax", default_TP, 2.0, NaN, 0.0)

# ╔═╡ 13e2f05a-baa8-48b2-ab51-6d84974c6a8b
md"""
### Graph plots
"""

# ╔═╡ b1297129-9c2c-4fc1-bd95-193ac035919f
GasChromatographySystems.plot_graph(GCxGC_FM_simp(30.0, 0.25, 0.25, "Rxi17SilMS", default_TP, 1.0, 2.0, 0.25, 0.25, "Wax", default_TP, 2.0, NaN, 0.0))

# ╔═╡ 988315f5-bd90-4658-a1d2-613d3515e9bf
GasChromatographySystems.plot_graph_with_flow(example_GCxGC_FM_simp(),0)

# ╔═╡ 442f58e6-44ac-43a2-8086-af19084ab617
md"""
### Simulation
"""

# ╔═╡ 3789a920-3bd7-4b3f-a020-a919e77c5fd5
sys_GCxGC_FM = example_GCxGC_FM_simp()

# ╔═╡ eea80749-c264-4110-896d-b458dc4555ea
solutes_GCxGC_FM = GasChromatographySystems.common_solutes(db, sys_GCxGC_FM).Name

# ╔═╡ 13b4ead0-754a-46e3-ad02-d81d53702921
par_GCxGC_FM = GasChromatographySystems.graph_to_parameters(sys_GCxGC_FM, db, solutes_GCxGC_FM)

# ╔═╡ 8bbce5dc-a0ba-497d-a25f-977fbd70c8ec
paths_GCxGC_FM = GasChromatographySystems.all_paths(sys_GCxGC_FM.g, 5)[2]

# ╔═╡ c2ac38e3-f0a1-4c2e-b416-f189599d3aa9
sim_GCxGC_FM = GasChromatographySystems.simulate_along_paths(sys_GCxGC_FM, paths_GCxGC_FM, par_GCxGC_FM; refocus=trues(ne(sys_GCxGC_FM.g)), τ₀_focus=0.1.*ones(length(par_GCxGC_FM[1].sub)))

# ╔═╡ 01b34209-dc48-43d9-8424-94d630e9bead
pl_GCxGC_FM = GasChromatographySystems.peaklist_GCxGC(sim_GCxGC_FM[2][1][1], sim_GCxGC_FM[2][1][2])

# ╔═╡ bf442d42-303d-4367-9eb5-16dce4422182
p_GCxGC_FM = GasChromatographySystems.plot_GCxGC(pl_GCxGC_FM, sys_GCxGC_FM)

# ╔═╡ c7889497-f66c-4611-ae92-1da339c92df2
p_GCxGC_FM

# ╔═╡ 96c540c6-b868-4056-931c-5ec8fe04bc20
md"""
# End
"""

# ╔═╡ Cell order:
# ╠═c774e732-ee3f-11ed-3a6e-35b02ecb1068
# ╟─e80f4b40-79af-4a8b-a720-0f6cd0262cf7
# ╠═0dd386aa-1673-4eb9-965a-adeed1c7f4f8
# ╠═1eacdd92-1ca2-419e-94a1-f0619198c6c1
# ╠═d02056cf-1116-4396-9432-caa5011e67da
# ╠═29921324-acf5-471d-96d0-9887646b9387
# ╠═19fefff3-29d4-4946-b443-c20487b9c742
# ╠═0d3776fe-0fc8-4262-9d11-e2b9fed1ad1e
# ╠═28f5c420-0f0d-471f-aa5a-f7c4b42a9322
# ╠═cb23a0dd-35f6-4a11-8ba0-67fd8257ef89
# ╠═9a26363e-9fd7-4cae-8b8b-8104e861dd3d
# ╠═71fcaf11-b474-424d-8076-d7882183588a
# ╠═06c64864-0d52-4317-a1b3-414cba3b3e9b
# ╠═c0b83ac2-27dd-4119-84e2-8bd3b930dbf8
# ╠═4d922d8f-8262-49de-ab5b-11d2fc0e14af
# ╠═a4de2a51-8fd7-4d9c-aef0-c300ed8952b0
# ╠═0fb1e8c2-2c84-430b-ad00-b8b8a2134310
# ╠═d54f6427-d941-4dab-bca2-61551868c396
# ╠═ff0014cb-cea0-4a41-9028-d28150f45ffd
# ╠═e2021079-b4ba-46b1-99d7-8dfebe9eeb34
# ╠═0120987d-757a-404f-a00a-6c6bd26ac115
# ╠═9069f693-8a0d-4593-a87b-6be4c532560f
# ╠═66ed5c2a-d5f6-452b-bb69-5dd138ad388b
# ╠═099aa77a-39fd-484c-a8a9-654b6033a8b7
# ╠═0ca9155f-16e3-4a56-a128-c31c1ef43118
# ╠═c9f120e0-94b9-43a5-b60d-3d84873b5c3d
# ╠═2ba6124a-0e10-4340-8023-5d5e7dec1899
# ╠═182a76d2-d2fe-47d0-9824-704b15103cfe
# ╠═f4af960a-6c7f-4348-8faf-e5b205484057
# ╠═cca1fcbc-2988-45cf-b6d3-4076a0aaeda1
# ╠═0820d89c-70da-4e2b-bfc9-1099e277bc1d
# ╠═ee8ba306-eac0-4ff2-b509-6411e5d138d8
# ╠═d2936bc3-08b1-4fb9-8eb5-38d63d0c02a9
# ╠═3a48baff-0cec-495d-8888-41161c4b2d7e
# ╠═e63232a0-1ed0-4a48-8069-069d2b737526
# ╠═7ca5d83f-7f9d-40ee-b5dd-f2a227f8fa27
# ╠═5a9b2a5c-0ca1-4d95-989f-7ad43390c84e
# ╠═82122a6e-8504-44da-80c2-6002106c0fa3
# ╠═e2d02c1c-4659-44ad-b7ce-47acdba2454f
# ╠═4a8e2f78-80f1-4fc1-863e-6d66d1d19b30
# ╠═286da99e-7610-483b-8fd1-434f59a16de4
# ╠═00fae38e-572a-4ebe-992a-a43c0a7ab1b4
# ╠═3e599175-5caf-44ba-be79-89957568bcd6
# ╠═076f8fd5-3c15-479b-90a0-bf5e8693845a
# ╠═4bc17bb0-ba03-48ee-81dd-13fc52bc28bb
# ╠═fd59e29d-4648-4310-b4d9-a0a66ac28fe3
# ╠═7002a7ec-447d-4937-9e17-d6adbfa91b28
# ╠═92af9dc0-c997-4819-8617-532d81f6606f
# ╠═6a3adb74-fcb6-4ec2-9a91-4d595cc29050
# ╠═124ab4c6-7a1d-463b-98e9-8ebfe45dd45e
# ╠═c7889497-f66c-4611-ae92-1da339c92df2
# ╠═5ab9c8d2-3ff1-47c7-82ff-a8d6ebe82437
# ╠═f06ec257-ec20-4f90-a0a2-488ac648e491
# ╠═2cd02dce-9e24-4342-ba4f-f1b583e7d407
# ╠═13e2f05a-baa8-48b2-ab51-6d84974c6a8b
# ╠═b1297129-9c2c-4fc1-bd95-193ac035919f
# ╠═988315f5-bd90-4658-a1d2-613d3515e9bf
# ╠═442f58e6-44ac-43a2-8086-af19084ab617
# ╠═3789a920-3bd7-4b3f-a020-a919e77c5fd5
# ╠═eea80749-c264-4110-896d-b458dc4555ea
# ╠═13b4ead0-754a-46e3-ad02-d81d53702921
# ╠═8bbce5dc-a0ba-497d-a25f-977fbd70c8ec
# ╠═c2ac38e3-f0a1-4c2e-b416-f189599d3aa9
# ╠═01b34209-dc48-43d9-8424-94d630e9bead
# ╠═bf442d42-303d-4367-9eb5-16dce4422182
# ╠═96c540c6-b868-4056-931c-5ec8fe04bc20
