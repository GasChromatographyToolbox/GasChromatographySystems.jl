### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ e9d64f0e-95b7-11ed-2fee-ab579cd02ebc
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

# ╔═╡ 139bfe64-17ac-455d-84a2-09c9564156cd
md"""
# Test of different Systems
"""

# ╔═╡ f6c6443d-0e39-4e82-b2c4-49f2954e5b49
md"""
## Definition of systems
"""

# ╔═╡ 719dc69d-2692-4d41-868e-14b664019d57
md"""
### Temperature programs
"""

# ╔═╡ f484559b-95d5-4155-aad6-fdeed3d238bc
default_TP = GasChromatographySystems.TemperatureProgram([0.0, 60.0, 1800.0, 60.0], [40.0, 40.0, 280.0, 280.0])

# ╔═╡ 2ff47732-7276-41e0-a412-4f2573e78fc3
begin
	a_gf = [[50.0, 50.0] [0.0, 0.0] [1.0, 1.0] [-5.0, -5.0]]
	gf(x) = GasChromatographySimulator.gradient(x, a_gf)
	gradient_TP = GasChromatographySystems.TemperatureProgram([0.0, 1800.0], [40.0, 240.0], gf, a_gf)
end

# ╔═╡ 513461ee-f3e9-4005-b3fa-915ac90adbef
md"""
### 1. One Column
"""

# ╔═╡ ec4674a5-5793-4b7b-98c5-da0a3d82108b
begin
	g1 = Graphs.SimpleDiGraph(2)
	Graphs.add_edge!(g1, 1, 2) # Inj -> GC column -> Det
	# pressure points:
	pp1 = Array{GasChromatographySystems.PressurePoint}(undef, Graphs.nv(g1))
	pp1[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 1800.0], [200000.0, 200000.0]) # inlet 
	pp1[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 1800.0], [101300.0, 101300.0]) # outlet
	# modules
	modules1 = Array{GasChromatographySystems.AbstractModule}(undef, Graphs.ne(g1))
	modules1[1] = GasChromatographySystems.ModuleColumn("GC column", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	# system
	sys1 = GasChromatographySystems.update_system(GasChromatographySystems.System(g1, pp1, modules1, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ 878efab2-1a84-4afa-94c5-d500857abc6d
GasChromatographySystems.plot_graph(sys1; color=:green)

# ╔═╡ 1fb81739-4d52-40ca-bb5e-7ad276990003
md"""
### 2. Series of Columns
"""

# ╔═╡ 42f628e1-2121-4cb1-9b0a-17c62f5346d0
begin
	g2 = SimpleDiGraph(3)
	add_edge!(g2, 1, 2) # Inj -> TL column -> GC column inlet
	add_edge!(g2, 2, 3) # GC column inlet -> GC column -> 2nd TL column inlet
	#add_edge!(g2, 3, 4) # 2nd TL column inlet -> TL column -> Det
	# pressure points
	pp2 = Array{GasChromatographySystems.PressurePoint}(undef, nv(g2))
	pp2[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 600.0, 1200.0], [200000.0, 200000.0, 200000.0]) # inlet 
	pp2[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 1800.0], [NaN, NaN]) 
	pp2[3] = GasChromatographySystems.PressurePoint("p₃", [0.0, 1800.0], [101300.0, 101300.0]) 
	#pp2[4] = GasChromatographySystems.PressurePoint("p₄", [0.0, 1800.0], [101300.0, 101300.0]) # outlet
	# modules
	modules2 = Array{GasChromatographySystems.AbstractModule}(undef, ne(g2))
	modules2[1] = GasChromatographySystems.ModuleColumn("GC column 1", 10.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	modules2[2] = GasChromatographySystems.ModuleColumn("GC column 2", 20.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	#modules2[3] = GasChromatographySystems.ModuleColumn("GC column 3", 10.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	# system
	sys2 = GasChromatographySystems.update_system(GasChromatographySystems.System(g2, pp2, modules2, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ 545fe3c2-9a5a-4cc2-936b-426786867029
GasChromatographySystems.plot_graph(sys2; color=:green)

# ╔═╡ 252bd934-7823-45c2-a5ec-0700561453d5
md"""
### 3. Split
"""

# ╔═╡ 0a0fff3d-1744-4a2a-977f-68a784e1ccc9
begin
	g3 = SimpleDiGraph(4)
	add_edge!(g3, 1, 2) # Inj -> GC column -> Split point
	add_edge!(g3, 2, 3) # Split point -> TL column -> Det 1
	add_edge!(g3, 2, 4) # Split point -> TL column -> Det 2
	# pressure points
	pp3 = Array{GasChromatographySystems.PressurePoint}(undef, nv(g3))
	pp3[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 1800.0], [200000.0, 250000.0]) # inlet 
	pp3[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 1800.0], [NaN, NaN]) # 
	pp3[3] = GasChromatographySystems.PressurePoint("p₃", [0.0, 1800.0], [eps(Float64), eps(Float64)]) # outlet 1 
	pp3[4] = GasChromatographySystems.PressurePoint("p₄", [0.0, 1800.0], [101300.0, 101300.0]) # outlet 2
	# modules
	modules3 = Array{GasChromatographySystems.AbstractModule}(undef, ne(g3))
	modules3[1] = GasChromatographySystems.ModuleColumn("GC column", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	modules3[2] = GasChromatographySystems.ModuleColumn("TL column 1", 0.5, 0.1e-3, 0.0e-6, "", 300.0)
	modules3[3] = GasChromatographySystems.ModuleColumn("TL column 2", 1.0, 0.15e-3, 0.15e-6, "SLB5ms", 300.0)
	# system
	sys3 = GasChromatographySystems.update_system(GasChromatographySystems.System(g3, pp3, modules3, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ ea5d611b-6b21-4185-94fc-ce3a6293c0dd
GasChromatographySystems.plot_graph(sys3; color=:green)

# ╔═╡ 6e185a7c-9bb3-4238-a0c1-efaf8a76df88
GasChromatographySystems.plot_graph(sys3.g, ["V₁ => V₂", "V₂ => V₃", "V₂ => V₄"], ["V₁", "V₂", "V₃", "V₄"]; color=:lightblue, node_size=70, nlabels_fontsize=26, elabels_fontsize=26)

# ╔═╡ 8c4137af-9be7-4224-a00a-d60d4dbd7fcf
md"""
### 4. Series of Columns with Temperature Gradient
"""

# ╔═╡ 9d776592-b23c-4859-b5ce-071a6b1cfbdd
begin
	g4 = SimpleDiGraph(4)
	add_edge!(g4, 1, 2) # TL column
	add_edge!(g4, 2, 3) # GC column
	add_edge!(g4, 3, 4) # 2nd GC column with grsdient program
	# pressure points
	pp4 = Array{GasChromatographySystems.PressurePoint}(undef, nv(g4))
	pp4[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 1800.0], [200000.0, 200000.0]) # inlet 
	pp4[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 1800.0], [NaN, NaN]) 
	pp4[3] = GasChromatographySystems.PressurePoint("p₃", [0.0, 1800.0], [NaN, NaN]) 
	pp4[4] = GasChromatographySystems.PressurePoint("p₄", [0.0, 1800.0], [0.0, 0.0]) # outlet
	# modules
	modules4 = Array{GasChromatographySystems.AbstractModule}(undef, ne(g4))
	modules4[1] = GasChromatographySystems.ModuleColumn("TL column", 0.2, 0.25e-3, 0.25e-6, "SLB5ms", 320.0)
	modules4[2] = GasChromatographySystems.ModuleColumn("GC column", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	modules4[3] = GasChromatographySystems.ModuleColumn("TGGC column", 1.0, 0.1e-3, 0.1e-6, "Wax", gradient_TP)
	# system
	sys4 = GasChromatographySystems.update_system(GasChromatographySystems.System(g4, pp4, modules4, GasChromatographySystems.Options()))
end

# ╔═╡ dbce409d-6c17-478a-8fd0-4f9444f85202
GasChromatographySystems.plot_graph(sys4; color=:orange)

# ╔═╡ 06f73bb4-3263-47f7-a8de-a4f56478ebfe
md"""
### 5. Combination of 2. and 3.
"""

# ╔═╡ 8655ab62-f1a0-4abb-8856-2eb6ca9736f9
begin
	# setup of the Graph
	g5 = SimpleDiGraph(9)
	add_edge!(g5, 1, 2) # Inj->GC1->pressure regulator 
	add_edge!(g5, 2, 3) # pressure regulator -> GC2 -> split
	add_edge!(g5, 3, 4) # Split -> TL1 -> Det1
	add_edge!(g5, 3, 5) # Split -> TL2 -> Modulator inlet
	add_edge!(g5, 5, 6) # Modulator inlet -> Modulator -> GC3 inlet
	add_edge!(g5, 6, 7) # GC3 inlet -> GC3 -> Split
	add_edge!(g5, 7, 8) # Split -> TL3 -> Det2
	add_edge!(g5, 7, 9) # Split -> TL4 -> Det3
	#edge_labels = ["GC₁", "GC₂", "TL₁", "TL₂", "Modulator", "GC₃", "TL₃", "TL₄"]
	#node_labels = ["p₁", "p₂", "p₃", "p₄", "p₅", "p₆", "p₇", "p₈", "p₉"]
	# pressure points
	pp5 = Array{GasChromatographySystems.PressurePoint}(undef, nv(g5))
	pp5[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 1800.0], [200000.0, 200000.0]) # inlet 
	pp5[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 1800.0], [150000.0, 150000.0]) # 
	pp5[3] = GasChromatographySystems.PressurePoint("p₃", [0.0, 1800.0], [NaN, NaN]) #  
	pp5[4] = GasChromatographySystems.PressurePoint("p₄", [0.0, 1800.0], [101300.0, 101300.0]) # outlet 
	pp5[5] = GasChromatographySystems.PressurePoint("p₅", [0.0, 1800.0], [NaN, NaN]) #
	pp5[6] = GasChromatographySystems.PressurePoint("p₆", [0.0, 1800.0], [NaN, NaN]) # 
	pp5[7] = GasChromatographySystems.PressurePoint("p₇", [0.0, 1800.0], [NaN, NaN]) #
	pp5[8] = GasChromatographySystems.PressurePoint("p₈", [0.0, 1800.0], [eps(Float64), eps(Float64)])#[0.0, 0.0]) # outlet 2
	pp5[9] = GasChromatographySystems.PressurePoint("p₉", [0.0, 1800.0], [101300.0, 101300.0]) # outlet 3
	# modules
	modules5 = Array{GasChromatographySystems.AbstractModule}(undef, ne(g5))
	modules5[1] = GasChromatographySystems.ModuleColumn("GC column 1", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	modules5[2] = GasChromatographySystems.ModuleColumn("GC column 2", 30.0, 0.25e-3, 0.25e-6, "SPB50", default_TP)
	modules5[3] = GasChromatographySystems.ModuleColumn("TL column 1", 3.0, 0.25e-3, 0.25e-6, "SPB50", 300.0)
	modules5[4] = GasChromatographySystems.ModuleColumn("TL column 2", 1.0, 0.25e-3, 0.25e-6, "SPB50", 300.0)
	modules5[5] = GasChromatographySystems.ModuleColumn("Modulator", 0.3, 0.1e-3, 0.1e-6, "Wax", 300.0) # replace later with ModuleModulator
	modules5[6] = GasChromatographySystems.ModuleColumn("GC column 3", 4.0, 0.1e-3, 0.1e-6, "Wax", default_TP)
	modules5[7] = GasChromatographySystems.ModuleColumn("TL column 3", 2.5, 0.1e-3, 0.1e-6, "", 300.0)
	modules5[8] = GasChromatographySystems.ModuleColumn("TL column 4", 1.5, 0.1e-3, 0.1e-6, "", 300.0)
	# system
	sys5 = GasChromatographySystems.update_system(GasChromatographySystems.System(g5, pp5, modules5, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ bc10e868-d48b-492b-ad65-170f38e4e340
GasChromatographySystems.plot_graph(sys5; lay=Stress())

# ╔═╡ 7cc0cc2b-db2f-40aa-8f76-d99b1958e817
md"""
### 6. Split and Merge
"""

# ╔═╡ e9607708-458d-4a85-8ba6-0f404e35d441
begin
	# setup of the Graph -> ATTENTION to the right order of the modules
	# look at `collect(edges(sys_loop.g))`
	g6 = SimpleDiGraph(6)
	add_edge!(g6, 1, 2) # Inj->TL1->Split 
	add_edge!(g6, 2, 3) # Split -> TL2
	add_edge!(g6, 2, 5) # GC2
	add_edge!(g6, 3, 4) # GC1
	add_edge!(g6, 4, 5) # TL3
	add_edge!(g6, 5, 6) # TL4 -> Det
	# pressure points
	pp6 = Array{GasChromatographySystems.PressurePoint}(undef, nv(g6))
	pp6[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 1800.0], [600000.0, 600000.0]) # inlet 
	pp6[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 1800.0], [NaN, NaN]) # 
	pp6[3] = GasChromatographySystems.PressurePoint("p₃", [0.0, 1800.0], [NaN, NaN]) #  
	pp6[4] = GasChromatographySystems.PressurePoint("p₄", [0.0, 1800.0], [NaN, NaN]) #
	pp6[5] = GasChromatographySystems.PressurePoint("p₅", [0.0, 1800.0], [NaN, NaN]) #
	pp6[6] = GasChromatographySystems.PressurePoint("p₆", [0.0, 1800.0], [101300.0, 101300.0]) # outlet
	# modules
	modules6 = Array{GasChromatographySystems.AbstractModule}(undef, ne(g6))
	modules6[1] = GasChromatographySystems.ModuleColumn("1 TL1", 0.5, 0.25e-3, 0.0e-6, "", default_TP)
	modules6[2] = GasChromatographySystems.ModuleColumn("2 TL2", 1.0, 0.25e-3, 0.0e-6, "", default_TP)
	modules6[3] = GasChromatographySystems.ModuleColumn("5 GC2", 2.0, 0.1e-3, 0.1e-6, "SLB5ms", default_TP) 
	modules6[4] = GasChromatographySystems.ModuleColumn("3 GC1", 10.0, 0.25e-3, 0.25e-6, "Wax", default_TP)
	modules6[5] = GasChromatographySystems.ModuleColumn("4 TL3", 1.0, 0.1e-3, 0.0e-6, "", default_TP)
	
	modules6[6] = GasChromatographySystems.ModuleColumn("6 TL4", 1.5, 0.1e-3, 0.0e-6, "", 300.0)
	# system
	sys6 = GasChromatographySystems.update_system(GasChromatographySystems.System(g6, pp6, modules6, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ 56c482e8-078c-4440-9dbd-51386abf1b4b
GasChromatographySystems.plot_graph(sys6; node_size=64, nlabels_fontsize=30, elabels_fontsize=26)

# ╔═╡ 208b87ea-35ab-4db0-9758-c7003d3326b8
md"""
### 7. Combination of all
"""

# ╔═╡ 8717ee7c-3f68-4c51-857a-a105ec9cd1bc
begin
	# setup of the Graph -> ATTENTION to the right order of the modules
	# look at `collect(edges(sys_loop.g))`
	g7 = SimpleDiGraph(7)
	add_edge!(g7, 1, 2)
	add_edge!(g7, 2, 3)
	add_edge!(g7, 2, 4)
	add_edge!(g7, 4, 5)
	add_edge!(g7, 4, 6)
	add_edge!(g7, 5, 6)
	add_edge!(g7, 6, 7)
	# pressure points
	pp7 = Array{GasChromatographySystems.PressurePoint}(undef, nv(g7))
	pp7[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 1800.0], [300000.0, 300000.0]) # inlet 
	pp7[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 1800.0], [NaN, NaN]) # 
	pp7[3] = GasChromatographySystems.PressurePoint("p₃", [0.0, 1800.0], [101300, 101300]) #  
	pp7[4] = GasChromatographySystems.PressurePoint("p₄", [0.0, 1800.0], [NaN, NaN]) #
	pp7[5] = GasChromatographySystems.PressurePoint("p₅", [0.0, 1800.0], [NaN, NaN]) #
	pp7[6] = GasChromatographySystems.PressurePoint("p₆", [0.0, 1800.0], [NaN, NaN]) # outlet
	pp7[7] = GasChromatographySystems.PressurePoint("p₇", [0.0, 1800.0], [eps(Float64), eps(Float64)])
	# modules
	modules7 = Array{GasChromatographySystems.AbstractModule}(undef, ne(g7))
	modules7[1] = GasChromatographySystems.ModuleColumn("1 TL1", 5.0, 0.25e-3, 0.25e-6, "Wax", 300.0)
	modules7[2] = GasChromatographySystems.ModuleColumn("2 TGGC", 1.0, 0.25e-3, 0.25e-6, "Wax", gradient_TP)
	modules7[3] = GasChromatographySystems.ModuleColumn("3 TL2", 0.1, 0.25e-3, 0.0, "", default_TP) 
	modules7[4] = GasChromatographySystems.ModuleColumn("4 GC", 10.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	modules7[5] = GasChromatographySystems.ModuleColumn("5 TL3", 0.1, 0.15e-3, 0.0, "", 300.0)
	modules7[6] = GasChromatographySystems.ModuleColumn("6 TL4", 0.1, 0.15e-3, 0.0, "", 300.0)
	modules7[7] = GasChromatographySystems.ModuleColumn("7 TL5", 0.2, 0.1e-3, 0.0, "", 300.0)
	# system
	sys7 = GasChromatographySystems.update_system(GasChromatographySystems.System(g7, pp7, modules7, GasChromatographySystems.Options()))
end

# ╔═╡ 895fd137-0f02-450b-9a72-bf0dc505b214
GasChromatographySystems.plot_graph(sys7; color=:pink)

# ╔═╡ 03190a64-09cf-4a9b-9b86-e2b50392481c
md"""
### 8. GCxGC
"""

# ╔═╡ 65bbadca-a950-4f42-9502-31f99332b0e4
260/5*60

# ╔═╡ c9941593-0311-473a-9c0c-f764d41ac7e7
GCxGC_TP = GasChromatographySystems.TemperatureProgram([0.0, 120.0, 3120.0, 300.0], [40.0, 40.0, 300.0, 300.0])

# ╔═╡ d571712f-981e-4d1e-95ce-b4dbf7de50c6
begin
	g8 = SimpleDiGraph(3)
	add_edge!(g8, 1, 2) # 1st-D GC
	add_edge!(g8, 2, 3) # 2nd-D GC
	#add_edge!(g2, 3, 4) # 2nd TL column inlet -> TL column -> Det
	# pressure points
	pp8 = Array{GasChromatographySystems.PressurePoint}(undef, nv(g8))
	pp8[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 120.0, 3120.0, 300.0], [159550.0, 159550.0, 265600.0, 265600.0]) # inlet 
	pp8[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 120.0, 3120.0, 300.0], [NaN, NaN, NaN, NaN]) 
	pp8[3] = GasChromatographySystems.PressurePoint("p₃", [0.0, 120.0, 3120.0, 300.0], [eps(Float64), eps(Float64), eps(Float64), eps(Float64)]) 
	#pp2[4] = GasChromatographySystems.PressurePoint("p₄", [0.0, 1800.0], [101300.0, 101300.0]) # outlet
	# modules
	modules8 = Array{GasChromatographySystems.AbstractModule}(undef, ne(g8))
	modules8[1] = GasChromatographySystems.ModuleColumn("GC column 1", 30.0, 0.25e-3, 0.25e-6, "Rxi5ms", GCxGC_TP)
	modules8[2] = GasChromatographySystems.ModuleColumn("GC column 2", 2.0, 0.25e-3, 0.25e-6, "Rxi17SilMS", GCxGC_TP)
	#modules2[3] = GasChromatographySystems.ModuleColumn("GC column 3", 10.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	# system
	sys8 = GasChromatographySystems.update_system(GasChromatographySystems.System(g8, pp8, modules8, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ 4b379226-1646-4b13-b64a-7b55b66442f5
md"""
## Pressure and flow calculation
"""

# ╔═╡ 7783216f-6a73-44aa-98f7-e0ca628ba145
function plot_flow_over_time(sys)
	#plotly()
	F_func = GasChromatographySystems.flow_functions(sys)
	com_timesteps = GasChromatographySystems.common_timesteps(sys)
	trange = 0:sum(com_timesteps)/100:sum(com_timesteps)
	p_flow = Plots.plot(xlabel="time in s", ylabel="flow in mL/min")
	for i=1:ne(sys.g)
		Plots.plot!(p_flow, trange, F_func[i].(trange)*60*1e6, label="F_$(sys.modules[i].name)")
	end
	return p_flow
end

# ╔═╡ 2263a0a6-c25a-4564-b47d-e097d3e60feb
function plot_pressure_over_time(sys)
	#plotly()
	p_func = GasChromatographySystems.pressure_functions(sys)
	com_timesteps = GasChromatographySystems.common_timesteps(sys)
	trange = 0:sum(com_timesteps)/100:sum(com_timesteps)
	p_pres = Plots.plot(xlabel="time in s", ylabel="pressure in Pa(a)", legend=:topleft)
	for i=1:nv(sys.g)
		Plots.plot!(p_pres, trange, p_func[i].(trange), label="$(sys.pressurepoints[i].name)")
	end
	return p_pres
end

# ╔═╡ 245b6182-1e0e-4ae2-9b40-1b62835edec3
md"""
### 1. One Column
"""

# ╔═╡ 63a6a5a0-a4d3-430e-af78-ba446c100db5
GasChromatographySystems.flow_balance(sys1)

# ╔═╡ 40c2806b-ed48-4948-8e8c-2e7c635b76f9
GasChromatographySystems.solve_balance(sys1)

# ╔═╡ 0ac6c03b-5be1-41ea-ad99-380a670eee62
GasChromatographySystems.plot_graph_with_flow(sys1, 0)

# ╔═╡ 4413a361-eb58-42f8-9c06-b49c2f0d7c04
plot_flow_over_time(sys1)

# ╔═╡ d7e335a0-7b81-40d0-86c1-77f139928939
plot_pressure_over_time(sys1)

# ╔═╡ 041406a3-799c-4ee2-b8f6-0a139d3ff310
md"""
### 2. Series of Columns
"""

# ╔═╡ daae9b7d-a0c3-4014-bf75-4234169c8b03
GasChromatographySystems.flow_balance(sys2)

# ╔═╡ 3aaf0983-3732-4e6b-900d-8614bed7f2b2
GasChromatographySystems.solve_balance(sys2)

# ╔═╡ ca7c13ae-6f9e-469a-85c5-de1031b56cb7
GasChromatographySystems.plot_graph_with_flow(sys2, 0)

# ╔═╡ cd0fca4d-4e20-4acc-a1d4-77277d97e2ed
plot_flow_over_time(sys2)

# ╔═╡ c8f88634-8846-4ef9-a720-103b3f788c80
plot_pressure_over_time(sys2)

# ╔═╡ ce562339-cb62-4967-ac79-38ef605167c0
md"""
### 3. Split
"""

# ╔═╡ 6c4e4dfd-3e02-4122-91ae-7b37d40113a9
GasChromatographySystems.flow_balance(sys3)

# ╔═╡ e231a9fa-3e9b-406b-886b-c8358c3b85a5
GasChromatographySystems.solve_balance(sys3)

# ╔═╡ 156e822b-75df-4841-ae8c-9b007d3373a2
GasChromatographySystems.plot_graph_with_flow(sys3, 0)

# ╔═╡ 11edd09e-6ce8-416b-ae82-e62e3c62eaf4
plot_flow_over_time(sys3)

# ╔═╡ c61ea89d-6d83-4caa-b30d-ebc09292bc26
p_sys3 = plot_pressure_over_time(sys3)

# ╔═╡ 4a3f20d2-9005-4273-8652-3e86a030226c
itp_p3 = GasChromatographySystems.interpolate_pressure_functions(sys3)

# ╔═╡ cfe00cb8-b64d-4499-b08f-d56711bf28fe
begin
	plotly()
	com_timesteps = GasChromatographySystems.common_timesteps(sys3)
	trange = 0:sum(com_timesteps)/100:sum(com_timesteps)
	Plots.scatter!(p_sys3, trange, itp_p3[2](trange))
	p_sys3
end

# ╔═╡ ab96fd01-3763-4a1e-b781-fe61c9e9f396
md"""
### 4. Series of Columns with Temperature Gradient
"""

# ╔═╡ b6296b1f-8458-470a-a47d-1835c3116ead
GasChromatographySystems.flow_balance(sys4)

# ╔═╡ fa9896da-74ec-46ba-8350-3690089293cb
GasChromatographySystems.solve_balance(sys4)

# ╔═╡ 45c1e624-2065-4e58-838b-bb4098133505
GasChromatographySystems.plot_graph_with_flow(sys4, 0)

# ╔═╡ 02d1f4db-c606-410d-aa48-8ba21cdf9303
plot_flow_over_time(sys4)

# ╔═╡ 7b5df074-d8a9-420f-b82c-dc2449ee997e
plot_pressure_over_time(sys4)

# ╔═╡ bea439b0-c428-426f-b7eb-8ab15d64f1e2
md"""
### 5. Combination of 2. and 3.
"""

# ╔═╡ 9359fc04-ea17-4551-9923-52c9d7509f2c
GasChromatographySystems.solve_balance(sys5)

# ╔═╡ 459ddf82-89d5-47a0-a535-699c015202bf
GasChromatographySystems.plot_graph_with_flow(sys5, 0, lay = Stress(), arrow_shift=0.6, nlabels_fontsize=16, elabels_fontsize=16)

# ╔═╡ e0fad228-0ba2-427e-8702-2d4640456e53
plot_flow_over_time(sys5)

# ╔═╡ e461d166-a7da-40fc-a07d-2c2d1165f03f
plot_pressure_over_time(sys5)

# ╔═╡ 655f92eb-4065-4cee-8c54-11b8e1db52b0
md"""
### 6. Split and Merge
"""

# ╔═╡ bfd7ed6a-f293-448f-b827-a914396f101b
GasChromatographySystems.flow_balance(sys6)

# ╔═╡ 3946d9e6-78e3-41df-b204-4fd780721242
GasChromatographySystems.solve_balance(sys6)

# ╔═╡ 0aef6b6a-c58d-4bc8-812a-e958973598b9
GasChromatographySystems.plot_graph_with_flow(sys6, 0)

# ╔═╡ a8b54992-6dec-42f2-b29b-7fa4f2128809
begin
	gr()
	p_F_sys6 = plot_flow_over_time(sys6)
	Plots.plot!(p_F_sys6, legend=:right)
	p_F_sys6
end

# ╔═╡ e842ff4e-b999-4794-9b73-d6754bb01d55
#savefig(p_F_sys6, "Flow_loop.svg")

# ╔═╡ 306a6211-5709-43ae-9154-2efed03c62e0
begin
	gr()
	p_p_sys6 = plot_pressure_over_time(sys6)
	Plots.plot!(p_p_sys6, legend=:right)
	Plots.ylims!(4.9e5, 6.05e5)
	p_p_sys6
end

# ╔═╡ 2128c82c-e497-4be3-9164-9cd689c19e6e
#savefig(p_p_sys6, "pressure_loop.svg")

# ╔═╡ 90c15fa2-5b79-41fc-804e-90eb70705dac
md"""
### 7. Combination of all
"""

# ╔═╡ 5252ef6b-1a09-4421-a18e-1fd351064392
GasChromatographySystems.flow_balance(sys7)

# ╔═╡ 12275b56-2711-4751-a298-d8dcf26747d3
GasChromatographySystems.solve_balance(sys7)

# ╔═╡ a441981d-3bfe-46be-ba13-3a69409e9e94
GasChromatographySystems.plot_graph_with_flow(sys7, 0)

# ╔═╡ 9972fb95-3509-43a6-8503-0105bcddd220
plot_flow_over_time(sys7)

# ╔═╡ 65e29f97-dc7d-4c1a-82be-e2e048052a0b
plot_pressure_over_time(sys7)

# ╔═╡ 975eed72-5f77-4e44-bf33-6be175638ebd
md"""
### 8. GCxGC
"""

# ╔═╡ 06b9ba3f-4dff-4fba-afce-923b95fae9c3
GasChromatographySystems.plot_graph_with_flow(sys8, 0)

# ╔═╡ af30bbcd-d6e4-4bd7-97f0-3243c514194a
plot_flow_over_time(sys8)

# ╔═╡ e9943b1c-2397-4cde-8c23-90f5f4a3f2c7
plot_pressure_over_time(sys8)

# ╔═╡ Cell order:
# ╠═e9d64f0e-95b7-11ed-2fee-ab579cd02ebc
# ╠═139bfe64-17ac-455d-84a2-09c9564156cd
# ╠═f6c6443d-0e39-4e82-b2c4-49f2954e5b49
# ╠═719dc69d-2692-4d41-868e-14b664019d57
# ╠═f484559b-95d5-4155-aad6-fdeed3d238bc
# ╠═2ff47732-7276-41e0-a412-4f2573e78fc3
# ╠═513461ee-f3e9-4005-b3fa-915ac90adbef
# ╠═ec4674a5-5793-4b7b-98c5-da0a3d82108b
# ╠═878efab2-1a84-4afa-94c5-d500857abc6d
# ╠═1fb81739-4d52-40ca-bb5e-7ad276990003
# ╠═42f628e1-2121-4cb1-9b0a-17c62f5346d0
# ╠═545fe3c2-9a5a-4cc2-936b-426786867029
# ╠═252bd934-7823-45c2-a5ec-0700561453d5
# ╠═0a0fff3d-1744-4a2a-977f-68a784e1ccc9
# ╠═ea5d611b-6b21-4185-94fc-ce3a6293c0dd
# ╠═6e185a7c-9bb3-4238-a0c1-efaf8a76df88
# ╠═8c4137af-9be7-4224-a00a-d60d4dbd7fcf
# ╠═9d776592-b23c-4859-b5ce-071a6b1cfbdd
# ╠═dbce409d-6c17-478a-8fd0-4f9444f85202
# ╠═06f73bb4-3263-47f7-a8de-a4f56478ebfe
# ╠═8655ab62-f1a0-4abb-8856-2eb6ca9736f9
# ╠═bc10e868-d48b-492b-ad65-170f38e4e340
# ╠═7cc0cc2b-db2f-40aa-8f76-d99b1958e817
# ╠═e9607708-458d-4a85-8ba6-0f404e35d441
# ╠═56c482e8-078c-4440-9dbd-51386abf1b4b
# ╠═208b87ea-35ab-4db0-9758-c7003d3326b8
# ╠═8717ee7c-3f68-4c51-857a-a105ec9cd1bc
# ╠═895fd137-0f02-450b-9a72-bf0dc505b214
# ╠═03190a64-09cf-4a9b-9b86-e2b50392481c
# ╠═65bbadca-a950-4f42-9502-31f99332b0e4
# ╠═c9941593-0311-473a-9c0c-f764d41ac7e7
# ╠═d571712f-981e-4d1e-95ce-b4dbf7de50c6
# ╠═4b379226-1646-4b13-b64a-7b55b66442f5
# ╠═7783216f-6a73-44aa-98f7-e0ca628ba145
# ╠═2263a0a6-c25a-4564-b47d-e097d3e60feb
# ╠═245b6182-1e0e-4ae2-9b40-1b62835edec3
# ╠═63a6a5a0-a4d3-430e-af78-ba446c100db5
# ╠═40c2806b-ed48-4948-8e8c-2e7c635b76f9
# ╠═0ac6c03b-5be1-41ea-ad99-380a670eee62
# ╠═4413a361-eb58-42f8-9c06-b49c2f0d7c04
# ╠═d7e335a0-7b81-40d0-86c1-77f139928939
# ╠═041406a3-799c-4ee2-b8f6-0a139d3ff310
# ╠═daae9b7d-a0c3-4014-bf75-4234169c8b03
# ╠═3aaf0983-3732-4e6b-900d-8614bed7f2b2
# ╠═ca7c13ae-6f9e-469a-85c5-de1031b56cb7
# ╠═cd0fca4d-4e20-4acc-a1d4-77277d97e2ed
# ╠═c8f88634-8846-4ef9-a720-103b3f788c80
# ╠═ce562339-cb62-4967-ac79-38ef605167c0
# ╠═6c4e4dfd-3e02-4122-91ae-7b37d40113a9
# ╠═e231a9fa-3e9b-406b-886b-c8358c3b85a5
# ╠═156e822b-75df-4841-ae8c-9b007d3373a2
# ╠═11edd09e-6ce8-416b-ae82-e62e3c62eaf4
# ╠═c61ea89d-6d83-4caa-b30d-ebc09292bc26
# ╠═4a3f20d2-9005-4273-8652-3e86a030226c
# ╠═cfe00cb8-b64d-4499-b08f-d56711bf28fe
# ╠═ab96fd01-3763-4a1e-b781-fe61c9e9f396
# ╠═b6296b1f-8458-470a-a47d-1835c3116ead
# ╠═fa9896da-74ec-46ba-8350-3690089293cb
# ╠═45c1e624-2065-4e58-838b-bb4098133505
# ╠═02d1f4db-c606-410d-aa48-8ba21cdf9303
# ╠═7b5df074-d8a9-420f-b82c-dc2449ee997e
# ╠═bea439b0-c428-426f-b7eb-8ab15d64f1e2
# ╠═9359fc04-ea17-4551-9923-52c9d7509f2c
# ╠═459ddf82-89d5-47a0-a535-699c015202bf
# ╠═e0fad228-0ba2-427e-8702-2d4640456e53
# ╠═e461d166-a7da-40fc-a07d-2c2d1165f03f
# ╠═655f92eb-4065-4cee-8c54-11b8e1db52b0
# ╠═bfd7ed6a-f293-448f-b827-a914396f101b
# ╠═3946d9e6-78e3-41df-b204-4fd780721242
# ╠═0aef6b6a-c58d-4bc8-812a-e958973598b9
# ╠═a8b54992-6dec-42f2-b29b-7fa4f2128809
# ╠═e842ff4e-b999-4794-9b73-d6754bb01d55
# ╠═306a6211-5709-43ae-9154-2efed03c62e0
# ╠═2128c82c-e497-4be3-9164-9cd689c19e6e
# ╠═90c15fa2-5b79-41fc-804e-90eb70705dac
# ╠═5252ef6b-1a09-4421-a18e-1fd351064392
# ╠═12275b56-2711-4751-a298-d8dcf26747d3
# ╠═a441981d-3bfe-46be-ba13-3a69409e9e94
# ╠═9972fb95-3509-43a6-8503-0105bcddd220
# ╠═65e29f97-dc7d-4c1a-82be-e2e048052a0b
# ╠═975eed72-5f77-4e44-bf33-6be175638ebd
# ╠═06b9ba3f-4dff-4fba-afce-923b95fae9c3
# ╠═af30bbcd-d6e4-4bd7-97f0-3243c514194a
# ╠═e9943b1c-2397-4cde-8c23-90f5f4a3f2c7
