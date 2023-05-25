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
	using Interpolations
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
default_TP = GasChromatographySystems.TemperatureProgram([0.0, 1800.0], [40.0, 240.0])

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
	pp3[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 1800.0], [NaN, NaN]) # inlet 
	pp3[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 1800.0], [NaN, NaN]) # 
	pp3[3] = GasChromatographySystems.PressurePoint("p₃", [0.0, 1800.0], [eps(Float64), eps(Float64)]) # outlet 1 
	pp3[4] = GasChromatographySystems.PressurePoint("p₄", [0.0, 1800.0], [101300.0, 101300.0]) # outlet 2
	# modules
	modules3 = Array{GasChromatographySystems.AbstractModule}(undef, ne(g3))
	modules3[1] = GasChromatographySystems.ModuleColumn("GC column", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP, 1.0/60e6)
	modules3[2] = GasChromatographySystems.ModuleColumn("TL column 1", 0.5, 0.1e-3, 0.0e-6, "", 300.0)
	modules3[3] = GasChromatographySystems.ModuleColumn("TL column 2", 1.0, 0.15e-3, 0.15e-6, "SLB5ms", 300.0)
	# system
	sys3 = GasChromatographySystems.update_system(GasChromatographySystems.System(g3, pp3, modules3, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ ea5d611b-6b21-4185-94fc-ce3a6293c0dd
GasChromatographySystems.plot_graph(sys3; color=:green)

# ╔═╡ dd4b227d-6424-4a42-8f21-b3a405649138
GasChromatographySystems.flow_balance(sys3)

# ╔═╡ 31b01b3e-9f16-445c-bd22-1db767ef90dc
GasChromatographySystems.substitute_unknown_flows(sys3)

# ╔═╡ 2bf42ff7-3a22-4a54-bd8c-ba0dabb4b8ff
GasChromatographySystems.solve_balance(sys3)

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
	pp5[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 1800.0], [NaN, NaN]) # 
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
GasChromatographySystems.plot_graph(sys5; color=:orange)

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
	modules6[1] = GasChromatographySystems.ModuleColumn("1 TL1", 0.5, 0.25e-3, 0.25e-6, "Wax", 300.0)
	modules6[2] = GasChromatographySystems.ModuleColumn("2 TL2", 1.0, 0.25e-3, 0.25e-6, "Wax", 300.0)
	modules6[3] = GasChromatographySystems.ModuleColumn("5 GC2", 2.0, 0.1e-3, 0.1e-6, "SLB5ms", default_TP) 
	modules6[4] = GasChromatographySystems.ModuleColumn("3 GC1", 10.0, 0.25e-3, 0.25e-6, "Wax", default_TP)
	modules6[5] = GasChromatographySystems.ModuleColumn("4 TL3", 1.0, 0.1e-3, 0.1e-6, "Wax", default_TP)
	
	modules6[6] = GasChromatographySystems.ModuleColumn("6 TL4", 1.5, 0.1e-3, 0.1e-6, "", 300.0)
	# system
	sys6 = GasChromatographySystems.update_system(GasChromatographySystems.System(g6, pp6, modules6, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ 56c482e8-078c-4440-9dbd-51386abf1b4b
GasChromatographySystems.plot_graph(sys6; color=:orange)

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

# ╔═╡ a65c868c-fa26-4a7f-be94-2d86ee9518fc
md"""
## Selection of solutes
"""

# ╔═╡ b424f5f1-aca8-40e5-aa0d-32fc7c4d4ee6
db = DataFrame(CSV.File(joinpath(dirname(pwd()), "data", "database_Kcentric.csv"), header=1, silencewarnings=true))

# ╔═╡ 2eff2573-137d-4a85-9bc9-e43d081b8d77
GasChromatographySystems.common_solutes(db, sys7)

# ╔═╡ 79bc683a-8d51-4d81-a27c-e42663720917
md"""
## Simulation
"""

# ╔═╡ a220280c-ced8-4eb9-a880-0c03524d44f1
md"""
### 1. One Column
"""

# ╔═╡ 457df2b2-8bfa-4e04-b4f2-ae1ca8902a4e
paths1 = GasChromatographySystems.all_paths(sys1.g, 1)[2]

# ╔═╡ f106b5fd-e17a-4876-8562-e82670a1ab25
sim_res1 = GasChromatographySystems.simulate_along_paths(sys1, paths1, db, GasChromatographySystems.common_solutes(db, sys1).Name)

# ╔═╡ 2220bbf8-8fb4-4430-9ecc-e7f2b734ccf9
GasChromatographySystems.positive_flow(sys1)

# ╔═╡ bc88c459-5eb5-4ecc-ae1d-a28b6252c9f3
F_func = GasChromatographySystems.flow_functions(sys1)

# ╔═╡ bda71b84-91e3-4928-bb68-cb69de2382cc
F_func[1](0)

# ╔═╡ 10dfc2a9-e0d1-4086-bd78-1de721220d1a
GasChromatographySystems.common_solutes(db, sys1).Name[1:17:36]

# ╔═╡ cd920789-0a68-4fcc-b79c-820a0960e9be
md"""
### 2. Series of Columns
"""

# ╔═╡ 18235f5a-b0d7-4547-9754-9bb18f80270e
paths2 = GasChromatographySystems.all_paths(sys2.g, 1)[2]

# ╔═╡ b054af17-0532-4cdd-b9fd-84d4b084201a
# 2 Problems:
# 1. slow simulation -> propably because of p(t) function -> try to approximate with Interpolations.jl
# 2. change of order of solutes -> check the change of initial values, t0 and τ0 need to be in the correct order for the solutes

# ╔═╡ f8890221-36bb-48a3-a19a-fb5ed3b8df9a
function change_initial(par::GasChromatographySimulator.Parameters, pl)
	# copys the parameters `par` and changes the values of par.sub[i].τ₀ and par.sub[i].t₀ to pl.tR[], pl.τR[]
	newsub = Array{GasChromatographySimulator.Substance}(undef, length(par.sub))
	for i=1:length(par.sub)
		ii = findfirst(par.sub[i].name.==pl.Name)
		newsub[i] = GasChromatographySimulator.Substance(par.sub[i].name, par.sub[i].CAS, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀, par.sub[i].ann, par.sub[i].Cag, pl.tR[ii], pl.τR[ii])
	end
	newpar = GasChromatographySimulator.Parameters(par.col, par.prog, newsub, par.opt)
	return newpar
end

# ╔═╡ ca10ef72-67fc-4f47-bbda-6bf51a9f9991
sim_res2 = GasChromatographySystems.simulate_along_paths(sys2, paths2, db, GasChromatographySystems.common_solutes(db, sys2).Name[1:5])

# ╔═╡ 3bc00d5e-efba-4eb5-9ced-3bb93a06d786
sim_res2[2][1][1]

# ╔═╡ 832485b2-5030-4b8b-bc04-0cefba6e5a6d
sim_res2[4][2].sub

# ╔═╡ 507c0fbc-6bac-4c3b-bda0-b621a75b9fd5
change_initial(sim_res2[4][2], sim_res2[2][1][1]).sub

# ╔═╡ 82cd23de-2d46-4153-990e-2e335ab30b63
p_func = GasChromatographySystems.pressure_functions(sys2)

# ╔═╡ fd32cd81-e9d5-4a56-ab68-1dd9442bd2cf
p_func[1]

# ╔═╡ 58ded2ac-794e-430b-b28b-73da2146440d
p_itp = LinearInterpolation((0:18:1800, ), p_func[2].(0:18:1800), extrapolation_bc=Flat())

# ╔═╡ 800c0d24-b08a-4664-add6-b616b6a7b7a0
begin
	Plots.plot(0:18:1800, p_func[2].(0:18:1800))
	Plots.plot!(0:1:1800, p_itp.(0:1:1800))
end

# ╔═╡ 22135b8e-2de0-4c14-beac-8ae0d54925d4
p_func[3]

# ╔═╡ f35d3c3a-4d08-49a6-af71-94c313d59c14
com_times = sum(GasChromatographySystems.common_timesteps(sys1))

# ╔═╡ e847daa4-dae1-4f20-906d-45bc71f9a977
function interpolate_pressure_functions(sys; dt=1e-3)
	tend = sum(GasChromatographySystems.common_timesteps(sys))
	trange = 0:dt:tend
	p_func = GasChromatographySystems.pressure_functions(sys)
	p_itp = Array{Any}(undef, length(p_func))
	for i=1:length(p_func)
		p_itp[i] = LinearInterpolation((trange, ), p_func[i].(trange), extrapolation_bc=Flat())
	end
	return p_itp
end

# ╔═╡ 1a9ee707-000f-4d7e-be20-eba8014cba42
p_interp = interpolate_pressure_functions(sys1)[1]

# ╔═╡ fd5fa781-2564-4873-b866-0db495b7e471
@time p_interp(563.456)

# ╔═╡ 48dfacf2-58b0-4ba8-b553-5b6e8df25dc4
@time p_func[1](563.456)

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
# ╠═dd4b227d-6424-4a42-8f21-b3a405649138
# ╠═31b01b3e-9f16-445c-bd22-1db767ef90dc
# ╠═2bf42ff7-3a22-4a54-bd8c-ba0dabb4b8ff
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
# ╠═a65c868c-fa26-4a7f-be94-2d86ee9518fc
# ╠═b424f5f1-aca8-40e5-aa0d-32fc7c4d4ee6
# ╠═2eff2573-137d-4a85-9bc9-e43d081b8d77
# ╠═79bc683a-8d51-4d81-a27c-e42663720917
# ╠═a220280c-ced8-4eb9-a880-0c03524d44f1
# ╠═457df2b2-8bfa-4e04-b4f2-ae1ca8902a4e
# ╠═f106b5fd-e17a-4876-8562-e82670a1ab25
# ╠═2220bbf8-8fb4-4430-9ecc-e7f2b734ccf9
# ╠═bc88c459-5eb5-4ecc-ae1d-a28b6252c9f3
# ╠═bda71b84-91e3-4928-bb68-cb69de2382cc
# ╠═10dfc2a9-e0d1-4086-bd78-1de721220d1a
# ╠═cd920789-0a68-4fcc-b79c-820a0960e9be
# ╠═18235f5a-b0d7-4547-9754-9bb18f80270e
# ╠═b054af17-0532-4cdd-b9fd-84d4b084201a
# ╠═3bc00d5e-efba-4eb5-9ced-3bb93a06d786
# ╠═f8890221-36bb-48a3-a19a-fb5ed3b8df9a
# ╠═832485b2-5030-4b8b-bc04-0cefba6e5a6d
# ╠═507c0fbc-6bac-4c3b-bda0-b621a75b9fd5
# ╠═ca10ef72-67fc-4f47-bbda-6bf51a9f9991
# ╠═82cd23de-2d46-4153-990e-2e335ab30b63
# ╠═fd32cd81-e9d5-4a56-ab68-1dd9442bd2cf
# ╠═800c0d24-b08a-4664-add6-b616b6a7b7a0
# ╠═58ded2ac-794e-430b-b28b-73da2146440d
# ╠═22135b8e-2de0-4c14-beac-8ae0d54925d4
# ╠═f35d3c3a-4d08-49a6-af71-94c313d59c14
# ╠═e847daa4-dae1-4f20-906d-45bc71f9a977
# ╠═1a9ee707-000f-4d7e-be20-eba8014cba42
# ╠═fd5fa781-2564-4873-b866-0db495b7e471
# ╠═48dfacf2-58b0-4ba8-b553-5b6e8df25dc4
