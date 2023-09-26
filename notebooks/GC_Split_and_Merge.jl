### A Pluto.jl notebook ###
# v0.19.26

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

# ╔═╡ 2578fc28-6e8c-4570-b45d-503ca6511863
html"""
	<style>
	  main {
		max-width: 900px;
	  }
	</style>
	"""

# ╔═╡ 139bfe64-17ac-455d-84a2-09c9564156cd
md"""
# GC System with two parallel columns

**Using spliting and merging of the flow.** 
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

# ╔═╡ 7cc0cc2b-db2f-40aa-8f76-d99b1958e817
md"""
### System
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
	modules6[3] = GasChromatographySystems.ModuleColumn("3 GC2", 2.0, 0.1e-3, 0.1e-6, "SLB5ms", default_TP) 
	modules6[4] = GasChromatographySystems.ModuleColumn("4 GC1", 10.0, 0.25e-3, 0.25e-6, "Wax", default_TP)
	modules6[5] = GasChromatographySystems.ModuleColumn("5 TL3", 1.0, 0.1e-3, 0.0e-6, "", default_TP)
	
	modules6[6] = GasChromatographySystems.ModuleColumn("6 TL4", 1.5, 0.1e-3, 0.0e-6, "", 300.0)
	# system
	sys6 = GasChromatographySystems.update_system(GasChromatographySystems.System(g6, pp6, modules6, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ 56c482e8-078c-4440-9dbd-51386abf1b4b
p_graph = GasChromatographySystems.plot_graph(sys6; node_size=64, arrow_size=30, nlabels_fontsize=30, elabels_fontsize=26)

# ╔═╡ ba4669d8-56f3-4e9d-a76e-e2da5f252ce5
#save("graph.svg", p_graph)

# ╔═╡ bebb2543-6f16-4563-8050-152381850110
pwd()

# ╔═╡ 4b379226-1646-4b13-b64a-7b55b66442f5
md"""
## Pressure and flow calculation
"""

# ╔═╡ bfd7ed6a-f293-448f-b827-a914396f101b
GasChromatographySystems.flow_balance(sys6)

# ╔═╡ 3946d9e6-78e3-41df-b204-4fd780721242
GasChromatographySystems.solve_balance(sys6)

# ╔═╡ 0aef6b6a-c58d-4bc8-812a-e958973598b9
p_graph_flow = GasChromatographySystems.plot_graph_with_flow(sys6, 0, arrow_size=30, node_size=90, nlabels_fontsize=20, elabels_fontsize=20, elabels_distance=28)

# ╔═╡ 81861823-5a64-41a5-aad0-117649ff2fc1
#save("graph_flow.svg", p_graph_flow)

# ╔═╡ a8b54992-6dec-42f2-b29b-7fa4f2128809
begin
	gr()
	p_F_sys6 = GasChromatographySystems.plot_flow_over_time(sys6)
	Plots.plot!(p_F_sys6, legend=:right)
	p_F_sys6
end

# ╔═╡ e842ff4e-b999-4794-9b73-d6754bb01d55
#savefig(p_F_sys6, "Flow_loop.svg")

# ╔═╡ 306a6211-5709-43ae-9154-2efed03c62e0
begin
	gr()
	p_p_sys6 = GasChromatographySystems.plot_pressure_over_time(sys6)
	Plots.plot!(p_p_sys6, legend=:right)
	Plots.ylims!(4.9e5, 6.05e5)
	p_p_sys6
end

# ╔═╡ 2128c82c-e497-4be3-9164-9cd689c19e6e
#savefig(p_p_sys6, "pressure_loop.svg")

# ╔═╡ 975eed72-5f77-4e44-bf33-6be175638ebd
md"""
## Simulation
"""

# ╔═╡ c057b9d4-feb7-414b-a1e5-7318c959825a
md"""
### Selection of solutes
"""

# ╔═╡ 83dbaf24-4764-426a-8f10-f6a3acef90d7
db = filter([:CAS] => x -> !ismissing(x), DataFrame(CSV.File(joinpath(dirname(pwd()), "data", "database_Blumberg2017.csv"), header=1, silencewarnings=true)))

# ╔═╡ 4f17cc8c-0004-44bd-b05d-97ee493235e0
unique(GasChromatographySystems.all_stationary_phases(sys6))

# ╔═╡ 0e5c952e-d0e4-442e-96db-f0a847831acd
selected_solutes = GasChromatographySystems.common_solutes(db, sys6).Name

# ╔═╡ 46414588-5b80-4995-80d1-51f058f6360b
md"""
### Graph to Parameters
"""

# ╔═╡ d2176d2e-9858-4a6b-be9e-be92768a2dbe
par = GasChromatographySystems.graph_to_parameters(sys6, db, selected_solutes; dt=60.0)

# ╔═╡ f821a4c0-4312-4d04-80f5-19f3b87ead97
md"""
### Paths
"""

# ╔═╡ b6f2e31c-9203-4101-aa7d-9aedee1de6b4
paths = GasChromatographySystems.all_paths(sys6.g, 2)[2]

# ╔═╡ c44a7a01-59df-4663-9280-96d2ed7de6e3
md"""
### Simulation
"""

# ╔═╡ 87757bee-79f3-48ff-92b7-c9d3b5e94b2d
res = GasChromatographySystems.simulate_along_paths(sys6, paths, db, selected_solutes, par; τ₀=0.5.*ones(length(selected_solutes)))

# ╔═╡ cd05457c-acf4-4bc3-9d47-2a0efdb2dd89
md"""
### Split at p₂

Split at p₂ to the two paths. Time of split is the same for all solutes. Split in the ration of the flows through edges `2 => 3` (`i=2`) and `2 => 5` (`i=3`) at this time.
"""

# ╔═╡ c5e3cb5f-11f2-476b-9e1c-afad3d024f19
t_split = res[2][1][1].tR[1]

# ╔═╡ 5b35524d-5f94-402c-a508-2f6d17e3bf42
F_func = GasChromatographySystems.flow_functions(sys6)

# ╔═╡ c0ba2b5d-199a-453f-8ca3-c7122239c279
collect(edges(sys6.g))

# ╔═╡ ab47be67-a524-4584-8925-02ae19a54298
F_func[1](t_split)*60e6

# ╔═╡ 5f0c8474-8d19-4e52-855c-23edb5478d41
F_func[2](t_split)*60e6

# ╔═╡ 22494320-33a0-4f8b-b739-038d6661afe5
F_func[3](t_split)*60e6

# ╔═╡ 722a5e35-1384-428c-bb48-719d90b06c1a
split_ratio_1 = F_func[2](t_split)/F_func[1](t_split)

# ╔═╡ aa28b189-2d36-49ba-b088-deab285aa479
split_ratio_2 = F_func[3](t_split)/F_func[1](t_split)

# ╔═╡ fba129ea-101b-41e6-9d75-ef493b9fd0a0
split_ratio_1 + split_ratio_2

# ╔═╡ 45b018aa-1e39-415c-8808-70fc54096212
begin
	p_chrom = GasChromatographySimulator.plot_chromatogram(res[2][1][end], (0, 1000.0))[1]
	GasChromatographySimulator.plot_chromatogram!(p_chrom, res[2][2][end], (0, 1000.0); mirror=true)
	p_chrom
end

# ╔═╡ 62764dd0-a127-4dc9-8cc5-b374464e231f
p1, t1, c1 = GasChromatographySimulator.plot_chromatogram(res[2][1][end], (0, 1000.0))

# ╔═╡ d4297e3b-75f9-4dd5-a73f-a00416e04a00
p2, t2, c2 = GasChromatographySimulator.plot_chromatogram(res[2][2][end], (0, 1000.0))

# ╔═╡ c1805860-59bf-4606-88d2-6e1af72a884c
begin
	p_chrom_ = Plots.plot(t1, split_ratio_1.*c1 .+ split_ratio_2.*c2, xlabel="time in s", label="detector signal")
	Plots.plot!(p_chrom_, t1, split_ratio_1.*c1.+ 0.4, label="GC1")
	Plots.plot!(p_chrom_, t1, split_ratio_2.*c2.+ 0.8, label="GC2")
	p_chrom_
end

# ╔═╡ 386a0b7b-f673-441e-83b8-ba9c533dda1b
#savefig(p_chrom_, "chrom_split_and_merge.svg")

# ╔═╡ 0e9d452e-5267-497c-ab80-c45c8551247a
function chromatogram(t, tR, τR; split_ratio=ones(length(tR)))
	g(t,tR,τR) = 1/sqrt(2*π*τR^2)*exp(-(t-tR)^2/(2*τR^2))
	chromatograms = Array{Array{Float64,1}}(undef, length(tR))
	for j=1:length(tR)
		chromatograms[j] = g.(t, tR[j], τR[j]).*split_ratio[j]
	end
	return sum(chromatograms)
end

# ╔═╡ 5938bed2-b209-45a8-8629-c5dff5ede982
paths

# ╔═╡ e073bc3d-4d70-4a67-8e94-8e63da808502
res[2][1][3]

# ╔═╡ 040f7251-1704-43d4-93fc-431d815ff548
begin
	t_TL1 = 0:0.01:1000.0
	c_TL1 = chromatogram(t_TL1, res[2][1][1].tR, res[2][1][1].τR)
	chrom_TL1 = Plots.plot(t_TL1, c_TL1, xlabel="time in s", legend=false, title="Simulation TL1", size=(600,300))
end

# ╔═╡ 1542e218-0b84-4b20-bde1-09b9e397d391
#savefig(chrom_TL1, "chrom_split_and_merge_TL1.svg")

# ╔═╡ 742b43f9-bfa5-4264-9555-6fba848112a1
begin
	t_TL2 = 0:0.01:1000.0
	c_TL2 = chromatogram(t_TL2, res[2][1][2].tR, res[2][1][2].τR, split_ratio=split_ratio_1.*ones(length(res[2][1][2].tR)))
	chrom_TL2 = Plots.plot(t_TL2, c_TL2, xlabel="time in s", legend=false, title="Simulation TL2", size=(600,300))
end

# ╔═╡ 24c490d9-ef9e-40fd-ad11-a5c9af11d6ad
#savefig(chrom_TL2, "chrom_split_and_merge_TL2.svg")

# ╔═╡ 4439fce1-fa7f-483c-abbb-f3b4951002a2
begin
	t_GC1 = 0:0.01:1000.0
	c_GC1 = chromatogram(t_GC1, res[2][1][3].tR, res[2][1][3].τR, split_ratio=split_ratio_1.*ones(length(res[2][1][2].tR)))
	chrom_GC1 = Plots.plot(t_GC1, c_GC1, xlabel="time in s", legend=false, title="Simulation GC1", size=(600,300))
end

# ╔═╡ b831681c-b6c3-4f1d-b1bb-a62ed16e5fa5
#savefig(chrom_GC1, "chrom_split_and_merge_GC1.svg")

# ╔═╡ 78354cd4-32f4-4e75-b38f-b7157cb19438
begin
	t_TL3 = 0:0.01:1000.0
	c_TL3 = chromatogram(t_TL3, res[2][1][4].tR, res[2][1][4].τR, split_ratio=split_ratio_1.*ones(length(res[2][1][2].tR)))
	chrom_TL3 = Plots.plot(t_TL3, c_TL3, xlabel="time in s", legend=false, title="Simulation TL3", size=(600,300))
end

# ╔═╡ 5354f9f4-95fb-40ab-9ddc-11e5519127f8
#savefig(chrom_TL3, "chrom_split_and_merge_TL3.svg")

# ╔═╡ d856bcf7-cfdc-4d09-8760-1fb9a1c8141b
begin
	t_TL4 = 0:0.01:1000.0
	c_TL4 = chromatogram(t_TL4, res[2][1][5].tR, res[2][1][5].τR, split_ratio=split_ratio_1.*ones(length(res[2][1][2].tR)))
	chrom_TL4 = Plots.plot(t_TL4, c_TL4, xlabel="time in s", legend=false, title="Simulation TL4", size=(600,300))
end

# ╔═╡ 9e79b73b-d174-46ce-8d9b-b15c6ab22110
#savefig(chrom_TL4, "chrom_split_and_merge_TL4a.svg")

# ╔═╡ e1ba1423-b301-4434-a33c-f7645e680b66
begin
	t_GC2 = 0:0.01:1000.0
	c_GC2 = chromatogram(t_GC2, res[2][2][2].tR, res[2][2][2].τR, split_ratio=split_ratio_2.*ones(length(res[2][1][2].tR)))
	chrom_GC2 = Plots.plot(t_GC2, c_GC2, xlabel="time in s", legend=false, title="Simulation GC2", size=(600,300))
end

# ╔═╡ db87a000-b2eb-484e-babc-a5dc8323ef87
#savefig(chrom_GC2, "chrom_split_and_merge_GC2.svg")

# ╔═╡ d764969f-2a1e-4b25-918f-c0345afd0ec3
begin
	t_TL4b = 0:0.01:1000.0
	c_TL4b = chromatogram(t_TL4b, res[2][2][3].tR, res[2][2][3].τR, split_ratio=split_ratio_2.*ones(length(res[2][1][2].tR)))
	chrom_TL4b = Plots.plot(t_TL4b, c_TL4b, xlabel="time in s", legend=false, title="Simulation TL4", size=(600,300))
end

# ╔═╡ c936a262-542d-43e5-a796-38e1993c40f3
#savefig(chrom_TL4b, "chrom_split_and_merge_TL4b.svg")

# ╔═╡ 376458a4-16b7-40a3-9c79-eefa580ba182
begin
	chrom_all = Plots.plot(t_TL4, c_TL4.+c_TL4b, xlabel="time in s", legend=false, title="Chromatogram at detector", size=(600,300))
end

# ╔═╡ 0d02887b-cb38-4c8f-8034-0e577cc78de6
#savefig(chrom_all, "chrom_split_and_merge_all.svg")

# ╔═╡ Cell order:
# ╠═e9d64f0e-95b7-11ed-2fee-ab579cd02ebc
# ╠═2578fc28-6e8c-4570-b45d-503ca6511863
# ╠═139bfe64-17ac-455d-84a2-09c9564156cd
# ╠═f6c6443d-0e39-4e82-b2c4-49f2954e5b49
# ╠═719dc69d-2692-4d41-868e-14b664019d57
# ╠═f484559b-95d5-4155-aad6-fdeed3d238bc
# ╠═7cc0cc2b-db2f-40aa-8f76-d99b1958e817
# ╠═e9607708-458d-4a85-8ba6-0f404e35d441
# ╠═56c482e8-078c-4440-9dbd-51386abf1b4b
# ╠═ba4669d8-56f3-4e9d-a76e-e2da5f252ce5
# ╠═bebb2543-6f16-4563-8050-152381850110
# ╠═4b379226-1646-4b13-b64a-7b55b66442f5
# ╠═bfd7ed6a-f293-448f-b827-a914396f101b
# ╠═3946d9e6-78e3-41df-b204-4fd780721242
# ╠═0aef6b6a-c58d-4bc8-812a-e958973598b9
# ╠═81861823-5a64-41a5-aad0-117649ff2fc1
# ╠═a8b54992-6dec-42f2-b29b-7fa4f2128809
# ╠═e842ff4e-b999-4794-9b73-d6754bb01d55
# ╠═306a6211-5709-43ae-9154-2efed03c62e0
# ╠═2128c82c-e497-4be3-9164-9cd689c19e6e
# ╠═975eed72-5f77-4e44-bf33-6be175638ebd
# ╠═c057b9d4-feb7-414b-a1e5-7318c959825a
# ╠═83dbaf24-4764-426a-8f10-f6a3acef90d7
# ╠═4f17cc8c-0004-44bd-b05d-97ee493235e0
# ╠═0e5c952e-d0e4-442e-96db-f0a847831acd
# ╠═46414588-5b80-4995-80d1-51f058f6360b
# ╠═d2176d2e-9858-4a6b-be9e-be92768a2dbe
# ╠═f821a4c0-4312-4d04-80f5-19f3b87ead97
# ╠═b6f2e31c-9203-4101-aa7d-9aedee1de6b4
# ╠═c44a7a01-59df-4663-9280-96d2ed7de6e3
# ╠═87757bee-79f3-48ff-92b7-c9d3b5e94b2d
# ╠═cd05457c-acf4-4bc3-9d47-2a0efdb2dd89
# ╠═c5e3cb5f-11f2-476b-9e1c-afad3d024f19
# ╠═5b35524d-5f94-402c-a508-2f6d17e3bf42
# ╠═c0ba2b5d-199a-453f-8ca3-c7122239c279
# ╠═ab47be67-a524-4584-8925-02ae19a54298
# ╠═5f0c8474-8d19-4e52-855c-23edb5478d41
# ╠═22494320-33a0-4f8b-b739-038d6661afe5
# ╠═722a5e35-1384-428c-bb48-719d90b06c1a
# ╠═aa28b189-2d36-49ba-b088-deab285aa479
# ╠═fba129ea-101b-41e6-9d75-ef493b9fd0a0
# ╠═45b018aa-1e39-415c-8808-70fc54096212
# ╠═62764dd0-a127-4dc9-8cc5-b374464e231f
# ╠═d4297e3b-75f9-4dd5-a73f-a00416e04a00
# ╠═c1805860-59bf-4606-88d2-6e1af72a884c
# ╠═386a0b7b-f673-441e-83b8-ba9c533dda1b
# ╠═0e9d452e-5267-497c-ab80-c45c8551247a
# ╠═5938bed2-b209-45a8-8629-c5dff5ede982
# ╠═e073bc3d-4d70-4a67-8e94-8e63da808502
# ╠═040f7251-1704-43d4-93fc-431d815ff548
# ╠═1542e218-0b84-4b20-bde1-09b9e397d391
# ╠═742b43f9-bfa5-4264-9555-6fba848112a1
# ╠═24c490d9-ef9e-40fd-ad11-a5c9af11d6ad
# ╠═4439fce1-fa7f-483c-abbb-f3b4951002a2
# ╠═b831681c-b6c3-4f1d-b1bb-a62ed16e5fa5
# ╠═78354cd4-32f4-4e75-b38f-b7157cb19438
# ╠═5354f9f4-95fb-40ab-9ddc-11e5519127f8
# ╠═d856bcf7-cfdc-4d09-8760-1fb9a1c8141b
# ╠═9e79b73b-d174-46ce-8d9b-b15c6ab22110
# ╠═e1ba1423-b301-4434-a33c-f7645e680b66
# ╠═db87a000-b2eb-484e-babc-a5dc8323ef87
# ╠═d764969f-2a1e-4b25-918f-c0345afd0ec3
# ╠═c936a262-542d-43e5-a796-38e1993c40f3
# ╠═376458a4-16b7-40a3-9c79-eefa580ba182
# ╠═0d02887b-cb38-4c8f-8034-0e577cc78de6
