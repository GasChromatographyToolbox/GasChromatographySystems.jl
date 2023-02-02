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
# Deans Switch
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
	g = SimpleDiGraph(13)
	add_edge!(g, 1, 2) # GC 
	add_edge!(g, 2, 3) # Deans
	add_edge!(g, 3, 4) # Deans
	add_edge!(g, 3, 7) # Deans
	add_edge!(g, 4, 5) # Deans
	add_edge!(g, 5, 6) # to Det1
	add_edge!(g, 7, 8) # Deans
	add_edge!(g, 8, 9) # to Det2
	add_edge!(g, 10, 4) # Deans
	add_edge!(g, 10, 12) # Deans, bypass
	add_edge!(g, 11, 10) # Deans, control in 1
	add_edge!(g, 12, 7) # Deans
	add_edge!(g, 13, 12) # Deans, control in 2
	# pressure points
	pp = Array{GasChromatographySystems.PressurePoint}(undef, nv(g))
	pp[1] = GasChromatographySystems.PressurePoint("p1", [0.0, 1800.0], [400000.0, 400000.0]) # inlet 
	pp[2] = GasChromatographySystems.PressurePoint("p2", [0.0, 1800.0], [NaN, NaN]) # 
	pp[3] = GasChromatographySystems.PressurePoint("p3", [0.0, 1800.0], [NaN, NaN]) # 
	pp[4] = GasChromatographySystems.PressurePoint("p4", [0.0, 1800.0], [NaN, NaN]) #
	pp[5] = GasChromatographySystems.PressurePoint("p5", [0.0, 1800.0], [NaN, NaN]) #
	pp[6] = GasChromatographySystems.PressurePoint("p6", [0.0, 1800.0], [101300.0, 101300.0]) # outlet Det1
	pp[7] = GasChromatographySystems.PressurePoint("p7", [0.0, 1800.0], [NaN, NaN]) #
	pp[8] = GasChromatographySystems.PressurePoint("p8", [0.0, 1800.0], [NaN, NaN]) #
	pp[9] = GasChromatographySystems.PressurePoint("p9", [0.0, 1800.0], [eps(Float64), eps(Float64)]) # outlet Det1
	pp[10] = GasChromatographySystems.PressurePoint("p10", [0.0, 1800.0], [NaN, NaN]) #
	pp[11] = GasChromatographySystems.PressurePoint("p11", [0.0, 1800.0], [150000.0, 150000.0]) # control in 1
	pp[12] = GasChromatographySystems.PressurePoint("p12", [0.0, 1800.0], [NaN, NaN]) #
	pp[13] = GasChromatographySystems.PressurePoint("p13", [0.0, 1800.0], [150000.0, 150000.0]) # control in 2
	# modules
	modules = Array{GasChromatographySystems.AbstractModule}(undef, ne(g))
	modules[1] = GasChromatographySystems.ModuleColumn("GC", 30, 0.25e-3, 0.25e-6, "Rxi17SilMS", default_TP)
	modules[2] = GasChromatographySystems.ModuleColumn("2->3", 4.5e-3, 0.1376e-3, 0.0e-6, "", default_TP)
	modules[3] = GasChromatographySystems.ModuleColumn("3->4", 5.0e-3, 0.1376e-3, 0.0e-6, "", default_TP)
	modules[4] = GasChromatographySystems.ModuleColumn("3->7", 5.0e-3, 0.1376e-3, 0.0e-6, "", default_TP)
	modules[5] = GasChromatographySystems.ModuleColumn("4->5", 2.8e-3, 0.1376e-3, 0.0e-6, "", default_TP)
	modules[6] = GasChromatographySystems.ModuleColumn("to Det1", 1.0, 0.18e-3, 0.0e-6, "", default_TP)
	modules[7] = GasChromatographySystems.ModuleColumn("7->8", 2.8e-3, 0.1376e-3, 0.0e-6, "", default_TP)
	modules[8] = GasChromatographySystems.ModuleColumn("to Det2", 1.0, 0.1e-3, 0.0e-6, "", default_TP)
	modules[9] = GasChromatographySystems.ModuleColumn("10->4", 5.0e-3, 0.1376e-3, 0.0e-6, "", default_TP)
	modules[10] = GasChromatographySystems.ModuleColumn("bypass", 20.0e-3, 0.0751e-3, 0.0e-6, "", default_TP)
	modules[11] = GasChromatographySystems.ModuleColumn("control in 1", 4.8e-3, 0.1376e-3, 0.0e-6, "", default_TP)
	modules[12] = GasChromatographySystems.ModuleColumn("12->7", 5.0e-3, 0.1376e-3, 0.0e-6, "", default_TP)
	modules[13] = GasChromatographySystems.ModuleColumn("control in 2", 4.8e-3, 0.1376e-3, 0.0e-6, "", default_TP)
	# system
	sys = GasChromatographySystems.update_system(GasChromatographySystems.System(g, pp, modules, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ 3dbd35e6-087d-4add-aa77-dc34ed429b46
collect(edges(sys.g))

# ╔═╡ 56c482e8-078c-4440-9dbd-51386abf1b4b
GasChromatographySystems.plot_graph(sys, lay=Stress())

# ╔═╡ 4b379226-1646-4b13-b64a-7b55b66442f5
md"""
## Pressure and flow calculation
"""

# ╔═╡ bfd7ed6a-f293-448f-b827-a914396f101b
GasChromatographySystems.flow_balance(sys)

# ╔═╡ 3946d9e6-78e3-41df-b204-4fd780721242
GasChromatographySystems.solve_balance(sys);

# ╔═╡ 0aef6b6a-c58d-4bc8-812a-e958973598b9
GasChromatographySystems.plot_graph_with_flow(sys, 0, lay=Stress())#, node_size=90, nlabels_fontsize=20, elabels_fontsize=20, elabels_distance=28)

# ╔═╡ a8b54992-6dec-42f2-b29b-7fa4f2128809
begin
	gr()
	p_F_sys = GasChromatographySystems.plot_flow_over_time(sys)
	Plots.plot!(p_F_sys, legend=:right)
	p_F_sys
end

# ╔═╡ e842ff4e-b999-4794-9b73-d6754bb01d55
#savefig(p_F_sys6, "Flow_loop.svg")

# ╔═╡ 306a6211-5709-43ae-9154-2efed03c62e0
begin
	plotly()
	p_p_sys = GasChromatographySystems.plot_pressure_over_time(sys)
	Plots.plot!(p_p_sys, legend=:right)
	#Plots.ylims!(4.9e5, 6.05e5)
	p_p_sys
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
db = filter([:CAS] => x -> !ismissing(x), DataFrame(CSV.File(joinpath(dirname(pwd()), "data", "database_Rxi17SilMS_Rxi5SilMS_.csv"), header=1, silencewarnings=true)))

# ╔═╡ 4f17cc8c-0004-44bd-b05d-97ee493235e0
unique(GasChromatographySystems.all_stationary_phases(sys))

# ╔═╡ 0e5c952e-d0e4-442e-96db-f0a847831acd
selected_solutes = GasChromatographySystems.common_solutes(db, sys).Name[11:20]

# ╔═╡ 46414588-5b80-4995-80d1-51f058f6360b
md"""
### Graph to Parameters
"""

# ╔═╡ d2176d2e-9858-4a6b-be9e-be92768a2dbe
par = GasChromatographySystems.graph_to_parameters(sys, db, selected_solutes; dt=60.0)

# ╔═╡ f821a4c0-4312-4d04-80f5-19f3b87ead97
md"""
### Paths
"""

# ╔═╡ b6f2e31c-9203-4101-aa7d-9aedee1de6b4
paths = GasChromatographySystems.all_paths(sys.g, 2)[2]

# ╔═╡ c44a7a01-59df-4663-9280-96d2ed7de6e3
md"""
### Simulation
"""

# ╔═╡ 87757bee-79f3-48ff-92b7-c9d3b5e94b2d
res = GasChromatographySystems.simulate_along_paths(sys, paths, db, selected_solutes, par; τ₀=0.5.*ones(length(selected_solutes)))

# ╔═╡ cd05457c-acf4-4bc3-9d47-2a0efdb2dd89
md"""
### Split at p3
"""

# ╔═╡ c5e3cb5f-11f2-476b-9e1c-afad3d024f19
t_split = res[2][2][1].tR

# ╔═╡ 5b35524d-5f94-402c-a508-2f6d17e3bf42
F_func = GasChromatographySystems.flow_functions(sys)

# ╔═╡ c0ba2b5d-199a-453f-8ca3-c7122239c279
collect(edges(sys.g))

# ╔═╡ ab47be67-a524-4584-8925-02ae19a54298
F_func[2].(t_split).*60e6

# ╔═╡ 5f0c8474-8d19-4e52-855c-23edb5478d41
F_func[3].(t_split).*60e6

# ╔═╡ 22494320-33a0-4f8b-b739-038d6661afe5
F_func[4].(t_split).*60e6

# ╔═╡ 722a5e35-1384-428c-bb48-719d90b06c1a
split_ratio_1 = F_func[3].(t_split)./F_func[2].(t_split)

# ╔═╡ aa28b189-2d36-49ba-b088-deab285aa479
split_ratio_2 = F_func[4].(t_split)./F_func[2].(t_split)

# ╔═╡ fba129ea-101b-41e6-9d75-ef493b9fd0a0
split_ratio_1 + split_ratio_2

# ╔═╡ a22e348f-353d-471c-8543-6e090adb6ecb
tend = sum(GasChromatographySystems.common_timesteps(sys))

# ╔═╡ 45b018aa-1e39-415c-8808-70fc54096212
begin # peak areas without split
	p_chrom = GasChromatographySimulator.plot_chromatogram(res[2][1][end], (0, tend))[1]
	GasChromatographySimulator.plot_chromatogram!(p_chrom, res[2][2][end], (0, tend); mirror=true)
	p_chrom
end

# ╔═╡ 7ff0bd0f-9dd4-4f87-98cb-178925745042
function chromatogram(t, tR, τR; split_ratio=ones(length(tR)))
	g(t,tR,τR) = 1/sqrt(2*π*τR^2)*exp(-(t-tR)^2/(2*τR^2))
	chromatograms = Array{Array{Float64,1}}(undef, length(tR))
	for j=1:length(tR)
		chromatograms[j] = g.(t, tR[j], τR[j]).*split_ratio[j]
	end
	return sum(chromatograms)
end

# ╔═╡ 880d893a-7359-47d4-8bab-d8560a9d7ec2
t = 0.0:0.01:tend

# ╔═╡ e77c0ad3-018e-4fba-beef-9b6626a436c4
c1 = chromatogram(t, res[2][1][end].tR, res[2][1][end].τR; split_ratio=split_ratio_1)

# ╔═╡ f0bd4d2d-0435-444a-b61e-eab74e995970
c2 = chromatogram(t, res[2][2][end].tR, res[2][2][end].τR; split_ratio=split_ratio_2)

# ╔═╡ ac2a4635-c57f-44ad-bc00-e396eb18e89a
begin
	p_chrom_ = Plots.plot(t, c1, xlabel="time in s", legend=false)
	Plots.plot!(p_chrom_, t, -1.0.*c2)
	p_chrom_
end

# ╔═╡ e0b310dd-074b-42fa-bece-ee4488fa0590
md"""
# End
"""

# ╔═╡ Cell order:
# ╠═e9d64f0e-95b7-11ed-2fee-ab579cd02ebc
# ╠═139bfe64-17ac-455d-84a2-09c9564156cd
# ╠═f6c6443d-0e39-4e82-b2c4-49f2954e5b49
# ╠═719dc69d-2692-4d41-868e-14b664019d57
# ╠═f484559b-95d5-4155-aad6-fdeed3d238bc
# ╠═7cc0cc2b-db2f-40aa-8f76-d99b1958e817
# ╠═e9607708-458d-4a85-8ba6-0f404e35d441
# ╠═3dbd35e6-087d-4add-aa77-dc34ed429b46
# ╠═56c482e8-078c-4440-9dbd-51386abf1b4b
# ╠═4b379226-1646-4b13-b64a-7b55b66442f5
# ╠═bfd7ed6a-f293-448f-b827-a914396f101b
# ╠═3946d9e6-78e3-41df-b204-4fd780721242
# ╠═0aef6b6a-c58d-4bc8-812a-e958973598b9
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
# ╠═a22e348f-353d-471c-8543-6e090adb6ecb
# ╠═45b018aa-1e39-415c-8808-70fc54096212
# ╠═7ff0bd0f-9dd4-4f87-98cb-178925745042
# ╠═880d893a-7359-47d4-8bab-d8560a9d7ec2
# ╠═e77c0ad3-018e-4fba-beef-9b6626a436c4
# ╠═f0bd4d2d-0435-444a-b61e-eab74e995970
# ╠═ac2a4635-c57f-44ad-bc00-e396eb18e89a
# ╠═e0b310dd-074b-42fa-bece-ee4488fa0590
