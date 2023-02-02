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
# Complex GC system

**With multiple columns and split points.**
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
TP1 = GasChromatographySystems.TemperatureProgram([0.0, 60.0, 1800.0, 60.0], [40.0, 40.0, 280.0, 280.0])

# ╔═╡ 0710b79d-8cf5-4f02-a042-237b86ef216f
(280-40)/(1440/60)

# ╔═╡ 06cfa2c9-4867-4c9d-96e5-eb021b8cdb43
420/60

# ╔═╡ b4bc0ec9-7184-4862-93d6-dc562de4aa81
TP2 = GasChromatographySystems.TemperatureProgram([0.0, 60.0, 1440.0, 420.0], [40.0, 40.0, 280.0, 280.0])

# ╔═╡ 06f73bb4-3263-47f7-a8de-a4f56478ebfe
md"""
### System
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
	pp5[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 60.0, 1800.0, 60.0], [450000.0, 450000.0, 600000.0, 600000.0]) # inlet 
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
	modules5[1] = GasChromatographySystems.ModuleColumn("GC column 1", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", TP1)
	modules5[2] = GasChromatographySystems.ModuleColumn("GC column 2", 30.0, 0.25e-3, 0.25e-6, "SPB50", TP1)
	modules5[3] = GasChromatographySystems.ModuleColumn("TL column 1", 5.0, 0.1e-3, 0.1e-6, "", 300.0)
	modules5[4] = GasChromatographySystems.ModuleColumn("TL column 2", 1.0, 0.25e-3, 0.0e-6, "", 300.0)
	modules5[5] = GasChromatographySystems.ModuleColumn("Modulator", 0.3, 0.1e-3, 0.1e-6, "Wax", TP2) # replace later with ModuleModulator
	modules5[6] = GasChromatographySystems.ModuleColumn("GC column 3", 2.0, 0.1e-3, 0.1e-6, "Wax", TP2)
	modules5[7] = GasChromatographySystems.ModuleColumn("TL column 3", 2.5, 0.1e-3, 0.0e-6, "", 300.0)
	modules5[8] = GasChromatographySystems.ModuleColumn("TL column 4", 1.5, 0.1e-3, 0.1e-6, "", 300.0)
	# system
	sys5 = GasChromatographySystems.update_system(GasChromatographySystems.System(g5, pp5, modules5, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ bc10e868-d48b-492b-ad65-170f38e4e340
GasChromatographySystems.plot_graph(sys5; lay=Stress())

# ╔═╡ 4b379226-1646-4b13-b64a-7b55b66442f5
md"""
## Pressure and flow calculation
"""

# ╔═╡ 3f9550a5-53c0-4630-aabd-177c1b2a4ccd
GasChromatographySystems.flow_balance(sys5)

# ╔═╡ 9359fc04-ea17-4551-9923-52c9d7509f2c
GasChromatographySystems.solve_balance(sys5)

# ╔═╡ 459ddf82-89d5-47a0-a535-699c015202bf
GasChromatographySystems.plot_graph_with_flow(sys5, 0, lay = Stress(), arrow_shift=0.6, nlabels_fontsize=16, elabels_fontsize=16)

# ╔═╡ e0fad228-0ba2-427e-8702-2d4640456e53
begin
	plotly()
	p_flow = GasChromatographySystems.plot_flow_over_time(sys5)
	Plots.ylims!(0.3,1.3)
	Plots.yticks!(0.3:0.1:1.3)
	p_flow
end

# ╔═╡ d30c1090-0c60-41ee-8c46-28f6ac013d49
#savefig(p_flow, "flow_complex.svg")

# ╔═╡ e461d166-a7da-40fc-a07d-2c2d1165f03f
begin
	p_pressure = GasChromatographySystems.plot_pressure_over_time(sys5)
	Plots.plot!(p_pressure, legend=:bottomright)
	Plots.ylims!(240000.0, 605000.0)
	p_pressure
end

# ╔═╡ 2c7657f7-35f2-483b-9a2a-716a60186b8a
#savefig(p_pressure, "pressure_complex.svg")

# ╔═╡ 975eed72-5f77-4e44-bf33-6be175638ebd
md"""
## Simulation
"""

# ╔═╡ 3570aceb-42ca-4696-b60f-79e78a5c02bc
md"""
### Selection of solutes
"""

# ╔═╡ f83611a0-319b-4c5c-8fd1-237c170dbf92
db = filter([:CAS] => x -> !ismissing(x), DataFrame(CSV.File(joinpath(dirname(pwd()), "data", "database_Blumberg2017.csv"), header=1, silencewarnings=true)))

# ╔═╡ b1a7cd90-9050-4d69-ac67-f8fb06508159
unique(GasChromatographySystems.all_stationary_phases(sys5))

# ╔═╡ 2554d565-f5fa-45be-a6c4-8f063d468845
selected_solutes = GasChromatographySystems.common_solutes(db, sys5).Name

# ╔═╡ 408164f8-0228-413c-8c8d-26b079d8aca1
md"""
### Graph to Parameters
"""

# ╔═╡ 38edeb8f-52c4-4246-84d6-82eb59e73773
par = GasChromatographySystems.graph_to_parameters(sys5, db, selected_solutes; dt=60.0)

# ╔═╡ 67c96cfc-ea58-4a2c-9b01-101fdacb6110
md"""
### Paths
"""

# ╔═╡ 2e682f84-1890-4ed9-90ce-2c46ced16df7
paths = GasChromatographySystems.all_paths(sys5.g, 3)[2]

# ╔═╡ 7a92fd8d-50e9-4551-b7ef-c6e72376ac4e
md"""
### Simulation
"""

# ╔═╡ 7d9ccb8f-240c-4bc3-a222-79c227b9bf69
begin
	ref = falses(ne(sys5.g))
	ref[6] = true
	ref
end

# ╔═╡ 8b9582a9-edc4-4252-83e2-30d1c76f58a5
res = GasChromatographySystems.simulate_along_paths(sys5, paths, db, selected_solutes, par; τ₀=0.5.*ones(length(selected_solutes)), τ₀_focus=0.1.*ones(length(selected_solutes)), refocus=ref)

# ╔═╡ 877b51b1-0bf1-4bb9-b0bb-dd25e96019ab
md"""
### Split at p₃

Time of split different for every solute at this point.
"""

# ╔═╡ 6ef0b010-1b0b-419e-a552-2286b5758af6
F_func = GasChromatographySystems.flow_functions(sys5)

# ╔═╡ e7db0062-02b8-412c-8a18-28bb39b69836
t_split_1 = res[2][1][2].tR

# ╔═╡ 44b2bb81-cc0d-4477-b48f-6e77f5ec76e0
collect(edges(sys5.g))

# ╔═╡ 5f9b52d8-a4f6-48b4-8fc2-7b7f0ca73a5c
paths

# ╔═╡ 5d17c71a-abeb-4b88-84c5-1bdcdddccdd6
F_func[2].(t_split_1).*60e6 # Flow comming from GC2

# ╔═╡ ac35bbc5-4325-4f9f-b21a-dbdd5a38ae3b
F_func[3].(t_split_1).*60e6 # Flow to 1st detector

# ╔═╡ 20eed3d5-fbdb-44e0-a83b-4f0990eebfed
F_func[4].(t_split_1).*60e6 # Flow to the 2nd and 3rd path

# ╔═╡ b9e8f414-b64b-4c7e-89ca-13c24efe6491
split_ratio_1 = F_func[3].(t_split_1)./F_func[2].(t_split_1)

# ╔═╡ 0d5708a4-75cc-4a8a-b736-062f3f1d4aba
split_ratio_23 = F_func[4].(t_split_1)./F_func[2].(t_split_1)

# ╔═╡ ea1d4291-8fa6-430f-9bf0-51ebca0250f8
md"""
### Split at p₇
"""

# ╔═╡ 641a6b04-d717-4621-9fe7-7cd8bcabfd83
t_split_2 = res[2][2][5].tR

# ╔═╡ 6f8d4b51-943b-4db2-892c-641547172646
F_func[6].(t_split_2).*60e6 # Flow comming from GC3

# ╔═╡ dac44d26-c00d-49a0-a774-556447e11c29
F_func[7].(t_split_2).*60e6 # Flow to 2nd detector (vacuum)

# ╔═╡ a6790f7d-4585-48c4-9501-c929d078a056
F_func[8].(t_split_2).*60e6 # Flow to 3rd detector (atmospheric)

# ╔═╡ 75fc18fe-8691-417a-8af9-07c7055ecae9
split_ratio_2 = F_func[7].(t_split_2)./F_func[6].(t_split_2)

# ╔═╡ 50aab675-7e7b-499d-833b-68a268a41e89
split_ratio_3 = F_func[8].(t_split_2)./F_func[6].(t_split_2)

# ╔═╡ d4d9f37c-a8b3-4a3a-8495-32c590e0a8ec
total_split_ratio_1 = split_ratio_1

# ╔═╡ 4881b2cb-bd0e-4330-9cfa-cce42ffb76f2
total_split_ratio_2 = split_ratio_23 .* split_ratio_2

# ╔═╡ 00ae0e5e-88bb-4154-9f67-f49843761d09
total_split_ratio_3 = split_ratio_23 .* split_ratio_3

# ╔═╡ 1f6d8da5-dd72-4b7d-98c5-108cb4527733
total_split_ratio_1 .+ total_split_ratio_2 .+ total_split_ratio_3

# ╔═╡ e9464022-6b07-439d-81cd-bca88cc8619a
md"""
### Chromatograms
"""

# ╔═╡ ae3f4430-7740-45c8-b49c-32fd350de7d2
md"""
#### After GC1
"""

# ╔═╡ 45293478-31cd-4d45-9791-1d5c355623d8
p_GC1, t_GC1, c_GC1 = GasChromatographySimulator.plot_chromatogram(res[2][1][1], (0, 1600.0))

# ╔═╡ e017eb1f-65f2-4c39-b90f-5371773cab26
md"""
#### After GC2
"""

# ╔═╡ 984e072b-52bf-4065-9f3e-58c7041ebd08
p_GC2, t_GC2, c_GC2 = GasChromatographySimulator.plot_chromatogram(res[2][1][2], (0, 1600.0))

# ╔═╡ cf089987-409d-4fec-bd19-c51263496572
md"""
#### At 1st detector (p₄)
"""

# ╔═╡ 4db43001-391b-4983-b10b-3de7cdc1d6a3
p_det1, t_det1, c_det1 = GasChromatographySimulator.plot_chromatogram(res[2][1][3], (0, 1600.0));

# ╔═╡ 6a0e8068-705e-4fcd-869a-07223253a2cf
function chromatogram(t, tR, τR; split_ratio=ones(length(tR)))
	g(t,tR,τR) = 1/sqrt(2*π*τR^2)*exp(-(t-tR)^2/(2*τR^2))
	chromatograms = Array{Array{Float64,1}}(undef, length(tR))
	for j=1:length(tR)
		chromatograms[j] = g.(t, tR[j], τR[j]).*split_ratio[j]
	end
	return sum(chromatograms)
end

# ╔═╡ 1bb0626d-2869-498f-b9c1-ce2dd14be382
c_det1_split = chromatogram(t_det1, res[2][1][3].tR, res[2][1][3].τR, split_ratio=total_split_ratio_1)

# ╔═╡ b3ce2b29-e811-4491-af0a-962d20e04148
begin
	gr()
	p_Det1 = Plots.plot(t_det1, c_det1_split, xlabel="time in s", label="1st detector", legend=false)
end

# ╔═╡ d63ee52f-4a6b-43fb-92da-51895b377109
#savefig(p_Det1, "chrom_complex_Det1.svg")

# ╔═╡ a8b5f914-5d7b-4650-be71-07857d5be4e4
md"""
#### After GC3 (p₇)
"""

# ╔═╡ 9afa3e82-2c39-4c1a-b88d-fc07d4298dc8
t_GC3 = 0.0:0.001:1600.0

# ╔═╡ 9fff7679-8370-4ef1-ad1d-d04d9c0513bc
c_GC3_split = chromatogram(t_GC3, res[2][2][5].tR, res[2][2][5].τR, split_ratio=split_ratio_23)

# ╔═╡ 0fae19b2-f88f-4fbc-bc8f-4ffcc0f96b4c
Plots.plot(t_GC3, c_GC3_split, xlabel="time in s", label="after GC3", legend=:topleft)

# ╔═╡ 4df0ee58-51fd-448f-b7a5-b5336260292a
md"""
#### At 2nd detector (p₈)
"""

# ╔═╡ 8435ea41-e027-4cdc-912a-bc5103f156c6
t_det2 = 0.0:0.001:1600.0

# ╔═╡ d19e5d20-092e-4567-9980-6811e18ef61b
c_det2_split = chromatogram(t_det2, res[2][2][6].tR, res[2][2][6].τR, split_ratio=total_split_ratio_2)

# ╔═╡ 00ae08bd-2f29-485a-8b08-00f7e2b620b0
Plots.plot(t_det2, c_det2_split, xlabel="time in s", label="2nd detector", legend=:topleft)

# ╔═╡ 70906a07-1d8b-46bf-9af8-8c2daa9662cf
md"""
#### At 3rd detector (p₉)
"""

# ╔═╡ 44288d1a-78f0-4895-b409-77fe178a3581
t_det3 = 0.0:0.001:1600.0

# ╔═╡ 4aec32ae-f28c-4836-9222-8aa87f108a5b
c_det3_split = chromatogram(t_det3, res[2][3][6].tR, res[2][3][6].τR, split_ratio=total_split_ratio_3)

# ╔═╡ fae813eb-f8a4-4d27-a86b-cebf951e8b5b
Plots.plot(t_det3, c_det3_split, xlabel="time in s", label="3rd detector", legend=:topleft)

# ╔═╡ 98ebf65a-cbce-4e19-bd82-3d393cf9941d
md"""
#### Combine chromatograms at detectors
"""

# ╔═╡ 3a337817-8858-4e42-b0fc-96a1eda6967f
begin
	p_chrom = Plots.plot(t_det1, c_det1_split, xlabel="time in s", label="1st detector", legend=:topleft)
	Plots.plot!(p_chrom, t_det2, c_det2_split.+0.2, xlabel="time in s", label="2nd detector", legend=:topleft)
	Plots.plot!(p_chrom, t_det3, c_det3_split.+0.4, xlabel="time in s", label="3rd detector", legend=:topleft)
	p_chrom
end

# ╔═╡ 0a462311-33e4-4340-928f-46b78562283e
md"""
### 2D Chromatograms
"""

# ╔═╡ 5b8cb487-c57b-4d05-ab44-0ac8116163ad
md"""
#### At 2nd detector
"""

# ╔═╡ 046488dc-2031-4fc6-bd9f-1db781916949
function peaklist_GCxGC(pl_1, pl_2)
	CAS1 = GasChromatographySimulator.CAS_identification(pl_1.Name).CAS
	pl_GCxGC = DataFrame(Name = pl_1.Name, CAS = CAS1, tR1 = pl_1.tR, τR1 = pl_1.τR)
	tR2 = Array{Float64}(undef, length(pl_GCxGC.Name))
	τR2 = Array{Float64}(undef, length(pl_GCxGC.Name))
	Cat = Array{Array{String}}(undef, length(pl_GCxGC.Name))
	CAS2 = GasChromatographySimulator.CAS_identification(pl_2.Name).CAS
	for i=1:length(pl_GCxGC.Name)
		ii = findfirst(CAS1[i].==CAS2)
		tR2[i] = pl_2.tR[ii]
		τR2[i] = pl_2.τR[ii]
		i_db = findfirst(CAS1[i].==db.CAS)
		Cat[i] = unique([string(db.Cat_1[i_db]), string(db.Cat_2[i_db]), string(db.Cat_3[i_db]), string(db.Cat_4[i_db])])
	end
	pl_GCxGC[!, :tR2] = tR2
	pl_GCxGC[!, :τR2] = τR2
	pl_GCxGC[!, :Cat] = Cat
	return pl_GCxGC
end

# ╔═╡ 123c4821-a7d8-449b-86af-77dc55aae322
pl_GCxGC_2 = peaklist_GCxGC(res[2][2][4], res[2][2][6])

# ╔═╡ 0ea0b11e-c99d-4362-a304-ce021b406176
function plot_GCxGC(pl_GCxGC, sys; categories = String[])
	t1end = sum(GasChromatographySystems.common_timesteps(sys))
	t2end = ceil(maximum(pl_GCxGC.tR2.-pl_GCxGC.tR1))
	p_gcxgc = Plots.scatter(pl_GCxGC.tR1, pl_GCxGC.tR2.-pl_GCxGC.tR1, xlims=(0.0,t1end), ylims=(0.0, t2end), label="", xlabel=string(sys.modules[1].stationary_phase, ", tR1 in s"), ylabel=string(sys.modules[2].stationary_phase, ", tR2 in s"))
	if !isempty(categories)
		for i=1:length(categories) 
			pl_f = filter([:Cat] => x -> categories[i] in x, pl_GCxGC)
			Plots.scatter!(p_gcxgc, pl_f.tR1, pl_f.tR2.-pl_f.tR1, label=categories[i], legend=:bottomright)
		end
	end
	return p_gcxgc
end

# ╔═╡ 126f7ea0-67bd-4f67-b9f3-addc37abf15f
function plot_GCxGC_contour(pl_GCxGC; split_ratio=ones(length(pl_GCxGC.Name)))
	t¹ = 0.9*minimum(pl_GCxGC.tR1):(1.1*maximum(pl_GCxGC.tR1)-0.9*minimum(pl_GCxGC.tR1))/1000:1.1*maximum(pl_GCxGC.tR1)
	t² = 0.9*minimum(pl_GCxGC.tR2.-pl_GCxGC.tR1):(1.1*maximum(pl_GCxGC.tR2.-pl_GCxGC.tR1)-0.9*minimum(pl_GCxGC.tR2.-pl_GCxGC.tR1))/1000:1.1*maximum(pl_GCxGC.tR2.-pl_GCxGC.tR1)
	chrom2D = Array{Float64}(undef, length(t¹), length(t²), length(pl_GCxGC.Name))
	for j=1:length(t²)
		for i=1:length(t¹)
			for k=1:length(pl_GCxGC.Name)
				chrom2D[i,j,k] = split_ratio[k]*exp(-((t¹[i]-pl_GCxGC.tR1[k])^2/(2*pl_GCxGC.τR1[k]^2)+(t²[j]-(pl_GCxGC.tR2[k]-pl_GCxGC.tR1[k]))^2/(2*pl_GCxGC.τR2[k]^2)))
			end
		end
	end
	p_2D = Plots.contour(t¹, t², sum(chrom2D, dims=3)[:,:,1]', fill=false, legend=false, xlabel="tR1 in s", ylabel="tR2 in s")
	return p_2D
end

# ╔═╡ 5231a1a0-7835-4380-83a0-e85f9a6c24d7
begin
	p_Det2_2D = plot_GCxGC_contour(pl_GCxGC_2; split_ratio=total_split_ratio_2)
	Plots.ylims!(5.4, 6.4)
	p_Det2_2D
end

# ╔═╡ 10eac303-5e44-4785-8f66-a243894fe044
plot_GCxGC_contour(pl_GCxGC_2)

# ╔═╡ c6e0e23d-f156-4e34-be6a-38c8eb432819
plot_GCxGC(pl_GCxGC_2, sys5)

# ╔═╡ 6e9f2180-52bc-427d-bed7-9b8f89d6e4f4
md"""
#### At 3rd detector
"""

# ╔═╡ 69dcf445-d9c5-4a31-9779-c3da545b76f2
pl_GCxGC_3 = peaklist_GCxGC(res[2][3][4], res[2][3][6])

# ╔═╡ 2c6b1127-2f26-426e-a03f-b01916c08c06
begin
	p_Det3_2D = plot_GCxGC_contour(pl_GCxGC_3; split_ratio=total_split_ratio_3)
	#Plots.ylims!(5.4, 6.4)
	p_Det3_2D
end

# ╔═╡ e08e96bd-b937-44c6-91c0-908423b521ae
savefig(p_Det3_2D, "chrom_complex_Det3.svg")

# ╔═╡ 1ec22e27-f1af-49fc-894a-90a64805a7a1
plot_GCxGC(pl_GCxGC_3, sys5)

# ╔═╡ 4e25336a-06b5-4f5f-8610-b201cf96441c
md"""
# End
"""

# ╔═╡ Cell order:
# ╠═e9d64f0e-95b7-11ed-2fee-ab579cd02ebc
# ╠═139bfe64-17ac-455d-84a2-09c9564156cd
# ╠═f6c6443d-0e39-4e82-b2c4-49f2954e5b49
# ╠═719dc69d-2692-4d41-868e-14b664019d57
# ╠═f484559b-95d5-4155-aad6-fdeed3d238bc
# ╠═0710b79d-8cf5-4f02-a042-237b86ef216f
# ╠═06cfa2c9-4867-4c9d-96e5-eb021b8cdb43
# ╠═b4bc0ec9-7184-4862-93d6-dc562de4aa81
# ╠═06f73bb4-3263-47f7-a8de-a4f56478ebfe
# ╠═8655ab62-f1a0-4abb-8856-2eb6ca9736f9
# ╠═bc10e868-d48b-492b-ad65-170f38e4e340
# ╠═4b379226-1646-4b13-b64a-7b55b66442f5
# ╠═3f9550a5-53c0-4630-aabd-177c1b2a4ccd
# ╠═9359fc04-ea17-4551-9923-52c9d7509f2c
# ╠═459ddf82-89d5-47a0-a535-699c015202bf
# ╠═e0fad228-0ba2-427e-8702-2d4640456e53
# ╠═d30c1090-0c60-41ee-8c46-28f6ac013d49
# ╠═e461d166-a7da-40fc-a07d-2c2d1165f03f
# ╠═2c7657f7-35f2-483b-9a2a-716a60186b8a
# ╠═975eed72-5f77-4e44-bf33-6be175638ebd
# ╠═3570aceb-42ca-4696-b60f-79e78a5c02bc
# ╠═f83611a0-319b-4c5c-8fd1-237c170dbf92
# ╠═b1a7cd90-9050-4d69-ac67-f8fb06508159
# ╠═2554d565-f5fa-45be-a6c4-8f063d468845
# ╠═408164f8-0228-413c-8c8d-26b079d8aca1
# ╠═38edeb8f-52c4-4246-84d6-82eb59e73773
# ╠═67c96cfc-ea58-4a2c-9b01-101fdacb6110
# ╠═2e682f84-1890-4ed9-90ce-2c46ced16df7
# ╠═7a92fd8d-50e9-4551-b7ef-c6e72376ac4e
# ╠═7d9ccb8f-240c-4bc3-a222-79c227b9bf69
# ╠═8b9582a9-edc4-4252-83e2-30d1c76f58a5
# ╠═877b51b1-0bf1-4bb9-b0bb-dd25e96019ab
# ╠═6ef0b010-1b0b-419e-a552-2286b5758af6
# ╠═e7db0062-02b8-412c-8a18-28bb39b69836
# ╠═44b2bb81-cc0d-4477-b48f-6e77f5ec76e0
# ╠═5f9b52d8-a4f6-48b4-8fc2-7b7f0ca73a5c
# ╠═5d17c71a-abeb-4b88-84c5-1bdcdddccdd6
# ╠═ac35bbc5-4325-4f9f-b21a-dbdd5a38ae3b
# ╠═20eed3d5-fbdb-44e0-a83b-4f0990eebfed
# ╠═b9e8f414-b64b-4c7e-89ca-13c24efe6491
# ╠═0d5708a4-75cc-4a8a-b736-062f3f1d4aba
# ╠═ea1d4291-8fa6-430f-9bf0-51ebca0250f8
# ╠═641a6b04-d717-4621-9fe7-7cd8bcabfd83
# ╠═6f8d4b51-943b-4db2-892c-641547172646
# ╠═dac44d26-c00d-49a0-a774-556447e11c29
# ╠═a6790f7d-4585-48c4-9501-c929d078a056
# ╠═75fc18fe-8691-417a-8af9-07c7055ecae9
# ╠═50aab675-7e7b-499d-833b-68a268a41e89
# ╠═d4d9f37c-a8b3-4a3a-8495-32c590e0a8ec
# ╠═4881b2cb-bd0e-4330-9cfa-cce42ffb76f2
# ╠═00ae0e5e-88bb-4154-9f67-f49843761d09
# ╠═1f6d8da5-dd72-4b7d-98c5-108cb4527733
# ╠═e9464022-6b07-439d-81cd-bca88cc8619a
# ╠═ae3f4430-7740-45c8-b49c-32fd350de7d2
# ╠═45293478-31cd-4d45-9791-1d5c355623d8
# ╠═e017eb1f-65f2-4c39-b90f-5371773cab26
# ╠═984e072b-52bf-4065-9f3e-58c7041ebd08
# ╠═cf089987-409d-4fec-bd19-c51263496572
# ╠═4db43001-391b-4983-b10b-3de7cdc1d6a3
# ╠═6a0e8068-705e-4fcd-869a-07223253a2cf
# ╠═1bb0626d-2869-498f-b9c1-ce2dd14be382
# ╠═b3ce2b29-e811-4491-af0a-962d20e04148
# ╠═d63ee52f-4a6b-43fb-92da-51895b377109
# ╠═a8b5f914-5d7b-4650-be71-07857d5be4e4
# ╠═9afa3e82-2c39-4c1a-b88d-fc07d4298dc8
# ╠═9fff7679-8370-4ef1-ad1d-d04d9c0513bc
# ╠═0fae19b2-f88f-4fbc-bc8f-4ffcc0f96b4c
# ╠═4df0ee58-51fd-448f-b7a5-b5336260292a
# ╠═8435ea41-e027-4cdc-912a-bc5103f156c6
# ╠═d19e5d20-092e-4567-9980-6811e18ef61b
# ╠═00ae08bd-2f29-485a-8b08-00f7e2b620b0
# ╠═70906a07-1d8b-46bf-9af8-8c2daa9662cf
# ╠═44288d1a-78f0-4895-b409-77fe178a3581
# ╠═4aec32ae-f28c-4836-9222-8aa87f108a5b
# ╠═fae813eb-f8a4-4d27-a86b-cebf951e8b5b
# ╠═98ebf65a-cbce-4e19-bd82-3d393cf9941d
# ╠═3a337817-8858-4e42-b0fc-96a1eda6967f
# ╠═0a462311-33e4-4340-928f-46b78562283e
# ╠═5b8cb487-c57b-4d05-ab44-0ac8116163ad
# ╠═046488dc-2031-4fc6-bd9f-1db781916949
# ╠═123c4821-a7d8-449b-86af-77dc55aae322
# ╠═0ea0b11e-c99d-4362-a304-ce021b406176
# ╠═126f7ea0-67bd-4f67-b9f3-addc37abf15f
# ╠═5231a1a0-7835-4380-83a0-e85f9a6c24d7
# ╠═10eac303-5e44-4785-8f66-a243894fe044
# ╠═c6e0e23d-f156-4e34-be6a-38c8eb432819
# ╠═6e9f2180-52bc-427d-bed7-9b8f89d6e4f4
# ╠═69dcf445-d9c5-4a31-9779-c3da545b76f2
# ╠═2c6b1127-2f26-426e-a03f-b01916c08c06
# ╠═e08e96bd-b937-44c6-91c0-908423b521ae
# ╠═1ec22e27-f1af-49fc-894a-90a64805a7a1
# ╠═4e25336a-06b5-4f5f-8610-b201cf96441c
