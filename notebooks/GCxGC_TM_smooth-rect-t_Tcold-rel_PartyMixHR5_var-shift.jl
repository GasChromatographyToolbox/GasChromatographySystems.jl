### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 113c6e70-2168-11ee-3f7f-2775a67bab90
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
	using OrdinaryDiffEq
	using LsqFit
	using Interpolations
	using UrlDownload
	include(joinpath(dirname(pwd()), "src", "GasChromatographySystems.jl"))
	TableOfContents(depth=4)
end

# ╔═╡ 98217474-a16f-406a-83a7-17fee89c951a
html"""
	<style>
	  main {
		max-width: 1000px;
	  }
	</style>
"""

# ╔═╡ 22091a27-80e1-4e98-abe7-4b9652cb832c
md"""
# Test of the GCxGC with thermal modulator simulation

- smooth rectangle temperature function over time at modulation point (`spatial=false`, `tflank=20`, `algTM=Vern9` and `ng=true`)
- Tcold is relative (`Tcold_abs=false`)
- variation of the modulator shift and its influence on the retention times
"""

# ╔═╡ 1948baf0-1217-485a-8c8b-6e826bd5ef50
optTM = GasChromatographySystems.ModuleTMopt(Tcold_abs=false, ng=true, alg=Vern9(), tflank=20.0, sflank=Inf, dtinit=0.5e-7)

# ╔═╡ 71d12824-5fa7-4c51-94d0-46f6646df242
optCol = GasChromatographySystems.ModuleColumnOpt(ng=true)

# ╔═╡ 13d4415d-3dc4-4eed-aded-c3561b77a9eb
md"""
## Definition System
"""

# ╔═╡ ddde4e3d-4e72-4315-9555-5f8b3375b04c
md"""
### Temperature program
"""

# ╔═╡ 55194da2-9c5b-464b-bfe1-a65c48bc7801
begin
	tsteps_, Tsteps_ = GasChromatographySimulator.conventional_program([50.0, 2.0, 5.0, 200.0, 0.0, 5.0, 280.0, 10.0]) 
	GCxGC_TP = GasChromatographySystems.TemperatureProgram(tsteps_, Tsteps_)

	TP1 = GCxGC_TP
	TPM = GCxGC_TP
	TP2 = GCxGC_TP
	TTL = 280.0
end

# ╔═╡ 09438575-8370-475f-b2df-fda5155e8c82
md"""
### Columns
"""

# ╔═╡ d0379590-edf3-4017-b8b8-07bd6080f757
begin
	# 1st D 
	L1 = 29.5 # m
	d1 = 0.25 # mm
	df1 = 0.25 # µm
	sp1 = "Rxi5SilMS"
	
	# modulator
	LM = [0.30, 0.005, 0.9, 0.005, 0.3]#[0.3, 0.005, 0.9, 0.005, 0.3]
	dM = 0.25
	dfM = 0.25
	spM = "Rxi17SilMS"

	# 2nd D
	L2 = 0.245#0.56 # m
	d2 = 0.25 # mm
	df2 = 0.25 # µm
	sp2 = "Rxi17SilMS"

	# TL
	LTL = 0.235
	dTL = 0.25
	dfTL = 0.25
	spTL = "Rxi17SilMS"
end

# ╔═╡ 23b71d48-26bd-496d-9c10-e2d3ce2bc13c
sum(LM) + L2 + LTL

# ╔═╡ 0e349acc-46b4-4734-abb7-668ec1225c53
md"""
### Modulator
"""

# ╔═╡ de542d7e-bd0d-47a4-8977-5ab1e64a26f4
begin
	PM = 4.0
	shift = 0.0:0.1:4.0#3.35#(PM-0.35)#0.0#-0.35*2
	thot = 0.35
	ratio = (PM-thot)/PM
	Thot = 25.0
	Tcold = -130.0
end

# ╔═╡ daf7484b-ac12-4515-8e26-adbc8e214cbc
md"""
### Flows and pressures
"""

# ╔═╡ 83d10e79-3ab6-418c-80ba-653a50fa5d21
begin
	F = 0.8
	pin = NaN
	pout = 0.0
end

# ╔═╡ 5dd6d2b7-1d6f-49ef-b684-a453a50a9e96
md"""
### System
"""

# ╔═╡ 5883701b-96f6-431e-b29b-b5cc0be940e0
begin
	sys = Array{GasChromatographySystems.System}(undef, length(shift))
	Threads.@threads for i=1:length(shift)
		sys[i] = GasChromatographySystems.GCxGC_TM(L1, d1, df1, sp1, TP1, L2, d2, df2, sp2, TP2, LTL, dTL, dfTL, spTL, TTL, LM, dM, dfM, spM, shift[i], PM, ratio, Thot, Tcold, TPM, F, pin, pout; optTM=optTM, optCol=optCol)
	end
end

# ╔═╡ 6afa596b-70dd-4b01-899b-5d9ef1349491
GasChromatographySystems.plot_flow_over_time(sys[1]; dt=30.0)

# ╔═╡ 7ddfe44e-fde9-422f-98de-27b80e8517be
GasChromatographySystems.plot_pressure_over_time(sys[1]; dt=30.0)

# ╔═╡ 03b1cae0-cf56-4e0f-b0ad-fdfa736e0fa8
md"""
### Substances
"""

# ╔═╡ 397d7248-e069-4958-849d-6761cb20283d
begin
	db = DataFrame(urldownload("https://raw.githubusercontent.com/JanLeppert/RetentionData/main/Databases/GCSim_database_nonflag.csv"))
	#db.Tchar = db.Tchar .- 273.15
	insertcols!(db, 1, :No => collect(1:length(db.Name)))
	filter!([:phi0] => x -> x == 0.001, db)
	db
end

# ╔═╡ 6b2c6634-c13a-4b71-93ff-d66ab82854c9
unique(db.Phase)

# ╔═╡ 35bbe26a-f82c-4f60-a2d0-74b052037ae7
selected_solutes_ = GasChromatographySystems.common_solutes(db, sys[1]).Name

# ╔═╡ 510bf59c-a99e-4b30-a82e-4bbc1b336b7b
# 85 -> multiple entries in database for one of the stat. phases

# ╔═╡ 1964abf5-a627-4014-8348-e81145647104
md"""
## Simulation
"""

# ╔═╡ b284207c-70b5-4da5-94e9-0fe70680b0e3
paths = GasChromatographySystems.all_paths(sys[1].g, 1)[2]

# ╔═╡ bbb75e0b-ccab-4c32-9e33-0d8915dba95d
md"""
## Chromatograms
"""

# ╔═╡ 8572809a-0ea8-4520-8768-6d7f3288d06f
i_select = 35

# ╔═╡ b586403d-12ac-4b38-b4de-f974fe71b5dc
md"""
### 1D Chromatograms
"""

# ╔═╡ 7e945bad-0a59-4a5e-8c63-cef21ee847d0
#collect_chrom(sim[2][1], sys; markings=true)

# ╔═╡ f1f1ffdb-bc7e-48fd-b688-21f375ab21eb
#GasChromatographySystems.chrom_marked(sim[2][1][4], PM, ratio, shift1; nτ=6)[1]

# ╔═╡ 7739a83c-d001-47cd-aee2-36917106c2ad
md"""
### 2D Chromatogram
"""

# ╔═╡ f268df70-af5a-482f-85e6-84216681ecf4
md"""
## Peaklist
"""

# ╔═╡ 271bb7a6-9a98-4336-9571-c131cc125486
md"""
## Measurement
"""

# ╔═╡ b7c2070b-16bf-41f5-b021-160e54a53a21
meas = DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/GCxGC/Simulation/PartyMixHR5_RT.csv", header=1, silencewarnings=true, stringtype=String))

# ╔═╡ 7580ebee-b061-4c11-aaad-0e8310e5277d
begin
	b = [meas.CAS[x] in filter([:Phase] => y -> y == sp1 || y == sp2, db).CAS for x in 1:length(meas.CAS)]
	b[findall(ismissing.(b))] .= false
	possible_solutes_meas = meas.Name[findall(b)] 
end

# ╔═╡ b35d871c-4f74-4290-95a4-24ab7115c591
b

# ╔═╡ 262fba99-1b75-4803-8f34-a0619de559a9
selected_solutes = possible_solutes_meas[[4, 5, 8, 9, 11, 12]]#selected_solutes_[36:38]#[[36, 37, 38, 82, 83, 84, 129, 130, 131, 132]]#[Not([1, 29, 30])]#[[14,16,17]]

# ╔═╡ 5f3c1ce1-50bb-44ba-abac-d3ccc82fdb51
begin
	par = Array{Array{GasChromatographySimulator.Parameters,1}}(undef, length(shift))
	Threads.@threads for i=1:length(shift)
		par[i] = GasChromatographySystems.graph_to_parameters(sys[i], db, selected_solutes)
	end
end

# ╔═╡ e2e1ee2a-f98d-47fd-8c22-cc441f2d12c2
begin
	Plots.plot(cumsum(tsteps_), par[1][1].prog.T_itp(0.0, cumsum(tsteps_)).-273.15, label="")
end

# ╔═╡ a2d062a4-bd59-4ab7-9316-dc769b2c3c09
begin
	TM_itp = par[i_select][3].prog.T_itp
	Plots.plot(0.0:0.01:sys[i_select].modules[3].PM, TM_itp.(sys[i_select].modules[3].length/2, 0.0:0.01:sys[i_select].modules[3].PM).-273.115, title="temperature at modulator point, modulation shift = $(shift[i_select])s", xlabel="time in s", ylabel="temperature in °C", label="")
	Plots.plot!([-shift[i_select], -shift[i_select]], [TM_itp(sys[i_select].modules[3].length/2,-shift[i_select]), TM_itp(sys[i_select].modules[3].length/2,(1+ratio)/2*PM-shift[i_select])].-273.15, c=:black, label="")
	Plots.plot!([ratio*PM-shift[i_select], ratio*PM-shift[i_select]], [TM_itp(sys[i_select].modules[3].length/2,ratio*PM-shift[i_select]-eps()), TM_itp(sys[i_select].modules[3].length/2,(1+ratio)/2*PM-shift[i_select])].-273.15, c=:black, linestyle=:dash, label="")
	Plots.plot!([PM-shift[i_select], PM-shift[i_select]], [TM_itp(sys[i_select].modules[3].length/2,PM-shift[i_select]), TM_itp(sys[i_select].modules[3].length/2,(1+ratio)/2*PM-shift[i_select])].-273.15, c=:black, label="")

end

# ╔═╡ b8a6cdec-457b-4fb0-b6f8-60f713f4a015
begin
	Plots.plot(xlabel="position x in m", ylabel="temperature in °C", title="temperature at modulator point, modulation shift = $(shift[i_select])s", legend=false)
	Plots.plot!(0.0:sys[i_select].modules[3].length/100.0:sys[i_select].modules[3].length, TM_itp.(0.0:sys[i_select].modules[3].length/100.0:sys[i_select].modules[3].length, sys[i_select].modules[3].ratio*sys[i_select].modules[3].PM-sys[i_select].modules[3].shift).-273.15, label="t=$(sys[i_select].modules[3].ratio*sys[i_select].modules[3].PM-sys[i_select].modules[3].shift)s")
	Plots.plot!(0.0:sys[i_select].modules[3].length/100.0:sys[i_select].modules[3].length, TM_itp.(0.0:sys[i_select].modules[3].length/100.0:sys[i_select].modules[3].length, (1+3*sys[i_select].modules[3].ratio)/4*sys[i_select].modules[3].PM-sys[i_select].modules[3].shift).-273.15, label="t=$((1+3*sys[i_select].modules[3].ratio)/4*sys[i_select].modules[3].PM-sys[i_select].modules[3].shift)s")
	Plots.plot!(0.0:sys[i_select].modules[3].length/100.0:sys[i_select].modules[3].length, TM_itp.(0.0:sys[i_select].modules[3].length/100.0:sys[i_select].modules[3].length, (1+sys[i_select].modules[3].ratio)/2*sys[i_select].modules[3].PM-sys[i_select].modules[3].shift).-273.15, label="t=$((1+sys[i_select].modules[3].ratio)/2*sys[i_select].modules[3].PM-sys[i_select].modules[3].shift)s")
	Plots.plot!(0.0:sys[i_select].modules[3].length/100.0:sys[i_select].modules[3].length, TM_itp.(0.0:sys[i_select].modules[3].length/100.0:sys[i_select].modules[3].length, 3.95-sys[i_select].modules[3].shift).-273.15, label="t=$(3.95-sys[i_select].modules[3].shift)s")
	Plots.plot!(0.0:sys[i_select].modules[3].length/100.0:sys[i_select].modules[3].length, TM_itp.(0.0:sys[i_select].modules[3].length/100.0:sys[i_select].modules[3].length, sys[i_select].modules[3].PM-sys[i_select].modules[3].shift).-273.15, label="t=$(sys[i_select].modules[3].PM-sys[i_select].modules[3].shift)s")
end

# ╔═╡ e72a6c41-2b14-45cb-81e7-263610c3c8ae
begin
	Plots.plot(0.0:0.01:(sys[i_select].modules[3].PM*100.0), TM_itp.(sys[i_select].modules[3].length/2, 0.0:0.01:(sys[i_select].modules[3].PM*100.0)).-273.115, title="temperature at modulator point, modulation shift = $(shift[i_select])s", xlabel="time in s", ylabel="temperature in °C", label="")
	Plots.plot!(0.0:0.01:(sys[i_select].modules[3].PM*100.0), par[i_select][1].prog.T_itp(0.0, 0.0:0.01:(sys[i_select].modules[3].PM*100.0)).-273.15, label="")
end

# ╔═╡ e0f0bb8f-1be0-4d84-9e0d-d86d3df65d6b
par[1]

# ╔═╡ 2b41aa34-d994-440b-a8e1-dd6337008a48
begin
	sim = Array{Tuple{Vector{String}, Vector{Vector{DataFrame}}, Vector{Vector{Any}}, Vector{GasChromatographySimulator.Parameters}}}(undef, length(shift))
	Threads.@threads for i=1:length(shift)
		sim[i] = GasChromatographySystems.simulate_along_paths(sys[i], paths, par[i])
	end
end

# ╔═╡ d912d940-5b39-4bb9-a136-1f5e67ad7d82
GasChromatographySystems.chrom(sim[i_select][2][1][end])[1]

# ╔═╡ 582499db-7b5a-4054-92c2-d0e2a97cd91f
slice_mat, t_D1, t_D2, c_slices, t_, chrom_sliced_sum  = GasChromatographySystems.chrom2d(sim[i_select][2][1][end], sys[i_select], PM)

# ╔═╡ e5d4fb3a-f05f-4533-9666-08be32ef19ef
begin
	pl_GCxGC = Array{DataFrame}(undef, length(shift))
	Threads.@threads for i=1:length(shift)
		pl_GCxGC[i] = GasChromatographySystems.peaklist_GCxGC(sim[i][2][1][8], PM)
	end
	pl_GCxGC
end

# ╔═╡ 5f67ab8d-3b3d-4bbd-81c8-04eec158ef6e
begin
	plotly()
	x1 = floor(minimum(pl_GCxGC[i_select].tR1*0.99))
	x2 = ceil(maximum(pl_GCxGC[i_select].tR1*1.01))
	y1 = floor(minimum(pl_GCxGC[i_select].tR2*0.99))
	y2 = ceil(maximum(pl_GCxGC[i_select].tR2*1.01))
	p_contour = Plots.contour(t_D1[1:end-1], t_D2[1], slice_mat[1:end-1,1:end]', levels=10, xlabel="tR1 in s", ylabel="tR2 in s")#, c=:jet1)
	Plots.scatter!(p_contour, pl_GCxGC[i_select].tR1, pl_GCxGC[i_select].tR2, msize=3, shape=:diamond, c=:black, markeralpha=1)
	Plots.plot!(p_contour, xlims=(x1, x2), ylims=(y1, y2))

	p_heatmap = Plots.heatmap(t_D1[1:end-1], t_D2[1], slice_mat[1:end-1,1:end]', levels=10, xlabel="tR1 in s", ylabel="tR2 in s", c=:jet1)
	Plots.scatter!(p_heatmap, pl_GCxGC[i_select].tR1, pl_GCxGC[i_select].tR2, msize=3, shape=:diamond, c=:black, markeralpha=1)
	Plots.plot!(p_heatmap, xlims=(x1, x2), ylims=(y1, y2))

	Plots.plot(p_contour, p_heatmap, size=(900,800), legend=false)
end

# ╔═╡ 189034b6-09e5-4e6b-80b3-e442ad30705c
pl_GCxGC

# ╔═╡ 3ed4ad21-c8e4-428a-b6d3-b37cfdbd0778
begin
	compare = Array{DataFrame}(undef, length(shift))
	Threads.@threads for i=1:length(shift)
		compare[i] = GasChromatographySystems.comparison_meas_sim(meas, pl_GCxGC[i])
		filter!([:tR1_meas, :tR1_sim] => (x, y) -> !ismissing(x) && !ismissing(y), compare[i])
	end
	compare
end

# ╔═╡ aabb8d13-85a8-4252-9c71-99b78dd9e2cf
md"""
## Change of retention times with modulation shift
"""

# ╔═╡ 76236614-3c78-4c60-9eda-d906337362e7
begin
	ΔtR1s = Array{Float64}(undef, length(shift), length(compare[1].Name))
	ΔtR2s = Array{Float64}(undef, length(shift), length(compare[1].Name))
	Names = Array{String}(undef, length(compare[1].Name))
	p_ΔtR1s = Plots.plot(xlabel="modulation shift in s", ylabel="ΔtR1 in s")
	p_ΔtR2s = Plots.plot(xlabel="modulation shift in s", ylabel="ΔtR2 in s")
	for j=1:length(compare[1].Name)
		for i=1:length(shift)
			ΔtR1s[i,j] = compare[i].ΔtR1[j]
			ΔtR2s[i,j] = compare[i].ΔtR2[j]
			
		end
		Names[j] = compare[1].Name[j]
		Plots.plot!(p_ΔtR1s, shift, ΔtR1s[:,j], shape=:diamond, label=Names[j])
		Plots.plot!(p_ΔtR2s, shift, ΔtR2s[:,j], shape=:diamond, label=Names[j])
	end
	Plots.plot(p_ΔtR1s, p_ΔtR2s, legend=false, size=(900,800))
end

# ╔═╡ dba04cae-048e-40e4-91bd-f7f962b083e7
findfirst(shift.==3.4)

# ╔═╡ fb7236b5-98fd-4aac-86ed-5e7118b19ed6
meas_chrom = sort!(DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/GCxGC/Simulation/PartyMixHR5.csv", header=1, silencewarnings=true)), :RT)

# ╔═╡ 8260f156-3548-44a6-bba7-ce8860b4f7c8
meas_chrom.RT[118660:118670]

# ╔═╡ 2ff72067-0026-4f63-8640-59a911a968bc
ii = 118661#findall(meas_chrom.RT.≈20.6308)

# ╔═╡ 6d16fcfa-d922-4559-80cb-de166525cc5d
function massspectra(meas_chrom, i)
	sig = Array(meas_chrom[i, :])[3:end]
	mz = 0.0:1.0:(length(sig)-1)
	i_sig = findall(sig.>0.0)
	return DataFrame(mz=mz[i_sig], abundance=sig[i_sig])
end

# ╔═╡ 24053468-95d2-4fb2-8254-7c4aebcc7750
mz = massspectra(meas_chrom, ii)

# ╔═╡ f7660784-a782-4253-bcd6-15bed8838487
Plots.bar(mz.mz, mz.abundance)

# ╔═╡ d8c9f126-9e65-43cd-a378-d14286da6f13
function mz_filter(meas_chrom, mz)
	RT = meas_chrom.RT
	mz_trace = meas_chrom[!, mz+3]
	return DataFrame(RT=RT, mz_trace=mz_trace)
end

# ╔═╡ 02235d58-147e-4fab-92b3-60c6b1e5b736
begin
	mz105 = mz_filter(meas_chrom, 105)
	mz120 = mz_filter(meas_chrom, 120)
	mz162 = mz_filter(meas_chrom, 162)
	Plots.plot(mz105.RT.*60.0, mz105.mz_trace, label="mz=105")
	Plots.plot!(mz120.RT.*60.0, mz120.mz_trace, label="mz=120")
	Plots.plot!(mz162.RT.*60.0, mz162.mz_trace, label="mz=162")
end

# ╔═╡ a22cb9e0-5fd8-4172-8b75-4a13bbc481e0
Symbol("mz$(ii)_trace")

# ╔═╡ 42901987-35b2-4690-a8ec-6b865ac6166e
begin
	p_chrom_meas = Plots.plot(meas_chrom.RT.*60, meas_chrom.TIC)

	t_sum, c_sum = GasChromatographySystems.chrom(sim[i_select][2][1][end])[2:3]
	Plots.plot!(p_chrom_meas, t_sum, c_sum.*1e4, label="sim")

	# mark selected point for mass spectra
	Plots.scatter!(p_chrom_meas, (meas_chrom.RT[ii].*60, meas_chrom.TIC[ii]), shape=:diamond, label="mz")
end

# ╔═╡ c8809c64-33d3-4371-b5b5-6ee41baf2a3f
shift

# ╔═╡ b28b299c-1c2f-4638-a91f-a970c92e4f5a
begin
	c_slices_m, t_D1_m, t_D2_m = GasChromatographySystems.chrom_slicing(meas_chrom.RT.*60, meas_chrom.TIC, PM)

	# interpolation in tR2 for the same time basis
	itps = Array{Interpolations.Extrapolation}(undef, length(c_slices_m))
	for i=1:length(c_slices_m)
		itps[i] = LinearInterpolation((t_D2_m[i],), c_slices_m[i], extrapolation_bc=Flat())
	end
	chrom_mat = Array{Float64}(undef, length(c_slices_m), length(0.0:0.01:4.0))
	for i=1:length(c_slices_m)
		chrom_mat[i,:] = itps[i].(0.0:0.01:4.0)
	end
	chrom_mat
end;

# ╔═╡ 847d5b23-967f-4316-941f-fde283200b36
begin
	plotly()
	chrom_mat_mod = Array{Float64}(undef, size(chrom_mat)[1], size(chrom_mat)[2])
	for j=1:size(chrom_mat)[2]
		for i=1:size(chrom_mat)[1]
			if chrom_mat[i,j] > 25000.0
				chrom_mat_mod[i,j] = 25000.0
			else
				chrom_mat_mod[i,j] = chrom_mat[i,j]
			end
		end
	end
	#gr()
	p_meas = Plots.heatmap(0.0:4.0:(size(chrom_mat)[1]-1)*4.0, 0.0:0.01:4.0, chrom_mat_mod', c=:jet1, colorbar=false, legend=:bottomright);
	# add markers for the measured RTs
#	Plots.scatter!(meas.tR1[1:5], meas.tR2[1:5], markersize=4, c=:red, label="alcohols");# alcohols
#	Plots.scatter!(meas.tR1[6:22], meas.tR2[6:22], markersize=4, c=:orange, label="terpenes") # terpenes
#	Plots.scatter!(meas.tR1[23:28], meas.tR2[23:28], markersize=4, c=:yellow, label="phenones"); # phenones
#	Plots.scatter!(meas.tR1[29:35], meas.tR2[29:35], markersize=4, c=:lawngreen, label="ketones"); # ketones

	# simulation
	# spot marker
	Plots.scatter!(pl_GCxGC[i_select].tR1, pl_GCxGC[i_select].tR2, c=:lightblue, m=:diamond, markeralpha=1, msize=3, label="simulated", xlims=(0.0, 3500.0), ylims=(0.0, 4.0), xlabel="tR1 in s", ylabel="tR2 in s");
	# line marker
	RT1 = Tuple{Float64, Float64}[]
	RT2 = Tuple{Float64, Float64}[]
	for i=1:length(compare[i_select].tR1_meas)
		if ismissing(compare[i_select].tR1_meas[i])==false && ismissing(compare[i_select].tR1_sim[i])==false
			push!(RT1, (compare[i_select].tR1_meas[i], compare[i_select].tR1_sim[i]))
			push!(RT2, (compare[i_select].tR2_meas[i], compare[i_select].tR2_sim[i]))
		end
	end
	Plots.plot!(RT1, RT2, c=:lightblue, label="")

	p_meas
end

# ╔═╡ 954597a2-b847-4cf0-8c1c-ff4b4b4b670b
md"""
# End
"""

# ╔═╡ Cell order:
# ╠═113c6e70-2168-11ee-3f7f-2775a67bab90
# ╠═98217474-a16f-406a-83a7-17fee89c951a
# ╠═22091a27-80e1-4e98-abe7-4b9652cb832c
# ╠═1948baf0-1217-485a-8c8b-6e826bd5ef50
# ╠═71d12824-5fa7-4c51-94d0-46f6646df242
# ╟─13d4415d-3dc4-4eed-aded-c3561b77a9eb
# ╟─ddde4e3d-4e72-4315-9555-5f8b3375b04c
# ╠═55194da2-9c5b-464b-bfe1-a65c48bc7801
# ╠═e2e1ee2a-f98d-47fd-8c22-cc441f2d12c2
# ╟─09438575-8370-475f-b2df-fda5155e8c82
# ╠═d0379590-edf3-4017-b8b8-07bd6080f757
# ╠═23b71d48-26bd-496d-9c10-e2d3ce2bc13c
# ╟─0e349acc-46b4-4734-abb7-668ec1225c53
# ╠═de542d7e-bd0d-47a4-8977-5ab1e64a26f4
# ╠═a2d062a4-bd59-4ab7-9316-dc769b2c3c09
# ╟─b8a6cdec-457b-4fb0-b6f8-60f713f4a015
# ╠═e72a6c41-2b14-45cb-81e7-263610c3c8ae
# ╟─daf7484b-ac12-4515-8e26-adbc8e214cbc
# ╠═83d10e79-3ab6-418c-80ba-653a50fa5d21
# ╠═6afa596b-70dd-4b01-899b-5d9ef1349491
# ╠═7ddfe44e-fde9-422f-98de-27b80e8517be
# ╟─5dd6d2b7-1d6f-49ef-b684-a453a50a9e96
# ╠═5883701b-96f6-431e-b29b-b5cc0be940e0
# ╟─03b1cae0-cf56-4e0f-b0ad-fdfa736e0fa8
# ╠═397d7248-e069-4958-849d-6761cb20283d
# ╠═6b2c6634-c13a-4b71-93ff-d66ab82854c9
# ╠═35bbe26a-f82c-4f60-a2d0-74b052037ae7
# ╠═7580ebee-b061-4c11-aaad-0e8310e5277d
# ╠═b35d871c-4f74-4290-95a4-24ab7115c591
# ╠═262fba99-1b75-4803-8f34-a0619de559a9
# ╠═510bf59c-a99e-4b30-a82e-4bbc1b336b7b
# ╟─1964abf5-a627-4014-8348-e81145647104
# ╠═5f3c1ce1-50bb-44ba-abac-d3ccc82fdb51
# ╠═e0f0bb8f-1be0-4d84-9e0d-d86d3df65d6b
# ╠═b284207c-70b5-4da5-94e9-0fe70680b0e3
# ╠═2b41aa34-d994-440b-a8e1-dd6337008a48
# ╟─bbb75e0b-ccab-4c32-9e33-0d8915dba95d
# ╠═8572809a-0ea8-4520-8768-6d7f3288d06f
# ╟─b586403d-12ac-4b38-b4de-f974fe71b5dc
# ╠═7e945bad-0a59-4a5e-8c63-cef21ee847d0
# ╠═d912d940-5b39-4bb9-a136-1f5e67ad7d82
# ╠═f1f1ffdb-bc7e-48fd-b688-21f375ab21eb
# ╟─7739a83c-d001-47cd-aee2-36917106c2ad
# ╠═582499db-7b5a-4054-92c2-d0e2a97cd91f
# ╠═5f67ab8d-3b3d-4bbd-81c8-04eec158ef6e
# ╟─f268df70-af5a-482f-85e6-84216681ecf4
# ╠═e5d4fb3a-f05f-4533-9666-08be32ef19ef
# ╟─271bb7a6-9a98-4336-9571-c131cc125486
# ╠═b7c2070b-16bf-41f5-b021-160e54a53a21
# ╠═189034b6-09e5-4e6b-80b3-e442ad30705c
# ╠═3ed4ad21-c8e4-428a-b6d3-b37cfdbd0778
# ╟─aabb8d13-85a8-4252-9c71-99b78dd9e2cf
# ╠═76236614-3c78-4c60-9eda-d906337362e7
# ╠═dba04cae-048e-40e4-91bd-f7f962b083e7
# ╠═fb7236b5-98fd-4aac-86ed-5e7118b19ed6
# ╠═8260f156-3548-44a6-bba7-ce8860b4f7c8
# ╠═2ff72067-0026-4f63-8640-59a911a968bc
# ╠═24053468-95d2-4fb2-8254-7c4aebcc7750
# ╠═f7660784-a782-4253-bcd6-15bed8838487
# ╠═6d16fcfa-d922-4559-80cb-de166525cc5d
# ╠═d8c9f126-9e65-43cd-a378-d14286da6f13
# ╠═02235d58-147e-4fab-92b3-60c6b1e5b736
# ╠═a22cb9e0-5fd8-4172-8b75-4a13bbc481e0
# ╠═42901987-35b2-4690-a8ec-6b865ac6166e
# ╠═c8809c64-33d3-4371-b5b5-6ee41baf2a3f
# ╠═b28b299c-1c2f-4638-a91f-a970c92e4f5a
# ╠═847d5b23-967f-4316-941f-fde283200b36
# ╠═954597a2-b847-4cf0-8c1c-ff4b4b4b670b
