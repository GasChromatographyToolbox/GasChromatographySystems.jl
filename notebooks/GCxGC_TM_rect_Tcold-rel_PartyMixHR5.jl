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

- rectangle modulation point (`spatial=false`, `tflank=Inf`, `algTM="simplifiedTM"` and `ng=true`)
- Tcold is relative (`Tcold_abs=false`)
"""

# ╔═╡ 1948baf0-1217-485a-8c8b-6e826bd5ef50
optTM = GasChromatographySystems.ModuleTMOptions(Tcold_abs=false, ng=true, alg="simplifiedTM", tflank=Inf, sflank=Inf)

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
	shift = 3.5#(PM-0.35)#0.0#-0.35*2
	thot = 0.35
	ratio = (PM-thot)/PM
	Thot = 25.0
	Tcold = -50.0
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
	sys = GasChromatographySystems.GCxGC_TM(L1, d1, df1, sp1, TP1, L2, d2, df2, sp2, TP2, LTL, dTL, dfTL, spTL, TTL, LM, dM, dfM, spM, shift, PM, ratio, Thot, Tcold, TPM, F, pin, pout; optTM=optTM)
end

# ╔═╡ 6afa596b-70dd-4b01-899b-5d9ef1349491
GasChromatographySystems.plot_flow_over_time(sys; dt=30.0)

# ╔═╡ 7ddfe44e-fde9-422f-98de-27b80e8517be
GasChromatographySystems.plot_pressure_over_time(sys; dt=30.0)

# ╔═╡ d1eca540-5621-4928-9d74-ab0e0e3ca22a
optTM

# ╔═╡ 8e365f1d-708e-4e4d-abb4-b186f901f169
TP1

# ╔═╡ 13b47494-1183-4c8f-b66f-9c4aaa0f2a3b
1800/60

# ╔═╡ c193bbda-82bb-459a-8677-d907d5cdd75f
([0.0, 120.0, 1800.0, 0.0, 960.0, 600.0], [143.43, 143.43, 396.33, 396.33, 593.86, 593.86])

# ╔═╡ 51ac4361-66b7-4f3b-88fb-f218f98e9073


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
selected_solutes_ = GasChromatographySystems.common_solutes(db, sys).Name

# ╔═╡ 262fba99-1b75-4803-8f34-a0619de559a9
selected_solutes = selected_solutes_[[33, 34, 35, 36, 37, 38, 80, 82, 83, 84, 86, 127, 129, 130, 131, 132]]#[Not([1, 29, 30])]#[[14,16,17]]

# ╔═╡ 510bf59c-a99e-4b30-a82e-4bbc1b336b7b
# 85 -> multiple entries in database for one of the stat. phases

# ╔═╡ 1964abf5-a627-4014-8348-e81145647104
md"""
## Simulation
"""

# ╔═╡ 5f3c1ce1-50bb-44ba-abac-d3ccc82fdb51
par = GasChromatographySystems.graph_to_parameters(sys, db, selected_solutes)

# ╔═╡ e2e1ee2a-f98d-47fd-8c22-cc441f2d12c2
begin
	Plots.plot(cumsum(tsteps_), par[1].prog.T_itp(0.0, cumsum(tsteps_)).-273.15, label="")
end

# ╔═╡ a2d062a4-bd59-4ab7-9316-dc769b2c3c09
begin
	TM_itp = par[3].prog.T_itp
	Plots.plot(0.0:0.01:sys.modules[3].PM, TM_itp.(sys.modules[3].L/2, 0.0:0.01:sys.modules[3].PM).-273.115, title="temperature at modulator point", xlabel="time in s", ylabel="temperature in °C", label="")
	Plots.plot!([-shift, -shift], [TM_itp(sys.modules[3].L/2,-shift), TM_itp(sys.modules[3].L/2,(1+ratio)/2*PM-shift)].-273.15, c=:black, label="")
	Plots.plot!([ratio*PM-shift, ratio*PM-shift], [TM_itp(sys.modules[3].L/2,ratio*PM-shift-eps()), TM_itp(sys.modules[3].L/2,(1+ratio)/2*PM-shift)].-273.15, c=:black, linestyle=:dash, label="")
	Plots.plot!([PM-shift, PM-shift], [TM_itp(sys.modules[3].L/2,PM-shift), TM_itp(sys.modules[3].L/2,(1+ratio)/2*PM-shift)].-273.15, c=:black, label="")
end

# ╔═╡ b8a6cdec-457b-4fb0-b6f8-60f713f4a015
begin
	Plots.plot(xlabel="position x in m", ylabel="temperature in °C", title="temperature at modulator point", legend=false)
	Plots.plot!(0.0:sys.modules[3].L/100.0:sys.modules[3].L, TM_itp.(0.0:sys.modules[3].L/100.0:sys.modules[3].L, sys.modules[3].ratio*sys.modules[3].PM-sys.modules[3].shift).-273.15, label="t=$(sys.modules[3].ratio*sys.modules[3].PM-sys.modules[3].shift)s")
	Plots.plot!(0.0:sys.modules[3].L/100.0:sys.modules[3].L, TM_itp.(0.0:sys.modules[3].L/100.0:sys.modules[3].L, (1+3*sys.modules[3].ratio)/4*sys.modules[3].PM-sys.modules[3].shift).-273.15, label="t=$((1+3*sys.modules[3].ratio)/4*sys.modules[3].PM-sys.modules[3].shift)s")
	Plots.plot!(0.0:sys.modules[3].L/100.0:sys.modules[3].L, TM_itp.(0.0:sys.modules[3].L/100.0:sys.modules[3].L, (1+sys.modules[3].ratio)/2*sys.modules[3].PM-sys.modules[3].shift).-273.15, label="t=$((1+sys.modules[3].ratio)/2*sys.modules[3].PM-sys.modules[3].shift)s")
	Plots.plot!(0.0:sys.modules[3].L/100.0:sys.modules[3].L, TM_itp.(0.0:sys.modules[3].L/100.0:sys.modules[3].L, 3.95-sys.modules[3].shift).-273.15, label="t=$(3.95-sys.modules[3].shift)s")
	Plots.plot!(0.0:sys.modules[3].L/100.0:sys.modules[3].L, TM_itp.(0.0:sys.modules[3].L/100.0:sys.modules[3].L, sys.modules[3].PM-sys.modules[3].shift).-273.15, label="t=$(sys.modules[3].PM-sys.modules[3].shift)s")
end

# ╔═╡ e72a6c41-2b14-45cb-81e7-263610c3c8ae
begin
	Plots.plot(0.0:0.01:(sys.modules[3].PM*100.0), TM_itp.(sys.modules[3].L/2, 0.0:0.01:(sys.modules[3].PM*100.0)).-273.115, title="temperature at modulator point", xlabel="time in s", ylabel="temperature in °C", label="")
	Plots.plot!(0.0:0.01:(sys.modules[3].PM*100.0), par[1].prog.T_itp(0.0, 0.0:0.01:(sys.modules[3].PM*100.0)).-273.15, label="")
end

# ╔═╡ b284207c-70b5-4da5-94e9-0fe70680b0e3
paths = GasChromatographySystems.all_paths(sys.g, 1)[2]

# ╔═╡ 2b41aa34-d994-440b-a8e1-dd6337008a48
sim = GasChromatographySystems.simulate_along_paths(sys, paths, par)

# ╔═╡ 6b7f9f70-2391-4797-ada2-9da8f22c411e
sim[2][1][2]

# ╔═╡ 122f3335-653f-4abd-b9a9-81d57fa12d64
new_par_sys, df_A = GasChromatographySystems.slicing(sim[2][1][2], sys.modules[3].PM, sys.modules[3].ratio, sys.modules[3].shift, par[3]; abstol=sys.modules[3].opt.abstol, reltol=sys.modules[3].opt.reltol, alg=sys.modules[3].opt.alg)

# ╔═╡ 44be951b-d837-413a-970e-c295a89a87e8
df_A

# ╔═╡ 6d91a65e-204b-47eb-94c6-2ca933fa94d1
sim[2][1][2].tR

# ╔═╡ 1fc5c501-b9fb-43b6-a757-9502fa2a55b5
t0 = fld.(df_A.t0 .+ shift, PM).*PM .- shift

# ╔═╡ 8dc288a9-1f2e-4e9e-a833-d165b55c9843
tR = ceil.(t0 .+ PM*ratio, digits=3)

# ╔═╡ c463a51f-873e-4661-842c-ece71bfc4bd2
TR = new_par_sys.prog.T_itp.(new_par_sys.col.L, tR).-273.15

# ╔═╡ 44bc71d6-7049-47fa-a05a-5357c4f7e7ab
tR[10] 

# ╔═╡ f1595846-2492-4e08-a673-eb8e0403be94
tmod = ceil(mod(tR[25]+0.35+shift, PM), digits=4)

# ╔═╡ 11dfd391-4377-4012-b56a-f5ccdcdc28f8
tmod < PM*ratio

# ╔═╡ 2d487515-a2b3-4c4a-8698-e63b21eb70e9
ifelse(tmod < PM*ratio, Tcold, Thot)

# ╔═╡ bb382d48-f6cb-42ee-87fa-9330e68179aa
mod(tR[9]+shift, PM)

# ╔═╡ 54525c7e-4cb9-4fd5-8df8-1445150c0c09
mod(tR[10]+shift, PM)

# ╔═╡ b2fa6bd9-f381-4dfb-ba34-d8a426a29914
tR[9]+shift - PM*fld(tR[9]+shift, PM)

# ╔═╡ ac7ddeab-58a5-419c-869e-d23c4d9c7b46
tR[10]+shift - PM*fld(tR[10]+shift, PM)

# ╔═╡ 605bab67-6413-4d7a-8888-9f52b91c3fc0
rem(shift,PM,RoundUp)

# ╔═╡ b6e656e2-33a7-47e1-96fa-acff3e8374e0
mod(tR[9], PM) +shift

# ╔═╡ b1908b5a-0cff-44a0-8605-87029e7e198f
mod(tR[10], PM) +shift

# ╔═╡ 16e68147-0159-418f-8c2d-c3bc4e437de5
tR[10]-548.3

# ╔═╡ 88988bb8-d9e7-4750-95a7-d99f930a5d86
t0[10] + PM*ratio

# ╔═╡ 0bc2c8e4-8fe1-414a-8517-ae51b2b74f6b
ceil(tR[7], digits=2)

# ╔═╡ 9a3a7bcf-27ce-4914-aed9-59949837cedb
sort(df_A, :t0).Name

# ╔═╡ 339e4504-a0f7-4d26-b69e-c1ab874bc3ae
new_par_sys.sub

# ╔═╡ 623bf815-957f-4c29-8bfb-aa2a811c31e7
function rect(x, x0, min, max, width)
	if abs(x-x0) > width/2
		min
	elseif abs(x-x0) == width/2
		(min+max)/2
	elseif abs(x-x0) < width/2
		max
	end
end

# ╔═╡ 3eb664ff-6b09-4b9e-86e5-6979e2629b81
function therm_mod(t, shift, PM, ratio, Tcold, Thot; flank=20) 
	# add warning, if flank value is to low -> jumps in the function
	width = (1-ratio)*PM
	tmod = mod(t+shift, PM)
	tstart = ratio*PM
	if flank == Inf # rectangle function
		# ceil() to assure, that the rounding errors do not affect the rectangle function 
		return rect(ceil(tmod, digits=4), (1+ratio)*PM/2, Tcold, Thot, (1-ratio)*PM)
	else # smoothed rectangle
		return smooth_rectangle.(tmod, tstart, width, Tcold, Thot; flank=flank) 
	end
end

# ╔═╡ d69f07bd-83c6-4364-9a7e-11e1afc1e729
therm_mod(tR[23], shift, PM, ratio, Tcold, Thot; flank=Inf)

# ╔═╡ 627bfbc2-97f8-4875-bdc6-53384730f00e
Plots.plot(0.0:0.001:20.0, therm_mod.(0.0:0.001:20.0, shift, PM, ratio, Tcold, Thot; flank=Inf))

# ╔═╡ 8e765549-8c8b-4a21-a747-6c17a86a298f
(Tcold+Thot)/2

# ╔═╡ 3b885c67-d52d-4395-8885-c4aebd893b1b
Plots.plot(-4.0:0.001:4.1, rect.((-4.0:0.001:4.1), (1+ratio)*PM/2, Tcold, Thot, (1-ratio)*PM))

# ╔═╡ d5e0d4ee-5f59-42ec-979d-c6d4213a14ac
(1+ratio)*PM/2

# ╔═╡ da553041-fc23-44fd-b917-9522b626a10e
(3-ratio)*PM/2

# ╔═╡ f92d4d80-b1cb-4e3f-b658-a53786e28583
begin
	sort_df_A = sort(df_A, :t0)
	Name = sort_df_A.Name
	kR = Array{Float64}(undef, length(tR))
	uR = Array{Float64}(undef, length(tR))
	τ₀ = Array{Float64}(undef, length(tR))
	for i=1:length(tR)
		i_sub = findfirst(Name[i] .== df_A.Name)
		kR[i] = GasChromatographySimulator.retention_factor(new_par_sys.col.L, tR[i], new_par_sys.prog.T_itp, new_par_sys.col.d, new_par_sys.col.df, new_par_sys.sub[i_sub].Tchar, new_par_sys.sub[i_sub].θchar, new_par_sys.sub[i_sub].ΔCp, new_par_sys.sub[i_sub].φ₀)

		rM = GasChromatographySimulator.mobile_phase_residency(new_par_sys.col.L, tR[i], new_par_sys.prog.T_itp, new_par_sys.prog.Fpin_itp, new_par_sys.prog.pout_itp, new_par_sys.col.L, new_par_sys.col.d, new_par_sys.col.gas; ng=new_par_sys.opt.ng, vis=new_par_sys.opt.vis, control=new_par_sys.opt.control)
		uR[i] = 1/(rM*(1+kR[i]))

		τ₀[i] = new_par_sys.sub[i_sub].τ₀
	end
end

# ╔═╡ b98d9e89-d3a6-4386-9115-2facf2134638
kR

# ╔═╡ 391808f8-1e03-4407-b3e4-d7031bf2f7f3
uR

# ╔═╡ 878eaf01-87ca-4c7a-9265-92953459965c
τ₀

# ╔═╡ 933a6d8d-c4a9-4d18-9c6e-55771491da9d
τR = new_par_sys.col.L./uR

# ╔═╡ 98df4471-8cdd-4be0-9116-20f6ac8ce1fe
begin
	Plots.plot(540.0:0.01:560.0, TM_itp.(sys.modules[3].L/2, 540.0:0.01:560.0).-273.115, title="temperature at modulator point", xlabel="time in s", ylabel="temperature in °C", label="")
	Plots.plot!(540.0:0.01:560.0, par[1].prog.T_itp(0.0, 540.0:0.01:560.0).-273.15, label="")
end

# ╔═╡ bbb75e0b-ccab-4c32-9e33-0d8915dba95d
md"""
## Chromatograms
"""

# ╔═╡ b586403d-12ac-4b38-b4de-f974fe71b5dc
md"""
### 1D Chromatograms
"""

# ╔═╡ 7e945bad-0a59-4a5e-8c63-cef21ee847d0
#collect_chrom(sim[2][1], sys; markings=true)

# ╔═╡ d912d940-5b39-4bb9-a136-1f5e67ad7d82
GasChromatographySystems.chrom(sim[2][1][end])[1]

# ╔═╡ f1f1ffdb-bc7e-48fd-b688-21f375ab21eb
#GasChromatographySystems.chrom_marked(sim[2][1][4], PM, ratio, shift1; nτ=6)[1]

# ╔═╡ 7739a83c-d001-47cd-aee2-36917106c2ad
md"""
### 2D Chromatogram
"""

# ╔═╡ 582499db-7b5a-4054-92c2-d0e2a97cd91f
slice_mat, t_D1, t_D2, c_slices, t_, chrom_sliced_sum  = GasChromatographySystems.chrom2d(sim[2][1][end], sys, PM)

# ╔═╡ f268df70-af5a-482f-85e6-84216681ecf4
md"""
## Peaklist
"""

# ╔═╡ e5d4fb3a-f05f-4533-9666-08be32ef19ef
pl_GCxGC = GasChromatographySystems.peaklist_GCxGC(sim[2][1][8], PM)

# ╔═╡ 5f67ab8d-3b3d-4bbd-81c8-04eec158ef6e
begin
	plotly()
	x1 = floor(minimum(pl_GCxGC.tR1*0.99))
	x2 = ceil(maximum(pl_GCxGC.tR1*1.01))
	y1 = floor(minimum(pl_GCxGC.tR2*0.99))
	y2 = ceil(maximum(pl_GCxGC.tR2*1.01))
	p_contour = Plots.contour(t_D1[1:end-1], t_D2[1], slice_mat[1:end-1,1:end]', levels=10, xlabel="tR1 in s", ylabel="tR2 in s")#, c=:jet1)
	Plots.scatter!(p_contour, pl_GCxGC.tR1, pl_GCxGC.tR2, msize=3, shape=:diamond, c=:black, markeralpha=1)
	Plots.plot!(p_contour, xlims=(x1, x2), ylims=(y1, y2))

	p_heatmap = Plots.heatmap(t_D1[1:end-1], t_D2[1], slice_mat[1:end-1,1:end]', levels=10, xlabel="tR1 in s", ylabel="tR2 in s", c=:jet1)
	Plots.scatter!(p_heatmap, pl_GCxGC.tR1, pl_GCxGC.tR2, msize=3, shape=:diamond, c=:black, markeralpha=1)
	Plots.plot!(p_heatmap, xlims=(x1, x2), ylims=(y1, y2))

	Plots.plot(p_contour, p_heatmap, size=(900,800), legend=false)
end

# ╔═╡ ab0c7a5a-7b6b-4b01-a283-dc4b157092ef
GasChromatographySystems.peaklist_GCxGC(sim[2][1][8], sim[2][1][2], PM,)

# ╔═╡ 15547cdd-7c1c-49a7-9c4a-3ecca725c231
Int(4.0)

# ╔═╡ a2a3c1a8-ccf2-449a-8ee2-265a4270d3f6
Int(4.1)

# ╔═╡ c126a17e-8ee3-4e8b-b43d-036a07768795
round(Int, 3.6)

# ╔═╡ 93e2068d-760e-403e-a6fa-984a237ad6c3
md"""
## Checks
"""

# ╔═╡ 9bc833db-7538-4520-9a80-ab49f9bbeb81
GasChromatographySystems.check_peakwidths(sim[2][1][5])

# ╔═╡ b872f155-d36f-4a2d-9741-342a108d6b5f
GasChromatographySystems.check_duration_modulation(sim[2][1], sim[4], PM, ratio)

# ╔═╡ daed9ba6-a44c-48a0-999b-5862838253c2
GasChromatographySystems.check_area(sim[2][1][8])

# ╔═╡ 27645404-48fd-488c-8b2e-66300f4b3303
md"""
## Traces 1st Modulation
"""

# ╔═╡ c20db2f4-b5b8-4436-9d4b-df1333659d25
#=begin
	plotly()
	p_xt_check = Plots.plot(xlabel="time in s", ylabel="x in m", legend=true)
	p_τt_check = Plots.plot(xlabel="time in s", ylabel="τ in s", legend=true)
	p_Tt_check = Plots.plot(xlabel="time in s", ylabel="T in °C", legend=true)
	p_lnkt_check = Plots.plot(xlabel="time in s", ylabel="lnk", legend=true)
	first_step_TM1_z = Array{Float64}(undef, length(sim[3][1][3]))
	first_step_TM1_t = Array{Float64}(undef, length(sim[3][1][3]))
	for i=1:length(sim[3][1][3])
		trace = GasChromatographySystems.traces(sim[3][1][3], sim[4][3], i)
		lbl = string(split(sim[4][3].sub[i].ann, "_")[1], "_", sim[4][3].sub[i].name)
		Plots.plot!(p_xt_check, trace.t, trace.z, label=lbl, marker=:circle)
		Plots.plot!(p_τt_check, trace.t, sqrt.(trace.τ²), label=lbl, marker=:circle)
		Plots.plot!(p_Tt_check, trace.t, trace.T.-273.15, label=lbl, marker=:circle)
		Plots.plot!(p_lnkt_check, trace.t, log.(trace.k), label=lbl, marker=:circle)
		first_step_TM1_z[i] = trace.z[2] - trace.z[1]
		first_step_TM1_t[i] = trace.t[2] - trace.t[1]
	end
	Plots.plot(p_xt_check, p_τt_check, p_Tt_check, p_lnkt_check, size=(1000,800))
end=#

# ╔═╡ d4c0a1c4-68ef-408a-9919-1bd737b6c229
md"""
## Traces 2nd Modulation
"""

# ╔═╡ 2737c1e7-6239-41d7-8868-242bb881f436
#=begin
	plotly()
	p_xt_check_ = Plots.plot(xlabel="time in s", ylabel="x in m", legend=true)
	p_τt_check_ = Plots.plot(xlabel="time in s", ylabel="τ in s", legend=true)
	p_Tt_check_ = Plots.plot(xlabel="time in s", ylabel="T in °C", legend=true)
	p_lnkt_check_ = Plots.plot(xlabel="time in s", ylabel="lnk", legend=true)
	first_step_TM2_z = Array{Float64}(undef, length(sim[3][1][5]))
	first_step_TM2_t = Array{Float64}(undef, length(sim[3][1][5]))
	for i=1:length(sim[3][1][3])
		trace = GasChromatographySystems.traces(sim[3][1][5], sim[4][5], i)
		lbl = string(join(split(sim[4][5].sub[i].ann, "_")[1:2], "_"), "_", sim[4][5].sub[i].name)
		Plots.plot!(p_xt_check_, trace.t, trace.z, label=lbl, marker=:circle)
		Plots.plot!(p_τt_check_, trace.t, sqrt.(trace.τ²), label=lbl, marker=:circle)
		Plots.plot!(p_Tt_check_, trace.t, trace.T.-273.15, label=lbl, marker=:circle)
		Plots.plot!(p_lnkt_check_, trace.t, log.(trace.k), label=lbl, marker=:circle)
		first_step_TM2_z[i] = trace.z[2] - trace.z[1]
		first_step_TM2_t[i] = trace.t[2] - trace.t[1]
	end
	Plots.plot(p_xt_check_, p_τt_check_, p_Tt_check_, p_lnkt_check_, size=(1000,800))
end=#

# ╔═╡ 271bb7a6-9a98-4336-9571-c131cc125486
md"""
## Measurement
"""

# ╔═╡ b7c2070b-16bf-41f5-b021-160e54a53a21
meas = DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/GCxGC/Simulation/PartyMixHR5_RT.csv", header=1, silencewarnings=true, stringtype=String))

# ╔═╡ 189034b6-09e5-4e6b-80b3-e442ad30705c
pl_GCxGC

# ╔═╡ 2143e6c5-b407-47b9-a4c4-674ac51384fd
pl_GCxGC.tR1./60

# ╔═╡ 3ed4ad21-c8e4-428a-b6d3-b37cfdbd0778
compare = GasChromatographySystems.comparison_meas_sim(meas, pl_GCxGC)

# ╔═╡ e17c92b0-2db3-4e29-b56f-c7b5650c9d95
avreltR1 = sum(abs.(collect(skipmissing(compare.relΔtR1_percent))))/length(collect(skipmissing(compare.relΔtR1_percent)))

# ╔═╡ c190dd0f-6f09-4d18-9aa5-6b40beddfbf3
avreltR2 = sum(abs.(collect(skipmissing(compare.relΔtR2_percent))))/length(collect(skipmissing(compare.relΔtR2_percent)))

# ╔═╡ fb7236b5-98fd-4aac-86ed-5e7118b19ed6
meas_chrom = sort!(DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/GCxGC/Simulation/PartyMixHR5.csv", header=1, silencewarnings=true)), :RT)

# ╔═╡ 5bbe0038-c95b-4ca3-a114-7134ae5beb23
1237.85/60

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
	Plots.scatter!(p_chrom_meas, 468.0:4.0:500.0, 1000.0.*ones(length(468.0:4.0:500.0)), label="no shift")
	Plots.scatter!(p_chrom_meas, (468.0:4.0:500.0).+shift, 1000.0.*ones(length(468.0:4.0:500.0)), label="+shift")
	Plots.scatter!(p_chrom_meas, (468.0:4.0:500.0).-shift, 1000.0.*ones(length(468.0:4.0:500.0)), label="-shift")

	t_sum, c_sum = GasChromatographySystems.chrom(sim[2][1][end])[2:3]
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
	Plots.scatter!(pl_GCxGC.tR1, pl_GCxGC.tR2, c=:lightblue, m=:diamond, markeralpha=1, msize=3, label="simulated", xlims=(0.0, 3500.0), ylims=(0.0, 4.0), xlabel="tR1 in s", ylabel="tR2 in s");
	# line marker
	RT1 = Tuple{Float64, Float64}[]
	RT2 = Tuple{Float64, Float64}[]
	for i=1:length(compare.tR1_meas)
		if ismissing(compare.tR1_meas[i])==false && ismissing(compare.tR1_sim[i])==false
			push!(RT1, (compare.tR1_meas[i], compare.tR1_sim[i]))
			push!(RT2, (compare.tR2_meas[i], compare.tR2_sim[i]))
		end
	end
	Plots.plot!(RT1, RT2, c=:black, linewidth=2, label="")

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
# ╠═b8a6cdec-457b-4fb0-b6f8-60f713f4a015
# ╠═e72a6c41-2b14-45cb-81e7-263610c3c8ae
# ╟─daf7484b-ac12-4515-8e26-adbc8e214cbc
# ╠═83d10e79-3ab6-418c-80ba-653a50fa5d21
# ╠═6afa596b-70dd-4b01-899b-5d9ef1349491
# ╠═7ddfe44e-fde9-422f-98de-27b80e8517be
# ╟─5dd6d2b7-1d6f-49ef-b684-a453a50a9e96
# ╠═5883701b-96f6-431e-b29b-b5cc0be940e0
# ╠═d1eca540-5621-4928-9d74-ab0e0e3ca22a
# ╠═8e365f1d-708e-4e4d-abb4-b186f901f169
# ╠═13b47494-1183-4c8f-b66f-9c4aaa0f2a3b
# ╠═c193bbda-82bb-459a-8677-d907d5cdd75f
# ╠═51ac4361-66b7-4f3b-88fb-f218f98e9073
# ╟─03b1cae0-cf56-4e0f-b0ad-fdfa736e0fa8
# ╠═397d7248-e069-4958-849d-6761cb20283d
# ╠═6b2c6634-c13a-4b71-93ff-d66ab82854c9
# ╠═35bbe26a-f82c-4f60-a2d0-74b052037ae7
# ╠═262fba99-1b75-4803-8f34-a0619de559a9
# ╠═510bf59c-a99e-4b30-a82e-4bbc1b336b7b
# ╟─1964abf5-a627-4014-8348-e81145647104
# ╠═5f3c1ce1-50bb-44ba-abac-d3ccc82fdb51
# ╠═b284207c-70b5-4da5-94e9-0fe70680b0e3
# ╠═2b41aa34-d994-440b-a8e1-dd6337008a48
# ╠═6b7f9f70-2391-4797-ada2-9da8f22c411e
# ╠═122f3335-653f-4abd-b9a9-81d57fa12d64
# ╠═44be951b-d837-413a-970e-c295a89a87e8
# ╠═6d91a65e-204b-47eb-94c6-2ca933fa94d1
# ╠═1fc5c501-b9fb-43b6-a757-9502fa2a55b5
# ╠═8dc288a9-1f2e-4e9e-a833-d165b55c9843
# ╠═c463a51f-873e-4661-842c-ece71bfc4bd2
# ╠═44bc71d6-7049-47fa-a05a-5357c4f7e7ab
# ╠═f1595846-2492-4e08-a673-eb8e0403be94
# ╠═11dfd391-4377-4012-b56a-f5ccdcdc28f8
# ╠═2d487515-a2b3-4c4a-8698-e63b21eb70e9
# ╠═d69f07bd-83c6-4364-9a7e-11e1afc1e729
# ╠═bb382d48-f6cb-42ee-87fa-9330e68179aa
# ╠═54525c7e-4cb9-4fd5-8df8-1445150c0c09
# ╠═b2fa6bd9-f381-4dfb-ba34-d8a426a29914
# ╠═ac7ddeab-58a5-419c-869e-d23c4d9c7b46
# ╠═605bab67-6413-4d7a-8888-9f52b91c3fc0
# ╠═b6e656e2-33a7-47e1-96fa-acff3e8374e0
# ╠═b1908b5a-0cff-44a0-8605-87029e7e198f
# ╠═3eb664ff-6b09-4b9e-86e5-6979e2629b81
# ╠═627bfbc2-97f8-4875-bdc6-53384730f00e
# ╠═16e68147-0159-418f-8c2d-c3bc4e437de5
# ╠═88988bb8-d9e7-4750-95a7-d99f930a5d86
# ╠═0bc2c8e4-8fe1-414a-8517-ae51b2b74f6b
# ╠═9a3a7bcf-27ce-4914-aed9-59949837cedb
# ╠═339e4504-a0f7-4d26-b69e-c1ab874bc3ae
# ╠═623bf815-957f-4c29-8bfb-aa2a811c31e7
# ╠═8e765549-8c8b-4a21-a747-6c17a86a298f
# ╠═3b885c67-d52d-4395-8885-c4aebd893b1b
# ╠═d5e0d4ee-5f59-42ec-979d-c6d4213a14ac
# ╠═da553041-fc23-44fd-b917-9522b626a10e
# ╠═f92d4d80-b1cb-4e3f-b658-a53786e28583
# ╠═b98d9e89-d3a6-4386-9115-2facf2134638
# ╠═391808f8-1e03-4407-b3e4-d7031bf2f7f3
# ╠═878eaf01-87ca-4c7a-9265-92953459965c
# ╠═933a6d8d-c4a9-4d18-9c6e-55771491da9d
# ╠═98df4471-8cdd-4be0-9116-20f6ac8ce1fe
# ╟─bbb75e0b-ccab-4c32-9e33-0d8915dba95d
# ╟─b586403d-12ac-4b38-b4de-f974fe71b5dc
# ╠═7e945bad-0a59-4a5e-8c63-cef21ee847d0
# ╠═d912d940-5b39-4bb9-a136-1f5e67ad7d82
# ╠═f1f1ffdb-bc7e-48fd-b688-21f375ab21eb
# ╟─7739a83c-d001-47cd-aee2-36917106c2ad
# ╠═582499db-7b5a-4054-92c2-d0e2a97cd91f
# ╠═5f67ab8d-3b3d-4bbd-81c8-04eec158ef6e
# ╟─f268df70-af5a-482f-85e6-84216681ecf4
# ╠═e5d4fb3a-f05f-4533-9666-08be32ef19ef
# ╠═ab0c7a5a-7b6b-4b01-a283-dc4b157092ef
# ╠═15547cdd-7c1c-49a7-9c4a-3ecca725c231
# ╠═a2a3c1a8-ccf2-449a-8ee2-265a4270d3f6
# ╠═c126a17e-8ee3-4e8b-b43d-036a07768795
# ╟─93e2068d-760e-403e-a6fa-984a237ad6c3
# ╠═9bc833db-7538-4520-9a80-ab49f9bbeb81
# ╠═b872f155-d36f-4a2d-9741-342a108d6b5f
# ╠═daed9ba6-a44c-48a0-999b-5862838253c2
# ╟─27645404-48fd-488c-8b2e-66300f4b3303
# ╠═c20db2f4-b5b8-4436-9d4b-df1333659d25
# ╟─d4c0a1c4-68ef-408a-9919-1bd737b6c229
# ╠═2737c1e7-6239-41d7-8868-242bb881f436
# ╟─271bb7a6-9a98-4336-9571-c131cc125486
# ╠═b7c2070b-16bf-41f5-b021-160e54a53a21
# ╠═189034b6-09e5-4e6b-80b3-e442ad30705c
# ╠═2143e6c5-b407-47b9-a4c4-674ac51384fd
# ╟─3ed4ad21-c8e4-428a-b6d3-b37cfdbd0778
# ╠═e17c92b0-2db3-4e29-b56f-c7b5650c9d95
# ╠═c190dd0f-6f09-4d18-9aa5-6b40beddfbf3
# ╠═fb7236b5-98fd-4aac-86ed-5e7118b19ed6
# ╠═5bbe0038-c95b-4ca3-a114-7134ae5beb23
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
