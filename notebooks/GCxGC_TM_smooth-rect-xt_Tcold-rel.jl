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

- smooth rectangle modulation point in time and space (`spatial=true`, `tflank=20`, `algTM=Vern9()` and `ng=false`)
- Tcold is relative (`Tcold_abs=false`)
"""

# ╔═╡ 470c8573-a0f5-4ab2-9385-c256491f648d
spatial = true

# ╔═╡ b7b00dda-40d3-40d8-86e5-377d48b189fc
flank = 20

# ╔═╡ 14a1a970-3cfe-4086-94d2-d0c7ed6dab15
algTM = Vern9()

# ╔═╡ b5efd58e-666b-43c0-9b98-31085c8c36fe
ng = false

# ╔═╡ 9297b636-6daa-40e9-b312-e05f5e5da5fd
Tcold_abs = false

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
	tsteps_, Tsteps_ = GasChromatographySimulator.conventional_program([50.0, 1.0, 3.0, 200.0, 0.0, 10.0, 225.0, 10.0]) 
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
	L1 = 30.6 # m
	d1 = 0.25 # mm
	df1 = 0.25 # µm
	sp1 = "ZB1ms"
	
	# modulator
	LM = [0.3, 0.005, 0.9, 0.005, 0.3]#[0.3, 0.005, 0.9, 0.005, 0.3]
	dM = 0.1
	dfM = 0.1
	spM = "Stabilwax"

	# 2nd D
	L2 = 0.72#0.56 # m
	d2 = 0.1 # mm
	df2 = 0.1 # µm
	sp2 = "Stabilwax"

	# TL
	LTL = 0.235
	dTL = 0.1
	dfTL = 0.1
	spTL = "Stabilwax"
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
	shift1 = 0.0
	shift2 = 0.0
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
	sys = GasChromatographySystems.GCxGC_TM(L1, d1, df1, sp1, TP1, L2, d2, df2, sp2, TP2, LTL, dTL, dfTL, spTL, TTL, LM, dM, dfM, spM, shift1, shift2, PM, ratio, Thot, Tcold, TPM, F, pin, pout)
end

# ╔═╡ 03b1cae0-cf56-4e0f-b0ad-fdfa736e0fa8
md"""
### Substances
"""

# ╔═╡ 397d7248-e069-4958-849d-6761cb20283d
begin
	db = DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/Estimator/GCsim_d_renamed.csv", header=1, silencewarnings=true))
	db.Tchar = db.Tchar .- 273.15
	insertcols!(db, 1, :No => collect(1:length(db.Name)))
	db
end

# ╔═╡ 35bbe26a-f82c-4f60-a2d0-74b052037ae7
selected_solutes_ = GasChromatographySystems.common_solutes(db, sys).Name

# ╔═╡ 262fba99-1b75-4803-8f34-a0619de559a9
selected_solutes = selected_solutes_[[14,16,17]]

# ╔═╡ 1964abf5-a627-4014-8348-e81145647104
md"""
## Simulation
"""

# ╔═╡ 5f3c1ce1-50bb-44ba-abac-d3ccc82fdb51
par = GasChromatographySystems.graph_to_parameters(sys, db, selected_solutes; Tcold_abs=Tcold_abs, spatial=spatial, tflank=flank, ng=ng)

# ╔═╡ a2d062a4-bd59-4ab7-9316-dc769b2c3c09
begin
	TM_itp = par[3].prog.T_itp
	Plots.plot(0.0:0.01:sys.modules[3].PM, TM_itp.(sys.modules[3].length/2, 0.0:0.01:sys.modules[3].PM).-273.115, title="temperature at modulator point", xlabel="time in s", ylabel="temperature in °C", label="")
	Plots.plot!([0.0, 0.0], [TM_itp(sys.modules[3].length/2,0.0), TM_itp(sys.modules[3].length/2,(1+ratio)/2*PM)].-273.15, c=:black, label="")
	Plots.plot!([ratio*PM, ratio*PM], [TM_itp(sys.modules[3].length/2,0.0), TM_itp(sys.modules[3].length/2,(1+ratio)/2*PM)].-273.15, c=:black, linestyle=:dash, label="")
	Plots.plot!([PM, PM], [TM_itp(sys.modules[3].length/2,0.0), TM_itp(sys.modules[3].length/2,(1+ratio)/2*PM)].-273.15, c=:black, label="")
end

# ╔═╡ b8a6cdec-457b-4fb0-b6f8-60f713f4a015
begin
	Plots.plot(xlabel="position x in m", ylabel="temperature in °C", title="temperature at modulator point", legend=false)
	Plots.plot!(0.0:sys.modules[3].length/100.0:sys.modules[3].length, TM_itp.(0.0:sys.modules[3].length/100.0:sys.modules[3].length, sys.modules[3].ratio*sys.modules[3].PM).-273.15, label="t=$(sys.modules[3].ratio*sys.modules[3].PM)s")
	Plots.plot!(0.0:sys.modules[3].length/100.0:sys.modules[3].length, TM_itp.(0.0:sys.modules[3].length/100.0:sys.modules[3].length, (1+3*sys.modules[3].ratio)/4*sys.modules[3].PM).-273.15, label="t=$((1+3*sys.modules[3].ratio)/4*sys.modules[3].PM)s")
	Plots.plot!(0.0:sys.modules[3].length/100.0:sys.modules[3].length, TM_itp.(0.0:sys.modules[3].length/100.0:sys.modules[3].length, (1+sys.modules[3].ratio)/2*sys.modules[3].PM).-273.15, label="t=$((1+sys.modules[3].ratio)/2*sys.modules[3].PM)s")
	Plots.plot!(0.0:sys.modules[3].length/100.0:sys.modules[3].length, TM_itp.(0.0:sys.modules[3].length/100.0:sys.modules[3].length, 3.95).-273.15, label="t=3.95s")
	Plots.plot!(0.0:sys.modules[3].length/100.0:sys.modules[3].length, TM_itp.(0.0:sys.modules[3].length/100.0:sys.modules[3].length, sys.modules[3].PM).-273.15, label="t=$(sys.modules[3].PM)s")
end

# ╔═╡ e72a6c41-2b14-45cb-81e7-263610c3c8ae
begin
	Plots.plot(0.0:0.01:(sys.modules[3].PM*100.0), TM_itp.(sys.modules[3].length/2, 0.0:0.01:(sys.modules[3].PM*100.0)).-273.115, title="temperature at modulator point", xlabel="time in s", ylabel="temperature in °C", label="")
	Plots.plot!(0.0:0.01:(sys.modules[3].PM*100.0), par[1].prog.T_itp(0.0, 0.0:0.01:(sys.modules[3].PM*100.0)).-273.15, label="")
end

# ╔═╡ b284207c-70b5-4da5-94e9-0fe70680b0e3
paths = GasChromatographySystems.all_paths(sys.g, 1)[2]

# ╔═╡ 2b41aa34-d994-440b-a8e1-dd6337008a48
sim = GasChromatographySystems.simulate_along_paths(sys, paths, par; algTM=algTM, dt=sys.modules[3].length*1e-7)

# ╔═╡ baf56852-3559-4f11-a688-8d2724ed7fe5
md"""
## Determine tR1 and tR2
"""

# ╔═╡ 22d1d757-ef2f-40a5-945a-c153de565781
fit_envel = GasChromatographySystems.fit_envelope(sim[2][1][8], sim[2][1][2])

# ╔═╡ dff26d84-0a03-448f-94e1-cabd3ff7cc44
fit_D1 = GasChromatographySystems.fit_gauss_D1(sim[2][1][8], sim[2][1][2], PM)

# ╔═╡ a70a7f0a-53a8-47ec-ace2-ae59f1aa3fbd
fit_D2 = GasChromatographySystems.fit_gauss_D2(sim[2][1][8], PM)

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
GasChromatographySystems.chrom(sim[2][1][3])[1]

# ╔═╡ f1f1ffdb-bc7e-48fd-b688-21f375ab21eb
GasChromatographySystems.chrom_marked(sim[2][1][4], PM, ratio, shift1; nτ=6)[1]

# ╔═╡ 363533f9-0bc4-4ea9-82ad-b8b22a42e323
md"""
### Enveloping curve
"""

# ╔═╡ 7c2c251b-5cc7-4dc6-b385-d6a408628acc
tall, call, ts, cs = GasChromatographySystems.chrom(sim[2][1][end]; nτ=6)[2:end]

# ╔═╡ bd810db6-ed39-482e-8107-af3169ff09be
begin
	gr()
	@. model_g(x,p) = p[3]/sqrt(2*π*p[2]^2)*exp(-(x-p[1])^2/(2*p[2]^2))
	p_envel = GasChromatographySystems.collect_chrom(sim[2][1], sys; markings=true)[end]
	t1, t2 = p_envel[1][1][:x_extrema]
	t = t1:(t2-t1)/1000:t2
	ymax = Array{Float64}(undef, length(fit_envel.Name))
	for i=1:length(fit_envel.Name)
		ymax[i] = model_g(fit_envel.fits[i].param[1], fit_envel.fits[i].param)
		Plots.scatter!(p_envel, fit_envel.tRs[i], fit_envel.heights[i], msize=2, c=i+1)
		Plots.plot!(p_envel, t, model_g(t, fit_envel.fits[i].param), c=i+1, linestyle=:dash)
		Plots.plot!(p_envel, [fit_envel.fits[i].param[1], fit_envel.fits[i].param[1]], [0.0, ymax[i]], c=i+1, linestyle=:dot)
		Plots.scatter!(p_envel, (fit_envel.fits[i].param[1], ymax[i]), c=i+1, shape=:diamond)
	end
	Plots.plot!(p_envel, xticks=0:4:3000, legend=false, ylims=(-0.01*maximum(ymax), 1.1*maximum(ymax)), title="Enveloping curve")
	p_envel
end

# ╔═╡ 7a81859f-7c0c-4115-9cca-367753898427
md"""
### Projection on 1st dimension
"""

# ╔═╡ 8cbc706f-1d3f-49ac-8f1c-d281b2318dc1
begin
	plotly()
	p_proj1stD = Plots.plot(legend=false, xlabel="time in s")

	mins_ = minimum.(fit_D1.tRs)
	maxs_ = maximum.(fit_D1.tRs)
	for i=1:length(fit_envel.Name)
		ymax = model_g(fit_D1.fits[i].param[1], fit_D1.fits[i].param)
		for j=1:length(fit_D1.tRs[i]) 
			Plots.plot!(p_proj1stD, [fit_D1.tRs[i][j], fit_D1.tRs[i][j]], [0.0, fit_envel.heights[i][j]], c=i+1, linewidth=10, linealpha=0.5)
		end
		Plots.scatter!(p_proj1stD, fit_D1.tRs[i], fit_D1.heights[i], msize=2, c=i+1)

		tt = mins_[i]:(maxs_[i]-mins_[i])/1000.0:maxs_[i]
		Plots.plot!(p_proj1stD, tt, model_g(tt, fit_D1.fits[i].param), c=i+1, linestyle=:dash)

		Plots.plot!(p_proj1stD, [fit_D1.fits[i].param[1], fit_D1.fits[i].param[1]], [0.0, ymax], c=i+1, linestyle=:dot)
		Plots.scatter!(p_proj1stD, (fit_D1.fits[i].param[1], ymax), c=i+1, shape=:diamond)
	end
	Plots.plot!(p_proj1stD, xticks=0:4:3000, legend=false, title="Projection 1st D")
	p_proj1stD
end

# ╔═╡ 2491363c-4d1a-4baf-8c47-f85f550d59f8
md"""
### Projection on 2nd dimension
"""

# ╔═╡ fd5310bd-005a-4bf5-b6b7-79b137a101f7
begin
	plotly()
	p_proj2ndD = Plots.plot(legend=false, xlabel="time in s", title="Projection 2nd D")

	for i=1:length(cs)
		Plots.plot!(p_proj2ndD, ts[i].-fld.(ts[i], PM).*PM, cs[i], c=findfirst(sim[2][1][end].Name[i].==selected_solutes)+1)
	end
	mins = minimum.(fit_D2.tRs)
	maxs = maximum.(fit_D2.tRs)
	for i=1:length(fit_D2.Name)
		ymax = model_g(fit_D2.fits[i].param[1], fit_D2.fits[i].param)
		col = findfirst(fit_D2.Name[i].==selected_solutes)+1
		Plots.scatter!(p_proj2ndD, fit_D2.tRs[i], fit_D2.heights[i], msize=2, c=col)

		tt = 0.9*mins[i]:(1.1*maxs[i]-0.9*mins[i])/1000.0:1.1*maxs[i]
		Plots.plot!(p_proj2ndD, tt, model_g(tt, fit_D2.fits[i].param), c=col, linestyle=:dash)

		Plots.plot!(p_proj2ndD, [fit_D2.fits[i].param[1], fit_D2.fits[i].param[1]], [0.0, ymax], c=col, linestyle=:dot)
		Plots.scatter!(p_proj2ndD, (fit_D2.fits[i].param[1], ymax), c=col, shape=:diamond)
	end
	p_proj2ndD
end

# ╔═╡ 7739a83c-d001-47cd-aee2-36917106c2ad
md"""
### 2D Chromatogram
"""

# ╔═╡ 582499db-7b5a-4054-92c2-d0e2a97cd91f
slice_mat, t_D1, t_D2, c_slices, t_, chrom_sliced_sum, chrom_sliced  = GasChromatographySystems.chrom2d(sim[2][1][end], sys)

# ╔═╡ f268df70-af5a-482f-85e6-84216681ecf4
md"""
## Peaklist
"""

# ╔═╡ e5d4fb3a-f05f-4533-9666-08be32ef19ef
pl_GCxGC = GasChromatographySystems.peaklist_GCxGC(sim[2][1][8], sim[2][1][2], PM)

# ╔═╡ 5f67ab8d-3b3d-4bbd-81c8-04eec158ef6e
begin
	plotly()
	x1 = floor(minimum(pl_GCxGC.tR1*0.99))
	x2 = ceil(maximum(pl_GCxGC.tR1*1.01))
	y1 = floor(minimum(pl_GCxGC.tR2*0.99))
	y2 = ceil(maximum(pl_GCxGC.tR2*1.01))
	p_contour = Plots.contour(t_D1[1:end-1], t_D2[1], slice_mat', levels=10, xlabel="tR1 in s", ylabel="tR2 in s")#, c=:jet1)
	Plots.scatter!(p_contour, pl_GCxGC.tR1, pl_GCxGC.tR2, msize=3, shape=:diamond, c=:black, markeralpha=1)
	Plots.plot!(p_contour, xlims=(x1, x2), ylims=(y1, y2))

	p_heatmap = Plots.heatmap(t_D1[1:end-1], t_D2[1], slice_mat', levels=10, xlabel="tR1 in s", ylabel="tR2 in s", c=:jet1)
	Plots.scatter!(p_heatmap, pl_GCxGC.tR1, pl_GCxGC.tR2, msize=3, shape=:diamond, c=:black, markeralpha=1)
	Plots.plot!(p_heatmap, xlims=(x1, x2), ylims=(y1, y2))

	Plots.plot(p_contour, p_heatmap, size=(900,800), legend=false)
end

# ╔═╡ 93e2068d-760e-403e-a6fa-984a237ad6c3
md"""
## Checks
"""

# ╔═╡ 9bc833db-7538-4520-9a80-ab49f9bbeb81
GasChromatographySystems.check_peakwidths(sim[2][1][5])

# ╔═╡ ab04f869-20f0-4234-b281-f6c32ac5a36e
function duration_in_module(pl_array, par)
	Δts = Array{DataFrame}(undef, length(pl_array))
	for i=1:length(pl_array)
		CAS_par = [par[i].sub[x].CAS for x in 1:length(par[i].sub)]
		ann_par = [par[i].sub[x].ann for x in 1:length(par[i].sub)]

		Δt = Array{Float64}(undef, length(CAS_par))
		name = Array{String}(undef, length(CAS_par))
		for j=1:length(CAS_par)
			jj = GasChromatographySystems.common_index(pl_array[i], CAS_par[j], ann_par[j])
			Δt[j] = pl_array[i].tR[jj] - par[i].sub[j].t₀
			name[j] = pl_array[i].Name[jj]
		end
		df_Δt = DataFrame(Name=name, CAS=CAS_par, Δt=Δt, Annotations=ann_par)
		Δts[i] = df_Δt
	end
	return Δts
end

# ╔═╡ b872f155-d36f-4a2d-9741-342a108d6b5f
GasChromatographySystems.check_duration_modulation(sim[2][1], sim[4], PM, ratio)

# ╔═╡ daed9ba6-a44c-48a0-999b-5862838253c2
GasChromatographySystems.check_area(sim[2][1][8])

# ╔═╡ 27645404-48fd-488c-8b2e-66300f4b3303
md"""
## Traces 1st Modulation
"""

# ╔═╡ c20db2f4-b5b8-4436-9d4b-df1333659d25
begin
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
end

# ╔═╡ d4c0a1c4-68ef-408a-9919-1bd737b6c229
md"""
## Traces 2nd Modulation
"""

# ╔═╡ 2737c1e7-6239-41d7-8868-242bb881f436
begin
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
end

# ╔═╡ 271bb7a6-9a98-4336-9571-c131cc125486
md"""
## Measurement
"""

# ╔═╡ b7c2070b-16bf-41f5-b021-160e54a53a21
meas = DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/Estimator/ZB1xWax/meas_GCxGC.csv", header=1, silencewarnings=true, stringtype=String))

# ╔═╡ 189034b6-09e5-4e6b-80b3-e442ad30705c
pl_GCxGC

# ╔═╡ 3ed4ad21-c8e4-428a-b6d3-b37cfdbd0778
compare = GasChromatographySystems.comparison_meas_sim(meas, pl_GCxGC)

# ╔═╡ a8de50cb-e4a6-41dc-9f7b-e88d742c0efa
sum(collect(skipmissing(compare.relΔtR2_percent)))/length(collect(skipmissing(compare.relΔtR2_percent)))

# ╔═╡ c190dd0f-6f09-4d18-9aa5-6b40beddfbf3
sum(collect(skipmissing(compare.relΔtR2_percent)))/length(collect(skipmissing(compare.relΔtR2_percent)))

# ╔═╡ fb7236b5-98fd-4aac-86ed-5e7118b19ed6
meas_chrom = sort!(DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/Estimator/ZB1xWax/CSV/KetAlkPhenHR3Mod1.csv", header=1, silencewarnings=true)), :RT)

# ╔═╡ a2cc7971-bd69-4623-86a1-7e7d43c6b05c
c_slices_m, t_D1_m, t_D2_m = GasChromatographySystems.chrom_slicing(meas_chrom.RT.*60, meas_chrom.TIC, PM)

# ╔═╡ 6cf28744-7748-4d88-bfd7-f8e428732228
begin
	itps = Array{Interpolations.Extrapolation}(undef, length(c_slices_m))
	for i=1:length(c_slices_m)
		itps[i] = LinearInterpolation((t_D2_m[i],), c_slices_m[i], extrapolation_bc=Flat())
	end
end

# ╔═╡ e82034ef-57fa-4f78-b663-7465a97ff127
sim[4]

# ╔═╡ b28b299c-1c2f-4638-a91f-a970c92e4f5a
begin
	chrom_mat = Array{Float64}(undef, length(c_slices_m), length(0.0:0.001:4.0))
	for i=1:length(c_slices_m)
		chrom_mat[i,:] = itps[i].(0.0:0.001:4.0)
	end
	chrom_mat
end

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
	p_meas = Plots.heatmap(0.0:4.0:(size(chrom_mat)[1]-1)*4.0, 0.0:0.001:4.0, chrom_mat_mod', c=:jet1, colorbar=false);
	# add markers for the measured RTs
	Plots.scatter!(meas.tR1[1:5], meas.tR2[1:5], markersize=4, c=:red, label="alcohols");# alcohols
	Plots.scatter!(meas.tR1[6:22], meas.tR2[6:22], markersize=4, c=:orange, label="terpenes") # terpenes
	Plots.scatter!(meas.tR1[23:28], meas.tR2[23:28], markersize=4, c=:yellow, label="phenones"); # phenones
	Plots.scatter!(meas.tR1[29:35], meas.tR2[29:35], markersize=4, c=:lawngreen, label="ketones"); # ketones

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
	Plots.plot!(RT1, RT2, c=:lightblue, label="")

	p_meas
end

# ╔═╡ 954597a2-b847-4cf0-8c1c-ff4b4b4b670b
md"""
# End
"""

# ╔═╡ Cell order:
# ╟─113c6e70-2168-11ee-3f7f-2775a67bab90
# ╠═98217474-a16f-406a-83a7-17fee89c951a
# ╠═22091a27-80e1-4e98-abe7-4b9652cb832c
# ╠═470c8573-a0f5-4ab2-9385-c256491f648d
# ╠═b7b00dda-40d3-40d8-86e5-377d48b189fc
# ╠═14a1a970-3cfe-4086-94d2-d0c7ed6dab15
# ╠═b5efd58e-666b-43c0-9b98-31085c8c36fe
# ╠═9297b636-6daa-40e9-b312-e05f5e5da5fd
# ╟─13d4415d-3dc4-4eed-aded-c3561b77a9eb
# ╟─ddde4e3d-4e72-4315-9555-5f8b3375b04c
# ╠═55194da2-9c5b-464b-bfe1-a65c48bc7801
# ╟─09438575-8370-475f-b2df-fda5155e8c82
# ╠═d0379590-edf3-4017-b8b8-07bd6080f757
# ╠═23b71d48-26bd-496d-9c10-e2d3ce2bc13c
# ╠═a8de50cb-e4a6-41dc-9f7b-e88d742c0efa
# ╟─0e349acc-46b4-4734-abb7-668ec1225c53
# ╠═de542d7e-bd0d-47a4-8977-5ab1e64a26f4
# ╟─a2d062a4-bd59-4ab7-9316-dc769b2c3c09
# ╟─b8a6cdec-457b-4fb0-b6f8-60f713f4a015
# ╟─e72a6c41-2b14-45cb-81e7-263610c3c8ae
# ╟─daf7484b-ac12-4515-8e26-adbc8e214cbc
# ╠═83d10e79-3ab6-418c-80ba-653a50fa5d21
# ╟─5dd6d2b7-1d6f-49ef-b684-a453a50a9e96
# ╠═5883701b-96f6-431e-b29b-b5cc0be940e0
# ╟─03b1cae0-cf56-4e0f-b0ad-fdfa736e0fa8
# ╠═397d7248-e069-4958-849d-6761cb20283d
# ╠═35bbe26a-f82c-4f60-a2d0-74b052037ae7
# ╠═262fba99-1b75-4803-8f34-a0619de559a9
# ╟─1964abf5-a627-4014-8348-e81145647104
# ╠═5f3c1ce1-50bb-44ba-abac-d3ccc82fdb51
# ╠═b284207c-70b5-4da5-94e9-0fe70680b0e3
# ╠═2b41aa34-d994-440b-a8e1-dd6337008a48
# ╟─baf56852-3559-4f11-a688-8d2724ed7fe5
# ╠═22d1d757-ef2f-40a5-945a-c153de565781
# ╠═dff26d84-0a03-448f-94e1-cabd3ff7cc44
# ╠═a70a7f0a-53a8-47ec-ace2-ae59f1aa3fbd
# ╟─bbb75e0b-ccab-4c32-9e33-0d8915dba95d
# ╟─b586403d-12ac-4b38-b4de-f974fe71b5dc
# ╠═7e945bad-0a59-4a5e-8c63-cef21ee847d0
# ╠═d912d940-5b39-4bb9-a136-1f5e67ad7d82
# ╠═f1f1ffdb-bc7e-48fd-b688-21f375ab21eb
# ╟─363533f9-0bc4-4ea9-82ad-b8b22a42e323
# ╠═7c2c251b-5cc7-4dc6-b385-d6a408628acc
# ╟─bd810db6-ed39-482e-8107-af3169ff09be
# ╟─7a81859f-7c0c-4115-9cca-367753898427
# ╟─8cbc706f-1d3f-49ac-8f1c-d281b2318dc1
# ╟─2491363c-4d1a-4baf-8c47-f85f550d59f8
# ╟─fd5310bd-005a-4bf5-b6b7-79b137a101f7
# ╟─7739a83c-d001-47cd-aee2-36917106c2ad
# ╠═582499db-7b5a-4054-92c2-d0e2a97cd91f
# ╟─5f67ab8d-3b3d-4bbd-81c8-04eec158ef6e
# ╟─f268df70-af5a-482f-85e6-84216681ecf4
# ╠═e5d4fb3a-f05f-4533-9666-08be32ef19ef
# ╟─93e2068d-760e-403e-a6fa-984a237ad6c3
# ╠═9bc833db-7538-4520-9a80-ab49f9bbeb81
# ╟─ab04f869-20f0-4234-b281-f6c32ac5a36e
# ╠═b872f155-d36f-4a2d-9741-342a108d6b5f
# ╠═daed9ba6-a44c-48a0-999b-5862838253c2
# ╟─27645404-48fd-488c-8b2e-66300f4b3303
# ╟─c20db2f4-b5b8-4436-9d4b-df1333659d25
# ╟─d4c0a1c4-68ef-408a-9919-1bd737b6c229
# ╟─2737c1e7-6239-41d7-8868-242bb881f436
# ╟─271bb7a6-9a98-4336-9571-c131cc125486
# ╟─b7c2070b-16bf-41f5-b021-160e54a53a21
# ╠═189034b6-09e5-4e6b-80b3-e442ad30705c
# ╟─3ed4ad21-c8e4-428a-b6d3-b37cfdbd0778
# ╠═c190dd0f-6f09-4d18-9aa5-6b40beddfbf3
# ╟─fb7236b5-98fd-4aac-86ed-5e7118b19ed6
# ╠═a2cc7971-bd69-4623-86a1-7e7d43c6b05c
# ╟─6cf28744-7748-4d88-bfd7-f8e428732228
# ╠═e82034ef-57fa-4f78-b663-7465a97ff127
# ╟─b28b299c-1c2f-4638-a91f-a970c92e4f5a
# ╠═847d5b23-967f-4316-941f-fde283200b36
# ╠═954597a2-b847-4cf0-8c1c-ff4b4b4b670b
