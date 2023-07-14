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

# ╔═╡ 22091a27-80e1-4e98-abe7-4b9652cb832c
md"""
# Test of the GCxGC with thermal modulator simulation
"""

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
	L1 = 30.0 # m
	d1 = 0.25 # mm
	df1 = 0.25 # µm
	sp1 = "ZB1ms"
	
	# modulator
	LM = [0.3, 0.005, 0.9, 0.005, 0.3]
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
	Tcold = -70#-100.0
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
	db[!, :No] = collect(1:length(db.Name))
	db
end

# ╔═╡ 35bbe26a-f82c-4f60-a2d0-74b052037ae7
selected_solutes_ = GasChromatographySystems.common_solutes(db, sys).Name

# ╔═╡ 262fba99-1b75-4803-8f34-a0619de559a9
#selected_solutes = selected_solutes_#[10:20]#[[14,16,17]]

# ╔═╡ 09a051d7-9017-4df4-b643-9d0c5f04078f
#selected_solutes = ["4-Methyl-2-pentanone", "2-Hexanone", "2-Decanone", "2-Undecanone", "2-Dodecanone", "2-Pentadecanone"]

# ╔═╡ e322553a-23a0-4b5e-96f1-33732c2ad55a
selected_solutes = selected_solutes_

# ╔═╡ cf547169-b5b1-4fef-8ce3-414db6fa8a2b
# if I select the ketones separatly the GCxGC simulation is successfull
# if I select all substances, than the ketones are not simulated successfully (peak width after 1st Modulator point to big)

# ╔═╡ 1964abf5-a627-4014-8348-e81145647104
md"""
## Simulation
"""

# ╔═╡ 5f3c1ce1-50bb-44ba-abac-d3ccc82fdb51
par = GasChromatographySystems.graph_to_parameters(sys, db, selected_solutes; Tcold_abs=true)

# ╔═╡ b284207c-70b5-4da5-94e9-0fe70680b0e3
paths = GasChromatographySystems.all_paths(sys.g, 1)[2]

# ╔═╡ 2b41aa34-d994-440b-a8e1-dd6337008a48
sim = GasChromatographySystems.simulate_along_paths(sys, paths, par; abstolTM=1e-10, reltolTM=1e-8, algTM=Vern9())

# ╔═╡ baf56852-3559-4f11-a688-8d2724ed7fe5
md"""
## Determine tR1 and tR2
"""

# ╔═╡ 31a51f15-f021-4934-bc55-cb3bba7794c6
function heights_of_peaks(pl)
	heights = Array{Float64}(undef, length(pl.tR))
 	for i=1:length(pl.tR)
		heights[i] = (GasChromatographySimulator.chromatogram([pl.tR[i]], [pl.tR[i]], [pl.τR[i]])[1].*pl.A[i])
	end
	#Plots.plot(t_, chrom_sliced)
	heights
end

# ╔═╡ 7e7c73cc-1142-40d2-be46-00b9f760bbf4
function fit_gauss_D1(pl_D2, pl_D1)
	heights = heights_of_peaks(pl_D2)
	tR = pl_D2.tR
	Name = unique(pl_D2.Name)
	nsub = length(Name)

	@. model_g(x,p) = p[3]/sqrt(2*π*p[2]^2)*exp(-(x-p[1])^2/(2*p[2]^2))
	
	sort_heights = Array{Array{Float64}}(undef, nsub)
	sort_tR = Array{Array{Float64}}(undef, nsub)
	fits = Array{LsqFit.LsqFitResult}(undef, nsub)
	for i=1:nsub
		ii_name = findall(Name[i].==pl_D2.Name)
		ii = findfirst(Name[i].==pl_D1.Name)
		sort_heights[i] = heights[ii_name]
		sort_tR[i] = tR[ii_name]
		fits[i] = curve_fit(model_g, sort_tR[i], sort_heights[i], [pl_D1.tR[ii], pl_D1.τR[ii], 1.0])
	end
	return DataFrame(Name=Name, tRs=sort_tR, heights=sort_heights, fits=fits)
end

# ╔═╡ fd8f4663-d5f3-445f-b3c4-4b902a8ae2e4
function fit_gauss_D2(pl_D2, PM) # shift ?
	heights = heights_of_peaks(pl_D2)
	tR = pl_D2.tR .- fld.(pl_D2.tR, PM).*PM
	Name = unique(pl_D2.Name)
	nsub = length(Name)

	@. model_g(x,p) = p[3]/sqrt(2*π*p[2]^2)*exp(-(x-p[1])^2/(2*p[2]^2))
	
	sort_heights = Array{Array{Float64}}(undef, nsub)
	sort_tR = Array{Array{Float64}}(undef, nsub)
	fits = Array{LsqFit.LsqFitResult}(undef, nsub)
	for i=1:nsub
		ii_name = findall(Name[i].==pl_D2.Name)
		#ii = findfirst(Name[i].==pl_D1.Name)
		sort_heights[i] = heights[ii_name]
		sort_tR[i] = tR[ii_name]
		mean_tR = sum(sort_tR[i])/length(sort_tR[i])
		fits[i] = curve_fit(model_g, sort_tR[i], sort_heights[i], [mean_tR, mean_tR/10, 1.0])
	end
	return DataFrame(Name=Name, tRs=sort_tR, heights=sort_heights, fits=fits)
end

# ╔═╡ 22d1d757-ef2f-40a5-945a-c153de565781
fit_D1 = fit_gauss_D1(sim[2][1][8], sim[2][1][2])

# ╔═╡ a70a7f0a-53a8-47ec-ace2-ae59f1aa3fbd
fit_D2 = fit_gauss_D2(sim[2][1][8], PM)

# ╔═╡ 05699339-f326-4b45-9e1d-6b8898dc47c8
function peaklist_GCxGC(pl_end, pl_1D, PM)
	fit_D1 = fit_gauss_D1(pl_end, pl_1D)
	fit_D2 = fit_gauss_D2(pl_end, PM)
	tR1 = Array{Float64}(undef, length(fit_D1.Name))
	tR2 = Array{Float64}(undef, length(fit_D1.Name))
	for i=1:length(fit_D1.Name)
		tR1[i] = fit_D1.fits[i].param[1]
		ii = findfirst(fit_D1.Name[i].==fit_D2.Name)
		tR2[i] = fit_D2.fits[ii].param[1]
	end
	pl_GCxGC = DataFrame(Name=fit_D1.Name, tR1=tR1, tR2=tR2)
end

# ╔═╡ 6c787874-7b2e-4b9a-a3d5-5a1bc5abf41b
pl_GCxGC = peaklist_GCxGC(sim[2][1][8], sim[2][1][2], PM)

# ╔═╡ bbb75e0b-ccab-4c32-9e33-0d8915dba95d
md"""
## Chromatograms
"""

# ╔═╡ 0e86e900-c0b6-4dd0-9121-8d428411febe
function chrom(pl; nτ=6)
	tstart = Array{Float64}(undef, length(pl.Name))
	tend = Array{Float64}(undef, length(pl.Name))
	t = Array{Array{Float64}}(undef, length(pl.Name))
	c = Array{Array{Float64}}(undef, length(pl.Name))
	for i=1:length(pl.Name)
		tstart[i] = pl.tR[i] - nτ * pl.τR[i]
		tend[i] = pl.tR[i] + nτ * pl.τR[i]
		t[i] = collect(tstart[i]:(2*nτ*pl.τR[i]/100):tend[i])
		c[i] = GasChromatographySimulator.chromatogram(t[i], [pl.tR[i]], [pl.τR[i]])*pl.A[i]
	end
	t1 = minimum(tstart)
	t2 = maximum(tend)
	dt = 2*nτ*minimum(pl.τR)/100
	t_sum = collect(t1:dt:t2)
	c_sum = fill(0.0, length(t_sum))
	for i=1:length(pl.Name)
		c_ = GasChromatographySimulator.chromatogram(t_sum, [pl.tR[i]], [pl.τR[i]])*pl.A[i]
		c_sum = c_sum .+ c_
	end

	names = unique(pl.Name)
	
	p_chrom = Plots.plot(t_sum, c_sum, xlabel="time in s", label="Chromatogram")
	for i=1:length(pl.Name)
		i_names = findfirst(pl.Name[i].==names)
		#if i > length(names)
		#	lbl = ""
		#else
			lbl = pl.Name[i]
		#end
		Plots.plot!(p_chrom, t[i], c[i], label=lbl, color=i_names+1)
	end
	Plots.plot!(p_chrom, ylims=(-0.02*maximum(c_sum), 1.02*maximum(c_sum)), xlims=(minimum(t_sum), maximum(t_sum)))
	return p_chrom, t_sum, c_sum, t, c 
end

# ╔═╡ 4ab0a12d-bdce-4c79-a2d2-14405be66c8b
function chrom_marked(pl, PM, ratio, shift; nτ=6)
	p_chrom, t_sum, c_sum, t, c = chrom(pl; nτ=nτ)
	max_y = maximum(c_sum)
	# add modulation period
	n = unique(fld.(t_sum, PM + shift))
	for i=1:length(n)
		Plots.plot!(p_chrom, [n[i]*PM+shift, n[i]*PM+shift], [0.0, 1.1*maximum(c_sum)], c=:black, label="")
		Plots.plot!(p_chrom, [n[i]*PM+shift+PM*ratio, n[i]*PM+shift+PM*ratio], [0.0, 1.1*maximum(c_sum)], c=:black, linestyle=:dash, label="")
	end
	Plots.plot!(p_chrom, ylims=(-0.02*maximum(c_sum), 1.02*maximum(c_sum)), xlims=(minimum(t_sum), maximum(t_sum)))
	return p_chrom, t_sum, c_sum, t, c
end

# ╔═╡ 3b459167-a48b-47ab-891c-5de71a304c3d
begin
	gr()
	ptest, t_sum, c_sum, t, c  = chrom(sim[2][1][8])
	ptest
end

# ╔═╡ 67778138-8ba8-4417-95a7-971de0977734
sim[2][1][8]

# ╔═╡ 93e2068d-760e-403e-a6fa-984a237ad6c3
md"""
### Checks
"""

# ╔═╡ 13025994-1118-491e-90b0-f12a9e1a0b7f
filter([:τR] => x -> x > 0.5, sim[2][1][3]) # solutions with a to broad peak width after TM1

# ╔═╡ de0b3bf4-c244-4852-98da-44b0cc8dd51b
filter([:τR] => x -> x > 0.5, sim[2][1][8]) # solutions with a to broad peak width at the end

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

# ╔═╡ 6caa8902-e681-4f6b-b6da-ecdfc76d25de
Δts = duration_in_module(sim[2][1], sim[4])

# ╔═╡ d5694843-c834-42c8-b9d1-d2fdb4b79364
Δts[3]

# ╔═╡ 61fecca4-0cd8-4655-bd85-de54865925be
ratio*PM

# ╔═╡ e1644c22-1daf-41ef-ab2d-ac664dac21d1
function check_duration_modulation(pl_array, par, PM, ratio)
	Δts = duration_in_module(pl_array, par)
	ok_TM1 = Array{Bool}(undef, length(Δts[3].Δt))
	for i=1:length(Δts[3].Δt)
		if (Δts[3].Δt[i] > ratio*PM - 2/3*(1-ratio)*PM) && (Δts[3].Δt[i] < PM)
			ok_TM1[i] = true
		else
			ok_TM1[i] = false
		end
	end
	ok_TM2 = Array{Bool}(undef, length(Δts[5].Δt))
	for i=1:length(Δts[5].Δt)
		if (Δts[5].Δt[i] > ratio*PM - 2/3*(1-ratio)*PM) && (Δts[5].Δt[i] < PM)
			ok_TM2[i] = true
		else
			ok_TM2[i] = false
		end
	end
	index_TM1 = findall(ok_TM1.==false)
	index_TM2 = findall(ok_TM2.==false)
	TM1 = DataFrame(Index=index_TM1, Name=Δts[3].Name[index_TM1], CAS=Δts[3].CAS[index_TM1], tR=pl_array[3].tR[index_TM1], τR=pl_array[3].τR[index_TM1], Δt=Δts[3].Δt[index_TM1])
	TM2 = DataFrame(Index=index_TM2, Name=Δts[5].Name[index_TM2], CAS=Δts[5].CAS[index_TM2], tR=pl_array[5].tR[index_TM2], τR=pl_array[5].τR[index_TM2], Δt=Δts[5].Δt[index_TM2])
	return TM1, TM2
end

# ╔═╡ b872f155-d36f-4a2d-9741-342a108d6b5f
check_duration_modulation(sim[2][1], sim[4], PM, ratio)

# ╔═╡ e1097a60-ef24-4142-8aad-44d96fd27162
function check_peakwidths(pl; τ_threshold = 0.5)
	pl_f = filter([:τR] => x -> x > τ_threshold, pl)
	names = unique(pl_f.Name)
	return names
end

# ╔═╡ 8a524e35-88c8-475b-a33c-50d80b25414e
check_peakwidths(sim[2][1][8])

# ╔═╡ 9bc833db-7538-4520-9a80-ab49f9bbeb81
check_peakwidths(sim[2][1][5])

# ╔═╡ dccab549-442c-4c59-82a0-a1e756896c37
function check_area(pl)
	names = unique(pl.Name)
	area = Array{Float64}(undef, length(names))
	for i=1:length(names)
		area[i] = sum(filter([:Name] => x -> x == names[i], pl).A)
	end
	return DataFrame(Name=names, sum_A = area)
end

# ╔═╡ daed9ba6-a44c-48a0-999b-5862838253c2
check_area(sim[2][1][8])

# ╔═╡ a968c5e4-43e9-444f-af10-776275b034a6
# add test for the duration in the modulator point, should be < PM

# ╔═╡ 2264cdc9-722a-4fe3-8833-c3f7c7b5968f
function traces(sol, par, i_select)
	z = sol[i_select].t
	
		tt = Array{Float64}(undef, length(sol[i_select].t))
		ττ = Array{Float64}(undef, length(sol[i_select].t))
		TT = Array{Float64}(undef, length(sol[i_select].t))
		kk = Array{Float64}(undef, length(sol[i_select].t))
		for j=1:length(sol[i_select].t)
			tt[j] = sol[i_select].u[j][1]
			ττ[j] = sol[i_select].u[j][2]
			TT[j] = par.prog.T_itp(sol[i_select].t[j], sol[i_select].u[j][1])
			kk[j] = GasChromatographySimulator.retention_factor(sol[i_select].t[j], sol[i_select].u[j][1], par.col, par.prog, par.sub[i_select], par.opt)
		end
	return trace = DataFrame(z=z, t=tt, τ²=ττ, T=TT, k=kk)
end

# ╔═╡ 81026230-e998-4cd1-b04b-b967995df5f9
sim[2][1][2]

# ╔═╡ 030cef9c-13b2-44b0-8451-a5c1caee8e4c
sim[2][1][3]

# ╔═╡ 029c457d-af54-442f-96f0-9a486e73f31d
sim[4]

# ╔═╡ 5ce4212c-67cc-4306-8e2b-fb7e8927abdd
check_duration_modulation(sim[2][1], sim[4], PM, ratio)

# ╔═╡ d95ac504-4e34-432f-837e-8c12dd8ce6f8
ratio*PM - 2/3*(1-ratio)*PM

# ╔═╡ 2f8c4181-4912-4603-bd19-6b433111a77b
sim[2][1][5]

# ╔═╡ c20db2f4-b5b8-4436-9d4b-df1333659d25
begin
	plotly()
	p_xt_check = Plots.plot(xlabel="time in s", ylabel="x in m", legend=true)
	p_τt_check = Plots.plot(xlabel="time in s", ylabel="τ in s", legend=true)
	p_Tt_check = Plots.plot(xlabel="time in s", ylabel="T in °C", legend=true)
	p_lnkt_check = Plots.plot(xlabel="time in s", ylabel="lnk", legend=true)
	for i=35:41
		trace = traces(sim[3][1][3], sim[4][3], i)
		lbl = string(split(sim[4][3].sub[i].ann, "_")[1], "_", sim[4][3].sub[i].name)
		Plots.plot!(p_xt_check, trace.t, trace.z, label=lbl, marker=:circle)
		Plots.plot!(p_τt_check, trace.t, sqrt.(trace.τ²), label=lbl, marker=:circle)
		Plots.plot!(p_Tt_check, trace.t, trace.T.-273.15, label=lbl, marker=:circle)
		Plots.plot!(p_lnkt_check, trace.t, log.(trace.k), label=lbl, marker=:circle)
	end
	Plots.plot(p_xt_check, p_τt_check, p_Tt_check, p_lnkt_check, size=(800,800))
end

# ╔═╡ 2737c1e7-6239-41d7-8868-242bb881f436
begin
	plotly()
	p_xt_check_ = Plots.plot(xlabel="time in s", ylabel="x in m", legend=true)
	p_τt_check_ = Plots.plot(xlabel="time in s", ylabel="τ in s", legend=true)
	p_Tt_check_ = Plots.plot(xlabel="time in s", ylabel="T in °C", legend=true)
	p_lnkt_check_ = Plots.plot(xlabel="time in s", ylabel="lnk", legend=true)
	for i=139:151
		trace = traces(sim[3][1][5], sim[4][5], i)
		lbl = string(join(split(sim[4][5].sub[i].ann, "_")[1:2], "_"), "_", sim[4][5].sub[i].name)
		Plots.plot!(p_xt_check_, trace.t, trace.z, label=lbl, marker=:circle)
		Plots.plot!(p_τt_check_, trace.t, sqrt.(trace.τ²), label=lbl, marker=:circle)
		Plots.plot!(p_Tt_check_, trace.t, trace.T.-273.15, label=lbl, marker=:circle)
		Plots.plot!(p_lnkt_check_, trace.t, log.(trace.k), label=lbl, marker=:circle)
	end
	Plots.plot(p_xt_check_, p_τt_check_, p_Tt_check_, p_lnkt_check_, size=(800,800))
end

# ╔═╡ b982fae3-da16-4aba-8e93-e8dbff571413
begin
	plotly()
	p_τt = Plots.plot(xlabel="time in s", ylabel="τ in s", legend=false)
	for i=1:length(sim[3][1][3])
		trace = traces(sim[3][1][3], sim[4][3], i)
		Plots.plot!(p_τt, trace.t, sqrt.(trace.τ²))
	end
	p_τt
end

# ╔═╡ 04d0312b-ecaf-4f16-acfa-65a8811188a8
begin
	plotly()
	p_Tt = Plots.plot(xlabel="time in s", ylabel="temperature in °C", legend=false)
	for i=1:length(sim[3][1][3])
		trace = traces(sim[3][1][3], sim[4][3], i)
		Plots.plot!(p_Tt, trace.t, trace.T.-273.15)
	end
	p_Tt
end

# ╔═╡ 727273fd-12ad-4fe5-b336-b3dcf1a13773
begin
	plotly()
	p_kt = Plots.plot(xlabel="time in s", ylabel="k", legend=false)
	for i=1:length(sim[3][1][3])
		trace = traces(sim[3][1][3], sim[4][3], i)
		Plots.plot!(p_kt, trace.t, trace.k)
	end
	p_kt
end

# ╔═╡ 7f97133b-52a2-4af8-a340-41cc9aae8b04
begin
	plotly()
	p_lnkt = Plots.plot(xlabel="time in s", ylabel="lnk", legend=false)
	for i=1:length(sim[3][1][3])
		trace = traces(sim[3][1][3], sim[4][3], i)
		Plots.plot!(p_lnkt, trace.t, log.(trace.k))
	end
	p_lnkt
end

# ╔═╡ b18f7ee6-19fb-465c-ba43-5a158b21821e
chrom_marked(sim[2][1][4], PM, ratio, shift2)[1]

# ╔═╡ adab4c1c-d861-4df7-9787-768177dfc380
function collect_chrom(pl_array, sys; markings=true)
	# collect chromatograms for all segments
	# chromatograms of segments, which are followed by a ModuleTM should be marked with the modulation period (coldjet, hotjet)
	# add titles (module name)
	p = Array{Any}(undef, length(pl_array))
	for i=1:length(pl_array)
		if markings == true && i < length(pl_array)
			if typeof(sys.modules[i+1]) == GasChromatographySystems.ModuleTM
				p[i] = chrom_marked(pl_array[i], sys.modules[i+1].PM, sys.modules[i+1].ratio, sys.modules[i+1].shift)[1]
				
			else
				p[i] = chrom(pl_array[i])[1]
			end
		else
			p[i] = chrom(pl_array[i])[1]
		end
		Plots.plot!(p[i], title=sys.modules[i].name)
	end
	return p
end

# ╔═╡ 76998e3f-7e04-480b-b36e-2eacc621cedd
#Plots.plot(collect_chrom(sim[2][1], sys; markings=false)..., legend=false)

# ╔═╡ 8a91adc9-c3a9-4b77-bda8-b982c050401d
function chrom_slicing(t, c, PM, shift) # correctly account for the shift!!!
	n = Int.(fld.(collect(t), PM+shift)) # number of the slices
	slices = Array{Array{Float64, 1}, 1}(undef, length(unique(n)))
	t_D2 = Array{Array{Float64, 1}, 1}(undef, length(unique(n)))
	for i=1:length(unique(n))
		i1 = findfirst(unique(n)[i].==n)
		i2 = findlast(unique(n)[i].==n)
		slices[i] = c[i1:i2]
		t_D2[i] = t[i1:i2] .- unique(n)[i] * PM
	end
	t_D1 = 0.0:PM:t[end] # shift?
	return slices, t_D1, t_D2
end

# ╔═╡ 355df747-6062-41f1-b1a9-55f9ef2f5faf
pl_GCxGC

# ╔═╡ 271bb7a6-9a98-4336-9571-c131cc125486
md"""
## Measurement
"""

# ╔═╡ b7c2070b-16bf-41f5-b021-160e54a53a21
meas = DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/Estimator/ZB1xWax/meas_GCxGC.csv", header=1, silencewarnings=true, stringtype=String))

# ╔═╡ 189034b6-09e5-4e6b-80b3-e442ad30705c
pl_GCxGC

# ╔═╡ 75f44ad0-ee52-49d0-a59f-17b87f516122
[findfirst(meas.Name[x].==pl_GCxGC.Name) for x in 1:length(meas.Name)]

# ╔═╡ 07fed619-53d1-48d2-a5f9-f44d0e107afa
function comparison_meas_sim(meas, pl_sim)
	Name = meas.Name
	tR1_meas = meas.tR1
	tR2_meas = meas.tR2
	index = [findfirst(meas.Name[x].==pl_sim.Name) for x in 1:length(meas.Name)]
	tR1_sim = Array{Union{Missing,Float64}}(undef, length(meas.Name))
	tR2_sim = Array{Union{Missing,Float64}}(undef, length(meas.Name))
	for i=1:length(meas.Name)
		if isnothing(index[i])
			tR1_sim[i] = missing
			tR2_sim[i] = missing
		else
			tR1_sim[i] = pl_sim.tR1[index[i]]
			tR2_sim[i] = pl_sim.tR2[index[i]]
		end
	end
	comp = DataFrame(Name=Name, tR1_meas=tR1_meas, tR1_sim=tR1_sim, ΔtR1=tR1_meas.-tR1_sim, relΔtR1_percent=(tR1_meas.-tR1_sim)./tR1_meas.*100.0, tR2_meas=tR2_meas, tR2_sim=tR2_sim, ΔtR2=tR2_meas.-tR2_sim, relΔtR2_percent=(tR2_meas.-tR2_sim)./tR2_meas.*100.0)
	#for i=1:length(comp.Name)
	#	ii = findfirst(comp.Name[i].==meas.Name)
	#	if ismissing(meas.tR1[ii])
	#		comp[i, :tR1_meas] = NaN
	#		comp[i, :tR2_meas] = NaN
	#		comp[i, :ΔtR1] = NaN
	#		comp[i, :ΔtR2] = NaN
	#	else
	#		comp[i, :tR1_meas] = meas.tR1[ii]
	#		comp[i, :tR2_meas] = meas.tR2[ii]
	#		comp[i, :ΔtR1] = meas.tR1[ii] - comp[i, :tR1_sim]
	#		comp[i, :ΔtR2] = meas.tR2[ii] - comp[i, :tR2_sim]
	#	end
	#end
	return comp
end

# ╔═╡ 3ed4ad21-c8e4-428a-b6d3-b37cfdbd0778
compare = comparison_meas_sim(meas, pl_GCxGC)

# ╔═╡ a8de50cb-e4a6-41dc-9f7b-e88d742c0efa
sum(collect(skipmissing(compare.relΔtR2_percent)))/length(collect(skipmissing(compare.relΔtR2_percent)))

# ╔═╡ c190dd0f-6f09-4d18-9aa5-6b40beddfbf3
sum(collect(skipmissing(compare.relΔtR2_percent)))/length(collect(skipmissing(compare.relΔtR2_percent)))

# ╔═╡ fb7236b5-98fd-4aac-86ed-5e7118b19ed6
meas_chrom = sort!(DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/Estimator/ZB1xWax/CSV/KetAlkPhenHR3Mod1.csv", header=1, silencewarnings=true)), :RT)

# ╔═╡ 64b6cb30-93b5-4a62-98c0-70ef421b29ba
function chrom_slicing(t1, c, PM)
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

# ╔═╡ 1b5831af-a02f-46d4-9101-5ed9d20cce42
# correct?
# creating the final chromatogram and slicing it up every modulation period
function chrom2d(pl_final, sys)
	t_ = 0.0:0.01:sum(sys.modules[1].temperature.timesteps)
	chrom_sliced = Array{Array{Float64}}(undef, length(pl_final.tR))
	for i=1:length(pl_final.tR)
		chrom_sliced[i] = GasChromatographySimulator.chromatogram(collect(t_), [pl_final.tR[i]], [pl_final.τR[i]]).*pl_final.A[i]
	end
	chrom_sliced_sum = chrom_sliced[1]
	for i=2:length(chrom_sliced)
		chrom_sliced_sum = chrom_sliced_sum .+ chrom_sliced[i]
	end
	#Plots.plot(t_, chrom_sliced_sum)
	# determin the index of the ModuleTM
	c_slices, t_D1, t_D2 = chrom_slicing(t_, chrom_sliced_sum, sys.modules[5].PM, sys.modules[5].shift)
	slice_mat = Array{Float64}(undef, length(c_slices)-1, length(t_D2[1]))
	for j=1:length(t_D2[1])
		for i=1:(length(c_slices)-1)
			slice_mat[i,j] = c_slices[i][j]
		end
	end
	return slice_mat, t_D1, t_D2, c_slices, t_, chrom_sliced_sum, chrom_sliced 
end

# ╔═╡ 46e6fe39-93c3-4131-80bf-a99f365036ff
slice_mat, t_D1, t_D2, c_slices, t_, chrom_sliced_sum, chrom_sliced  = chrom2d(sim[2][1][end], sys)

# ╔═╡ 343e745b-283e-4896-ba06-b3a3d3d943e4
begin
	gr()
	Plots.contour(t_D1[1:end-1], t_D2[1], slice_mat', levels=10, xlabel="tR1 in s", ylabel="tR2 in s")#, c=:jet1)
	Plots.scatter!(pl_GCxGC.tR1, pl_GCxGC.tR2, msize=3, shape=:diamond, c=:gray, markeralpha=0.5)
	#Plots.plot!(xlims=(2500,3100), ylims=(0,1))
end

# ╔═╡ a2cc7971-bd69-4623-86a1-7e7d43c6b05c
c_slices_m, t_D1_m, t_D2_m = chrom_slicing(meas_chrom.RT.*60, meas_chrom.TIC, PM)

# ╔═╡ 6cf28744-7748-4d88-bfd7-f8e428732228
begin
	itps = Array{Interpolations.Extrapolation}(undef, length(c_slices_m))
	for i=1:length(c_slices_m)
		itps[i] = LinearInterpolation((t_D2_m[i],), c_slices_m[i], extrapolation_bc=Flat())
	end
end

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
	gr()
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
	# shift by 15*PM in the 1st dimension
	# why?
	# - measured retention times from different single measurement, where the retention times match the 2D chromatogram
	# - these retention times are shifted by 15 PM respectively to the 2D chromatogram of the measurement of all samples together
	offset_D1 = 0.0#15.0
	Plots.scatter!(meas.tR1[1:5].-offset_D1*PM, meas.tR2[1:5], markersize=4, c=:red, label="alcohols");# alcohols
	Plots.scatter!(meas.tR1[6:22].-offset_D1*PM, meas.tR2[6:22], markersize=4, c=:orange, label="terpenes") # terpenes
	Plots.scatter!(meas.tR1[23:28].-offset_D1*PM, meas.tR2[23:28], markersize=4, c=:yellow, label="phenones"); # phenones
	Plots.scatter!(meas.tR1[29:35].-offset_D1*PM, meas.tR2[29:35], markersize=4, c=:lawngreen, label="ketones"); # ketones

	# simulation
	offset_D2 = 0.0
	Plots.scatter!(pl_GCxGC.tR1.-offset_D1*PM, pl_GCxGC.tR2.+offset_D2, c=:lightblue, m=:diamond, markeralpha=1, msize=3, label="simulated", xlims=(0.0, 3500.0), ylims=(0.0, 4.0), xlabel="tR1 in s", ylabel="tR2 in s");
	#Plots.plot!(p_meas, xlims=(0.0, 3500.0), ylims=(0.0, 4.0), xlabel="tR1 in s", ylabel="tR2 in s")
	p_meas
end

# ╔═╡ 954597a2-b847-4cf0-8c1c-ff4b4b4b670b
md"""
# End
"""

# ╔═╡ Cell order:
# ╠═113c6e70-2168-11ee-3f7f-2775a67bab90
# ╠═22091a27-80e1-4e98-abe7-4b9652cb832c
# ╠═13d4415d-3dc4-4eed-aded-c3561b77a9eb
# ╠═ddde4e3d-4e72-4315-9555-5f8b3375b04c
# ╠═55194da2-9c5b-464b-bfe1-a65c48bc7801
# ╠═09438575-8370-475f-b2df-fda5155e8c82
# ╠═d0379590-edf3-4017-b8b8-07bd6080f757
# ╠═23b71d48-26bd-496d-9c10-e2d3ce2bc13c
# ╠═a8de50cb-e4a6-41dc-9f7b-e88d742c0efa
# ╠═0e349acc-46b4-4734-abb7-668ec1225c53
# ╠═de542d7e-bd0d-47a4-8977-5ab1e64a26f4
# ╠═daf7484b-ac12-4515-8e26-adbc8e214cbc
# ╠═83d10e79-3ab6-418c-80ba-653a50fa5d21
# ╠═5dd6d2b7-1d6f-49ef-b684-a453a50a9e96
# ╠═5883701b-96f6-431e-b29b-b5cc0be940e0
# ╠═03b1cae0-cf56-4e0f-b0ad-fdfa736e0fa8
# ╠═397d7248-e069-4958-849d-6761cb20283d
# ╠═35bbe26a-f82c-4f60-a2d0-74b052037ae7
# ╠═262fba99-1b75-4803-8f34-a0619de559a9
# ╠═09a051d7-9017-4df4-b643-9d0c5f04078f
# ╠═e322553a-23a0-4b5e-96f1-33732c2ad55a
# ╠═cf547169-b5b1-4fef-8ce3-414db6fa8a2b
# ╠═1964abf5-a627-4014-8348-e81145647104
# ╠═5f3c1ce1-50bb-44ba-abac-d3ccc82fdb51
# ╠═b284207c-70b5-4da5-94e9-0fe70680b0e3
# ╠═2b41aa34-d994-440b-a8e1-dd6337008a48
# ╠═8a524e35-88c8-475b-a33c-50d80b25414e
# ╠═baf56852-3559-4f11-a688-8d2724ed7fe5
# ╠═31a51f15-f021-4934-bc55-cb3bba7794c6
# ╠═7e7c73cc-1142-40d2-be46-00b9f760bbf4
# ╠═fd8f4663-d5f3-445f-b3c4-4b902a8ae2e4
# ╠═22d1d757-ef2f-40a5-945a-c153de565781
# ╠═a70a7f0a-53a8-47ec-ace2-ae59f1aa3fbd
# ╠═05699339-f326-4b45-9e1d-6b8898dc47c8
# ╠═6c787874-7b2e-4b9a-a3d5-5a1bc5abf41b
# ╠═bbb75e0b-ccab-4c32-9e33-0d8915dba95d
# ╠═0e86e900-c0b6-4dd0-9121-8d428411febe
# ╠═4ab0a12d-bdce-4c79-a2d2-14405be66c8b
# ╠═3b459167-a48b-47ab-891c-5de71a304c3d
# ╠═67778138-8ba8-4417-95a7-971de0977734
# ╠═93e2068d-760e-403e-a6fa-984a237ad6c3
# ╠═9bc833db-7538-4520-9a80-ab49f9bbeb81
# ╠═13025994-1118-491e-90b0-f12a9e1a0b7f
# ╠═de0b3bf4-c244-4852-98da-44b0cc8dd51b
# ╠═ab04f869-20f0-4234-b281-f6c32ac5a36e
# ╠═6caa8902-e681-4f6b-b6da-ecdfc76d25de
# ╠═d5694843-c834-42c8-b9d1-d2fdb4b79364
# ╠═61fecca4-0cd8-4655-bd85-de54865925be
# ╠═b872f155-d36f-4a2d-9741-342a108d6b5f
# ╠═e1644c22-1daf-41ef-ab2d-ac664dac21d1
# ╠═e1097a60-ef24-4142-8aad-44d96fd27162
# ╠═dccab549-442c-4c59-82a0-a1e756896c37
# ╠═daed9ba6-a44c-48a0-999b-5862838253c2
# ╠═a968c5e4-43e9-444f-af10-776275b034a6
# ╠═2264cdc9-722a-4fe3-8833-c3f7c7b5968f
# ╠═81026230-e998-4cd1-b04b-b967995df5f9
# ╠═030cef9c-13b2-44b0-8451-a5c1caee8e4c
# ╠═029c457d-af54-442f-96f0-9a486e73f31d
# ╠═5ce4212c-67cc-4306-8e2b-fb7e8927abdd
# ╠═d95ac504-4e34-432f-837e-8c12dd8ce6f8
# ╠═2f8c4181-4912-4603-bd19-6b433111a77b
# ╠═c20db2f4-b5b8-4436-9d4b-df1333659d25
# ╠═2737c1e7-6239-41d7-8868-242bb881f436
# ╠═b982fae3-da16-4aba-8e93-e8dbff571413
# ╠═04d0312b-ecaf-4f16-acfa-65a8811188a8
# ╠═727273fd-12ad-4fe5-b336-b3dcf1a13773
# ╠═7f97133b-52a2-4af8-a340-41cc9aae8b04
# ╠═b18f7ee6-19fb-465c-ba43-5a158b21821e
# ╠═adab4c1c-d861-4df7-9787-768177dfc380
# ╠═76998e3f-7e04-480b-b36e-2eacc621cedd
# ╠═8a91adc9-c3a9-4b77-bda8-b982c050401d
# ╠═1b5831af-a02f-46d4-9101-5ed9d20cce42
# ╠═46e6fe39-93c3-4131-80bf-a99f365036ff
# ╠═343e745b-283e-4896-ba06-b3a3d3d943e4
# ╠═355df747-6062-41f1-b1a9-55f9ef2f5faf
# ╠═271bb7a6-9a98-4336-9571-c131cc125486
# ╠═b7c2070b-16bf-41f5-b021-160e54a53a21
# ╠═189034b6-09e5-4e6b-80b3-e442ad30705c
# ╠═75f44ad0-ee52-49d0-a59f-17b87f516122
# ╠═07fed619-53d1-48d2-a5f9-f44d0e107afa
# ╠═3ed4ad21-c8e4-428a-b6d3-b37cfdbd0778
# ╠═c190dd0f-6f09-4d18-9aa5-6b40beddfbf3
# ╠═fb7236b5-98fd-4aac-86ed-5e7118b19ed6
# ╠═a2cc7971-bd69-4623-86a1-7e7d43c6b05c
# ╠═64b6cb30-93b5-4a62-98c0-70ef421b29ba
# ╠═6cf28744-7748-4d88-bfd7-f8e428732228
# ╠═b28b299c-1c2f-4638-a91f-a970c92e4f5a
# ╠═847d5b23-967f-4316-941f-fde283200b36
# ╠═954597a2-b847-4cf0-8c1c-ff4b4b4b670b
