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

- smooth rectangle modulation point in time (`spatial=false`, `tflank=20`, `algTM=Vern9()` ang `ng=true`)
- Tcold is absolute (`Tcold_abs=true`)
"""

# ╔═╡ 470c8573-a0f5-4ab2-9385-c256491f648d
spatial = false

# ╔═╡ b7b00dda-40d3-40d8-86e5-377d48b189fc
flank = 20

# ╔═╡ 14a1a970-3cfe-4086-94d2-d0c7ed6dab15
algTM = Vern9()

# ╔═╡ b5efd58e-666b-43c0-9b98-31085c8c36fe
ng = true

# ╔═╡ 9297b636-6daa-40e9-b312-e05f5e5da5fd
Tcold_abs = true

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
	Tcold = -70.0
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
selected_solutes = selected_solutes_[Not([22, 28, 29])]#[[14,16,17]]

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

# ╔═╡ 2301a74c-2290-4de6-992d-bc887b66269e
sys.modules[3].length*1e-6

# ╔═╡ 80cd5d72-7629-45dc-a10e-53b3db29f21b
# old version
function simulate_along_paths_og(sys, paths, par_sys; t₀=zeros(length(par_sys[1].sub)), τ₀=zeros(length(par_sys[1].sub)), nτ=6, abstolTM=1e-10, reltolTM=1e-8, algTM=Vern9(), refocus=falses(ne(sys.g)), τ₀_focus=zeros(length(par_sys[1].sub)), kwargsTM...)
	
	E = collect(edges(sys.g))
	peaklists = Array{Array{DataFrame,1}}(undef, length(paths))
	solutions = Array{Array{Any,1}}(undef, length(paths))
	path_pos = Array{String}(undef, length(paths))
	new_par_sys = Array{GasChromatographySimulator.Parameters}(undef, length(par_sys))
	visited_E = falses(length(E))

	# i -> path number
	# j -> segment/module number
	for i=1:length(paths)
		i_par = GasChromatographySystems.index_parameter(sys.g, paths[i])
		if GasChromatographySystems.path_possible(sys, paths[i]) == true
			path_pos[i] = "path is possible"
			peaklists_ = Array{DataFrame}(undef, length(i_par))
			solutions_ = Array{Any}(undef, length(i_par))
			#As_ = Array{DataFrame}(undef, length(i_par))
			for j=1:length(i_par)
				if (i>1) && (all(visited_E[1:i_par[j]].==true))
					# was the segment already simulated in a previous simulated path?
					# look in all previous paths for the correct result -> the simulation correlated to the same edge and where this edge is connected to only previouse visited edges
					i_path = 0
					i_edge = 0
					for k=1:i-1 # previous paths
						i_par_previous = GasChromatographySystems.index_parameter(sys.g, paths[k])
						if length(i_par_previous) < j
							j0 = length(i_par_previous)
						else
							j0 = j
						end
						if all(x->x in i_par_previous[1:j0], i_par[1:j]) == true # all edges up to j are the same between the two paths
							i_path = k
							i_edge = findfirst(i_par[j].==i_par_previous)
						end
					end
					# re-use the results
					peaklists_[j] = peaklists[i_path][i_edge]
					solutions_[j] = solutions[i_path][i_edge]
				else # new simulated segments
					if j == 1 # first segment, directly after injection, it is assumed to be a segment of type `ModuleColumn`
						new_par_sys[i_par[j]] = GasChromatographySystems.change_initial(par_sys[i_par[j]], t₀, τ₀)
						peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(new_par_sys[i_par[j]])
						peaklists_[j][!,:A] = ones(length(peaklists_[j].Name)) # add a relativ area factor, splitting of `A` at split points is not accounted for yet, is only used for the slicing of peaks at modulators
					elseif typeof(sys.modules[i_par[j]]) == GasChromatographySystems.ModuleTM
						if refocus[i_par[j]] == true
							τ₀=τ₀_focus
						else
							τ₀=peaklists_[j-1].τR # PM?
						end
						
						new_par_sys[i_par[j]], df_A = GasChromatographySystems.slicing(peaklists_[j-1], sys.modules[i_par[j]].PM, sys.modules[i_par[j]].ratio, sys.modules[i_par[j]].shift, par_sys[i_par[j]]; nτ=nτ, τ₀=τ₀, abstol=abstolTM, reltol=reltolTM, alg=algTM)
						if algTM == "simplifiedTM"
							peaklists_[j], solutions_[j] = approximate_modulator(new_par_sys[i_par[j]], df_A, sys.modules[i_par[j]].PM, sys.modules[i_par[j]].ratio, sys.modules[i_par[j]].shift)
						else
							#kwargsTM = (dtmax = dtmaxTM, )
							peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(new_par_sys[i_par[j]]; kwargsTM...)
							GasChromatographySystems.add_A_to_pl!(peaklists_[j], df_A)
						end
						if maximum(peaklists_[j].τR) > sys.modules[i_par[j]].PM
							return @warn "Peak width of focussed peaks > modulation period. Simulation is aborted."
						end
					else # ModuleColumn
						new_par_sys[i_par[j]] = GasChromatographySystems.change_initial(par_sys[i_par[j]], peaklists_[j-1])
						#else
						#	new_par_sys[i_par[j]] = GasChromatographySystems.change_initial(par_sys[i_par[j]], peaklists_[j-1])
						#end
						peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(new_par_sys[i_par[j]])
						GasChromatographySystems.add_A_to_pl!(peaklists_[j], peaklists_[j-1])
					end
				end
			end
			visited_E[i_par] .= true
			peaklists[i] = peaklists_
			solutions[i] = solutions_
		else # path not possible
			neg_flow_modules = sys.modules[findall(paths[i][findall(GasChromatographySystems.positive_flow(sys)[i_par].==false)].==edges(sys.g))]
			str_neg_flow = "path not possible: "
			for j=1:length(neg_flow_modules)
				str_neg_flow = str_neg_flow*"Flow in module >$(neg_flow_modules[j].name)< becomes negative during the program. "
			end
			path_pos[i] = str_neg_flow
		end
	end
	return path_pos, peaklists, solutions, new_par_sys
end

# ╔═╡ 2b41aa34-d994-440b-a8e1-dd6337008a48
sim = simulate_along_paths_og(sys, paths, par; algTM=algTM, dt=sys.modules[3].length*1e-7)#, dtmax=sys.modules[3].length/100, dtmin=eps())

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
#GasChromatographySystems.chrom(sim[2][1][3])[1]

# ╔═╡ f1f1ffdb-bc7e-48fd-b688-21f375ab21eb
#GasChromatographySystems.chrom_marked(sim[2][1][4], PM, ratio, shift1; nτ=6)[1]

# ╔═╡ 363533f9-0bc4-4ea9-82ad-b8b22a42e323
md"""
### Enveloping curve
"""

# ╔═╡ 7c2c251b-5cc7-4dc6-b385-d6a408628acc
tall, call, ts, cs = GasChromatographySystems.chrom(sim[2][1][end]; nτ=6)[2:end]

# ╔═╡ 6d1ff697-e2d5-4582-b932-7a47c0c0256c
@. model_g(x,p) = p[3]/sqrt(2*π*p[2]^2)*exp(-(x-p[1])^2/(2*p[2]^2))

# ╔═╡ bd810db6-ed39-482e-8107-af3169ff09be
#=begin
	gr()
	#@. model_g(x,p) = p[3]/sqrt(2*π*p[2]^2)*exp(-(x-p[1])^2/(2*p[2]^2))
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
end=#

# ╔═╡ 7a81859f-7c0c-4115-9cca-367753898427
md"""
### Projection on 1st dimension
"""

# ╔═╡ 8cbc706f-1d3f-49ac-8f1c-d281b2318dc1
#=begin
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
end=#

# ╔═╡ 2491363c-4d1a-4baf-8c47-f85f550d59f8
md"""
### Projection on 2nd dimension
"""

# ╔═╡ fd5310bd-005a-4bf5-b6b7-79b137a101f7
#=begin
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
end=#

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

# ╔═╡ b872f155-d36f-4a2d-9741-342a108d6b5f
GasChromatographySystems.check_duration_modulation(sim[2][1], sim[4], PM, ratio)

# ╔═╡ daed9ba6-a44c-48a0-999b-5862838253c2
GasChromatographySystems.check_area(sim[2][1][8])

# ╔═╡ 27645404-48fd-488c-8b2e-66300f4b3303
md"""
## Traces 1st Modulation
"""

# ╔═╡ 09b5b754-9218-4784-977d-48759600037b
sim[3][1][1]

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
	for i=21:25#length(sim[3][1][3])
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

# ╔═╡ 8c156d25-1cd4-493e-ba6e-3371d9cdbd3f
sum(abs.(collect(skipmissing(compare.relΔtR1_percent))))/length(collect(skipmissing(compare.relΔtR1_percent)))

# ╔═╡ c190dd0f-6f09-4d18-9aa5-6b40beddfbf3
sum(abs.(collect(skipmissing(compare.relΔtR2_percent))))/length(collect(skipmissing(compare.relΔtR2_percent)))

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
	p_meas = Plots.heatmap(0.0:4.0:(size(chrom_mat)[1]-1)*4.0, 0.0:0.001:4.0, chrom_mat_mod', c=:jet1, colorbar=false, legend=:bottomright);
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

# ╔═╡ 5486bee0-cb8a-406f-ad00-937a8b985486
md"""
# Mod - rerun sim if result wrong
"""

# ╔═╡ 881bb34f-7061-408d-a64e-a4bbb672162c
function check_duration_modulation(pl_array, par, PM)
	#Δts = duration_in_module(pl_array, par)
	i_TM1 = 3
	i_TM2 = 5
	ok_TM1 = Array{Bool}(undef, length(par[i_TM1].sub))
	for i=1:length(par[i_TM1].sub)
		ii = GasChromatographySystems.common_index(pl_array[i_TM1], par[i_TM1].sub[i].CAS, join(split(par[i_TM1].sub[i].ann, ", ")[1:end-1], ", "))
		if fld(pl_array[i_TM1].tR[ii], PM) == fld(par[i_TM1].sub[i].t₀, PM)
			ok_TM1[i] = true
		else
			ok_TM1[i] = false
		end
	end
	#ok_TM2 = Array{Bool}(undef, length(Δts[5].Δt))
	#for i=1:length(Δts[5].Δt)
	#	if (Δts[5].Δt[i] > ratio*PM) && (Δts[5].Δt[i] < PM) # change this condition
	#		ok_TM2[i] = true
	#	else
	#		ok_TM2[i] = false
	#	end
	#end
	#index_TM1 = findall(ok_TM1.==false)
	#index_TM2 = findall(ok_TM2.==false)
	#TM1 = DataFrame(Index=index_TM1, Name=Δts[3].Name[index_TM1], CAS=Δts[3].CAS[index_TM1], tR=pl_array[3].tR[index_TM1], τR=pl_array[3].τR[index_TM1], Δt=Δts[3].Δt[index_TM1], Annotations=Δts[3].Annotations[index_TM1])
	#TM2 = DataFrame(Index=index_TM2, Name=Δts[5].Name[index_TM2], CAS=Δts[5].CAS[index_TM2], tR=pl_array[5].tR[index_TM2], τR=pl_array[5].τR[index_TM2], Δt=Δts[5].Δt[index_TM2], Annotations=Δts[5].Annotations[index_TM2])
	return ok_TM1#, TM2
end

# ╔═╡ 846a8397-f6d0-490c-857b-3d4cf691dc65
chck_dur = GasChromatographySystems.check_duration_modulation(sim[2][1], sim[4], PM, ratio)

# ╔═╡ c20db2f4-b5b8-4436-9d4b-df1333659d25
begin
	plotly()
	p_xt_check = Plots.plot(xlabel="time in s", ylabel="x in m", legend=true)
	p_τt_check = Plots.plot(xlabel="time in s", ylabel="τ in s", legend=true)
	p_Tt_check = Plots.plot(xlabel="time in s", ylabel="T in °C", legend=true)
	p_lnkt_check = Plots.plot(xlabel="time in s", ylabel="lnk", legend=true)
	first_step_TM1_z = Array{Float64}(undef, length(sim[3][1][3]))
	first_step_TM1_t = Array{Float64}(undef, length(sim[3][1][3]))
	for i=chck_dur[1].Index#length(sim[3][1][3])
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

# ╔═╡ b77246af-6887-4c1b-a6cc-17fbd9bc6e6d
chck_dur[1].Index # solutions on the 1st modulator spot, which elute to late (in a following modulation periode)

# ╔═╡ ee1afda3-4e74-403f-a42c-67a0363cf095
i = chck_dur[1].Index[end]

# ╔═╡ 84b9987c-746b-4576-b3a0-c4a4c2aa12ae
sim[4][3].sub[i].t₀

# ╔═╡ 19a6f833-ba99-40dd-9b80-3de7ae4fa330
sol_1 = GasChromatographySimulator.solving_odesystem_r(sim[4][3].col, sim[4][3].prog, sim[4][3].sub[i], sim[4][3].opt; dt=sys.modules[3].length*1e-7) # as in the original simulation

# ╔═╡ 031f3e45-9e38-4c8c-9eee-9393b520098d
sol_1.u[end]

# ╔═╡ ec6ab0ff-f4f1-45ff-b0ae-71ef0a4436e5
sol_2 = GasChromatographySimulator.solving_odesystem_r(sim[4][3].col, sim[4][3].prog, sim[4][3].sub[i], sim[4][3].opt; dt=sys.modules[3].length*1e-8) # decreased initial stepwidth

# ╔═╡ a6c7ed29-d369-42b7-a71f-4771fc0492f2
sol_2.u[end]

# ╔═╡ eddfdb96-0071-416b-9737-1a363be9dd42
typeof(sol_2)

# ╔═╡ fb3a3e22-3c5c-4079-9498-a5f5c2c119db
sol_2.prob

# ╔═╡ b79811dd-b0ec-4830-9a78-ab2a298ca8a4
sol_2.stats

# ╔═╡ 09272842-7f7c-4db0-b8a8-299b8c412119
sol_2.alg

# ╔═╡ 7a1e23c8-b431-4a01-b639-468e494bb4e5
sol_3 = GasChromatographySimulator.solving_odesystem_r(sim[4][3].col, sim[4][3].prog, sim[4][3].sub[i], sim[4][3].opt; dt=sys.modules[3].length*1e-10)

# ╔═╡ f1d9b432-e39c-4451-8d52-2e7ee5647115
sol_3.u[end]

# ╔═╡ dcb53d85-710e-4601-87d6-81c329d6ec75
sol_3.stats

# ╔═╡ d962563c-dccc-4997-9305-2d4fbf48840d
begin
	Plots.plot(sol_1)
	Plots.plot!(sol_2)
	Plots.plot!(sol_3)
end

# ╔═╡ 49dc357d-5c4c-4dd1-a92e-5a22f414feff
GasChromatographySimulator.solving_odesystem_r(sim[4][3].col, sim[4][3].prog, sim[4][3].sub[9], sim[4][3].opt; dt=sys.modules[3].length*1e-7, dtmax=sys.modules[3].length*1e-5) # as in the original simulation

# ╔═╡ 3c0dcfce-2d5c-4879-8da6-206bed7a15eb
md"""
## New `simulate_along_paths`-function
"""

# ╔═╡ 15e4e972-c2b5-4d93-96f8-4ab664793a93
function simulate_along_paths(sys, paths, par_sys; t₀=zeros(length(par_sys[1].sub)), τ₀=zeros(length(par_sys[1].sub)), nτ=6, abstolTM=1e-10, reltolTM=1e-8, algTM=Vern9(), refocus=falses(ne(sys.g)), τ₀_focus=zeros(length(par_sys[1].sub)), kwargsTM...)
	
	E = collect(edges(sys.g))
	peaklists = Array{Array{DataFrame,1}}(undef, length(paths))
	solutions = Array{Array{Any,1}}(undef, length(paths))
	path_pos = Array{String}(undef, length(paths))
	new_par_sys = Array{GasChromatographySimulator.Parameters}(undef, length(par_sys))
	visited_E = falses(length(E))

	# i -> path number
	# j -> segment/module number
	for i=1:length(paths)
		i_par = GasChromatographySystems.index_parameter(sys.g, paths[i])
		if GasChromatographySystems.path_possible(sys, paths[i]) == true
			path_pos[i] = "path is possible"
			peaklists_ = Array{DataFrame}(undef, length(i_par))
			solutions_ = Array{Any}(undef, length(i_par))
			#As_ = Array{DataFrame}(undef, length(i_par))
			for j=1:length(i_par)
				if (i>1) && (all(visited_E[1:i_par[j]].==true))
					# was the segment already simulated in a previous simulated path?
					# look in all previous paths for the correct result -> the simulation correlated to the same edge and where this edge is connected to only previouse visited edges
					i_path = 0
					i_edge = 0
					for k=1:i-1 # previous paths
						i_par_previous = GasChromatographySystems.index_parameter(sys.g, paths[k])
						if length(i_par_previous) < j
							j0 = length(i_par_previous)
						else
							j0 = j
						end
						if all(x->x in i_par_previous[1:j0], i_par[1:j]) == true # all edges up to j are the same between the two paths
							i_path = k
							i_edge = findfirst(i_par[j].==i_par_previous)
						end
					end
					# re-use the results
					peaklists_[j] = peaklists[i_path][i_edge]
					solutions_[j] = solutions[i_path][i_edge]
				else # new simulated segments
					if j == 1 # first segment, directly after injection, it is assumed to be a segment of type `ModuleColumn`
						new_par_sys[i_par[j]] = GasChromatographySystems.change_initial(par_sys[i_par[j]], t₀, τ₀)
						peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(new_par_sys[i_par[j]])
						peaklists_[j][!,:A] = ones(length(peaklists_[j].Name)) # add a relativ area factor, splitting of `A` at split points is not accounted for yet, is only used for the slicing of peaks at modulators
					elseif typeof(sys.modules[i_par[j]]) == GasChromatographySystems.ModuleTM
						# put this in a separate function
						if refocus[i_par[j]] == true
							τ₀=τ₀_focus
						else
							τ₀=peaklists_[j-1].τR # PM?
						end
						new_par_sys[i_par[j]], df_A = GasChromatographySystems.slicing(peaklists_[j-1], sys.modules[i_par[j]].PM, sys.modules[i_par[j]].ratio, sys.modules[i_par[j]].shift, par_sys[i_par[j]]; nτ=nτ, τ₀=τ₀, abstol=abstolTM, reltol=reltolTM, alg=algTM)
						if algTM == "simplifiedTM"
							peaklists_[j], solutions_[j] = approximate_modulator(new_par_sys[i_par[j]], df_A, sys.modules[i_par[j]].PM, sys.modules[i_par[j]].ratio, sys.modules[i_par[j]].shift)
						else
							sol = Array{Any}(undef, length(new_par_sys[i_par[j]].sub))
							for i_sub=1:length(new_par_sys[i_par[j]].sub)
								dt = new_par_sys[i_par[j]].col.L/1e6
								sol[i_sub] = GasChromatographySimulator.solving_odesystem_r(new_par_sys[i_par[j]].col, new_par_sys[i_par[j]].prog, new_par_sys[i_par[j]].sub[i_sub], new_par_sys[i_par[j]].opt; kwargsTM..., dt=dt)
								tR = sol[i_sub].u[end][1]
								while (fld(tR, sys.modules[i_par[j]].PM) > fld(new_par_sys[i_par[j]].sub[i_sub].t₀, sys.modules[i_par[j]].PM) && dt > eps()) || (sol[i_sub].retcode != ReturnCode.Success && dt > eps()) # solute elutes not in the same modulation periode or the solving failed
									dt = dt/10 # reduce initial step-width
									@warn "Retention time $(tR) surpasses modulation period, t₀ = $(new_par_sys[i_par[j]].sub[i_sub].t₀), PM = $(sys.modules[i_par[j]].PM). Initial step-width dt is decreased ($(dt))."
									sol[i_sub] = GasChromatographySimulator.solving_odesystem_r(new_par_sys[i_par[j]].col, new_par_sys[i_par[j]].prog, new_par_sys[i_par[j]].sub[i_sub], new_par_sys[i_par[j]].opt; kwargsTM..., dt=dt)
									tR = sol[i_sub].u[end][1]
								end
							end
							peaklists_[j] = GasChromatographySimulator.peaklist(sol, new_par_sys[i_par[j]])
							GasChromatographySystems.add_A_to_pl!(peaklists_[j], df_A)
							solutions_[j] = sol
						end
						if maximum(peaklists_[j].τR) > sys.modules[i_par[j]].PM
							return @warn "Peak width of focussed peaks > modulation period. Simulation is aborted."
						end
					else # ModuleColumn
						new_par_sys[i_par[j]] = GasChromatographySystems.change_initial(par_sys[i_par[j]], peaklists_[j-1])
						peaklists_[j], solutions_[j] = GasChromatographySimulator.simulate(new_par_sys[i_par[j]])
						GasChromatographySystems.add_A_to_pl!(peaklists_[j], peaklists_[j-1])
					end
				end
			end
			visited_E[i_par] .= true
			peaklists[i] = peaklists_
			solutions[i] = solutions_
		else # path not possible
			neg_flow_modules = sys.modules[findall(paths[i][findall(GasChromatographySystems.positive_flow(sys)[i_par].==false)].==edges(sys.g))]
			str_neg_flow = "path not possible: "
			for j=1:length(neg_flow_modules)
				str_neg_flow = str_neg_flow*"Flow in module >$(neg_flow_modules[j].name)< becomes negative during the program. "
			end
			path_pos[i] = str_neg_flow
		end
	end
	return path_pos, peaklists, solutions, new_par_sys
end

# ╔═╡ 441d1e8a-9acf-44c4-90eb-aad70e252607
[paths[1][1:4]]

# ╔═╡ 792b9be7-16af-4374-b5ad-5b2f84c7727a
sim__ = simulate_along_paths(sys, [paths[1][1:end]], par; algTM=algTM)

# ╔═╡ 0e8c5b92-2f28-423a-a6cb-739852d689ea
begin
	plotly()
	p_xt_check__ = Plots.plot(xlabel="time in s", ylabel="x in m", legend=true)
	p_τt_check__ = Plots.plot(xlabel="time in s", ylabel="τ in s", legend=true)
	p_Tt_check__ = Plots.plot(xlabel="time in s", ylabel="T in °C", legend=true)
	p_lnkt_check__ = Plots.plot(xlabel="time in s", ylabel="lnk", legend=true)
	for i=25:38#chck_dur[1].Index#length(sim[3][1][3])
		trace = GasChromatographySystems.traces(sim__[3][1][3], sim__[4][3], i)
		lbl = string(split(sim__[4][3].sub[i].ann, "_")[1], "_", sim__[4][3].sub[i].name)
		Plots.plot!(p_xt_check__, trace.t, trace.z, label=lbl, marker=:circle)
		Plots.plot!(p_τt_check__, trace.t, sqrt.(trace.τ²), label=lbl, marker=:circle)
		Plots.plot!(p_Tt_check__, trace.t, trace.T.-273.15, label=lbl, marker=:circle)
		Plots.plot!(p_lnkt_check__, trace.t, log.(trace.k), label=lbl, marker=:circle)
	end
	Plots.plot(p_xt_check__, p_τt_check__, p_Tt_check__, p_lnkt_check__, size=(1000,800))
end

# ╔═╡ 8f04614f-f019-4d02-8fc5-b98732274515
sim__[2][1][3]

# ╔═╡ 57f78673-a6ad-43f1-963a-bd93aaa1ae6d
Plots.plot(sim__[3][1][3][205])

# ╔═╡ 9ed995f9-9919-45f3-bc32-c41ec61cec65
sim__[3][1][3][205].stats

# ╔═╡ b5275005-aaf2-4b93-b3ed-c970e5351b5d
sim__[3][1][3][205].retcode

# ╔═╡ 0f417041-f3a1-43bd-bbb2-b4185f03df21
sim__[3][1][3][205].retcode != ReturnCode.Success

# ╔═╡ 6f347e4c-fd73-4b72-b59b-bde0f39c7a46
sim__[3][1][3][205]

# ╔═╡ 9597a2b3-e137-4c7b-9aef-304ab7fa4a4c
GasChromatographySystems.check_duration_modulation(sim__[2][1], sim__[4], PM, ratio)

# ╔═╡ edbfc946-ac0b-40f6-8c8b-669cb47a6d40
pl_GCxGC__ = GasChromatographySystems.peaklist_GCxGC(sim__[2][1][8], sim__[2][1][2], PM)

# ╔═╡ 14c496bf-867f-4e9a-b2b8-7732f6b69d7e
slice_mat__, t_D1__, t_D2__, c_slices__, t___, chrom_sliced_sum__, chrom_sliced__  = GasChromatographySystems.chrom2d(sim__[2][1][end], sys)

# ╔═╡ c2cb92d1-463a-4572-b327-102e74b7afb7
begin
	plotly()
	x1__ = floor(minimum(pl_GCxGC__.tR1*0.99))
	x2__ = ceil(maximum(pl_GCxGC__.tR1*1.01))
	y1__ = floor(minimum(pl_GCxGC__.tR2*0.99))
	y2__ = ceil(maximum(pl_GCxGC__.tR2*1.01))
	p_contour__ = Plots.contour(t_D1__[1:end-1], t_D2__[1], slice_mat__', levels=10, xlabel="tR1 in s", ylabel="tR2 in s")#, c=:jet1)
	Plots.scatter!(p_contour__, pl_GCxGC__.tR1, pl_GCxGC__.tR2, msize=3, shape=:diamond, c=:black, markeralpha=1)
	Plots.plot!(p_contour__, xlims=(x1__, x2__), ylims=(y1__, y2__))

	p_heatmap__ = Plots.heatmap(t_D1__[1:end-1], t_D2__[1], slice_mat__', levels=10, xlabel="tR1 in s", ylabel="tR2 in s", c=:jet1)
	Plots.scatter!(p_heatmap__, pl_GCxGC__.tR1, pl_GCxGC__.tR2, msize=3, shape=:diamond, c=:black, markeralpha=1)
	Plots.plot!(p_heatmap__, xlims=(x1__, x2__), ylims=(y1__, y2__))

	Plots.plot(p_contour__, p_heatmap__, size=(900,800), legend=false)
end

# ╔═╡ 8359f4af-4eb2-4f6e-b67b-502858041286
begin
	fit_envel__ = GasChromatographySystems.fit_envelope(sim__[2][1][8], sim__[2][1][2])
	fit_D1__ = GasChromatographySystems.fit_gauss_D1(sim__[2][1][8], sim__[2][1][2], PM)
	fit_D2__ = GasChromatographySystems.fit_gauss_D2(sim__[2][1][8], PM)
end

# ╔═╡ c9d900c5-2d02-4036-9c4d-d70ae0cfb2bf
fit_D1__

# ╔═╡ 92b36404-b1cc-48c6-8504-c53c0844cc7f
fit_envel__

# ╔═╡ 3841d718-5f12-4105-880b-0337bcbcc187
sim__[2][1][end]

# ╔═╡ 3788e180-41ac-49f4-876b-7deea5cf1e92
filter([:Name] => x -> x == "d-Valerolactam", sim__[2][1][4])

# ╔═╡ b906c42e-bbb0-480c-910e-adc88ad503c6
filter([:Name] => x -> x == "d-Valerolactam", sim__[2][1][end])

# ╔═╡ 8bb60845-9bb9-47c6-8a86-2ca955da93fe
sum(filter([:Name] => x -> x == "d-Valerolactam", sim__[2][1][end]).A)

# ╔═╡ 823c6edd-1abe-4e07-85ab-f7f4acff55b1
# substance 7 (d-Valerolactam) -> no peak in heatmap? -> gets splitted a second time in the 2nd modulator spot, the handling of this case seems to make no sense -> d-Valerolactam peaks are a magnitude broader!!!

# ╔═╡ c8a23d3e-8633-456a-ae34-372246008131
minimum(fit_D1__.tRs[7])

# ╔═╡ ac5bc5aa-ecb3-49ce-9cb5-c7879f976c66
begin
	plotly()
	#@. model_g(x,p) = p[3]/sqrt(2*π*p[2]^2)*exp(-(x-p[1])^2/(2*p[2]^2))
	p_envel_7 = GasChromatographySystems.collect_chrom(sim__[2][1], sys; markings=true)[end]
	t1__, t2__ = p_envel_7[1][1][:x_extrema]
	t__ = t1__:(t2__-t1__)/1000:t2__
	ymax__ = Array{Float64}(undef, length(fit_envel__.Name))
	for i=7:7#1:length(fit_envel__.Name)
		ymax__[i] = model_g(fit_envel__.fits[i].param[1], fit_envel__.fits[i].param)
		Plots.scatter!(p_envel_7, fit_envel__.tRs[i], fit_envel__.heights[i], msize=2, c=i+1)
		Plots.plot!(p_envel_7, t__, model_g(t__, fit_envel__.fits[i].param), c=i+1, linestyle=:dash)
		Plots.plot!(p_envel_7, [fit_envel__.fits[i].param[1], fit_envel__.fits[i].param[1]], [0.0, ymax__[i]], c=i+1, linestyle=:dot)
		Plots.scatter!(p_envel_7, (fit_envel__.fits[i].param[1], ymax__[i]), c=i+1, shape=:diamond)
	end
	Plots.plot!(p_envel_7, xticks=0:4:3000, legend=false, ylims=(-0.01*maximum(ymax__), 1.1*maximum(ymax__)), title="Enveloping curve")
	p_envel_7
end

# ╔═╡ 84389904-9957-4135-b2e5-d4b6d98bd9a7
begin
	plotly()
	p_proj1stD_7 = Plots.plot(legend=false, xlabel="time in s")

	mins_ = minimum.(fit_D1__.tRs)
	maxs_ = maximum.(fit_D1__.tRs)
	for i=7:7#1:length(fit_D1__.Name)
		ymax = model_g(fit_D1__.fits[i].param[1], fit_D1__.fits[i].param)
		for j=1:length(fit_D1__.tRs[i]) 
			Plots.plot!(p_proj1stD_7, [fit_D1__.tRs[i][j], fit_D1__.tRs[i][j]], [0.0, fit_envel__.heights[i][j]], c=i+1, linewidth=10, linealpha=0.5)
		end
		Plots.scatter!(p_proj1stD_7, fit_D1__.tRs[i], fit_D1__.heights[i], msize=2, c=i+1)

		tt = mins_[i]:(maxs_[i]-mins_[i])/1000.0:maxs_[i]
		Plots.plot!(p_proj1stD_7, tt, model_g(tt, fit_D1__.fits[i].param), c=i+1, linestyle=:dash)

		Plots.plot!(p_proj1stD_7, [fit_D1__.fits[i].param[1], fit_D1__.fits[i].param[1]], [0.0, ymax], c=i+1, linestyle=:dot)
		Plots.scatter!(p_proj1stD_7, (fit_D1__.fits[i].param[1], ymax), c=i+1, shape=:diamond)
	end
	Plots.plot!(p_proj1stD_7, xticks=0:4:3000, legend=false, title="Projection 1st D")
	p_proj1stD_7
end

# ╔═╡ b2569687-0616-4f2f-bef3-7fdb20985b34
# substance 8 (Propiophenone) -> long tail in heatmap? -> artfact from substance 7 (d-Valerolactam) 

# ╔═╡ 5309cbf5-c4af-4a87-99fc-c5e0ae4bf970
begin
	plotly()
	p_proj1stD_8 = Plots.plot(legend=false, xlabel="time in s")

	mins_8 = minimum.(fit_D1__.tRs)
	maxs_8 = maximum.(fit_D1__.tRs)
	for i=8:8#1:length(fit_D1__.Name)
		ymax = model_g(fit_D1__.fits[i].param[1], fit_D1__.fits[i].param)
		for j=1:length(fit_D1__.tRs[i]) 
			Plots.plot!(p_proj1stD_8, [fit_D1__.tRs[i][j], fit_D1__.tRs[i][j]], [0.0, fit_envel__.heights[i][j]], c=i+1, linewidth=10, linealpha=0.5)
		end
		Plots.scatter!(p_proj1stD_8, fit_D1__.tRs[i], fit_D1__.heights[i], msize=2, c=i+1)

		tt = mins_8[i]:(maxs_8[i]-mins_8[i])/1000.0:maxs_8[i]
		Plots.plot!(p_proj1stD_8, tt, model_g(tt, fit_D1__.fits[i].param), c=i+1, linestyle=:dash)

		Plots.plot!(p_proj1stD_8, [fit_D1__.fits[i].param[1], fit_D1__.fits[i].param[1]], [0.0, ymax], c=i+1, linestyle=:dot)
		Plots.scatter!(p_proj1stD_8, (fit_D1__.fits[i].param[1], ymax), c=i+1, shape=:diamond)
	end
	Plots.plot!(p_proj1stD_8, xticks=0:4:3000, legend=false, title="Projection 1st D")
	p_proj1stD_8
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
# ╠═2301a74c-2290-4de6-992d-bc887b66269e
# ╠═2b41aa34-d994-440b-a8e1-dd6337008a48
# ╠═80cd5d72-7629-45dc-a10e-53b3db29f21b
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
# ╠═6d1ff697-e2d5-4582-b932-7a47c0c0256c
# ╠═bd810db6-ed39-482e-8107-af3169ff09be
# ╟─7a81859f-7c0c-4115-9cca-367753898427
# ╠═8cbc706f-1d3f-49ac-8f1c-d281b2318dc1
# ╟─2491363c-4d1a-4baf-8c47-f85f550d59f8
# ╠═fd5310bd-005a-4bf5-b6b7-79b137a101f7
# ╟─7739a83c-d001-47cd-aee2-36917106c2ad
# ╠═582499db-7b5a-4054-92c2-d0e2a97cd91f
# ╟─5f67ab8d-3b3d-4bbd-81c8-04eec158ef6e
# ╟─f268df70-af5a-482f-85e6-84216681ecf4
# ╠═e5d4fb3a-f05f-4533-9666-08be32ef19ef
# ╟─93e2068d-760e-403e-a6fa-984a237ad6c3
# ╠═b872f155-d36f-4a2d-9741-342a108d6b5f
# ╠═daed9ba6-a44c-48a0-999b-5862838253c2
# ╟─27645404-48fd-488c-8b2e-66300f4b3303
# ╠═c20db2f4-b5b8-4436-9d4b-df1333659d25
# ╠═09b5b754-9218-4784-977d-48759600037b
# ╟─d4c0a1c4-68ef-408a-9919-1bd737b6c229
# ╠═2737c1e7-6239-41d7-8868-242bb881f436
# ╟─271bb7a6-9a98-4336-9571-c131cc125486
# ╟─b7c2070b-16bf-41f5-b021-160e54a53a21
# ╠═189034b6-09e5-4e6b-80b3-e442ad30705c
# ╟─3ed4ad21-c8e4-428a-b6d3-b37cfdbd0778
# ╠═8c156d25-1cd4-493e-ba6e-3371d9cdbd3f
# ╠═c190dd0f-6f09-4d18-9aa5-6b40beddfbf3
# ╟─fb7236b5-98fd-4aac-86ed-5e7118b19ed6
# ╠═a2cc7971-bd69-4623-86a1-7e7d43c6b05c
# ╟─6cf28744-7748-4d88-bfd7-f8e428732228
# ╠═e82034ef-57fa-4f78-b663-7465a97ff127
# ╟─b28b299c-1c2f-4638-a91f-a970c92e4f5a
# ╠═847d5b23-967f-4316-941f-fde283200b36
# ╠═5486bee0-cb8a-406f-ad00-937a8b985486
# ╠═881bb34f-7061-408d-a64e-a4bbb672162c
# ╠═846a8397-f6d0-490c-857b-3d4cf691dc65
# ╠═b77246af-6887-4c1b-a6cc-17fbd9bc6e6d
# ╠═ee1afda3-4e74-403f-a42c-67a0363cf095
# ╠═84b9987c-746b-4576-b3a0-c4a4c2aa12ae
# ╠═19a6f833-ba99-40dd-9b80-3de7ae4fa330
# ╠═031f3e45-9e38-4c8c-9eee-9393b520098d
# ╠═ec6ab0ff-f4f1-45ff-b0ae-71ef0a4436e5
# ╠═a6c7ed29-d369-42b7-a71f-4771fc0492f2
# ╠═eddfdb96-0071-416b-9737-1a363be9dd42
# ╠═fb3a3e22-3c5c-4079-9498-a5f5c2c119db
# ╠═b79811dd-b0ec-4830-9a78-ab2a298ca8a4
# ╠═09272842-7f7c-4db0-b8a8-299b8c412119
# ╠═7a1e23c8-b431-4a01-b639-468e494bb4e5
# ╠═f1d9b432-e39c-4451-8d52-2e7ee5647115
# ╠═dcb53d85-710e-4601-87d6-81c329d6ec75
# ╠═d962563c-dccc-4997-9305-2d4fbf48840d
# ╠═49dc357d-5c4c-4dd1-a92e-5a22f414feff
# ╠═3c0dcfce-2d5c-4879-8da6-206bed7a15eb
# ╠═15e4e972-c2b5-4d93-96f8-4ab664793a93
# ╠═441d1e8a-9acf-44c4-90eb-aad70e252607
# ╠═792b9be7-16af-4374-b5ad-5b2f84c7727a
# ╠═0e8c5b92-2f28-423a-a6cb-739852d689ea
# ╠═8f04614f-f019-4d02-8fc5-b98732274515
# ╠═57f78673-a6ad-43f1-963a-bd93aaa1ae6d
# ╠═9ed995f9-9919-45f3-bc32-c41ec61cec65
# ╠═b5275005-aaf2-4b93-b3ed-c970e5351b5d
# ╠═0f417041-f3a1-43bd-bbb2-b4185f03df21
# ╠═6f347e4c-fd73-4b72-b59b-bde0f39c7a46
# ╠═9597a2b3-e137-4c7b-9aef-304ab7fa4a4c
# ╠═edbfc946-ac0b-40f6-8c8b-669cb47a6d40
# ╠═14c496bf-867f-4e9a-b2b8-7732f6b69d7e
# ╠═c2cb92d1-463a-4572-b327-102e74b7afb7
# ╠═8359f4af-4eb2-4f6e-b67b-502858041286
# ╠═c9d900c5-2d02-4036-9c4d-d70ae0cfb2bf
# ╠═92b36404-b1cc-48c6-8504-c53c0844cc7f
# ╠═3841d718-5f12-4105-880b-0337bcbcc187
# ╠═3788e180-41ac-49f4-876b-7deea5cf1e92
# ╠═b906c42e-bbb0-480c-910e-adc88ad503c6
# ╠═8bb60845-9bb9-47c6-8a86-2ca955da93fe
# ╠═823c6edd-1abe-4e07-85ab-f7f4acff55b1
# ╠═c8a23d3e-8633-456a-ae34-372246008131
# ╠═ac5bc5aa-ecb3-49ce-9cb5-c7879f976c66
# ╠═84389904-9957-4135-b2e5-d4b6d98bd9a7
# ╠═b2569687-0616-4f2f-bef3-7fdb20985b34
# ╠═5309cbf5-c4af-4a87-99fc-c5e0ae4bf970
# ╠═954597a2-b847-4cf0-8c1c-ff4b4b4b670b
