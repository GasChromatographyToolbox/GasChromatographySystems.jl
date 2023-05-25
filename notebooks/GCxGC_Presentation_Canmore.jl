### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ e1e9eba6-f3fa-11ed-0472-51879f762c01
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

# ╔═╡ 6c6ceef1-a56c-4d59-b8c8-46243d2cdafc
md"""
# Simplified GCxGC simulation
** for the presentation at the 20th GCxGC Symposium in Canmore, Canada, 2023.**

Using determined retention parameters with RetentionParameterEstimator.jl from temperature programmed measurements of two columns to simulate simplified GCxGC measurements and compare them to measurements.
"""

# ╔═╡ f2d311f8-8d4b-4261-96aa-4348e22c24c4
md"""
## Definition of the GCxGC system
"""

# ╔═╡ b76c3755-8191-4a83-8f5a-39c71cda85ad
begin # Placeholder, exact programs from Tillman Brehmer
	# also pressure values from the GC for the constant flow program
	tsteps_, Tsteps_ = GasChromatographySimulator.conventional_program([50.0, 1.0, 3.0, 200.0, 0.0, 10.0, 225.0, 10.0]) 
	TP1 = GasChromatographySystems.TemperatureProgram(tsteps_, Tsteps_)
end

# ╔═╡ 1a0489e2-7832-4d0a-85fa-187332af9c58
# exact lengths not known
# stat Phase 2nd dim placeholder
# for pressure program the system has to be defined manually
sys = GasChromatographySystems.GCxGC_TM_simp(25.0, 0.25, 0.25, "ZB1ms", TP1, 3.0, 0.1, 0.1, "Stabilwax", TP1, 0.8, NaN, 0.0; opt=GasChromatographySystems.Options(ng=true))

# ╔═╡ ce89ef57-0e36-46ae-90de-3c9e0da922f6
GasChromatographySystems.plot_graph_with_flow(sys, 0; lay=SquareGrid(cols=3), node_size=90, nlabels_fontsize=20, elabels_fontsize=20, elabels_distance=28)

# ╔═╡ 69b5feae-a0af-4c88-af86-54662640af8c
GasChromatographySystems.plot_graph_with_flow(sys, sum(tsteps_); lay=SquareGrid(cols=3), node_size=90, nlabels_fontsize=20, elabels_fontsize=20, elabels_distance=28)

# ╔═╡ c7b914ea-33b5-4cdc-bb6b-0442774cdf30
GasChromatographySystems.plot_pressure_over_time(sys)

# ╔═╡ f5eb0e48-fd3b-4cfe-a3e0-152c61ae0e13
md"""
## Retention data
"""

# ╔═╡ 8a361ebf-aead-46d1-aa97-26503c70e696
begin
	db_m1 = DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/Estimator/GCsim_renamed.csv", header=1, silencewarnings=true))
	db_m1.Tchar = db_m1.Tchar .- 273.15
	db_m1
end

# ╔═╡ 34fb18b0-2109-4f48-b086-b3be98dad24a
begin
	db_m2 = DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/Estimator/GCsim_d_renamed.csv", header=1, silencewarnings=true))
	db_m2.Tchar = db_m2.Tchar .- 273.15
	db_m2
end

# ╔═╡ b15c6562-77f7-42fd-bdd8-003a0a68ff50
solutes = GasChromatographySystems.common_solutes(db_m2, sys).Name

# ╔═╡ f02e43e1-bff9-4da5-bf49-72c340dc8465
md"""
## Parameters
"""

# ╔═╡ 142be678-b698-4f57-b35f-d9926be157f5
par_m1 = GasChromatographySystems.graph_to_parameters(sys, db_m1, solutes)

# ╔═╡ d5a64ddd-2503-4c3a-9d96-71c61b52897d
par_m2 = GasChromatographySystems.graph_to_parameters(sys, db_m2, solutes)

# ╔═╡ a7814d64-9559-4749-8c69-833fd263e1bb
md"""
## Simulations
"""

# ╔═╡ 754de3c1-e69b-4594-87cf-e6e08a0b6097
paths = GasChromatographySystems.all_paths(sys.g, 1)[2]

# ╔═╡ 96436272-484a-4288-b1ee-f8fa42bd257f
pl_m1 = GasChromatographySystems.simulate_along_paths(sys, paths, par_m1; refocus=trues(ne(sys.g)), τ₀_focus=0.1.*ones(length(par_m1[1].sub)))[2]

# ╔═╡ dda3ef3c-53d8-4f09-aa04-484ecd10f291
pl_m2 = GasChromatographySystems.simulate_along_paths(sys, paths, par_m2; refocus=trues(ne(sys.g)), τ₀_focus=0.1.*ones(length(par_m2[1].sub)))[2]

# ╔═╡ 88a5f050-1589-4a8b-8e2f-436d12ff5d97
md"""
## Measured Chromatogram
"""

# ╔═╡ 3cf4b576-75dc-4df9-b89b-7723a5b21dde
begin
	meas = DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/Estimator/meas_GCxGC.csv", header=1, silencewarnings=true))
end

# ╔═╡ 1e76c773-5d1d-4319-b093-fbba55a2c588
md"""
# End
"""

# ╔═╡ d32158cf-ac99-4ba3-8c44-98434a64a819
function plot_GCxGC(pl_GCxGC, sys; categories = String[])
	t1end = sum(GasChromatographySystems.common_timesteps(sys))
	t2end = ceil(maximum(pl_GCxGC.tR2.-pl_GCxGC.tR1))
	p_gcxgc = Plots.scatter(pl_GCxGC.tR1, pl_GCxGC.tR2.-pl_GCxGC.tR1, xlims=(0.0,t1end), ylims=(0.0, t2end), label="", xlabel=string(sys.modules[1].stationary_phase, ", t¹ in s"), ylabel=string(sys.modules[2].stationary_phase, ", t² in s"))
	if !isempty(categories)
		for i=1:length(categories) 
			pl_f = filter([:Cat] => x -> occursin(categories[i], x), pl_GCxGC)
			Plots.scatter!(p_gcxgc, pl_f.tR1, pl_f.tR2.-pl_f.tR1, label=categories[i], legend=:bottomright)
		end
	end
	return p_gcxgc
end

# ╔═╡ 6aa81fa6-155d-4d54-9e5d-530b35e1a8ba
function plot_GCxGC_contour_(pl_GCxGC; split_ratio=ones(length(pl_GCxGC.Name)), mod_τ1=1.0)
	t¹ = 0.9*minimum(pl_GCxGC.tR1):(1.1*maximum(pl_GCxGC.tR1)-0.9*minimum(pl_GCxGC.tR1))/1000:1.1*maximum(pl_GCxGC.tR1)
	t² = 0.9*minimum(pl_GCxGC.tR2.-pl_GCxGC.tR1):(1.1*maximum(pl_GCxGC.tR2.-pl_GCxGC.tR1)-0.9*minimum(pl_GCxGC.tR2.-pl_GCxGC.tR1))/1000:1.1*maximum(pl_GCxGC.tR2.-pl_GCxGC.tR1)
	chrom2D = Array{Float64}(undef, length(t¹), length(t²), length(pl_GCxGC.Name))
	for j=1:length(t²)
		for i=1:length(t¹)
			for k=1:length(pl_GCxGC.Name)
				chrom2D[i,j,k] = split_ratio[k]*exp(-((t¹[i]-pl_GCxGC.tR1[k])^2/(2*(mod_τ1*pl_GCxGC.τR1[k])^2)+(t²[j]-(pl_GCxGC.tR2[k]-pl_GCxGC.tR1[k]))^2/(2*pl_GCxGC.τR2[k]^2)))
			end
		end
	end
	#p_2D = Plots.contourf(t¹, t², sum(chrom2D, dims=3)[:,:,1]', fill=true, legend=false, xlabel="tR1 in s", ylabel="tR2 in s", c=:turbo, levels=40)
	return t¹, t², sum(chrom2D, dims=3)[:,:,1]'
end

# ╔═╡ 989ef538-976e-485d-8727-dbc9d6ea182e
function peaklist_GCxGC_with_categories(pl_D1, pl_D2, db)
	pl = GasChromatographySystems.peaklist_GCxGC(pl_D1, pl_D2)
	cat = Array{String}(undef, length(pl.Name))
	for i=1:length(pl.Name)
		ii = findfirst(pl.Name[i].==db.Name)
		cat[i] = ""
		for j=1:length(findall(occursin.("Cat", names(db))))
			if !ismissing(db[ii, Symbol("Cat_$(j)")])
				if j==1
					cat[i] = cat[i]*db[ii, Symbol("Cat_$(j)")]
				else
					cat[i] = cat[i]*", "*db[ii, Symbol("Cat_$(j)")]
				end
			end
		end
	end
	pl[!, :Cat] = cat
	return pl
end

# ╔═╡ 532c1405-3577-4de4-a88c-770fefb88925
pl_GCxGC_m1 = peaklist_GCxGC_with_categories(pl_m1[1][1], pl_m1[1][2], db_m1)

# ╔═╡ 535da8c3-0b16-42dd-8092-6665e6f7e552
begin
	plotly()
	p_GCxGC_m1 = plot_GCxGC(pl_GCxGC_m1, sys; categories = ["alcohol", "terpene", "phenone", "ketone"])
	#Plots.plot!(p_GCxGC, legend=:topleft)
	p_GCxGC_m1
end

# ╔═╡ 083a435a-e726-4817-a049-594a916a82d6
t¹_m1, t²_m1, chrom_m1 = plot_GCxGC_contour_(pl_GCxGC_m1; mod_τ1=3.0)

# ╔═╡ b44e068a-6850-42ed-8eca-4bb657d862e3
pl_GCxGC_m2 = peaklist_GCxGC_with_categories(pl_m2[1][1], pl_m2[1][2], db_m2)

# ╔═╡ 8cabade3-0e57-4809-abca-b3ca6b601ea0
begin
	plotly()
	p_GCxGC_m2 = plot_GCxGC(pl_GCxGC_m2, sys; categories = ["alcohol", "terpene", "phenone", "ketone"])
	#Plots.plot!(p_GCxGC, legend=:topleft)
	p_GCxGC_m2
end

# ╔═╡ b49415b6-934f-44f9-9b89-c8244ef0f43b
t¹_m2, t²_m2, chrom_m2 = plot_GCxGC_contour_(pl_GCxGC_m2; mod_τ1=3.0)

# ╔═╡ a5be250f-eb47-4c65-8d7d-3227186333d7
begin
	plotly()
	p_GCxGC_m2_ = Plots.heatmap(t¹_m2, t²_m2, chrom_m2.^(1//1), legend=false, xlabel="tR1 in s", ylabel="tR2 in s", c=:jet1, dpi=500)#, xlims=(0.0, 1700.0), ylims=(0.0, 5.0), grid=false)
	#gui()
	Plots.scatter!(p_GCxGC_m2_, pl_GCxGC_m1.tR1, pl_GCxGC_m1.tR2.-pl_GCxGC_m1.tR1, label="m1")
	p_GCxGC_m2_
end

# ╔═╡ f7f5f8f0-cc7d-4e52-8781-1ff4c70152fd
begin
	plotly()
	p_GCxGC_m1_ = Plots.heatmap(t¹_m1, t²_m1, chrom_m1.^(1//1), legend=false, xlabel="tR1 in s", ylabel="tR2 in s", c=:jet1, dpi=500)#, xlims=(0.0, 1700.0), ylims=(0.0, 5.0), grid=false)
	#gui()
	Plots.scatter!(p_GCxGC_m1_, pl_GCxGC_m2.tR1, pl_GCxGC_m2.tR2.-pl_GCxGC_m2.tR1, label="m2")
	p_GCxGC_m1_
end

# ╔═╡ 504093c0-7478-4636-8f38-45efa67467d6
begin
	# comparison measurment simulation
	
	comp = DataFrame(Name=pl_GCxGC_m2.Name, tR1_meas=pl_GCxGC_m2.tR1, tR1_sim=pl_GCxGC_m2.tR1, ΔtR1=pl_GCxGC_m2.tR1, tR2_meas=pl_GCxGC_m2.tR2.-pl_GCxGC_m2.tR1, tR2_sim=pl_GCxGC_m2.tR2.-pl_GCxGC_m2.tR1, ΔtR2=pl_GCxGC_m2.tR2.-pl_GCxGC_m2.tR1)
	for i=1:length(comp.Name)
		ii = findfirst(comp.Name[i].==meas.Name)
		if ismissing(meas.tR1[ii])
			comp[i, :tR1_meas] = NaN
			comp[i, :tR2_meas] = NaN
			comp[i, :ΔtR1] = NaN
			comp[i, :ΔtR2] = NaN
		else
			comp[i, :tR1_meas] = meas.tR1[ii]
			comp[i, :tR2_meas] = meas.tR2[ii]
			comp[i, :ΔtR1] = meas.tR1[ii] - comp[i, :tR1_sim]
			comp[i, :ΔtR2] = meas.tR2[ii] - comp[i, :tR2_sim]
		end
	end
	comp
end

# ╔═╡ Cell order:
# ╠═e1e9eba6-f3fa-11ed-0472-51879f762c01
# ╠═6c6ceef1-a56c-4d59-b8c8-46243d2cdafc
# ╠═f2d311f8-8d4b-4261-96aa-4348e22c24c4
# ╠═b76c3755-8191-4a83-8f5a-39c71cda85ad
# ╠═1a0489e2-7832-4d0a-85fa-187332af9c58
# ╠═ce89ef57-0e36-46ae-90de-3c9e0da922f6
# ╠═69b5feae-a0af-4c88-af86-54662640af8c
# ╠═c7b914ea-33b5-4cdc-bb6b-0442774cdf30
# ╠═f5eb0e48-fd3b-4cfe-a3e0-152c61ae0e13
# ╠═8a361ebf-aead-46d1-aa97-26503c70e696
# ╠═34fb18b0-2109-4f48-b086-b3be98dad24a
# ╠═b15c6562-77f7-42fd-bdd8-003a0a68ff50
# ╠═f02e43e1-bff9-4da5-bf49-72c340dc8465
# ╠═142be678-b698-4f57-b35f-d9926be157f5
# ╠═d5a64ddd-2503-4c3a-9d96-71c61b52897d
# ╠═a7814d64-9559-4749-8c69-833fd263e1bb
# ╠═754de3c1-e69b-4594-87cf-e6e08a0b6097
# ╠═96436272-484a-4288-b1ee-f8fa42bd257f
# ╠═dda3ef3c-53d8-4f09-aa04-484ecd10f291
# ╠═532c1405-3577-4de4-a88c-770fefb88925
# ╠═b44e068a-6850-42ed-8eca-4bb657d862e3
# ╠═535da8c3-0b16-42dd-8092-6665e6f7e552
# ╠═8cabade3-0e57-4809-abca-b3ca6b601ea0
# ╠═083a435a-e726-4817-a049-594a916a82d6
# ╠═b49415b6-934f-44f9-9b89-c8244ef0f43b
# ╠═f7f5f8f0-cc7d-4e52-8781-1ff4c70152fd
# ╠═a5be250f-eb47-4c65-8d7d-3227186333d7
# ╠═88a5f050-1589-4a8b-8e2f-436d12ff5d97
# ╠═3cf4b576-75dc-4df9-b89b-7723a5b21dde
# ╠═504093c0-7478-4636-8f38-45efa67467d6
# ╠═1e76c773-5d1d-4319-b093-fbba55a2c588
# ╠═d32158cf-ac99-4ba3-8c44-98434a64a819
# ╠═6aa81fa6-155d-4d54-9e5d-530b35e1a8ba
# ╠═989ef538-976e-485d-8727-dbc9d6ea182e
