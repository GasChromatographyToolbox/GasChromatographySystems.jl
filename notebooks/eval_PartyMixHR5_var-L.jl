### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 7656e75e-2e47-11ee-31b0-f9a794c245bb
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

# ╔═╡ dc716bc4-153c-4aef-9d00-55f99ab3f589
md"""
# Evaluation of variation of length 1st D column

Using the simulation of GCxGC with looped thermal modulation with a 
- smoothed rectangular temperature function over time at the modulation points
- uniform temperature over the spot
- the cooling temperature is relative (reduction of the programmed oven temperature by the value of the cooling temperature)
"""

# ╔═╡ c7b83c87-2a13-487c-bc45-6915317201d0
md"""
## Measurement
"""

# ╔═╡ 10353347-1f2c-4f00-9a7f-cf5518eec85d
meas = DataFrame(CSV.File("/Users/janleppert/Documents/sciebo/GCsim/GCxGC/Simulation/PartyMixHR5_RT.csv", header=1, silencewarnings=true, stringtype=String))

# ╔═╡ 4c15c7db-01dc-4de5-a60b-5fd56c503926
md"""
## Simulation results
"""

# ╔═╡ 67df7e29-48cc-4a9b-b752-c14b6aa0128e
md"""
### Shift = 3.5 s
"""

# ╔═╡ 41d3693b-bfb4-4732-9a59-d92640f1ab79
sim_shift35 = DataFrame(CSV.File("/Users/janleppert/Documents/GitHub/GasChromatographySystems/scripts/var-L1_PartyMixHR5_shift35.csv", header=1, silencewarnings=true, stringtype=String))

# ╔═╡ a474cf70-6b92-4757-bb40-e9e2e3dbfb32
md"""
### Shift = 0.0 s
"""

# ╔═╡ 9fd0b2b7-c423-49fc-8e89-6d8ec5c0e8ea
sim_shift0 = DataFrame(CSV.File("/Users/janleppert/Documents/GitHub/GasChromatographySystems/scripts/var-L1_PartyMixHR5_shift0.csv", header=1, silencewarnings=true, stringtype=String))

# ╔═╡ 0180cdbd-73d0-43a0-b583-766e18942edb
md"""
## Evaluation
"""

# ╔═╡ 35d87a01-a90d-4f9e-9a63-9c1393ae4492
function eval_varL1(sim)
	name = unique(sim.Name)
	L1 = unique(sim.L1)
	results = Array{DataFrame}(undef, length(name))
	p_tR1 = Array{Plots.Plot}(undef, length(name))
	p_tR2 = Array{Plots.Plot}(undef, length(name))
	p_tR1s = Plots.plot(xlabel="L1 in m", ylabel="1st dim retention time tR1 in s")
	p_tR2s = Plots.plot(xlabel="L1 in m", ylabel="2nd dim retention time tR2 in s")
	p_ΔtR1s = Plots.plot(xlabel="L1 in m", ylabel="1st dim difference meas-sim ΔtR1 in s")
	p_ΔtR2s = Plots.plot(xlabel="L1 in m", ylabel="2nd dim difference meas-sim ΔtR2 in s")
	L1_opt_tR1 = fill(NaN, length(name))
	L1_opt_tR2 = fill(NaN, length(name))
	for i=1:length(name)
		results[i] = sort(filter([:Name] => x -> x == name[i], sim), :L1)
		p_tR1[i] = Plots.plot(results[i].L1, results[i].tR1, xlabel="L1 in m", ylabel="1st dim retention time tR1 in s", title=results[i].Name[1], label="")
		p_tR2[i] = Plots.plot(results[i].L1, results[i].tR2, xlabel="L1 in m", ylabel="2nd dim retention time tR2 in s", title=results[i].Name[1], label="")
		Plots.plot!(p_tR1s, results[i].L1, results[i].tR1, label=results[i].Name[1])
		Plots.plot!(p_tR2s, results[i].L1, results[i].tR2, label=results[i].Name[1])
		if !ismissing(meas.tR1[findfirst(name[i].==meas.Name)]) && !isnothing(meas.tR1[findfirst(name[i].==meas.Name)])
			Plots.plot!(p_ΔtR1s, results[i].L1, meas.tR1[findfirst(name[i].==meas.Name)] .- results[i].tR1, label=results[i].Name[1])
			Plots.plot!(p_ΔtR2s, results[i].L1, meas.tR2[findfirst(name[i].==meas.Name)] .- results[i].tR2, label=results[i].Name[1])

			itp_tR1 = LinearInterpolation((-1.0.*(meas.tR1[findfirst(name[i].==meas.Name)] .- results[i].tR1), ), L1, extrapolation_bc=Flat())
			
			itp_tR2 = LinearInterpolation((meas.tR2[findfirst(name[i].==meas.Name)] .- results[i].tR2, ), L1, extrapolation_bc=Flat())

			L1_opt_tR1[i] = itp_tR1(0.0)
			L1_opt_tR2[i] = itp_tR2(0.0)
		end
	end
	return p_tR1, p_tR2, p_tR1s, p_tR2s, p_ΔtR1s, p_ΔtR2s, L1_opt_tR1, L1_opt_tR2
end

# ╔═╡ ef286a7a-246c-4e2c-a257-6ab301ad8590
res_shift35 = eval_varL1(sim_shift35)

# ╔═╡ 2fc80039-5174-4e93-b997-b5d5cdd8b3fa
res_shift0 = eval_varL1(sim_shift0)

# ╔═╡ ef1942c3-8076-4830-8c5c-aee353ca1add
res_shift35[5]

# ╔═╡ 6658863c-e278-48fd-9f5c-d8a7ac7eacfd
opt_tR1_shift35 = sum(res_shift35[7][isnan.(res_shift35[7]).==false])/length(res_shift35[7][isnan.(res_shift35[7]).==false]), maximum(res_shift35[7][isnan.(res_shift35[7]).==false]), minimum(res_shift35[7][isnan.(res_shift35[7]).==false])

# ╔═╡ fa47c6f3-806d-43cd-9983-b2ca27534be3
md"""
### Optima regarding tR1 - shift = 3.5s

The optima is where the difference between measured and simulated 1st dimension retention time is 0.

Mean value of the optima: $(opt_tR1_shift35[1]) m

Minima: $(opt_tR1_shift35[3]) m

Maxima: $(opt_tR1_shift35[2]) m

The mean value of 29.5 m is the same used to determine the optimal shift of 3.5s. 

"""

# ╔═╡ 5507f392-35b2-43cd-b1f8-0f4b8e92eb7f
res_shift0[5]

# ╔═╡ 40175e8d-606c-4139-a06e-b753005f61f9
opt_tR1_shift0 = sum(res_shift0[7][isnan.(res_shift0[7]).==false])/length(res_shift0[7][isnan.(res_shift0[7]).==false]), maximum(res_shift0[7][isnan.(res_shift0[7]).==false]), minimum(res_shift0[7][isnan.(res_shift0[7]).==false])

# ╔═╡ d18c7f46-262b-472b-bffe-91a63c72a7d4
md"""
### Optima regarding tR1 - shift = 0.0s

The optima is where the difference between measured and simulated 1st dimension retention time is 0.

Mean value of the optima: $(opt_tR1_shift0[1]) m

Minima: $(opt_tR1_shift0[3]) m

Maxima: $(opt_tR1_shift0[2]) m

The mean value of 29.5 m is the same used to determine the optimal shift of 3.5s. 

Plot of ΔtR1 over L1 looks similar to shift = 3.5s and values are similar.
"""

# ╔═╡ 7feb935b-aba6-4c7b-8138-a45777978c55
Plots.plot!(res_shift35[6], legend=:topleft)

# ╔═╡ e281594a-0faf-4f50-915c-5bca8abb6969
opt_tR2_shift35 = sum(res_shift35[8][isnan.(res_shift35[8]).==false])/length(res_shift35[8][isnan.(res_shift35[8]).==false]), maximum(res_shift35[8][isnan.(res_shift35[8]).==false]), minimum(res_shift35[8][isnan.(res_shift35[8]).==false])

# ╔═╡ ba11ac4b-08c8-44d1-a896-80fd9bf9e921
md"""
### Optima regarding tR2 - shift = 3.5 s

The the optima is where the difference between measured and simulated 2nd dimension retention time is 0.

Mean value of the optima: $(opt_tR2_shift35[1]) m

Minima: $(opt_tR2_shift35[3]) m

Maxima: $(opt_tR2_shift35[2]) m

For some solutes the optima is not inside the tested range.

The differences in tR2 is mostly ± 0.1 s.
"""

# ╔═╡ 369e25c7-46d6-4d89-836e-d2efbe61d4a8
md"""
### Optima regarding tR2 - shift = 0.0 s

The the optima is where the difference between measured and simulated 2nd dimension retention time is 0.

For all solutes the optima is not inside the tested range.

Plot of ΔtR2 over L1 looks similar to shift = 3.5s but with significant shifted values to larger differnces.
"""

# ╔═╡ 668d0819-67f6-4717-a8cf-411fcdc94178
Plots.plot!(res_shift0[6], legend=:topleft)

# ╔═╡ c67d6dbf-7160-4b76-9a7a-5bc6a14f437b
opt_tR2_shift0 = sum(res_shift0[8][isnan.(res_shift0[8]).==false])/length(res_shift0[8][isnan.(res_shift0[8]).==false]), maximum(res_shift0[8][isnan.(res_shift0[8]).==false]), minimum(res_shift0[8][isnan.(res_shift0[8]).==false])

# ╔═╡ Cell order:
# ╠═7656e75e-2e47-11ee-31b0-f9a794c245bb
# ╠═dc716bc4-153c-4aef-9d00-55f99ab3f589
# ╠═c7b83c87-2a13-487c-bc45-6915317201d0
# ╠═10353347-1f2c-4f00-9a7f-cf5518eec85d
# ╠═4c15c7db-01dc-4de5-a60b-5fd56c503926
# ╠═67df7e29-48cc-4a9b-b752-c14b6aa0128e
# ╠═41d3693b-bfb4-4732-9a59-d92640f1ab79
# ╠═a474cf70-6b92-4757-bb40-e9e2e3dbfb32
# ╠═9fd0b2b7-c423-49fc-8e89-6d8ec5c0e8ea
# ╟─0180cdbd-73d0-43a0-b583-766e18942edb
# ╠═35d87a01-a90d-4f9e-9a63-9c1393ae4492
# ╠═ef286a7a-246c-4e2c-a257-6ab301ad8590
# ╠═2fc80039-5174-4e93-b997-b5d5cdd8b3fa
# ╟─fa47c6f3-806d-43cd-9983-b2ca27534be3
# ╠═ef1942c3-8076-4830-8c5c-aee353ca1add
# ╠═6658863c-e278-48fd-9f5c-d8a7ac7eacfd
# ╠═d18c7f46-262b-472b-bffe-91a63c72a7d4
# ╠═5507f392-35b2-43cd-b1f8-0f4b8e92eb7f
# ╠═40175e8d-606c-4139-a06e-b753005f61f9
# ╠═ba11ac4b-08c8-44d1-a896-80fd9bf9e921
# ╠═7feb935b-aba6-4c7b-8138-a45777978c55
# ╠═e281594a-0faf-4f50-915c-5bca8abb6969
# ╠═369e25c7-46d6-4d89-836e-d2efbe61d4a8
# ╠═668d0819-67f6-4717-a8cf-411fcdc94178
# ╠═c67d6dbf-7160-4b76-9a7a-5bc6a14f437b
