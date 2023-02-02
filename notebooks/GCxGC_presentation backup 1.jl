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

# ╔═╡ 13e893d1-da6f-4081-a705-b71cfb4cb547
using Interpolations

# ╔═╡ 139bfe64-17ac-455d-84a2-09c9564156cd
md"""
# Simple GCxGC system
"""

# ╔═╡ f6c6443d-0e39-4e82-b2c4-49f2954e5b49
md"""
## Definition of system
"""

# ╔═╡ 719dc69d-2692-4d41-868e-14b664019d57
md"""
### Temperature programs
"""

# ╔═╡ c9941593-0311-473a-9c0c-f764d41ac7e7
GCxGC_TP = GasChromatographySystems.TemperatureProgram([0.0, 120.0, 3120.0, 300.0], [40.0, 40.0, 300.0, 300.0])

# ╔═╡ 03190a64-09cf-4a9b-9b86-e2b50392481c
md"""
### GCxGC
"""

# ╔═╡ d571712f-981e-4d1e-95ce-b4dbf7de50c6
begin
	g = SimpleDiGraph(3)
	add_edge!(g, 1, 2) # 1st-D GC
	add_edge!(g, 2, 3) # 2nd-D GC
	#add_edge!(g2, 3, 4) # 2nd TL column inlet -> TL column -> Det
	# pressure points
	pp = Array{GasChromatographySystems.PressurePoint}(undef, nv(g))
	pp[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 120.0, 3120.0, 300.0], [159550.0, 159550.0, 265600.0, 265600.0]) # inlet 
	pp[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 120.0, 3120.0, 300.0], [NaN, NaN, NaN, NaN]) 
	pp[3] = GasChromatographySystems.PressurePoint("p₃", [0.0, 120.0, 3120.0, 300.0], [eps(Float64), eps(Float64), eps(Float64), eps(Float64)]) 
	#pp2[4] = GasChromatographySystems.PressurePoint("p₄", [0.0, 1800.0], [101300.0, 101300.0]) # outlet
	# modules
	modules = Array{GasChromatographySystems.AbstractModule}(undef, ne(g))
	modules[1] = GasChromatographySystems.ModuleColumn("GC column 1", 30.0, 0.25e-3, 0.25e-6, "Rxi5ms", GCxGC_TP)
	modules[2] = GasChromatographySystems.ModuleColumn("GC column 2", 2.0, 0.25e-3, 0.25e-6, "Rxi17SilMS", GCxGC_TP)
	#modules2[3] = GasChromatographySystems.ModuleColumn("GC column 3", 10.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP)
	# system
	sys = GasChromatographySystems.update_system(GasChromatographySystems.System(g, pp, modules, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ 4b379226-1646-4b13-b64a-7b55b66442f5
md"""
## Pressure and flow calculation
"""

# ╔═╡ 7783216f-6a73-44aa-98f7-e0ca628ba145
function plot_flow_over_time(sys)
	plotly()
	F_func = GasChromatographySystems.flow_functions(sys)
	com_timesteps = GasChromatographySystems.common_timesteps(sys)
	trange = 0:sum(com_timesteps)/100:sum(com_timesteps)
	p_flow = Plots.plot(xlabel="time in s", ylabel="flow in mL/min")
	for i=1:ne(sys.g)
		Plots.plot!(p_flow, trange, F_func[i].(trange)*60*1e6, label="F_$(sys.modules[i].name)")
	end
	return p_flow
end

# ╔═╡ 06b9ba3f-4dff-4fba-afce-923b95fae9c3
GasChromatographySystems.plot_graph_with_flow(sys, 0)

# ╔═╡ af30bbcd-d6e4-4bd7-97f0-3243c514194a
plot_flow_over_time(sys)

# ╔═╡ b51877ef-a8a6-48d4-a521-5b1e01d2be40
function plot_pressure_over_time(sys)
	plotly()
	p_func = GasChromatographySystems.pressure_functions(sys)
	com_timesteps = GasChromatographySystems.common_timesteps(sys)
	trange = 0:sum(com_timesteps)/100:sum(com_timesteps)
	p_pres = Plots.plot(xlabel="time in s", ylabel="pressure in Pa(a)", legend=:topleft)
	for i=1:nv(sys.g)
		Plots.plot!(p_pres, trange, p_func[i].(trange), label="$(sys.pressurepoints[i].name)")
	end
	return p_pres
end

# ╔═╡ 9c9aed82-3015-419b-ae6a-aa6f97a9d753
plot_pressure_over_time(sys)

# ╔═╡ 3f50f014-6ab8-4f0b-8302-677e5fb6c6c3
function pressures_squared(sys)
	#p² = Array{Interpolations.Extrapolation}(undef, nv(sys.g))
	p² = Array{Any}(undef, nv(sys.g))
	for i=1:nv(sys.g)
		f = if isnan(sys.pressurepoints[i].pressure_steps[1])
			GasChromatographySimulator.steps_interpolation(sys.pressurepoints[i].timesteps, NaN.*ones(length(sys.pressurepoints[i].timesteps)))
		else
			GasChromatographySimulator.steps_interpolation(sys.pressurepoints[i].timesteps, identity.(sys.pressurepoints[i].pressure_steps).^2)
		end
		p²[i] = f
	end
	return p²
end

# ╔═╡ 24e4d1f4-62d6-49b5-bf0f-390ba9387f86
p2 = pressures_squared(sys)

# ╔═╡ d5e44df5-6874-45e1-a3ae-84fe1c7a1dd9
p2[1](0)

# ╔═╡ ab0f6b8a-9824-4d7a-ba60-5dd7d7b14f0a
pin = GasChromatographySimulator.steps_interpolation([0.0, 120.0, 3120.0, 300.0], [159550.0, 159550.0, 265600.0, 265600.0])

# ╔═╡ 9cf51a5f-8bcb-4288-9459-7199ab9af8fa
typeof(pin)

# ╔═╡ 7bb127c8-07ad-450a-afdb-aa1f5be8ecc7
typeof(pin) <: Interpolations.Extrapolation

# ╔═╡ d80e286c-bc33-45cb-b0e5-d572c5f3cf6e
pin(0)^2

# ╔═╡ c9883ccd-9d39-404e-b094-23127c2ba4af
com_timesteps = GasChromatographySystems.common_timesteps(sys)

# ╔═╡ 6c54d175-e861-45a1-b8b5-515bdcad3418
function pressures_squared_(sys)
	#p² = Array{Interpolations.Extrapolation}(undef, nv(sys.g))
	p² = Array{Any}(undef, nv(sys.g))
	for i=1:nv(sys.g)
		f = if isnan(sys.pressurepoints[i].pressure_steps[1])
			GasChromatographySimulator.steps_interpolation(sys.pressurepoints[i].timesteps, NaN.*ones(length(sys.pressurepoints[i].timesteps)))
		else
			GasChromatographySimulator.steps_interpolation(sys.pressurepoints[i].timesteps, identity.(sys.pressurepoints[i].pressure_steps))
		end
		g(t) = f(t).^2 
		p²[i] = g
	end
	return p²
end

# ╔═╡ c2d1cfce-1636-42ac-bb42-0ec102a9bb1a
p2_ = pressures_squared_(sys)

# ╔═╡ 17f6dfdc-ba95-4351-b9dd-9ffd284b4a74
p2_[1](0)

# ╔═╡ a65c868c-fa26-4a7f-be94-2d86ee9518fc
md"""
## Selection of solutes
"""

# ╔═╡ b424f5f1-aca8-40e5-aa0d-32fc7c4d4ee6
db = DataFrame(CSV.File(joinpath(dirname(pwd()), "data", "database_Kcentric.csv"), header=1, silencewarnings=true))

# ╔═╡ 0b0c0577-ec2e-450f-9f19-4a29da213efe
unique(GasChromatographySystems.all_stationary_phases(sys))

# ╔═╡ 2eff2573-137d-4a85-9bc9-e43d081b8d77
GasChromatographySystems.common_solutes(db, sys)

# ╔═╡ 5e824c58-4d19-47fa-a1f9-27b977cb0761
selected_solutes = GasChromatographySystems.common_solutes(db, sys).Name

# ╔═╡ e935ae05-c654-4365-a569-fed418f35fe1
md"""
## Graph to Parameters
"""

# ╔═╡ d060c9a3-d760-472f-b52b-3c481b39a85c
par = GasChromatographySystems.graph_to_parameters(sys, db, selected_solutes)

# ╔═╡ 5bc5797c-d846-4b63-8b99-9e0b638c3ea8
md"""
## Simulation by hand
"""

# ╔═╡ ac23297f-1640-4325-8494-62539bb4034e
par_1 = GasChromatographySystems.change_initial(par[1], zeros(length(selected_solutes)), zeros(length(selected_solutes)))

# ╔═╡ e232556d-37dd-4d2e-afc2-65b8da8c78e4
pl_1, sol_1 = GasChromatographySimulator.simulate(par_1)

# ╔═╡ 85adbc6b-a7fc-4abf-8b0b-f82f63456396
par_1.sub

# ╔═╡ deeab4df-4a41-4ab5-a388-ae224503a8db
pl_1

# ╔═╡ 8ef596fe-e22f-41e4-a711-0011fa16634d
par_2 = GasChromatographySystems.change_initial(par[2], pl_1)

# ╔═╡ 0bd81908-36f7-42b6-9538-ceabc5fbdb5b
pl_2, sol_2 = GasChromatographySimulator.simulate(par_2)

# ╔═╡ 423417fd-1cee-4803-9efe-b07815cc8678
begin
	pl_GCxGC = DataFrame(Name = pl_1.Name, tR1 = pl_1.tR, τR1 = pl_1.τR)
	tR2 = Array{Float64}(undef, length(pl_GCxGC.Name))
	τR2 = Array{Float64}(undef, length(pl_GCxGC.Name))
	CAS1 = GasChromatographySimulator.CAS_identification(pl_1.Name).CAS
	CAS2 = GasChromatographySimulator.CAS_identification(pl_2.Name).CAS
	for i=1:length(pl_GCxGC.Name)
		ii = findfirst(CAS1[i].==CAS2)
		tR2[i] = pl_2.tR[ii]
		τR2[i] = pl_2.τR[ii]
	end
	pl_GCxGC[!, :tR2] = tR2
	pl_GCxGC[!, :τR2] = τR2
	pl_GCxGC
end

# ╔═╡ 83008b1d-caa8-48d6-a164-7a326fe86007
Plots.scatter(pl_GCxGC.tR1, pl_GCxGC.tR2.-pl_GCxGC.tR1, xlabel="tR1 in s", ylabel="tR2 in s")

# ╔═╡ bab0c37b-36b8-4d65-8824-53ccbc3bd237
GasChromatographySystems.common_timesteps(sys)

# ╔═╡ d7162c08-64af-408c-974a-16574076f8bc
# create new parameters with interpolated pressure functions

# ╔═╡ b1c437b8-93c2-479b-ac36-407e25f2de55
itp_pf = GasChromatographySystems.interpolate_pressure_functions(sys)

# ╔═╡ 8a185a43-a3b3-4fdc-9839-6a2b4e5f8a5d
GasChromatographySystems.interpolate_pressure_functions(sys, dt=10)

# ╔═╡ e4304fb8-15d0-4723-a2f5-fa78a01a2aed
itp_prog = GasChromatographySimulator.Program(par_1.prog.time_steps, par_1.prog.temp_steps, par_1.prog.Fpin_steps, itp_pf[2].(cumsum(GasChromatographySystems.common_timesteps(sys))), par_1.prog.gf, par_1.prog.a_gf, par_1.prog.T_itp, itp_pf[1], itp_pf[2])

# ╔═╡ 2c097584-8b8d-42f4-8639-52931fb2876b
itp_par = GasChromatographySimulator.Parameters(par_1.col, itp_prog, par_1.sub, par_1.opt)

# ╔═╡ 193f16fe-bfa6-45c5-b0a2-e19c3baee385
itp_pl_1, itp_sol_1 = GasChromatographySimulator.simulate(itp_par)

# ╔═╡ d71f36d5-0840-4cb4-ae1f-3d8672079b31
pres, unk = GasChromatographySystems.solve_pressure(sys)

# ╔═╡ d494a040-3ab6-4450-8107-4bfd04f9d15f
@time pres[1](567.8)

# ╔═╡ 6f4691c6-4da1-4751-a5ef-9f2cd076b2db
@time itp_pf[2](567.8)

# ╔═╡ 9c269628-5742-424e-b7b3-dc36424ea03e
md"""
The solution function for the unknown pressure using Symbolics.jl and a substitution dictionary, whixh is a function of time, takes to much time to evaluate. Faster is the method using interpolation functions as used before. But for systems with split/merge-points, the pressure at these locations is not linear in time, even if all other pressure functions and temperature programs are. Here a small additional error is introduced (a similar error is allready in GasChromatographySimulator.jl for thermal gradients, the non-linear gradient is linear interpolated for stepsizes of 1e-3m).
"""

# ╔═╡ 79bc683a-8d51-4d81-a27c-e42663720917
md"""
## Simulation
"""

# ╔═╡ 457df2b2-8bfa-4e04-b4f2-ae1ca8902a4e
paths = GasChromatographySystems.all_paths(sys.g, 1)[2]

# ╔═╡ f106b5fd-e17a-4876-8562-e82670a1ab25
sim_res = GasChromatographySystems.simulate_along_paths(sys, paths, db, selected_solutes, par)

# ╔═╡ cb6840f8-5853-41ec-8de3-e6b793ccdf33
all(degree(sys.g).<3)

# ╔═╡ be6434d2-7151-435e-9567-a385015f1dc5
GasChromatographySystems.common_timesteps(sys)

# ╔═╡ 69a5ef12-f1c2-4a35-b184-33f8b8af5d00
diff(GasChromatographySystems.common_timesteps(sys))

# ╔═╡ c66e2dfd-5ae1-403a-bd04-164aa5cf6f46
lcm(120, 3000)

# ╔═╡ 34623543-b78f-494b-b83d-94b60ec34ad1
mod(120, 2820)

# ╔═╡ f56c2e49-2e71-4a0d-bdf1-e822b96c8891
mod(120, 320)

# ╔═╡ 97e010eb-1829-4e51-b347-25c91d3c9d16
3120/60

# ╔═╡ 49486c29-8cd1-4721-ab20-577778e84c89
300/60

# ╔═╡ 5373eaae-4372-4561-83cf-b62664bf28af
sim_res[2][1][1]

# ╔═╡ e57ab914-8db4-41ca-9950-1134e35e8c40
function peaklist_GCxGC(pl_1, pl_2)
	pl_GCxGC = DataFrame(Name = pl_1.Name, tR1 = pl_1.tR, τR1 = pl_1.τR)
	tR2 = Array{Float64}(undef, length(pl_GCxGC.Name))
	τR2 = Array{Float64}(undef, length(pl_GCxGC.Name))
	CAS1 = GasChromatographySimulator.CAS_identification(pl_1.Name).CAS
	CAS2 = GasChromatographySimulator.CAS_identification(pl_2.Name).CAS
	for i=1:length(pl_GCxGC.Name)
		ii = findfirst(CAS1[i].==CAS2)
		tR2[i] = pl_2.tR[ii]
		τR2[i] = pl_2.τR[ii]
	end
	pl_GCxGC[!, :tR2] = tR2
	pl_GCxGC[!, :τR2] = τR2
	return pl_GCxGC
end

# ╔═╡ 0bc7a04d-7534-4e0b-993c-a6ca65f8d363
peaklist_GCxGC(sim_res[2][1][1], sim_res[2][1][2])

# ╔═╡ Cell order:
# ╠═e9d64f0e-95b7-11ed-2fee-ab579cd02ebc
# ╠═139bfe64-17ac-455d-84a2-09c9564156cd
# ╠═f6c6443d-0e39-4e82-b2c4-49f2954e5b49
# ╠═719dc69d-2692-4d41-868e-14b664019d57
# ╠═c9941593-0311-473a-9c0c-f764d41ac7e7
# ╠═03190a64-09cf-4a9b-9b86-e2b50392481c
# ╠═d571712f-981e-4d1e-95ce-b4dbf7de50c6
# ╠═4b379226-1646-4b13-b64a-7b55b66442f5
# ╠═7783216f-6a73-44aa-98f7-e0ca628ba145
# ╠═06b9ba3f-4dff-4fba-afce-923b95fae9c3
# ╠═af30bbcd-d6e4-4bd7-97f0-3243c514194a
# ╠═9c9aed82-3015-419b-ae6a-aa6f97a9d753
# ╠═b51877ef-a8a6-48d4-a521-5b1e01d2be40
# ╠═3f50f014-6ab8-4f0b-8302-677e5fb6c6c3
# ╠═24e4d1f4-62d6-49b5-bf0f-390ba9387f86
# ╠═d5e44df5-6874-45e1-a3ae-84fe1c7a1dd9
# ╠═ab0f6b8a-9824-4d7a-ba60-5dd7d7b14f0a
# ╠═9cf51a5f-8bcb-4288-9459-7199ab9af8fa
# ╠═13e893d1-da6f-4081-a705-b71cfb4cb547
# ╠═7bb127c8-07ad-450a-afdb-aa1f5be8ecc7
# ╠═d80e286c-bc33-45cb-b0e5-d572c5f3cf6e
# ╠═c9883ccd-9d39-404e-b094-23127c2ba4af
# ╠═c4ef5c8d-f0bf-40e9-b94a-c0a3064a65e4
# ╠═6b883e45-8fd5-4a71-8a02-c8c3f066da6c
# ╠═6c54d175-e861-45a1-b8b5-515bdcad3418
# ╠═c2d1cfce-1636-42ac-bb42-0ec102a9bb1a
# ╠═17f6dfdc-ba95-4351-b9dd-9ffd284b4a74
# ╠═6c0ed1ab-c521-4918-b11f-c0da142747a1
# ╠═a65c868c-fa26-4a7f-be94-2d86ee9518fc
# ╠═b424f5f1-aca8-40e5-aa0d-32fc7c4d4ee6
# ╠═0b0c0577-ec2e-450f-9f19-4a29da213efe
# ╠═2eff2573-137d-4a85-9bc9-e43d081b8d77
# ╠═5e824c58-4d19-47fa-a1f9-27b977cb0761
# ╠═e935ae05-c654-4365-a569-fed418f35fe1
# ╠═d060c9a3-d760-472f-b52b-3c481b39a85c
# ╠═5bc5797c-d846-4b63-8b99-9e0b638c3ea8
# ╠═ac23297f-1640-4325-8494-62539bb4034e
# ╠═e232556d-37dd-4d2e-afc2-65b8da8c78e4
# ╠═85adbc6b-a7fc-4abf-8b0b-f82f63456396
# ╠═deeab4df-4a41-4ab5-a388-ae224503a8db
# ╠═8ef596fe-e22f-41e4-a711-0011fa16634d
# ╠═0bd81908-36f7-42b6-9538-ceabc5fbdb5b
# ╠═423417fd-1cee-4803-9efe-b07815cc8678
# ╠═83008b1d-caa8-48d6-a164-7a326fe86007
# ╠═bab0c37b-36b8-4d65-8824-53ccbc3bd237
# ╠═d7162c08-64af-408c-974a-16574076f8bc
# ╠═b1c437b8-93c2-479b-ac36-407e25f2de55
# ╠═8a185a43-a3b3-4fdc-9839-6a2b4e5f8a5d
# ╠═e4304fb8-15d0-4723-a2f5-fa78a01a2aed
# ╠═2c097584-8b8d-42f4-8639-52931fb2876b
# ╠═193f16fe-bfa6-45c5-b0a2-e19c3baee385
# ╠═d71f36d5-0840-4cb4-ae1f-3d8672079b31
# ╠═d494a040-3ab6-4450-8107-4bfd04f9d15f
# ╠═6f4691c6-4da1-4751-a5ef-9f2cd076b2db
# ╠═9c269628-5742-424e-b7b3-dc36424ea03e
# ╠═79bc683a-8d51-4d81-a27c-e42663720917
# ╠═457df2b2-8bfa-4e04-b4f2-ae1ca8902a4e
# ╠═f106b5fd-e17a-4876-8562-e82670a1ab25
# ╠═cb6840f8-5853-41ec-8de3-e6b793ccdf33
# ╠═be6434d2-7151-435e-9567-a385015f1dc5
# ╠═69a5ef12-f1c2-4a35-b184-33f8b8af5d00
# ╠═c66e2dfd-5ae1-403a-bd04-164aa5cf6f46
# ╠═34623543-b78f-494b-b83d-94b60ec34ad1
# ╠═f56c2e49-2e71-4a0d-bdf1-e822b96c8891
# ╠═97e010eb-1829-4e51-b347-25c91d3c9d16
# ╠═49486c29-8cd1-4721-ab20-577778e84c89
# ╠═5373eaae-4372-4561-83cf-b62664bf28af
# ╠═0bc7a04d-7534-4e0b-993c-a6ca65f8d363
# ╠═e57ab914-8db4-41ca-9950-1134e35e8c40
