### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ c423e03a-5d84-11f1-b415-c5b3781ea849
begin
	import Pkg
	# careful: this is _not_ a reproducible environment
	# activate the shared project environment
	Pkg.activate(Base.current_project())
	# instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()
	using DataFrames
	using GasChromatographySystems
	using GasChromatographySimulator
	using Plots
	using UrlDownload
	using PlutoUI
end

# ╔═╡ 28046cdf-eecb-4918-bac0-4b42bbaa8c66
using OrdinaryDiffEq

# ╔═╡ 2036a9e1-0c8e-42fb-9cc8-4910e58b779a
TableOfContents(depth=4)

# ╔═╡ ea3b2dd6-9a9f-4be7-9221-80a7467f3885
md"""
# Test flows and pressures in GCxGC system with pressure modulation
using a ModuleValve and a periodic valve program
"""

# ╔═╡ 3e0a9e10-0614-4723-9373-d4a572fb46f9
md"""
## Setup
"""

# ╔═╡ 0312d777-ee32-4aaa-9143-1bd7a85e7925
begin
	L1 = 3.0
    d1 = 0.1
    df1 = 0.1
    sp1 = "ZB1ms"
    TP1 = GasChromatographySystems.TemperatureProgram([30.0, 1.0, 20.0, 300.0, 1.0])
    L2 = 1.0
    d2 = 0.1
    df2 = 0.1
    sp2 = "Stabilwax"
    TP2 = GasChromatographySystems.TemperatureProgram([30.0, 1.0, 20.0, 300.0, 1.0])
    pin = 400000.0
    pout = 101300.0
    L_valve = 0.01
    d_open = 1.0 # mm
    d_closed = eps(Float64)
    TP_valve = 250.0
    VP = GasChromatographySystems.PeriodicValveProgram(1.0/3, 0.1, 1800.0; inverted=true)
    pmod = 390000.0
end

# ╔═╡ b149034e-830c-4de7-9450-9475bace775a
function GCxGC_DPM(L1, d1, df1, sp1, TP1, L2, d2, df2, sp2, TP2, pin, pout, L_valve, d_open, d_closed, TP_valve, VP, pmod; name="GCxGC_DPM", opt=GasChromatographySystems.Options(), opt_valve=GasChromatographySystems.ModuleValveOptions(ng=true, valve_initial_width=(:fixed, 0.0)), alg=Tsit5(), kwargs...)
	g = SimpleDiGraph(4)
	add_edge!(g, 1, 2) # Inj -> 1st GC column -> Mod point
	add_edge!(g, 2, 3) # Mod point -> 2nd GC column -> Det 
	add_edge!(g, 4, 2) # Pressure Valve line

	# pressure points
	pp = Array{GasChromatographySystems.PressurePoint}(undef, nv(g))
	#pins = pin*1000.0.*ones(length(com_timesteps))
	#nans = NaN.*ones(length(com_timesteps))
	if pout == 0.0
		pouts = eps(Float64)
	else 
		pouts = pout
	end
	pp[1] = GasChromatographySystems.PressurePoint("p₁", pin) # inlet 
	pp[2] = GasChromatographySystems.PressurePoint("p₂", NaN) # 
	pp[3] = GasChromatographySystems.PressurePoint("p₃", pouts) # outlet 1 
	pp[4] = GasChromatographySystems.PressurePoint("p₄", pmod) # pressure modulation
	
	# modules
	modules = Array{GasChromatographySystems.AbstractModule}(undef, ne(g))
	modules[1] = GasChromatographySystems.ModuleColumn("1 -> 2", L1, d1*1e-3, df1*1e-6, sp1, TP1, NaN; alg=alg, ng=true, abstol=1e-10, reltol=1e-8, kwargs...)
	modules[2] = GasChromatographySystems.ModuleColumn("2 -> 3", L2, d2*1e-3, df2*1e-6, sp2, TP2, NaN; alg=alg, ng=true, abstol=1e-10, reltol=1e-8, kwargs...)
	modules[3] = GasChromatographySystems.ModuleValve("4 -> 2", L_valve, d_open*1e-3, d_closed*1e-3, TP_valve, VP, opt_valve, kwargs...)
	# system
	sys_ = GasChromatographySystems.System(name, g, pp, modules, opt)
	sys = GasChromatographySystems.update_system(sys_)

	# add test for the defined pressures and flows
	return sys
end

# ╔═╡ e6dc8d08-8e7a-4cb5-b693-42d0bcbcfe3b
sys = GCxGC_DPM(
    L1, d1, df1, sp1, TP1, 
    L2, d2, df2, sp2, TP2, 
    pin, pout, 
    L_valve, d_open, d_closed, TP_valve, VP, 
    pmod;
	alg = Vern9()
)

# ╔═╡ 27a5bb62-d10f-465c-8a3c-fd4a62e5d3d9
GasChromatographySystems.plot_graph(sys; elabels_fontsize=10, nlabels_fontsize=13)

# ╔═╡ f613a73b-9373-42b0-8365-9495e9f7aa60
flow_bal = GasChromatographySystems.flow_balance(sys)

# ╔═╡ 3c457b28-5c5a-4e09-91f7-2c0e3beaec5e
subst_flow_bal = GasChromatographySystems.substitute_unknown_flows(sys; mode="λ") 

# ╔═╡ 7bc21000-16cb-4e00-9ef7-62e0ec6ce15b
solution = GasChromatographySystems.solve_balance(sys; mode="λ")

# ╔═╡ e4911037-5ae4-4404-9079-91fa3e016a82
p2fun = GasChromatographySystems.build_pressure_squared_functions(sys, solution; mode="λ")

# ╔═╡ bd2d5859-ed71-4821-9227-0c0fcb7714f6
md"""
## Flows and Pressures
"""

# ╔═╡ 4b0a943f-958d-4b90-8314-8415527d4454
flow_func = GasChromatographySystems.flow_functions(sys, p2fun; mode="λ")

# ╔═╡ 2c2a09dc-6f6e-4c00-be77-0d092498587d
flow_func[1](0.0)*60e6, flow_func[1](0.1)*60e6

# ╔═╡ 8273b7a6-02c8-4320-8608-307906d6f96a
flow_func[2](0.0)*60e6, flow_func[2](0.1)*60e6

# ╔═╡ df5b4ac0-45e6-4569-b311-d82fac3a316d
flow_func[3](0.0)*60e6, flow_func[3](0.1)*60e6

# ╔═╡ a91980e3-603a-4b99-8d08-074f8749d7f7
begin
	plot(xlabel="time in s", ylabel="flow in mL/min")
	for i=1:length(flow_func)
		plot!(0:0.01:10.0, flow_func[i].(0.0:0.01:10.0).*60e6, label=sys.modules[i].name, linewidth=4-i)
	end
	plot!(xlims=(0.0, 10.0))
end

# ╔═╡ 2ae83f86-9b41-4d66-9a0f-32356df597ff
p_func = GasChromatographySystems.pressure_functions(sys, p2fun, mode="λ")

# ╔═╡ 167f06ad-e9d5-44d4-95d7-ce755fc47dcb
p_func[2](0.0), p_func[2](0.1)

# ╔═╡ 468cbe01-ad73-470d-8da3-5bd5e2a1d770


# ╔═╡ 50f7f4da-0b8f-4e71-a4d8-91ba3a9342da
begin
	plot(xlabel="time in s", ylabel="pressure in Pa")
	for i=1:length(p_func)
		plot!(0:0.01:10.0, p_func[i].(0.0:0.01:10.0), label=sys.pressurepoints[i].name)
	end
	plot!(xlims=(0.0, 10.0))
end

# ╔═╡ 6e6fc931-de55-4529-a510-3f68adff4a39
lambdas = GasChromatographySystems.flow_permeabilities(sys)

# ╔═╡ ae429231-81eb-4387-9db9-8e745029a407
begin
	plot(xlabel="time in s", ylabel="flow permeability λ")
	for i=1:length(lambdas)
		plot!(0:0.01:10.0, lambdas[i].(0.0:0.01:10.0), label=sys.modules[i].name)
	end
	plot!(xlims=(0.0, 10.0), ylims=(0.0, 2e-14))
end

# ╔═╡ faf84638-a764-4048-9df0-b6c50210d7c5
begin
	kappas = GasChromatographySystems.flow_restrictions(sys)
	plot(xlabel="time in s", ylabel="flow restriction κ")
	for i=1:length(kappas)
		plot!(0:0.01:10.0, kappas[i].(0.0:0.01:10.0), label=sys.modules[i].name)
	end
	plot!(xlims=(0.0, 10.0))
end

# ╔═╡ 574697c2-525a-4eeb-a6b8-fb4fcf61c7b4
md"""
## Graph to Parameters
"""

# ╔═╡ 57a1fabc-335c-4e79-968d-8013cc5fa1ef
begin
	db = DataFrame(urldownload("https://raw.githubusercontent.com/GasChromatographyToolbox/GasChromatographySystems.jl/refs/heads/main/data/Database_GCxGC-TM.csv"))
	insertcols!(db, 1, :No => collect(1:length(db.Name)))
	db
end

# ╔═╡ a7e22b5d-efda-4711-80a8-7c395b14f710
unique(db.Phase)

# ╔═╡ f0540c92-bbfc-446c-b20d-759ad57a05ad
GasChromatographySystems.all_stationary_phases(sys)

# ╔═╡ 66435871-c559-4953-a473-3e8b5af3a0ad
selected_solutes = GasChromatographySystems.common_solutes(db, sys).Name[1:3]

# ╔═╡ 194e1d70-e178-4701-a681-91be450ef0c4
par = GasChromatographySystems.graph_to_parameters(sys, p2fun, db, selected_solutes; interp=true, dt=0.01, mode="λ")

# ╔═╡ cfbc7a41-b4b3-4904-9128-fb1743bfffb0
par_ = GasChromatographySystems.graph_to_parameters(sys, p2fun, db, selected_solutes; interp=false, dt=0.01, mode="λ")

# ╔═╡ 8b5ae010-94a5-4869-8636-094887105af1
par__ = GasChromatographySystems.graph_to_parameters(sys, p2fun, db, selected_solutes; interp=true, dt=0.001, mode="λ")

# ╔═╡ 995e012d-00c9-4b27-93fb-72e850c9fe93
md"""
### Check pressure in 'par'
"""

# ╔═╡ e586f105-e0ca-4c8a-96bc-71335b1a6992
md"""
`dt` in `GasChromatographySystems.graph_to_parameters(sys, p2fun, db, selected_solutes; interp=true, dt=0.01, mode="λ")` must be set to low values (0.01 s) to resolve for the pressure modulation.

But, for `dt < 1` the command `GasChromatographySystems.simulate_along_paths(sys, p2fun, edge_paths, par)` ends in an error: 

"InexactError: Int64(NaN)"
"""

# ╔═╡ 835feb54-79de-4f6f-a780-9f50762c9bed
par[2].prog.Fpin_itp(0.0), par_[2].prog.Fpin_itp(0.0)

# ╔═╡ 0ec29fcc-27d6-4754-a1b1-592ee35f6d6f
par[2].prog.Fpin_itp(0.05), par_[2].prog.Fpin_itp(0.05)

# ╔═╡ 0abd3a51-dc57-4e8b-b6d2-76bd86c5d866
par[2].prog.Fpin_itp(0.1), par_[2].prog.Fpin_itp(0.1)

# ╔═╡ 0ca0c247-20c9-49ce-b045-2bd73368e87d
par[2].prog.Fpin_itp(0.15)

# ╔═╡ 73afd13e-8bd8-4b0d-b830-df6cf2da8bca
par[2].prog.Fpin_itp(0.2)

# ╔═╡ 890ab4e7-3179-4c76-8f0f-330ca96350f5
par[2].prog.Fpin_itp(0.25)

# ╔═╡ a05a7ed4-859c-45e2-94dc-114c41606d91
par[2].prog.Fpin_itp(0.3)

# ╔═╡ 2cc55efe-cdf3-4623-ac2a-70737a7ceadd
par[2].prog.Fpin_itp(0.35)

# ╔═╡ a027a228-6f2e-4e25-ad95-0fcffb1328ca
let
	t_range = 0:0.001:2
	plot(t_range, par[1].prog.pout_itp.(t_range), label="par")
	plot!(t_range, par_[1].prog.pout_itp.(t_range), label="par_")
	plot!(t_range, par__[1].prog.pout_itp.(t_range), label="par__")
end

# ╔═╡ 3e2d2368-f10d-47a7-a7ee-bfd1b0425dbf
let
	t_range = 0:0.001:2
	plot(t_range, par_[1].prog.pout_itp.(t_range))
end

# ╔═╡ 857bf5b1-3fb9-48d4-af44-d39a89e871a1
md"""
## Paths
"""

# ╔═╡ 270df276-ea42-430b-a851-30ac38420aee
vertex_paths, edge_paths = GasChromatographySystems.all_paths(sys.g, sys.modules)

# ╔═╡ ef5ea6c2-50d9-4d9f-981b-739333a874cd
GasChromatographySystems._collect_random_vertex_paths(sys.g, 2)

# ╔═╡ 22330f52-df23-4a7c-b517-31838916fd8f
GasChromatographySystems.index_parameter(sys.g, edge_paths[1])

# ╔═╡ 1d135b53-d7f6-46e6-a8f9-6c415e079eb6
GasChromatographySystems.common_edges(edge_paths[1], edge_paths[1])

# ╔═╡ f48d6432-0eb7-47c1-b61a-16c3197b8666
GasChromatographySystems.positive_flow(sys, p2fun; mode="λ")

# ╔═╡ 4dff2122-a41f-4e5e-a6cc-09b506f4ba49
GasChromatographySystems.path_possible(sys, p2fun, edge_paths[1]; mode="λ")

# ╔═╡ c0f2de8d-9dc7-4103-8702-9ba263dbb4bf
md"""
## Simulation along path
**Check the results**
Functions to model the Spliting at the junction, where the valve is attached, where created by Cursor AI using the spliting at ModuleTM as template.
"""

# ╔═╡ f8740c05-4a9c-4120-aced-c941216a5da3
md"""
### sim
"""

# ╔═╡ 4384ed86-ef55-46ad-95fe-6e62a72c6765
sim = GasChromatographySystems.simulate_along_paths(sys, p2fun, edge_paths, par; nτ=2)

# ╔═╡ 885b27b9-5a1a-48ed-a77a-d26f7e14b594
p_chrom = GasChromatographySystems.chrom(sim[2][1][end]; nτ=6)[1]

# ╔═╡ 087223c2-7711-46cf-8c69-551c674ce9aa
let 
	p = plot(p_chrom, xlims=(18.0, 21.0), legend=:none)
	vline!(p, 18.0:1.0/3:21.0, linestyle=:dash, c=:black, opacity=0.3)
end

# ╔═╡ a5210818-30ac-4438-9f0f-c29334c29d43
let 
	p = plot(p_chrom, xlims=(88.0, 94.0), ylims=(-0.2, 2), legend=:none)
	vline!(p, 88.0:1.0/3:94.0, linestyle=:dash, c=:black, opacity=0.3)
end

# ╔═╡ c7c95758-dbe4-443d-be33-fe6ca2cd3955
let
	p = plot(p_chrom, xlims=(188.0, 195.0), ylims=(-0.1, 1), legend=:none)
	vline!(p, 188.0:1.0/3:195.0, linestyle=:dash, c=:black, opacity=0.3)
end

# ╔═╡ 338d8686-6c05-4ed3-bfa6-c17357437c96
sim[2][1][end]

# ╔═╡ 8b204df3-1668-4ef5-8670-b7b52b7f36f9
pl_GCxGC = GasChromatographySystems.peaklist_GCxGC(sim[2][1][end], 1.0/3; digits=6)

# ╔═╡ 9a65f7c8-73cd-4875-a9dd-c859b22288bc
chrom2d = GasChromatographySystems.chrom2d(sim[2][1][end], sys, 1.0/3)

# ╔═╡ 07b10ac1-264e-4fc4-8aa0-60aca9940347
contour(chrom2d[2][1:end-1], chrom2d[3][1], chrom2d[1]'); xlims!(0.0, 200.0)

# ╔═╡ d60e1340-6bc5-471c-a6f9-09dcf106c49c
begin
	plot(xlabel="1st D time in s", ylabel="2nd D time in s")
	scatter!(pl_GCxGC.tR1, pl_GCxGC.tR2, label="center")
	for i=1:length(pl_GCxGC.Name)
		plot!(pl_GCxGC.tR1s[i], pl_GCxGC.tR2s[i], label=pl_GCxGC.Name[i], markershape=:x)
	end
	plot!()
end

# ╔═╡ a838db9a-ad3d-4091-908f-b0c0a6e6cc3e
sim[3][1][1]

# ╔═╡ 33d5f3b9-f097-42aa-baf4-43ce834fd885
par

# ╔═╡ e8d54030-32cb-4cf0-8b41-4c4b90957b91
GasChromatographySimulator.local_plots("z", "t", sim[3][1][1], par[1]; uncertainty=true)

# ╔═╡ f5fc893d-e452-4e92-8ff1-ee26d9ff550e
GasChromatographySimulator.local_plots("t", "τ", sim[3][1][1], par[1]; uncertainty=true)

# ╔═╡ 1948763e-d6ce-4983-8a66-0bf1866c37a0
GasChromatographySimulator.local_plots("t", "u", sim[3][1][1], par[1]; uncertainty=true); xlims!(100.0, 102.0); ylims!(-0.001, 0.02)

# ╔═╡ 85d0e6c6-7d50-4bcc-9458-1634500ba7ba
sim[3][1][2]

# ╔═╡ 2f9f5782-f746-4356-bcdc-63aeec65a733
par[2]

# ╔═╡ 52fdf357-9e5b-48c6-83aa-f848da1a6dfb
sim[4][2]

# ╔═╡ af9ce04f-cb33-4786-a598-a0409ee4e50a
GasChromatographySimulator.local_plots("z", "t", sim[3][1][2], sim[4][2]; uncertainty=true)

# ╔═╡ 209d13a4-943a-4f52-89e1-b5ec153c3732
GasChromatographySimulator.local_plots("z", "t", sim[3][1][2], sim[4][2]; uncertainty=true); xlims!(0.99, 1.001); ylims!(88.0, 94.0)

# ╔═╡ b9eaaba0-ddf5-4ae5-867a-642a727ff466
sim[4][2]

# ╔═╡ 6264ee71-40df-4035-840c-fef15f08c5df
GasChromatographySimulator.local_plots("t", "u", sim[3][1][2], sim[4][2]; uncertainty=true); xlims!(88.0, 94.0); vline!([88.0:0.333:94.0]); ylims!(0.0, 1.0)

# ╔═╡ 9ad36658-f44c-421c-946c-f33a0a863056
GasChromatographySimulator.local_plots("t", "τ", sim[3][1][2], sim[4][2]; uncertainty=true); xlims!(88.0, 94.0); vline!([88.0:0.333:94.0])#; ylims!(88.0, 94.0)

# ╔═╡ 5086ba92-889f-499a-b514-b65b9575e3e4
sub2_result = filter(:Name => x -> x == "Methyl hexanoate", sim[2][1][2])

# ╔═╡ d33adb32-28a1-472d-adf6-78c3a3c03388
scatter(sub2_result.tR, sub2_result.uR)

# ╔═╡ 1605352d-53c7-4621-af2b-15d4081dc68c
scatter(diff(sub2_result.tR)); hline!([0.333])

# ╔═╡ fe1688f3-31a2-4bb0-a6db-187fc56f3f33
md"""
### sim_
"""

# ╔═╡ 55fbb841-ae2f-4043-befd-4244f133ccb5
sim_ = GasChromatographySystems.simulate_along_paths(sys, p2fun, edge_paths, par_)

# ╔═╡ 3e5c32b3-f292-4208-9b69-b8bff91f825e
p_chrom_ = GasChromatographySystems.chrom(sim_[2][1][end]; nτ=6)[1]

# ╔═╡ 0751d42e-f722-4e80-a795-ddb8629e36a1
let 
	p = plot(p_chrom_, xlims=(88.0, 94.0), ylims=(-0.2, 2), legend=:none)
	vline!(p, 88.0:1.0/3:94.0, linestyle=:dash, c=:black, opacity=0.3)
end

# ╔═╡ c0781512-c9a7-4aa8-8bc9-e8437213e5b6
let
	p = plot(p_chrom_, xlims=(188.0, 195.0), ylims=(-0.1, 1), legend=:none)
	vline!(p, 188.0:1.0/3:195.0, linestyle=:dash, c=:black, opacity=0.3)
end

# ╔═╡ 98eb2ae6-a5ee-4555-9dc9-73e7105195a9
sim_[2][1][end]

# ╔═╡ 7ab36787-c569-4264-93c6-2f313b7b7eec
pl_GCxGC_ = GasChromatographySystems.peaklist_GCxGC(sim_[2][1][end], 1.0/3; digits=6)

# ╔═╡ 12c43908-7b08-48af-a605-2f2cd549fe35
chrom2d_ = GasChromatographySystems.chrom2d(sim_[2][1][end], sys, 1.0/3)

# ╔═╡ bd363798-6e42-45e4-bd5e-180cc73f4ad7
contour(chrom2d_[2][1:end-1], chrom2d_[3][1], chrom2d_[1]'); xlims!(0.0, 200.0)

# ╔═╡ fc9f0821-9dc2-4a41-852c-c7647c5534fd
begin
	plot(xlabel="1st D time in s", ylabel="2nd D time in s")
	scatter!(pl_GCxGC_.tR1, pl_GCxGC_.tR2, label="center")
	for i=1:length(pl_GCxGC_.Name)
		plot!(pl_GCxGC_.tR1s[i], pl_GCxGC_.tR2s[i], label=pl_GCxGC_.Name[i], markershape=:x)
	end
	plot!()
end

# ╔═╡ 849ed98a-f1a8-40ca-a0d2-fcdf69ff42f0
let
	plot(xlabel="position on 1st column in m", ylabel="time in s")
	for i = 70:75
		plot!(sim_[3][1][2][i].t, [u[1] for u in sim_[3][1][2][i].u], label=sim_[2][1][2].Annotations[i])
	end
	#hline!([0.0:(1.0/3):sim__[2][1][1].tR[i]], linestyle=:dash, opacity=0.3, label="modulation")
	plot!()
end

# ╔═╡ cd3f5c30-8558-43bc-b2ca-0e8581381345
let
	plot(xlabel="time in s", ylabel="peak width in s")
	for i = 70:75
		plot!([u[1] for u in sim_[3][1][2][i].u], sqrt.([u[2] for u in sim_[3][1][2][i].u]), label=sim_[2][1][2].Annotations[i])
	end
	#hline!([0.0:(1.0/3):sim__[2][1][1].tR[i]], linestyle=:dash, opacity=0.3, label="modulation")
	plot!()
end

# ╔═╡ 1e40b3ae-ab28-4476-b683-39eb9c455a8f
md"""
### sim__
"""

# ╔═╡ 20e6b4e4-37a5-4554-86a0-56971edd6511
sim__ = GasChromatographySystems.simulate_along_paths(sys, p2fun, edge_paths, par__)

# ╔═╡ 17b76093-19df-477f-bde7-66575277c228
GasChromatographySimulator.plot_chromatogram(sim__[2][1][end], (0.0, 200.0))[1]

# ╔═╡ 9798b7ee-dd8c-4d15-b147-c86318e86879
p_chrom__ = GasChromatographySystems.chrom(sim__[2][1][end]; nτ=6)[1]

# ╔═╡ c0c7cf34-49f0-4f61-8967-cbfb7baeb57b
plot!(p_chrom__, xlims=(88.0, 94.0), ylims=(-0.2, 2), legend=:none)

# ╔═╡ 6929d1de-3fbb-41b6-8129-146123d017af
plot!(p_chrom__, xlims=(188.0, 195.0), ylims=(-0.1, 1), legend=:none)

# ╔═╡ 2628f824-b6b9-4800-895f-7fb1a43d0274
pl_GCxGC__ = GasChromatographySystems.peaklist_GCxGC(sim__[2][1][end], 1.0/3; digits=6)

# ╔═╡ 26e6200d-b513-4110-ae81-d5786dd7745a
chrom2d__ = GasChromatographySystems.chrom2d(sim__[2][1][end], sys, 1.0/3)

# ╔═╡ 1626182a-192c-4fe7-831f-8ddb6f27828e
contour(chrom2d__[2][1:end-1], chrom2d__[3][1], chrom2d__[1]'); xlims!(0.0, 200.0)

# ╔═╡ f82f3b25-853a-499a-83e7-18e481820aac
begin
	plot(xlabel="1st D time in s", ylabel="2nd D time in s")
	scatter!(pl_GCxGC__.tR1, pl_GCxGC__.tR2, label="center")
	for i=1:length(pl_GCxGC__.Name)
		plot!(pl_GCxGC__.tR1s[i], pl_GCxGC__.tR2s[i], label=pl_GCxGC__.Name[i], markershape=:x)
	end
	plot!()
end

# ╔═╡ 1ffc40a9-4cd7-41ed-a440-e828ea4659b3
# take look at the traces for certain selected splitted solutes, including a look at the velocities

# ╔═╡ 27f17ad5-a29e-468c-b919-844d3f4b7d40
let
	i = 1
	plot(sim__[3][1][1][i].t, [u[1] for u in sim__[3][1][1][i].u], label=sim__[2][1][1].Name[i])
	hline!([0.0:(1.0/3):sim__[2][1][1].tR[i]], linestyle=:dash, opacity=0.3, label="modulation")
	plot!(xlabel="position on 1st column in m", ylabel="time in s")
end

# ╔═╡ 1c4e882d-9a7f-471a-a81b-10e1495ecd9e
md"""
* not every modulation seems to be covered by the simulation in the 1st dimension (jump in t(x) due to lower flow in 1st dimension is missing for some modulations)
"""

# ╔═╡ 10e76ad5-1540-4cb2-a64e-dbd7ce44beb7
let
	i = 1
	plot([u[1] for u in sim__[3][1][1][i].u], sqrt.([u[2] for u in sim__[3][1][1][i].u]), label=sim__[2][1][1].Name[i])
	vline!([0.0:(1.0/3):sim__[2][1][1].tR[i]], linestyle=:dash, opacity=0.3, label="modulation")
	plot!(xlabel="time in s", ylabel="peak width in s")
end

# ╔═╡ 52ef20c0-b4dd-4924-ac96-d49e500bd527
let
	i = 2
	plot(sim__[3][1][1][i].t, [u[1] for u in sim__[3][1][1][i].u], label=sim__[2][1][1].Name[i])
	hline!([0.0:(1.0/3):sim__[2][1][1].tR[i]], linestyle=:dash, opacity=0.3, label="modulation")
	plot!(xlabel="position on 1st column in m", ylabel="time in s")
end

# ╔═╡ 934b850a-596d-497d-83e5-9aad6dac2913
let
	i = 2
	plot([u[1] for u in sim__[3][1][1][i].u], sqrt.([u[2] for u in sim__[3][1][1][i].u]), label=sim__[2][1][1].Name[i])
	vline!([0.0:(1.0/3):sim__[2][1][1].tR[i]], linestyle=:dash, opacity=0.3, label="modulation")
	plot!(xlabel="time in s", ylabel="peak width in s")
end

# ╔═╡ 297030a2-c51c-408b-99b9-92603e0f32ce
let
	i = 3
	plot(sim__[3][1][1][i].t, [u[1] for u in sim__[3][1][1][i].u], label=sim__[2][1][1].Name[i])
	hline!([0.0:(1.0/3):sim__[2][1][1].tR[i]], linestyle=:dash, opacity=0.3, label="modulation")
	plot!(xlabel="position on 1st column in m", ylabel="time in s")
end

# ╔═╡ 9ec0bf9c-0666-4d50-8376-b1c1ffcb7a5c
let
	i = 3
	plot([u[1] for u in sim__[3][1][1][i].u], sqrt.([u[2] for u in sim__[3][1][1][i].u]), label=sim__[2][1][1].Name[i])
	vline!([0.0:(1.0/3):sim__[2][1][1].tR[i]], linestyle=:dash, opacity=0.3, label="modulation")
	plot!(xlabel="time in s", ylabel="peak width in s")
end

# ╔═╡ 4fb0c622-8259-4b27-b6cc-ff31a646fa43
sim__[3][1][2]

# ╔═╡ 18174095-79a0-4784-98a1-654f1ca02906
sim__[2][1][2]

# ╔═╡ 7c6b12f5-f16a-406f-8ff5-3599264666e6
let
	i = 6
	j = 2
	plot(sim__[3][1][2][i].t, [u[1] for u in sim__[3][1][2][i].u], label=sim__[2][1][2].Name[i])

	hline!([sim__[2][1][1].tR[j]:(1.0/3):sim__[2][1][2].tR[i]], linestyle=:dash, opacity=0.3, label="modulation")
	plot!(xlabel="position on 1st column in m", ylabel="time in s")
end

# ╔═╡ 5db2d32a-4e46-468d-a86a-3d729b1533ec
let
	i = 6
	j = 2
	plot([u[1] for u in sim__[3][1][2][i].u], sqrt.([u[2] for u in sim__[3][1][2][i].u]), label=sim__[2][1][2].Name[i])
	vline!([sim__[2][1][1].tR[j]:(1.0/3):sim__[2][1][2].tR[i]], linestyle=:dash, opacity=0.3, label="modulation")
	plot!(xlabel="time in s", ylabel="peak width in s")
end

# ╔═╡ 5e13a0e0-a448-4495-a912-738032517e22
let
	i = 7
	j = 2
	plot(sim__[3][1][2][i].t, [u[1] for u in sim__[3][1][2][i].u], label=sim__[2][1][2].Name[i])

	hline!([sim__[2][1][1].tR[j]:(1.0/3):sim__[2][1][2].tR[i]], linestyle=:dash, opacity=0.3, label="modulation")
	plot!(xlabel="position on 1st column in m", ylabel="time in s")
end

# ╔═╡ 644ebf05-c432-4469-8b95-7d343766917d
let
	i = 7
	j = 2
	plot([u[1] for u in sim__[3][1][2][i].u], sqrt.([u[2] for u in sim__[3][1][2][i].u]), label=sim__[2][1][2].Name[i])
	vline!([sim__[2][1][1].tR[j]:(1.0/3):sim__[2][1][2].tR[i]], linestyle=:dash, opacity=0.3, label="modulation")
	plot!(xlabel="time in s", ylabel="peak width in s")
end

# ╔═╡ 472fb53c-ec32-4458-8730-7ec84dc6973f
md"""
## Simulation step-by-step
"""

# ╔═╡ bcb0eaf2-7e9d-4a9e-9207-864e67edcaeb
md"""
### 1st segment - ModuleColumn
"""

# ╔═╡ bb76a596-abb7-4296-a1b9-1f6561e9acf7
par[1]

# ╔═╡ 674154e2-18d8-4682-82c8-23f0d0fc760a
let
	t_range = 0.0:0.01:10.0
	plot(xlabel="time in s", ylabel="pressure in Pa")
	plot!(t_range, par[1].prog.Fpin_itp.(t_range), label="p1")
	plot!(t_range, par[1].prog.pout_itp.(t_range), label="p2")
end

# ╔═╡ 3816c160-e5db-4d4b-87e2-7215fa5d9de3
t₀ = zeros(length(par[1].sub))

# ╔═╡ 74c9f5bc-767a-4f38-96b0-53c5a45d51ab
τ₀ = zeros(length(par[1].sub))

# ╔═╡ cf1dc690-bf5c-4ca0-aeeb-3a637a1b8aa2
sim_1 = GasChromatographySystems.simulate_ModuleColumn(par[1], t₀, τ₀)

# ╔═╡ cc6823fb-e5ec-49ec-b438-49dce9eb30df
GasChromatographySimulator.plot_chromatogram(sim_1[2], (0.0, 200.0))[1]

# ╔═╡ dee6ccd9-85e7-4a1f-8f0c-6e27b3137805
sim_1_ = GasChromatographySystems.simulate_ModuleColumn(par_[1], t₀, τ₀)

# ╔═╡ a654f082-7d18-44c4-9605-7bd04adac250
GasChromatographySimulator.plot_chromatogram(sim_1_[2], (0.0, 200.0))[1]

# ╔═╡ 755a5b3d-b29d-4b96-8043-478499a410b3
sim_1__ = GasChromatographySystems.simulate_ModuleColumn(par__[1], t₀, τ₀)

# ╔═╡ 163c558e-6035-479e-9707-eb3a2f55b4db
GasChromatographySimulator.plot_chromatogram(sim_1__[2], (0.0, 200.0))[1]

# ╔═╡ d7b0b51d-3ef0-463a-96d2-fdd673609b5a
let
	i = 1
	plot(sim_1[3][i].t, [u[1] for u in sim_1[3][i].u], label=sim_1[2].Name[i])
	hline!([0.0:(1.0/3):sim_1[2].tR[i]], linestyle=:dash, opacity=0.3, label="modulation")
	plot!(xlabel="position on 1st column in m", ylabel="time in s")
end

# ╔═╡ 59a65655-69e2-4d67-8c04-31e7eb1731b6
let
	i = 1
	plot([u[1] for u in sim_1[3][i].u], sqrt.([u[2] for u in sim_1[3][i].u]), label=sim_1[2].Name[i])
	vline!([0.0:(1.0/3):sim_1[2].tR[i]], linestyle=:dash, opacity=0.3, label="modulation")
	plot!(xlabel="time in s", ylabel="peak width in s")
end

# ╔═╡ 47cb1b98-034f-40a9-8156-a18e64a79148
let
	i = 1
	plot(sim_1_[3][i].t, [u[1] for u in sim_1_[3][i].u], label=sim_1_[2].Name[i])
	hline!([0.0:(1.0/3):sim_1_[2].tR[i]], linestyle=:dash, opacity=0.3, label="modulation")
	plot!(xlabel="position on 1st column in m", ylabel="time in s")
end

# ╔═╡ 629fcc63-d440-4468-8d9e-cb5853c75f6b
let
	i = 1
	plot([u[1] for u in sim_1_[3][i].u], sqrt.([u[2] for u in sim_1_[3][i].u]), label=sim_1_[2].Name[i])
	vline!([0.0:(1.0/3):sim_1_[2].tR[i]], linestyle=:dash, opacity=0.3, label="modulation")
	plot!(xlabel="time in s", ylabel="peak width in s")
end

# ╔═╡ 9fce356d-b31f-4e88-a159-b0588418d21f
let
	i = 1
	plot(sim_1__[3][i].t, [u[1] for u in sim_1__[3][i].u], label=sim_1__[2].Name[i])
	hline!([0.0:(1.0/3):sim_1__[2].tR[i]], linestyle=:dash, opacity=0.3, label="modulation")
	plot!(xlabel="position on 1st column in m", ylabel="time in s")
end

# ╔═╡ 33666226-4639-4ca3-aa9f-e2df9b7d8c80
let
	i = 1
	plot([u[1] for u in sim_1__[3][i].u], sqrt.([u[2] for u in sim_1__[3][i].u]), label=sim_1__[2].Name[i])
	vline!([0.0:(1.0/3):sim_1__[2].tR[i]], linestyle=:dash, opacity=0.3, label="modulation")
	plot!(xlabel="time in s", ylabel="peak width in s")
end

# ╔═╡ fcfd2b6b-8a54-46bc-ae90-ae05548934c2
sqrt(sim_1__[3][1].u[end][2]), sqrt(sim_1__[3][2].u[end][2]), sqrt(sim_1__[3][3].u[end][2])

# ╔═╡ 599c5e89-be69-45b3-b987-440fbb6705b5
sqrt(sim_1[3][1].u[end][2]), sqrt(sim_1_[3][1].u[end][2]), sqrt(sim_1__[3][1].u[end][2])

# ╔═╡ 3e3661b1-5f72-4638-ab42-bf3f8c94fc6c
sqrt(sim_1[3][2].u[end][2]), sqrt(sim_1_[3][2].u[end][2]), sqrt(sim_1__[3][2].u[end][2])

# ╔═╡ d8b0baf8-e60d-4354-920f-744367ee3ac5
sqrt(sim_1[3][3].u[end][2]), sqrt(sim_1_[3][3].u[end][2]), sqrt(sim_1__[3][3].u[end][2])

# ╔═╡ bc514c91-ab3f-4bff-95ac-8d36140f904c
md"""
### 2nd segment
manually 'injecting' analytes onto the second column

* with modulated inlet pressure
* with constant inlet pressure
"""

# ╔═╡ Cell order:
# ╠═c423e03a-5d84-11f1-b415-c5b3781ea849
# ╠═2036a9e1-0c8e-42fb-9cc8-4910e58b779a
# ╠═ea3b2dd6-9a9f-4be7-9221-80a7467f3885
# ╠═3e0a9e10-0614-4723-9373-d4a572fb46f9
# ╠═0312d777-ee32-4aaa-9143-1bd7a85e7925
# ╠═b149034e-830c-4de7-9450-9475bace775a
# ╠═28046cdf-eecb-4918-bac0-4b42bbaa8c66
# ╠═e6dc8d08-8e7a-4cb5-b693-42d0bcbcfe3b
# ╠═27a5bb62-d10f-465c-8a3c-fd4a62e5d3d9
# ╠═f613a73b-9373-42b0-8365-9495e9f7aa60
# ╠═3c457b28-5c5a-4e09-91f7-2c0e3beaec5e
# ╠═7bc21000-16cb-4e00-9ef7-62e0ec6ce15b
# ╠═e4911037-5ae4-4404-9079-91fa3e016a82
# ╠═bd2d5859-ed71-4821-9227-0c0fcb7714f6
# ╠═4b0a943f-958d-4b90-8314-8415527d4454
# ╠═2c2a09dc-6f6e-4c00-be77-0d092498587d
# ╠═8273b7a6-02c8-4320-8608-307906d6f96a
# ╠═df5b4ac0-45e6-4569-b311-d82fac3a316d
# ╠═a91980e3-603a-4b99-8d08-074f8749d7f7
# ╠═2ae83f86-9b41-4d66-9a0f-32356df597ff
# ╠═167f06ad-e9d5-44d4-95d7-ce755fc47dcb
# ╠═468cbe01-ad73-470d-8da3-5bd5e2a1d770
# ╟─50f7f4da-0b8f-4e71-a4d8-91ba3a9342da
# ╠═6e6fc931-de55-4529-a510-3f68adff4a39
# ╠═ae429231-81eb-4387-9db9-8e745029a407
# ╠═faf84638-a764-4048-9df0-b6c50210d7c5
# ╠═574697c2-525a-4eeb-a6b8-fb4fcf61c7b4
# ╠═57a1fabc-335c-4e79-968d-8013cc5fa1ef
# ╠═a7e22b5d-efda-4711-80a8-7c395b14f710
# ╠═f0540c92-bbfc-446c-b20d-759ad57a05ad
# ╠═66435871-c559-4953-a473-3e8b5af3a0ad
# ╠═194e1d70-e178-4701-a681-91be450ef0c4
# ╠═cfbc7a41-b4b3-4904-9128-fb1743bfffb0
# ╠═8b5ae010-94a5-4869-8636-094887105af1
# ╠═995e012d-00c9-4b27-93fb-72e850c9fe93
# ╠═e586f105-e0ca-4c8a-96bc-71335b1a6992
# ╠═835feb54-79de-4f6f-a780-9f50762c9bed
# ╠═0ec29fcc-27d6-4754-a1b1-592ee35f6d6f
# ╠═0abd3a51-dc57-4e8b-b6d2-76bd86c5d866
# ╠═0ca0c247-20c9-49ce-b045-2bd73368e87d
# ╠═73afd13e-8bd8-4b0d-b830-df6cf2da8bca
# ╠═890ab4e7-3179-4c76-8f0f-330ca96350f5
# ╠═a05a7ed4-859c-45e2-94dc-114c41606d91
# ╠═2cc55efe-cdf3-4623-ac2a-70737a7ceadd
# ╠═a027a228-6f2e-4e25-ad95-0fcffb1328ca
# ╠═3e2d2368-f10d-47a7-a7ee-bfd1b0425dbf
# ╠═857bf5b1-3fb9-48d4-af44-d39a89e871a1
# ╠═270df276-ea42-430b-a851-30ac38420aee
# ╠═ef5ea6c2-50d9-4d9f-981b-739333a874cd
# ╠═22330f52-df23-4a7c-b517-31838916fd8f
# ╠═1d135b53-d7f6-46e6-a8f9-6c415e079eb6
# ╠═f48d6432-0eb7-47c1-b61a-16c3197b8666
# ╠═4dff2122-a41f-4e5e-a6cc-09b506f4ba49
# ╠═c0f2de8d-9dc7-4103-8702-9ba263dbb4bf
# ╠═f8740c05-4a9c-4120-aced-c941216a5da3
# ╠═4384ed86-ef55-46ad-95fe-6e62a72c6765
# ╠═885b27b9-5a1a-48ed-a77a-d26f7e14b594
# ╠═087223c2-7711-46cf-8c69-551c674ce9aa
# ╠═a5210818-30ac-4438-9f0f-c29334c29d43
# ╠═c7c95758-dbe4-443d-be33-fe6ca2cd3955
# ╠═338d8686-6c05-4ed3-bfa6-c17357437c96
# ╠═8b204df3-1668-4ef5-8670-b7b52b7f36f9
# ╠═9a65f7c8-73cd-4875-a9dd-c859b22288bc
# ╠═07b10ac1-264e-4fc4-8aa0-60aca9940347
# ╠═d60e1340-6bc5-471c-a6f9-09dcf106c49c
# ╠═a838db9a-ad3d-4091-908f-b0c0a6e6cc3e
# ╠═33d5f3b9-f097-42aa-baf4-43ce834fd885
# ╠═e8d54030-32cb-4cf0-8b41-4c4b90957b91
# ╠═f5fc893d-e452-4e92-8ff1-ee26d9ff550e
# ╠═1948763e-d6ce-4983-8a66-0bf1866c37a0
# ╠═85d0e6c6-7d50-4bcc-9458-1634500ba7ba
# ╠═2f9f5782-f746-4356-bcdc-63aeec65a733
# ╠═52fdf357-9e5b-48c6-83aa-f848da1a6dfb
# ╠═af9ce04f-cb33-4786-a598-a0409ee4e50a
# ╠═209d13a4-943a-4f52-89e1-b5ec153c3732
# ╠═b9eaaba0-ddf5-4ae5-867a-642a727ff466
# ╠═6264ee71-40df-4035-840c-fef15f08c5df
# ╠═9ad36658-f44c-421c-946c-f33a0a863056
# ╠═5086ba92-889f-499a-b514-b65b9575e3e4
# ╠═d33adb32-28a1-472d-adf6-78c3a3c03388
# ╠═1605352d-53c7-4621-af2b-15d4081dc68c
# ╠═fe1688f3-31a2-4bb0-a6db-187fc56f3f33
# ╠═55fbb841-ae2f-4043-befd-4244f133ccb5
# ╠═3e5c32b3-f292-4208-9b69-b8bff91f825e
# ╠═0751d42e-f722-4e80-a795-ddb8629e36a1
# ╠═c0781512-c9a7-4aa8-8bc9-e8437213e5b6
# ╠═98eb2ae6-a5ee-4555-9dc9-73e7105195a9
# ╠═7ab36787-c569-4264-93c6-2f313b7b7eec
# ╠═12c43908-7b08-48af-a605-2f2cd549fe35
# ╠═bd363798-6e42-45e4-bd5e-180cc73f4ad7
# ╠═fc9f0821-9dc2-4a41-852c-c7647c5534fd
# ╠═849ed98a-f1a8-40ca-a0d2-fcdf69ff42f0
# ╠═cd3f5c30-8558-43bc-b2ca-0e8581381345
# ╠═1e40b3ae-ab28-4476-b683-39eb9c455a8f
# ╠═20e6b4e4-37a5-4554-86a0-56971edd6511
# ╠═17b76093-19df-477f-bde7-66575277c228
# ╠═9798b7ee-dd8c-4d15-b147-c86318e86879
# ╠═c0c7cf34-49f0-4f61-8967-cbfb7baeb57b
# ╠═6929d1de-3fbb-41b6-8129-146123d017af
# ╠═2628f824-b6b9-4800-895f-7fb1a43d0274
# ╠═26e6200d-b513-4110-ae81-d5786dd7745a
# ╠═1626182a-192c-4fe7-831f-8ddb6f27828e
# ╠═f82f3b25-853a-499a-83e7-18e481820aac
# ╠═1ffc40a9-4cd7-41ed-a440-e828ea4659b3
# ╠═27f17ad5-a29e-468c-b919-844d3f4b7d40
# ╠═1c4e882d-9a7f-471a-a81b-10e1495ecd9e
# ╠═10e76ad5-1540-4cb2-a64e-dbd7ce44beb7
# ╠═52ef20c0-b4dd-4924-ac96-d49e500bd527
# ╠═934b850a-596d-497d-83e5-9aad6dac2913
# ╠═297030a2-c51c-408b-99b9-92603e0f32ce
# ╠═9ec0bf9c-0666-4d50-8376-b1c1ffcb7a5c
# ╠═4fb0c622-8259-4b27-b6cc-ff31a646fa43
# ╠═18174095-79a0-4784-98a1-654f1ca02906
# ╠═7c6b12f5-f16a-406f-8ff5-3599264666e6
# ╠═5db2d32a-4e46-468d-a86a-3d729b1533ec
# ╠═5e13a0e0-a448-4495-a912-738032517e22
# ╠═644ebf05-c432-4469-8b95-7d343766917d
# ╠═472fb53c-ec32-4458-8730-7ec84dc6973f
# ╠═bcb0eaf2-7e9d-4a9e-9207-864e67edcaeb
# ╠═bb76a596-abb7-4296-a1b9-1f6561e9acf7
# ╠═674154e2-18d8-4682-82c8-23f0d0fc760a
# ╠═3816c160-e5db-4d4b-87e2-7215fa5d9de3
# ╠═74c9f5bc-767a-4f38-96b0-53c5a45d51ab
# ╠═cf1dc690-bf5c-4ca0-aeeb-3a637a1b8aa2
# ╠═cc6823fb-e5ec-49ec-b438-49dce9eb30df
# ╠═dee6ccd9-85e7-4a1f-8f0c-6e27b3137805
# ╠═a654f082-7d18-44c4-9605-7bd04adac250
# ╠═755a5b3d-b29d-4b96-8043-478499a410b3
# ╠═163c558e-6035-479e-9707-eb3a2f55b4db
# ╠═d7b0b51d-3ef0-463a-96d2-fdd673609b5a
# ╠═59a65655-69e2-4d67-8c04-31e7eb1731b6
# ╠═47cb1b98-034f-40a9-8156-a18e64a79148
# ╠═629fcc63-d440-4468-8d9e-cb5853c75f6b
# ╠═9fce356d-b31f-4e88-a159-b0588418d21f
# ╠═33666226-4639-4ca3-aa9f-e2df9b7d8c80
# ╠═fcfd2b6b-8a54-46bc-ae90-ae05548934c2
# ╠═599c5e89-be69-45b3-b987-440fbb6705b5
# ╠═3e3661b1-5f72-4638-ab42-bf3f8c94fc6c
# ╠═d8b0baf8-e60d-4354-920f-744367ee3ac5
# ╠═bc514c91-ab3f-4bff-95ac-8d36140f904c
