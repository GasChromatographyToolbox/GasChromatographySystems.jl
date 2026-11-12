### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 1a2b3c4d-0001-4000-8000-000000000001
begin
	import Pkg
	# Prefer the GC-Workspace shared env (local GCSim + GCSys) when available
	ws = joinpath(@__DIR__, "..", "..", "GC-Workspace")
	if isfile(joinpath(ws, "Project.toml"))
		Pkg.activate(ws)
	else
		Pkg.activate(joinpath(@__DIR__, ".."))
	end
	using GasChromatographySystems
	using GasChromatographySimulator
	using GasChromatographySystems.Graphs
	using Plots
	using CSV, DataFrames
end

# ╔═╡ 1a2b3c4d-0002-4000-8000-000000000002
md"""
# Temporal `ModuleColumn` (reverse-capable)

Demo for opt-in ``x(t)`` migration on a Systems column edge (`domain=:temporal`):

1. **Forward regression** — finite-outlet series column, compare `:spatial` vs `:temporal`
2. **Trap-like reverse** — cold load (`pin > pout`), then hot reverse (`pin < pout`) with same-capillary carry of ``(x_0,\sigma_0)``

Coordinates: ``x=0`` = load inlet; reverse leaves toward ``x=0``.
"""

# ╔═╡ 1a2b3c4d-0003-4000-8000-000000000003
begin
	db_file = joinpath(@__DIR__, "..", "test", "data", "Database_test.csv")
	db = DataFrame(CSV.File(db_file, header=1, silencewarnings=true))
	insertcols!(db, 1, :No => collect(1:nrow(db)))
	solute = ["Undecane"]
end

# ╔═╡ 1a2b3c4d-0004-4000-8000-000000000004
md"""
## 1. Forward elution — `:spatial` vs `:temporal`
"""

# ╔═╡ 1a2b3c4d-0005-4000-8000-000000000005
begin
	TP = GasChromatographySystems.TemperatureProgram([0.0, 600.0], [40.0, 200.0])
	sys_pos = GasChromatographySystems.SeriesSystem(
		[10.0], [0.25], [0.25], ["SLB5ms"], [TP],
		NaN, 150000.0, 101300.0;
		name="forward_spatial", opt=GasChromatographySystems.Options(control="Pressure"),
	)
	sys_td = GasChromatographySystems.SeriesSystem(
		[10.0], [0.25], [0.25], ["SLB5ms"], [TP],
		NaN, 150000.0, 101300.0;
		name="forward_temporal", opt=GasChromatographySystems.Options(control="Pressure"),
		domain=:temporal,
	)
	paths = GasChromatographySystems.all_paths(sys_pos.g, 1)[2]

	p2_pos = GasChromatographySystems.build_pressure_squared_functions(sys_pos, GasChromatographySystems.solve_balance(sys_pos))
	par_pos = GasChromatographySystems.graph_to_parameters(sys_pos, p2_pos, db, solute)
	_, pls_pos, _, _ = GasChromatographySystems.simulate_along_paths(sys_pos, p2_pos, paths, par_pos)

	p2_td = GasChromatographySystems.build_pressure_squared_functions(sys_td, GasChromatographySystems.solve_balance(sys_td))
	par_td = GasChromatographySystems.graph_to_parameters(sys_td, p2_td, db, solute)
	_, pls_td, _, _ = GasChromatographySystems.simulate_along_paths(sys_td, p2_td, paths, par_td)

	pl_pos = pls_pos[1][1]
	pl_td = pls_td[1][1]
	(; tR_spatial=pl_pos.tR[1], tR_temporal=pl_td.tR[1], exit_end=pl_td.exit_end[1], τR_temporal=pl_td.τR[1])
end

# ╔═╡ 1a2b3c4d-0006-4000-8000-000000000006
md"""
## 2. Cold load → hot reverse (same capillary)

Short column as a stand-in trap. Load parks the band; reverse exits at ``x=0``.
"""

# ╔═╡ 1a2b3c4d-0007-4000-8000-000000000007
begin
	Ltrap = 0.5
	opt_td = GasChromatographySystems.ModuleColumnOptions(domain=:temporal)
	g = SimpleDiGraph(2)
	add_edge!(g, 1, 2)

	# --- load ---
	TP_load = GasChromatographySystems.TemperatureProgram([0.0, 30.0], [40.0, 40.0])
	pp_load = [
		GasChromatographySystems.PressurePoint("pin", 150000.0),
		GasChromatographySystems.PressurePoint("pout", 101300.0),
	]
	mod_load = GasChromatographySystems.ModuleColumn("trap", Ltrap, 0.25e-3, 0.25e-6, "SLB5ms", TP_load, opt_td)
	sys_load = GasChromatographySystems.update_system(
		GasChromatographySystems.System("trap_load", g, pp_load, [mod_load], GasChromatographySystems.Options(control="Pressure")),
	)
	p2_load = GasChromatographySystems.build_pressure_squared_functions(sys_load, GasChromatographySystems.solve_balance(sys_load))
	par_load = GasChromatographySystems.graph_to_parameters(sys_load, p2_load, db, solute)[1]
	_, pl_load, sol_load = GasChromatographySystems.simulate_ModuleColumn(par_load, [0.0], [0.0])

	x_park = sol_load[1].u[end][1]
	σ_park = sqrt(max(sol_load[1].u[end][2], 0.0))
	t_load = sol_load[1].t[end]

	# --- reverse ---
	TP_rel = GasChromatographySystems.TemperatureProgram([0.0, 120.0], [250.0, 250.0])
	pp_rel = [
		GasChromatographySystems.PressurePoint("pin", 101300.0),
		GasChromatographySystems.PressurePoint("pout", 150000.0),
	]
	mod_rel = GasChromatographySystems.ModuleColumn("trap", Ltrap, 0.25e-3, 0.25e-6, "SLB5ms", TP_rel, opt_td)
	sys_rel = GasChromatographySystems.update_system(
		GasChromatographySystems.System("trap_rel", g, pp_rel, [mod_rel], GasChromatographySystems.Options(control="Pressure")),
	)
	p2_rel = GasChromatographySystems.build_pressure_squared_functions(sys_rel, GasChromatographySystems.solve_balance(sys_rel))
	par_rel0 = GasChromatographySystems.graph_to_parameters(sys_rel, p2_rel, db, solute)[1]
	par_rel = GasChromatographySystems.change_initial_spatial(par_rel0, [0.0], [x_park], [σ_park])
	pl_rel, sol_rel = GasChromatographySimulator.simulate(par_rel)

	(; x_park, σ_park, t_load, tR_rel=pl_rel.tR[1], exit_end=pl_rel.exit_end[1], τR_rel=pl_rel.τR[1])
end

# ╔═╡ 61d16b0c-10f8-422b-9a9a-17c26e117177
pl_rel

# ╔═╡ 1a2b3c4d-0008-4000-8000-000000000008
begin
	ts = sol_rel[1].t
	xs = [u[1] for u in sol_rel[1].u]
	p = plot(ts, xs; xlabel="t / s (release-local)", ylabel="x / m", label="x(t)", title="Hot reverse flush", legend=:topright)
	hline!(p, [0.0, Ltrap]; label="ends", ls=:dash, c=:black)
	p
end

# ╔═╡ 1a2b3c4d-0009-4000-8000-000000000009
md"""
### Notes

- `ModuleColumnOptions(domain=:temporal)` → GasChromatographySimulator `Options(domain=:temporal)` → one `simulate` entry point.
- Negative flow on a `:temporal` edge is allowed by `path_possible`.
- Junction handoff to a next column still uses `tR`/`τR` when applicable.
- Open-end ``t_{50}^\mathrm{rel}`` percentiles stay in Experiments for now.
"""

# ╔═╡ Cell order:
# ╠═1a2b3c4d-0001-4000-8000-000000000001
# ╟─1a2b3c4d-0002-4000-8000-000000000002
# ╠═1a2b3c4d-0003-4000-8000-000000000003
# ╟─1a2b3c4d-0004-4000-8000-000000000004
# ╠═1a2b3c4d-0005-4000-8000-000000000005
# ╟─1a2b3c4d-0006-4000-8000-000000000006
# ╠═1a2b3c4d-0007-4000-8000-000000000007
# ╠═61d16b0c-10f8-422b-9a9a-17c26e117177
# ╠═1a2b3c4d-0008-4000-8000-000000000008
# ╟─1a2b3c4d-0009-4000-8000-000000000009
