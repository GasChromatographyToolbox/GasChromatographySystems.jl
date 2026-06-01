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
end

# ╔═╡ ea3b2dd6-9a9f-4be7-9221-80a7467f3885
md"""
# Test flows and pressures in GCxGC system with pressure modulation
using a ModuleValve and a periodic valve program
"""

# ╔═╡ 0312d777-ee32-4aaa-9143-1bd7a85e7925
begin
	L1 = 3.0
    d1 = 0.1
    df1 = 0.1
    sp1 = "Rxi5ms"
    TP1 = GasChromatographySystems.TemperatureProgram([30.0, 1.0, 20.0, 300.0, 1.0])
    L2 = 1.0
    d2 = 0.1
    df2 = 0.1
    sp2 = "Rxi17SimMS"
    TP2 = GasChromatographySystems.TemperatureProgram([30.0, 1.0, 20.0, 300.0, 1.0])
    pin = 400000.0
    pout = 101300.0
    L_valve = 0.01
    d_open = 1.0 # mm
    d_closed = eps(Float64)
    TP_valve = 250.0
    VP = GasChromatographySystems.PeriodicValveProgram(1.0/3, 0.1, 1800.0)
    pmod = 400000.0
end

# ╔═╡ e6dc8d08-8e7a-4cb5-b693-42d0bcbcfe3b
sys = GasChromatographySystems.GCxGC_DPM(
    L1, d1, df1, sp1, TP1, 
    L2, d2, df2, sp2, TP2, 
    pin, pout, 
    L_valve, d_open, d_closed, TP_valve, VP, 
    pmod
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

# ╔═╡ Cell order:
# ╠═c423e03a-5d84-11f1-b415-c5b3781ea849
# ╠═ea3b2dd6-9a9f-4be7-9221-80a7467f3885
# ╠═0312d777-ee32-4aaa-9143-1bd7a85e7925
# ╠═e6dc8d08-8e7a-4cb5-b693-42d0bcbcfe3b
# ╠═27a5bb62-d10f-465c-8a3c-fd4a62e5d3d9
# ╠═f613a73b-9373-42b0-8365-9495e9f7aa60
# ╠═3c457b28-5c5a-4e09-91f7-2c0e3beaec5e
# ╠═7bc21000-16cb-4e00-9ef7-62e0ec6ce15b
# ╠═e4911037-5ae4-4404-9079-91fa3e016a82
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
