### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ f27920a2-00d5-11ef-0209-79647fcea22b
begin
	using CSV, DataFrames
	using Plots#, CairoMakie, GraphMakie
	using Graphs#, NetworkLayout, 
	using Symbolics
	using GasChromatographySimulator
	using PlutoUI
	using GasChromatographySystems
	using NLsolve
	using HypertextLiteral
	TableOfContents(depth=4)
end

# ╔═╡ 71902736-dea3-4c85-b067-e5d8d7c9ccb2
p2fun_flow_dict = include(download("https://raw.githubusercontent.com/GasChromatographyToolbox/GasChromatographySystems.jl/main/notebooks/FlowCalcPaper/p2fun_DeansSwitchSystem_3Fin_%CE%BB.jl"))

# ╔═╡ 64c1c01a-8a09-40b1-bfbd-a7668cde9b0b
html"""<style>
main {
    max-width: 80%;
    margin-left: 1%;
    margin-right: 10% !important;
}
"""

# ╔═╡ 13569378-fba2-4551-8435-aaf2f13bfc66
md"""
# Deans Switch System

This notebook shows the example of a GC-System with a Deans Switch and two outlets at different pressures and using different heated zones. This system is described in section 4.2.2 of the corresponding paper **Generalized flow calculator for the gas flow in a network of capillaries used in gas chromatography** by _Jan Leppert_, _Tillman Brehmer_, _Peter Boeker_ and _Matthias Wüst_.

In the first part, "Free setup", all lengths and diameters of the capillaries in the system can be selected, as well as the temperature program of the GC oven, transferlines to atmospheric and vacuum detector and of the ingoing auxillary flow capillaries. The dimensions of the Deans Switch itself are fixed using the values from [1]. Th injector is assumed as an on-column injector with the same temperature as the GC oven. 

In the second part, "Estimate lengths", the lengths of the two outgoing capillaries from the Deans Switch to the detectors are estimated using a non-linear numeric solver to achieve the two conditions of split ratio ``f`` and same hold-up times between the two paths for a fixed temperature.

The third part, "Estimate control flows", calculates the needed control flows from the auxillary inlets to achieve the selected conditions of split ratio ``f``  and same hold-up times for the different temperatures of the GC oven temperature program, defined in the first part. 

[1] Boeker, P., Leppert, J., Mysliwietz, B., Lammers, P. S., **Comprehensive theory of the Deans’ switch as a variable flow splitter: fluid mechanics, mass balance, and system behavior.** _Analytical Chemistry_ 2013, 85, 9021–9030.
"""

# ╔═╡ a8bbc6cc-3bcb-43c8-9b00-f2b3e814c136
md"""
## Definition of the System
The function to define the Deans Switch system:
"""

# ╔═╡ a6a03024-d77d-499a-9617-5c142dad4db6
function DeansSwitch(Ls, ds, dfs, sps, Tinj, TP, TTL1, TTL2, Taux, Fin, Faux1, Faux2, pin, paux1, paux2, pout1, pout2; name="DeansSwitchSystem", opt=GasChromatographySystems.Options(), kwargs...)
	# defining the graph:
	# the graph consists of 21 vertices:
	g = SimpleDiGraph(21)
	# adding the edges one after another:
	add_edge!(g, 1, 2) # 1 Pre column in Injector
	add_edge!(g, 2, 3) # 2 Pre column in Oven
	add_edge!(g, 3, 4) # 3 GC column
	add_edge!(g, 4, 5) # 4 Deans inlet
	add_edge!(g, 5, 6) # 5 Deans to Out1 and Aux1
	add_edge!(g, 5, 7) # 6 Deans to Out2 and Aux2
	add_edge!(g, 6, 8) # 7 Deans outlet to Out1
	add_edge!(g, 7, 9) # 8 Deans outlet to Out2
	add_edge!(g, 8, 14) # 9 1st cap to Out1
	add_edge!(g, 9, 15) # 10 2nd cap to Out2
	add_edge!(g, 10, 6) # 11 Deans Aux1 from Bypass to Out1
	add_edge!(g, 10, 11) # 12 Deans Bypass
	add_edge!(g, 11, 7) # 13 Deans Aux2 from Bypass to Out2
	add_edge!(g, 12, 10) # 14 Deans Aux1 to Bypass
	add_edge!(g, 13, 11) # 15 Deans Aux2 to Bypass
	add_edge!(g, 14, 16) # 16 2nd cap to Out1
	add_edge!(g, 15, 17) # 17 2nd cap to Out2
	add_edge!(g, 16, 18) # 18 3rd cap to Out1
	add_edge!(g, 17, 19) # 19 3rd cap to Out2
	add_edge!(g, 20, 12) # 20 inlet Aux1
	add_edge!(g, 21, 13) # 21 inlet Aux2
	# defining the pressure points:
	pp = [GasChromatographySystems.PressurePoint("p$(i)", NaN) for i=1:nv(g)] # all other pressures are not known (NaN)
	pp[1] = GasChromatographySystems.PressurePoint("p1", pin) # inlet 
	pp[18] = GasChromatographySystems.PressurePoint("p18", pout1) # outlet 1
	pp[19] = GasChromatographySystems.PressurePoint("p19", pout2) # autlet 2
	pp[20] = GasChromatographySystems.PressurePoint("p20", paux1) # aux1
	pp[21] = GasChromatographySystems.PressurePoint("p21", paux2) # aux1
	# defining the modules:
	modules = Array{GasChromatographySystems.AbstractModule}(undef, ne(g))
	# Inj
	modules[1] = GasChromatographySystems.ModuleColumn("1->2", Ls[1], ds[1]*1e-3, dfs[1]*1e-6, sps[1], Tinj, Fin/60e6; ng=true, kwargs...)
	# pre column
	modules[2] = GasChromatographySystems.ModuleColumn("2->3", Ls[2], ds[2]*1e-3, dfs[2]*1e-6, sps[2], TP; ng=true, kwargs...)
	# GC column
	modules[3] = GasChromatographySystems.ModuleColumn("3->4", Ls[3], ds[3]*1e-3, dfs[3]*1e-6, sps[3], TP; ng=true, kwargs...)
	# 4 Deans inlet
	modules[4] = GasChromatographySystems.ModuleColumn("4->5", 4.5e-3, 0.1376e-3, 0.0e-6, "", TP; ng=true, kwargs...)
	# 5 Deans to Out1 and Aux1
	modules[5] = GasChromatographySystems.ModuleColumn("5->6", 5.0e-3, 0.1376e-3, 0.0e-6, "", TP; ng=true, kwargs...)
	# 6 Deans to Out2 and Aux2
	modules[6] = GasChromatographySystems.ModuleColumn("5->7", 5.0e-3, 0.1376e-3, 0.0e-6, "", TP; ng=true, kwargs...) 
	# 7 Deans outlet to Out1
	modules[7] = GasChromatographySystems.ModuleColumn("6->8", 2.8e-3, 0.1376e-3, 0.0e-6, "", TP; ng=true, kwargs...)
	# 8 Deans outlet to Out2
	modules[8] = GasChromatographySystems.ModuleColumn("7->9", 2.8e-3, 0.1376e-3, 0.0e-6, "", TP; ng=true, kwargs...)
	# 9 1st cap to Out1
	modules[9] = GasChromatographySystems.ModuleColumn("8->14, X1", Ls[4], ds[4]*1e-3, dfs[4]*1e-6, sps[4], TP; ng=true, kwargs...)
	# 10 2nd cap to Out2
	modules[10] = GasChromatographySystems.ModuleColumn("9->15, Y1", Ls[7], ds[7]*1e-3, dfs[7]*1e-6, sps[7], TP; ng=true, kwargs...)
	# 11 Deans Aux1 from Bypass to Out1
	modules[11] = GasChromatographySystems.ModuleColumn("10->6", 5.0e-3, 0.1376e-3, 0.0e-6, "", TP; ng=true, kwargs...) 
	# 12 Deans Bypass
	modules[12] = GasChromatographySystems.ModuleColumn("10->11", 0.02, 0.0751e-3, 0.0e-6, "", TP; ng=true, kwargs...)
	# 13 Deans Aux2 from Bypass to Out2
	modules[13] = GasChromatographySystems.ModuleColumn("11->7", 5.0e-3, 0.1376e-3, 0.0e-6, "", TP; ng=true, kwargs...)
	# 14 Deans Aux1 to Bypass
	modules[14] = GasChromatographySystems.ModuleColumn("12->10", 4.8e-3, 0.1376e-3, 0.0e-6, "", TP; ng=true, kwargs...)
	# 15 Deans Aux2 to Bypass
	modules[15] = GasChromatographySystems.ModuleColumn("13->11", 4.8e-3, 0.1376e-3, 0.0e-6, "", TP; ng=true, kwargs...)
	# 16 2nd cap to Out1
	modules[16] = GasChromatographySystems.ModuleColumn("14->16, X2", Ls[5], ds[5]*1e-3, dfs[5]*1e-6, sps[5], TP; ng=true, kwargs...)
	# 17 2nd cap to Out2
	modules[17] = GasChromatographySystems.ModuleColumn("15->17, Y2", Ls[8], ds[8]*1e-3, dfs[8]*1e-6, sps[8], TP; ng=true, kwargs...)
	# 18 3rd cap to Out1
	modules[18] = GasChromatographySystems.ModuleColumn("16->18, X3", Ls[6], ds[6]*1e-3, dfs[6]*1e-6, sps[6], TTL1; ng=true, kwargs...)
	# 19 3rd cap to Out2
	modules[19] = GasChromatographySystems.ModuleColumn("17->19, Y3", Ls[9], ds[9]*1e-3, dfs[9]*1e-6, sps[9], TTL2; ng=true, kwargs...)
	# 20 inlet Aux1
	modules[20] = GasChromatographySystems.ModuleColumn("20->12", Ls[10], ds[10]*1e-3, dfs[10]*1e-6, sps[10], Taux, Faux1/60e6; ng=true, kwargs...)
	# 21 inlet Aux2
	modules[21] = GasChromatographySystems.ModuleColumn("21->13", Ls[11], ds[11]*1e-3, dfs[11]*1e-6, sps[11], Taux, Faux2/60e6; ng=true, kwargs...)
	# construct the system:
	sys = GasChromatographySystems.update_system(GasChromatographySystems.System(name, g, pp, modules, GasChromatographySystems.Options()))
	return sys
end

# ╔═╡ ec81a4e5-8916-46f5-bf1d-6a4781963bbb
begin # manipulate positions for the graph plot
	posin = Dict(1 => (-4, 3),
			2 => (-4, 2), 
			3 => (-4, 0),
			4 => (-1, 0),
			5 => (0, 0),
			6 => (1, 1),
			7 => (1, -1),
			8 => (1, 2),
			9 => (1, -2),
			10 => (2, 1),
			11 => (2, -1),
			12 => (3, 1),
			13 => (3, -1),
			14 => (2, 3),
			15 => (2, -3),
			16 => (4, 3),
			17 => (4, -3),
           	18 => (5,3),
           	19 => (5, -3),
           	20 => (5, 1),
			21 => (5, -1))
	stressl = GasChromatographySystems.Stress(;pin=posin, reltols=0.0, abstolx=0.0, iterations=100)
end;

# ╔═╡ 7345afa9-27c7-472f-b8dc-a55b2fb34c78
begin # manipulate positions for the graph plot
	posin_ = Dict(1 => (-5, 0),
			2 => (-4, 0), 
			3 => (-3, 0),
			4 => (-1, 0),
			5 => (0, 0),
			6 => (1, 1),
			7 => (1, -1),
			8 => (1, 2),
			9 => (1, -2),
			10 => (2, 1),
			11 => (2, -1),
			12 => (3, 1),
			13 => (3, -1),
			14 => (2, 3),
			15 => (2, -3),
			16 => (4, 3),
			17 => (4, -3),
           	18 => (5,3),
           	19 => (5, -3),
           	20 => (5, 1),
			21 => (5, -1))
	stressl_ = GasChromatographySystems.Stress(;pin=posin_, reltols=0.0, abstolx=0.0, iterations=100)
end;

# ╔═╡ 38afe7e3-84cc-4304-bcc2-5aabdd71e27d
sys0 = DeansSwitch([0.02, 1.0, 30.0, 0.5, 0.5, 1.0, 1.0, 0.5, 0.25, 0.5, 0.5],
					[0.53, 0.53, 0.32, 0.15, 0.18, 0.18, 0.18, 0.15, 0.15, 1.0, 1.0], 
					[fill(0.0, 2); 0.32; fill(0.0, 8)], 
					[fill("", 2); "DB5"; fill("", 8)], 
					40.0, 40.0, 280.0, 320.0, 30.0, 
					1.0, 2.0, 3.0, 	
					NaN, NaN, NaN, 
					101300.0, 0.0; 
					name="DeansSwitch_Placeholder");

# ╔═╡ dc2e6dfb-7f5a-462f-8fe8-912afa76478f
md"""
The flow balance equations of the Deans Switch system:
"""

# ╔═╡ 5bcfa4d6-6c4b-44ac-9473-f12d410b0ac7
bal_eq = GasChromatographySystems.substitute_unknown_flows(sys0; mode="κ")

# ╔═╡ e2df094e-ab54-43a6-8a9d-778031a576c0
md"""
Loading the solutions for the squared pressures fom a file:
"""

# ╔═╡ b0451d2d-5a1b-4027-affb-f6531d11d1a5
md"""
Evaluate the loaded solutions to an array of Julia functions:
"""

# ╔═╡ 6c701a8d-d6eb-4303-b000-ce864f222658
p2fun_flow_eval = eval.(p2fun_flow_dict["p2fun"])

# ╔═╡ b0750420-37d2-4e7f-9c1e-584a4890db00
md"""
The graph of the Deans Switch system:
"""

# ╔═╡ bbeb3406-ab4d-45c4-9a87-714c6cdee3bd
GasChromatographySystems.plot_graph(sys0; lay=stressl_, elabels_fontsize=8, nlabels_fontsize=12)#; GasChromatographySystems.save("graph.png", gr)

# ╔═╡ be454bf6-e220-410f-b7d1-d7e892b7b641
md"""
## 1. Free setup
"""

# ╔═╡ fd6973a3-68df-49ea-9988-8fc2f26a25e7
md"""
### Temperature Program

Number of ramps: $(@bind n_ramp confirm(NumberField(0:1:100; default=3)))
"""

# ╔═╡ 2fdf479e-1c3b-46a6-9994-6c5ea465ff10
md"""
#### Resulting Flows
"""

# ╔═╡ 3f570097-406c-4d84-947f-419678aa279b
md"""
## 2. Estimate lengths

of capillaries X1 and Y1 for the following selected split of the flows at the Deans Switch ``f`` at the selcted oven temperature ``T``. All other parameters, including the diameters of capillaries X1 and Y1 are used from the selected values in section "Free setup".
"""

# ╔═╡ aa9810b9-3d2c-4c35-8267-c93998414b20
md"""
To estimate the length of the two capillaries for the two conditions of the same hold-up times and the fixed split between the flows at the split results in the system of two equations:

```math
t_{\text{M,O}}(L_\text{X1}, L_\text{Y1}) - t_{\text{M,MS}}(L_\text{X1}, L_\text{Y1}) = 0
```
```math
F_\text{E}(L_\text{X1}, L_\text{Y1}) - f F_\text{F}(L_\text{X1}, L_\text{Y1}) = 0
```

As the dependencies of the hold-up times on the lengths are quadratic and for the flows are proportional to the invers lengths, this is a non-linear system of equations in regard to ``L_\text{X1}`` and ``L_\text{Y1}``. This is done using the function `nlsolve(eq_sys, init, method=:trust_region)` of the `NLsolve.jl` package using the Newton method with trust regions. This is done in this notebook with the below defined function `estimate_LX1_LY1()`.
"""

# ╔═╡ 64055873-b815-43fe-adcc-fea47a82768d
begin 
	
	function ΔtM(Ls, ds, dfs, sps, Tinj, Toven, TTL1, TTL2, Taux, Fin, Faux1, Faux2, pout1, pout2, p2fun)
		sys_new = DeansSwitch(Ls, ds, dfs, sps, Tinj, Toven, TTL1, TTL2, Taux, Fin, Faux1, Faux2, NaN, NaN, NaN, pout1, pout2)

		tMp_Det1, tMp_Det2 = GasChromatographySystems.holdup_time_path(sys_new, p2fun, 2)
		#ΔtM(t) = tMp_Det1(t) - tMp_Det2(t)
		tM = GasChromatographySystems.holdup_time_functions(sys_new, p2fun)
		#tMX(t) = tM[9](t) + tM[16](t) + tM[18](t)
		#tMY(t) = tM[10](t) + tM[17](t) + tM[19](t)
		paths = GasChromatographySystems.all_paths(sys_new.g, 2)[2]
		return tMp_Det1, tMp_Det2, tM, paths, sys_new
	end

	function ΔtM(LY1, LY2, dY1, T, sys, p2fun)
		Ls = [sys.modules[1].L, sys.modules[2].L, sys.modules[3].L, sys.modules[9].L, sys.modules[16].L, sys.modules[18].L, LY1, LY2, sys.modules[19].L, sys.modules[20].L, sys.modules[21].L]
		ds = [sys.modules[1].d*1e3, sys.modules[2].d*1e3, sys.modules[3].d*1e3, sys.modules[9].d*1e3, sys.modules[16].d*1e3, sys.modules[18].d*1e3, dY1, sys.modules[17].d*1e3, sys.modules[19].d*1e3, sys.modules[20].d*1e3, sys.modules[21].d*1e3]
		dfs = [sys.modules[i].df*1e6 for i=[1, 2, 3, 9, 16, 18, 10, 17, 19, 20, 21]]
		sps = [sys.modules[i].sp for i=[1, 2, 3, 9, 16, 18, 10, 17, 19, 20, 21]]
		Tinj = sys.modules[1].T
		Toven = T#sys.modules[2].T
		TTL1 = sys.modules[18].T
		TTL2 = sys.modules[19].T
		Taux = sys.modules[20].T
		
		pin = sys.pressurepoints[1].P
		pout1 = sys.pressurepoints[18].P
		pout2 = sys.pressurepoints[19].P

		Fin = sys.modules[1].F*60e6
		Faux1 = sys.modules[20].F*60e6
		Faux2 = sys.modules[21].F*60e6
		
		tMp_Det1, tMp_Det2, tM, paths, sys_new = ΔtM(Ls, ds, dfs, sps, Tinj, Toven, TTL1, TTL2, Taux, Fin, Faux1, Faux2, pout1, pout2, p2fun)
		return tMp_Det1, tMp_Det2, tM, paths, sys_new
	end

	function ΔtM(LX1, LX2, LY1, LY2, dX1, dY1, T, sys, p2fun)
		Ls = [sys.modules[1].L, sys.modules[2].L, sys.modules[3].L, LX1, LX2, sys.modules[18].L, LY1, LY2, sys.modules[19].L, sys.modules[20].L, sys.modules[21].L]
		ds = [sys.modules[1].d*1e3, sys.modules[2].d*1e3, sys.modules[3].d*1e3, dX1, sys.modules[16].d*1e3, sys.modules[18].d*1e3, dY1, sys.modules[17].d*1e3, sys.modules[19].d*1e3, sys.modules[20].d*1e3, sys.modules[21].d*1e3]
		dfs = [sys.modules[i].df*1e6 for i=[1, 2, 3, 9, 16, 18, 10, 17, 19, 20, 21]]
		sps = [sys.modules[i].sp for i=[1, 2, 3, 9, 16, 18, 10, 17, 19, 20, 21]]
		Tinj = sys.modules[1].T
		Toven = T#sys.modules[2].T
		TTL1 = sys.modules[18].T
		TTL2 = sys.modules[19].T
		Taux = sys.modules[20].T
		
		pin = sys.pressurepoints[1].P
		pout1 = sys.pressurepoints[18].P
		pout2 = sys.pressurepoints[19].P

		Fin = sys.modules[1].F*60e6
		Faux1 = sys.modules[20].F*60e6
		Faux2 = sys.modules[21].F*60e6
		
		tMp_Det1, tMp_Det2, tM, paths, sys_new = ΔtM(Ls, ds, dfs, sps, Tinj, Toven, TTL1, TTL2, Taux, Fin, Faux1, Faux2, pout1, pout2, p2fun)
		return tMp_Det1, tMp_Det2, tM, paths, sys_new
	end

end

# ╔═╡ 32cb7766-20b2-45c4-9909-cf97153aa050
function estimate_LX1_LY1(LX2, LY2, dX1, dY1, f, T, sys, p2fun)

	# the system of nonlinear equations 
	function eq_sys(x)
		tMp_Det1, tMp_Det2, tMs, paths, sys_new = ΔtM(abs(x[1]), LX2, abs(x[2]), LY2, dX1, dY1, T, sys, p2fun)
		Fs = GasChromatographySystems.flow_functions(sys_new, p2fun)
		tM_Det1 = if tMs[5](0.0) < 0.0
			1e6#Inf
		else
			tMp_Det1(0.0)
		end
		tM_Det2 = if tMs[6](0.0) < 0.0
			1e6#Inf
		else
			tMp_Det2(0.0)
		end
		return [tM_Det1 - tM_Det2, Fs[5](0.0)*60e6 - f*Fs[6](0.0)*60e6]
	end
	sol = nlsolve(eq_sys, [1.0, 1.0], method=:trust_region)
	LX1, LY1 = if sol.f_converged == true
		abs(sol.zero[1]), abs(sol.zero[2])
	else
		NaN, NaN
	end
	return LX1, LY1, sol, eq_sys(sol.zero), eq_sys(abs.(sol.zero))
end

# ╔═╡ deceb6b7-08f7-4c45-9663-1062b58fefb5
md"""
#### Resulting lengths 
The result of the estimation:
"""

# ╔═╡ 835ea802-4256-479f-9267-0aaefd33d1b6
md"""
#### Resulting Flows
"""

# ╔═╡ 92b0763f-24cd-4bed-a582-76cea406079e
md"""
## 3. Estimate control flows
To achive a constant split ratio throughout the temperature program, the two auxillary flows/pressures have to be changed. The previously estimated lengths of capillaries X1 and Y1 are used. The flow split ratio can be changed here:
"""

# ╔═╡ 5e07ad07-f3de-45ac-b98c-91ed1ddce44c
md"""
To estimate two auxillary pressures the same approach as before with the capillary lengths is used. The two conditions of the same hold-up times and the fixed split between the flows at the split results in the system of two equations:

```math
t_{\text{M,O}}(F_\text{aux1}, F_\text{aux2}) - t_{\text{M,MS}}(F_\text{aux1}, F_\text{aux2}) = 0
```
```math
F_\text{E}(F_\text{aux1}, F_\text{aux2}) - f F_\text{F}(F_\text{aux1}, F_\text{aux2}) = 0
```

This system of equations is solved by using the function `nlsolve(eq_sys, init, method=:trust_region)` of the `NLsolve.jl` package using the Newton method with trust regions. This is done in this notebook with the below defined function `estimate_F()`.
"""

# ╔═╡ 00c438cc-58c5-4e30-bd07-41f928d88634
function ΔtM_F(Faux1, Faux2, Toven, sys, p2fun)
	# it is assumed, that Tinj=Toven (typical for on-column injection)
	index = [1, 2, 3, 9, 16, 18, 10, 17, 19, 20, 21]
	Ls = [sys.modules[i].L for i=index]
	ds = [sys.modules[i].d*1e3 for i=index]
	dfs = [sys.modules[i].df*1e6 for i=index]
	sps = [sys.modules[i].sp for i=index]
	
	TTL1 = sys.modules[18].T
	TTL2 = sys.modules[19].T
	Taux = sys.modules[20].T
	
	pin = sys.pressurepoints[1].P
	pout1 = sys.pressurepoints[18].P
	pout2 = sys.pressurepoints[19].P

	Fin = sys.modules[1].F*60e6
	
	tMp_Det1, tMp_Det2, tM, paths, sys_new = ΔtM(Ls, ds, dfs, sps, Toven, Toven, TTL1, TTL2, Taux, Fin, Faux1, Faux2, pout1, pout2, p2fun)
	return tMp_Det1, tMp_Det2, tM, paths, sys_new
end

# ╔═╡ 30e2fc65-642b-4210-bcac-2fdea6d1946f
function estimate_F(Toven, f, sys, p2fun)
	
	function eq_sys(x)
		tMp_Det1, tMp_Det2, tMs, paths, sys_new = ΔtM_F(x[1], x[2], Toven, sys, p2fun)
		Fs = GasChromatographySystems.flow_functions(sys_new, p2fun)
		tM_Det1 = if tMs[5](0.0) < 0.0
			1e6#Inf
		else
			tMp_Det1(0.0)
		end
		tM_Det2 = if tMs[6](0.0) < 0.0
			1e6#Inf
		else
			tMp_Det2(0.0)
		end
		return [tM_Det1 - tM_Det2, Fs[5](0.0)*60e6 - f*Fs[6](0.0)*60e6]
	end
	sol = nlsolve(eq_sys, [1.0, 1.0], method=:trust_region)
	Faux1, Faux2 = if sol.f_converged == true
		sol.zero[1], sol.zero[2]
	else
		NaN, NaN
	end
	return Faux1, Faux2, sol, eq_sys(sol.zero), eq_sys(abs.(sol.zero))
end

# ╔═╡ afb33442-9cd1-4378-ac25-f67533da8206
md"""
#### Resulting Flows
"""

# ╔═╡ ba9f11c7-7857-4008-bb03-614e488e731f
md"""
With the changing auxillary flows during the temperature program a constant flow split ratio during the whole temperature program is achieved (variations in the plot are due to rounding inaccuracies.
"""

# ╔═╡ d3759b64-ebb1-4daa-ba6a-302630d87f85
md"""
# End
"""

# ╔═╡ fbd5ee2f-092f-429b-8797-660d1107308f
function UI_TP(n_ramp)
	PlutoUI.combine() do Child
		@htl("""
		<table>
			<tr>
				<th>ramp [°C/min]</th>
				<th>T [°C]</th>
				<th>hold [min]</th>
			</tr>
			<tr>
				<td></td>
				<td><center>$(Child(NumberField(0.0:0.1:400.0; default=40.0)))</center></td>
				<td><center>$(Child(NumberField(0.0:0.1:100.0; default=1.0)))</center></td>
			</tr>
			$([
				@htl("
					<tr>
						<td><center>$(Child(NumberField(0.0:0.1:100.0; default=(i-1)*10.0 + 5.0)))</center></td>
						<td><center>$(Child(NumberField(0.0:0.1:400.0; default=((i+1)*100.0))))</center></td>
						<td><center>$(Child(NumberField(0.0:0.1:100.0; default=(i+1)*1.0)))</center></td>
					</tr>
				")
				for i=1:n_ramp
			])
		</table>
		""")
	end
end

# ╔═╡ 64146877-f6b4-4d1a-aeb9-147809b756ca
@bind Tprog_values confirm(UI_TP(n_ramp))

# ╔═╡ f0a02e52-a0c1-40e9-b479-b8162c63f445
begin
	ts, Ts = GasChromatographySimulator.conventional_program(Tprog_values)
	TP = GasChromatographySystems.TemperatureProgram(ts, Ts);
end;

# ╔═╡ e41477bc-247e-454e-864d-68a223ad59a2
md"""
Time slider (in s): $(@bind select_t Slider(0.0:1.0:sum(TP.time_steps); show_value=true))
"""

# ╔═╡ 0e476729-f41a-4d70-aeed-04632354dd31
md"""
Time slider (in s): $(@bind select_t_est Slider(0.0:1.0:sum(TP.time_steps); show_value=true))
"""

# ╔═╡ 5abb07b4-55d3-4632-8a7b-ca0233454dde
md"""
Time slider (in s): $(@bind select_t_aux Slider(0.0:1.0:sum(TP.time_steps); show_value=true))
"""

# ╔═╡ ed7ea63d-66f7-4fe2-be14-6f916e548d10
function UI_flowsplit(def, title)
	PlutoUI.combine() do Child
		@htl("""
		<table>
			<caption> $(title) </capiton>
			<tr>
				<td><i>f</i></td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[1])))</center></td>
			</tr>
		</table>
		""")
	end
end

# ╔═╡ 3b92560e-11d1-42a2-8220-5dac046025d6
@bind flowsplit confirm(UI_flowsplit((1.0), "Split ratio, est. Faux"))

# ╔═╡ 9eaf6a93-a5d4-4cda-a281-d9e95b8df697
function UI_flowsplit_T(def, title)
	PlutoUI.combine() do Child
		@htl("""
		<table>
			<caption> $(title) </capiton>
			<tr>
				<td><i>f</i></td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[1])))</center></td>
			</tr>
			<tr>
				<td>at <i>T</i> [°C]</td>
				<td><center>$(Child(NumberField(0.0:0.1:400.0; default=def[2])))</center></td>
			</tr>
		</table>
		""")
	end
end

# ╔═╡ dfb62248-29ba-4cf3-8c28-0f4ff47e5d46
@bind flowsplit_T_values confirm(UI_flowsplit_T((1.0, 340.0), "Split ratio, est. L"))

# ╔═╡ 3d49a31c-ecd9-4eb0-8790-f2543706ead5
md"""
As the auxillary inlet flows are constant (changes in the line are due to small rounding errors), the split ratio for other temperatures, than the selected temperature ($(flowsplit_T_values[2])°C) for the length estimation, is different from the selected ratio. 
"""

# ╔═╡ 385d3632-d091-4fde-8045-545733ce6fee
function UI_L_d_F(def, title)
	PlutoUI.combine() do Child
		@htl("""
		<table>
			<caption> $(title) </capiton>
			<tr>
				<td><i>L</i> [m]</td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[1])))</center></td>
			</tr>
			<tr>
				<td><i>d</i> [mm]</td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[2])))</center></td>
			</tr>
			<tr>
				<td><i>F</i> [mL/min]</td>
				<td><center>$(Child(NumberField(0.0:0.01:10.0; default=def[3])))</center></td>
			</tr>
		</table>
		""")
	end
end

# ╔═╡ cf3c98f1-3cb7-41ef-b7fc-eabded35f1e1
function UI_L_d_F_T(def, title)
	PlutoUI.combine() do Child
		@htl("""
		<table>
			<caption> $(title) </capiton>
			<tr>
				<td><i>L</i> [m]</td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[1])))</center></td>
			</tr>
			<tr>
				<td><i>d</i> [mm]</td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[2])))</center></td>
			</tr>
			<tr>
				<td><i>F</i> [mL/min]</td>
				<td><center>$(Child(NumberField(0.0:0.01:10.0; default=def[3])))</center></td>
			</tr>
			<tr>
				<td><i>T</i> [°C]</td>
				<td><center>$(Child(NumberField(0.0:0.1:400.0; default=def[4])))</center></td>
			</tr>
		</table>
		""")
	end
end

# ╔═╡ fbd3ff22-c286-4c9e-b60b-3f6373ca40e5
function UI_L_d_pin(def, title)
	PlutoUI.combine() do Child
		@htl("""
		<table>
			<caption> $(title) </capiton>
			<tr>
				<td><i>L</i> [m]</td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[1])))</center></td>
			</tr>
			<tr>
				<td><i>d</i> [mm]</td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[2])))</center></td>
			</tr>
			<tr>
				<td><i>p</i><sub>in</sub> [kPa]</td>
				<td><center>$(Child(NumberField(0.0:0.1:800.0; default=def[3])))</center></td>
			</tr>
		</table>
		""")
	end
end

# ╔═╡ 2138443d-5a25-4e4d-92ee-9f19f931ea18
function UI_L_d(def, title)
	PlutoUI.combine() do Child
		@htl("""
		<table>
			<caption> $(title) </capiton>
			<tr>
				<td><i>L</i> [m]</td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[1])))</center></td>
			</tr>
			<tr>
				<td><i>d</i> [mm]</td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[2])))</center></td>
			</tr>
		</table>
		""")
	end
end

# ╔═╡ ea46e236-dbf9-41f8-a661-ffe3c76affec
function UI_L_d_phase(def, title; sp=["Rxi5SilMS", "Rxi17SilMS", ""])
	PlutoUI.combine() do Child
		@htl("""
		<table>
			<caption> $(title) </capiton>
			<tr>
				<td><i>L</i> [m]</td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[1])))</center></td>
			</tr>
			<tr>
				<td><i>d</i> [mm]</td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[2])))</center></td>
			</tr>
			<tr>
				<td><i>d</i><sub>f</sub> [µm]</td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[3])))</center></td>
			</tr>
			<tr>
				<td>phase</td>
				<td><center>$(Child(Select(sp; default=sp[1])))</center></td>
			</tr>
		</table>
		""")
	end
end

# ╔═╡ 89499fa9-ca50-4fcf-9f06-7f3eca7c5194
function UI_L_d_phase_T(def, title; sp=["Rxi5SilMS", "Rxi17SilMS", ""])
	PlutoUI.combine() do Child
		@htl("""
		<table>
			<caption> $(title) </capiton>
			<tr>
				<td><i>L</i> [m]</td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[1])))</center></td>
			</tr>
			<tr>
				<td><i>d</i> [mm]</td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[2])))</center></td>
			</tr>
			<tr>
				<td><i>d</i><sub>f</sub> [µm]</td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[3])))</center></td>
			</tr>
			<tr>
				<td>phase</td>
				<td><center>$(Child(Select(sp; default=sp[1])))</center></td>
			</tr>
			<tr>
				<td><i>T</i> [°C]</td>
				<td><center>$(Child(NumberField(0.0:0.1:400.0; default=def[4])))</center></td>
			</tr>
		</table>
		""")
	end
end

# ╔═╡ 4d9d9e40-0e9b-4044-b764-3fabc18b0970
function UI_L_d_T(def, title)
	PlutoUI.combine() do Child
		@htl("""
		<table>
			<caption> $(title) </capiton>
			<tr>
				<td><i>L</i> [m]</td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[1])))</center></td>
			</tr>
			<tr>
				<td><i>d</i> [mm]</td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[2])))</center></td>
			</tr>
			<tr>
				<td><i>T</i> [°C]</td>
				<td><center>$(Child(NumberField(0.0:0.1:400.0; default=def[3])))</center></td>
			</tr>
		</table>
		""")
	end
end

# ╔═╡ 26064dc0-fe8b-4cbb-9f95-8a8ef3d69b7a
function UI_L_d_T_pout(def, title)
	PlutoUI.combine() do Child
		@htl("""
		<table>
			<caption> $(title) </capiton>
			<tr>
				<td><i>L</i> [m]</td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[1])))</center></td>
			</tr>
			<tr>
				<td><i>d</i> [mm]</td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[2])))</center></td>
			</tr>
			<tr>
				<td><i>T</i> [°C]</td>
				<td><center>$(Child(NumberField(0.0:0.1:400.0; default=def[3])))</center></td>
			</tr>
			<tr>
				<td><i>p</i><sub>out</sub></td>
				<td><center>$(Child(Select(["atm", "vac"]; default=def[4])))</center></td>
			</tr>
		</table>
		""")
	end
end

# ╔═╡ b9a2552d-d8e9-4d94-b0a5-097e4b03a843
function UI_system_flow()
	PlutoUI.combine() do Child
		@htl("""
		<table>
			<tr>
				<td></td>
				<td></td>
				<td></td>
				<td>$(Child(UI_L_d((1.0, 0.15), "X1")))</td>
				<td>$(Child(UI_L_d((0.5, 0.25), "X2")))</td>
				<td>$(Child(UI_L_d_T_pout((1.0, 0.25, 280.0, "atm"), "X3")))</td>
			</tr>
			<tr>
				<td></td>
				<td></td>
				<td></td>
				<td></td>
				<td></td>
				<td>$(Child(UI_L_d_F_T((0.5, 1.0, 0.5, 30.0), "aux1")))</td>
			</tr>
			<tr>
				<td>$(Child(UI_L_d_F((0.02, 0.53, 1.0), "Injector")))</td>
				<td>$(Child(UI_L_d((1.0, 0.53), "Precolumn")))</td>
				<td>$(Child(UI_L_d_phase((30.0, 0.32, 0.32), "GC column")))</td>
				<td></td>
				<td></td>
				<td></td>
			</tr>
			<tr>
				<td></td>
				<td></td>
				<td></td>
				<td></td>
				<td></td>
				<td>$(Child(UI_L_d_F_T((0.5, 1.0, 0.5, 30.0), "aux2")))</td>
			</tr>
			<tr>
				<td></td>
				<td></td>
				<td></td>
				<td>$(Child(UI_L_d((1.0, 0.25), "Y1")))</td>
				<td>$(Child(UI_L_d((1.0, 0.15), "Y2")))</td>
				<td>$(Child(UI_L_d_T_pout((0.25, 0.15, 320.0, "vac"), "Y3")))</td>
			</tr>
		</table>
		""")
	end
end

# ╔═╡ 39ea740d-a6de-4e68-aaa1-6517935c3a93
@bind values confirm(UI_system_flow())

# ╔═╡ ed57c1f5-76de-48a5-bd1e-3e7172b2ad48
function outpressure_select(pout_str)
	pout = if pout_str == "atm"
		101300.0
	else
		0.0
	end
	return pout
end

# ╔═╡ 2dc31d04-abdf-478d-a7cc-eea66b00f2ed
function DeansSwitch_select(values, TP)
	pX3 = outpressure_select(values[3][4])
	pY3 = outpressure_select(values[11][4])
		
	sys = DeansSwitch([values[5][1], values[6][1], values[7][1], values[1][1], values[2][1], values[3][1], values[9][1], values[10][1], values[11][1], values[4][1], values[8][1]],
					[values[5][2], values[6][2], values[7][2], values[1][2], values[2][2], values[3][2], values[9][2], values[10][2], values[11][2], values[4][2], values[8][2]], 
					[fill(0.0, 2); values[7][3]; fill(0.0, 8)], 
					[fill("", 2); values[7][4]; fill("", 8)], 
					TP, TP, values[3][3], values[11][3], values[4][4], 
					values[5][3], values[4][3], values[8][3], 	
					NaN, NaN, NaN, 
					pX3, pY3; 
					name="DeansSwitch_Flow")
	return sys
end

# ╔═╡ 2a413dc2-9c4f-4608-8f40-034fba017163
sys_flow = DeansSwitch_select(values, TP);

# ╔═╡ dc33f080-4ccf-42ad-b8a8-b192da84f1ec
Ffs = GasChromatographySystems.flow_functions(sys_flow, p2fun_flow_eval);

# ╔═╡ e99c055a-0cac-40a1-8969-3fb5e8ba5825
begin
	p_Fd = Plots.plot(cumsum(TP.time_steps)./60, Ffs[18].(cumsum(TP.time_steps))*60e6; label="F_O", c=1)
	Plots.plot!(p_Fd, cumsum(TP.time_steps)./60, Ffs[19].(cumsum(TP.time_steps))*60e6; label="F_MS", c=2, xlabel="time in min", ylabel="flow in mL/min", title="Flows to detectors")
	p_Faux = Plots.plot(cumsum(TP.time_steps)./60, Ffs[20].(cumsum(TP.time_steps))*60e6; label="F_aux1", c=3)
	Plots.plot!(p_Faux, cumsum(TP.time_steps)./60, Ffs[21].(cumsum(TP.time_steps))*60e6; label="F_aux2", c=4, xlabel="time in min", ylabel="flow in mL/min", title="Auxillary flows")
	Plots.plot(p_Fd, p_Faux, size=(800,400), legend=:top)
end

# ╔═╡ 795fee94-00d6-442b-aa32-19a296e9c182
begin
	p_fr = Plots.plot(cumsum(TP.time_steps)./60, Ffs[5].(cumsum(TP.time_steps))./Ffs[6].(cumsum(TP.time_steps)); label="", c=1, title="Flow split ratio f", size=(400,400))
end

# ╔═╡ 25687b22-dcc3-4571-af9e-e9f260133c3f
pfs = GasChromatographySystems.pressure_functions(sys_flow, p2fun_flow_eval);

# ╔═╡ a72bdfbf-4bc4-4d25-991b-ee55c29ae57f
begin
	p_pins = Plots.plot(cumsum(TP.time_steps)./60, pfs[1].(cumsum(TP.time_steps)); label="p_in", c=1)
	Plots.plot!(p_pins, cumsum(TP.time_steps)./60, pfs[20].(cumsum(TP.time_steps)); label="p_aux1", c=2, linewidth=2, xlabel="time in min", ylabel="pressure in Pa", title="Inlet pressures")
	Plots.plot!(p_pins, cumsum(TP.time_steps)./60, pfs[21].(cumsum(TP.time_steps)); label="p_aux2", c=3, size=(600,400), legend=:topleft)
end

# ╔═╡ 5d4d3407-0136-44ba-8130-774cd09aa58e
estL = estimate_LX1_LY1(values[2][1], values[10][1], values[1][2], values[9][2], flowsplit_T_values[1], flowsplit_T_values[2], sys_flow, p2fun_flow_eval)

# ╔═╡ e72d0302-b19a-4792-8aa0-3a2df88061a2
md"""
The estimated capillary lengths are:
* ``L_\text{X1} = `` $(round(estL[1]; sigdigits=3)) m
* ``L_\text{Y1} = `` $(round(estL[2]; sigdigits=3)) m
"""

# ╔═╡ 7cc3428c-6b3e-4acd-9c72-2a0b19270c4c
function DeansSwitch_update(sys, est)
	Ls = [sys.modules[1].L, sys.modules[2].L, sys.modules[3].L, est[1], sys.modules[16].L, sys.modules[18].L, est[2], sys.modules[17].L, sys.modules[19].L, sys.modules[20].L, sys.modules[21].L]
	ds = [sys.modules[i].d*1e3 for i=[1, 2, 3, 9, 16, 18, 10, 17, 19, 20, 21]]
	dfs = [sys.modules[i].df*1e6 for i=[1, 2, 3, 9, 16, 18, 10, 17, 19, 20, 21]]
	sps = [sys.modules[i].sp for i=[1, 2, 3, 9, 16, 18, 10, 17, 19, 20, 21]]
	Tinj = sys.modules[1].T
	Toven = sys.modules[2].T
	TTL1 = sys.modules[18].T
	TTL2 = sys.modules[19].T
	Taux = sys.modules[20].T
	
	pin = sys.pressurepoints[1].P
	pout1 = sys.pressurepoints[18].P
	pout2 = sys.pressurepoints[19].P

	Fin = sys.modules[1].F*60e6
	Faux1 = sys.modules[20].F*60e6
	Faux2 = sys.modules[21].F*60e6
	sys_up = DeansSwitch(Ls, ds, dfs, sps, Tinj, Toven, TTL1, TTL2, Taux, Fin, Faux1, Faux2, pin, NaN, NaN, pout1, pout2; name="DeansSwitch_estimate_X1_Y1")
	return sys_up
end

# ╔═╡ 171f9c3f-d986-479b-ad21-46f890187cd8
sys_est = DeansSwitch_update(sys_flow, estL);

# ╔═╡ 7b42338e-85ae-4418-aa39-34ddcdd80aaa
Ffs_est = GasChromatographySystems.flow_functions(sys_est, p2fun_flow_eval);

# ╔═╡ 8674b5b0-37b0-4447-b1bc-4d7bf2cd8c33
begin
	p_Fd_est = Plots.plot(cumsum(TP.time_steps)./60, Ffs_est[18].(cumsum(TP.time_steps))*60e6; label="F_O", c=1)
	Plots.plot!(p_Fd_est, cumsum(TP.time_steps)./60, Ffs_est[19].(cumsum(TP.time_steps))*60e6; label="F_MS", c=2, xlabel="time in min", ylabel="flow in mL/min", title="Flows to detectors")
	p_Faux_est = Plots.plot(cumsum(TP.time_steps)./60, Ffs_est[20].(cumsum(TP.time_steps))*60e6; label="F_aux1", c=3)
	Plots.plot!(p_Faux_est, cumsum(TP.time_steps)./60, Ffs_est[21].(cumsum(TP.time_steps))*60e6; label="F_aux2", c=4, xlabel="time in min", ylabel="flow in mL/min", title="Auxillary flows")
	Plots.plot(p_Fd_est, p_Faux_est, size=(800,400), legend=:top)
end

# ╔═╡ 05b6291f-d6f6-40bf-a090-97c531c10c49
begin
	p_fr_est = Plots.plot(cumsum(TP.time_steps)./60, Ffs_est[5].(cumsum(TP.time_steps))./Ffs_est[6].(cumsum(TP.time_steps)); label="", c=1, title="Flow split ratio f", size=(400,400))
end

# ╔═╡ 0c670989-8354-4fc6-8578-8032c79ac085
pfs_est = GasChromatographySystems.pressure_functions(sys_est, p2fun_flow_eval);

# ╔═╡ fafa9077-92ee-4d02-80e6-9c4db34388b8
begin
	p_pins_est = Plots.plot(cumsum(TP.time_steps)./60, pfs_est[1].(cumsum(TP.time_steps)); label="p_in", c=1)
	Plots.plot!(p_pins_est, cumsum(TP.time_steps)./60, pfs_est[20].(cumsum(TP.time_steps)); label="p_aux1", c=2, linewidth=2, xlabel="time in min", ylabel="pressure in Pa", title="Inlet pressures")
	Plots.plot!(p_pins_est, cumsum(TP.time_steps)./60, pfs_est[21].(cumsum(TP.time_steps)); label="p_aux2", c=3, size=(600,400), legend=:topleft)
end

# ╔═╡ e82d7276-2601-4ce9-b0f4-d15bbe8214b3
function DeansSwitch_update_Faux(sys, Faux1_, Faux2_)
	Ls = [sys.modules[i].L for i=[1, 2, 3, 9, 16, 18, 10, 17, 19, 20, 21]]
	ds = [sys.modules[i].d*1e3 for i=[1, 2, 3, 9, 16, 18, 10, 17, 19, 20, 21]]
	dfs = [sys.modules[i].df*1e6 for i=[1, 2, 3, 9, 16, 18, 10, 17, 19, 20, 21]]
	sps = [sys.modules[i].sp for i=[1, 2, 3, 9, 16, 18, 10, 17, 19, 20, 21]]
	Tinj = sys.modules[1].T
	Toven = sys.modules[2].T
	TTL1 = sys.modules[18].T
	TTL2 = sys.modules[19].T
	Taux = sys.modules[20].T
	
	pin = sys.pressurepoints[1].P
	pout1 = sys.pressurepoints[18].P
	pout2 = sys.pressurepoints[19].P

	Fin = sys.modules[1].F*60e6
	Faux1 = Faux1_
	Faux2 = Faux2_
	sys_up = DeansSwitch(Ls, ds, dfs, sps, Tinj, Toven, TTL1, TTL2, Taux, Fin, Faux1, Faux2, pin, NaN, NaN, pout1, pout2; name="DeansSwitch_estimate_X1_Y1")
	return sys_up
end

# ╔═╡ faaa3135-0135-47f2-ab1b-9f92cdcb309b
begin
	Faux(T) = estimate_F(T, values[9][1], sys_est, p2fun_flow_eval)[1:2]
	T(t) = GasChromatographySystems.module_temperature(sys_est.modules[1], sys_est)[5](0.0, t) - 273.15
	Faux_1(t) = Faux(T(t))[1]
	Faux_2(t) = Faux(T(t))[2]
	sys_est_Faux(t) = DeansSwitch_update_Faux(sys_est, Faux_1(t), Faux_2(t))
end;

# ╔═╡ 3c0680df-79b7-4619-8625-d5cfd76e5779
begin
	Ffs_aux_o = Array{Float64}(undef, length(cumsum(TP.time_steps)))
	Ffs_aux_ms = Array{Float64}(undef, length(cumsum(TP.time_steps)))
	Ffs_aux_5 = Array{Float64}(undef, length(cumsum(TP.time_steps)))
	Ffs_aux_6 = Array{Float64}(undef, length(cumsum(TP.time_steps)))
	pfs_aux_in = Array{Float64}(undef, length(cumsum(TP.time_steps)))
	pfs_aux1 = Array{Float64}(undef, length(cumsum(TP.time_steps)))
	pfs_aux2 = Array{Float64}(undef, length(cumsum(TP.time_steps)))
	for i=1:length(cumsum(TP.time_steps))
		Ffs_aux = GasChromatographySystems.flow_functions(sys_est_Faux(cumsum(TP.time_steps)[i]), p2fun_flow_eval)
		Ffs_aux_o[i] = Ffs_aux[18](cumsum(TP.time_steps)[i]) 
		Ffs_aux_ms[i] = Ffs_aux[19](cumsum(TP.time_steps)[i]) 
		Ffs_aux_5[i] = Ffs_aux[5](cumsum(TP.time_steps)[i]) 
		Ffs_aux_6[i] = Ffs_aux[6](cumsum(TP.time_steps)[i]) 
		pfs_aux = GasChromatographySystems.pressure_functions(sys_est_Faux(cumsum(TP.time_steps)[i]), p2fun_flow_eval)
		pfs_aux_in[i] = pfs_aux[1](cumsum(TP.time_steps)[i]) 
		pfs_aux1[i] = pfs_aux[20](cumsum(TP.time_steps)[i]) 
		pfs_aux2[i] = pfs_aux[21](cumsum(TP.time_steps)[i]) 
	end
end;

# ╔═╡ 014b9528-a883-4599-a5dc-f51db29938ed
begin
	p_fr_aux = Plots.plot(cumsum(TP.time_steps)./60, Ffs_aux_5./Ffs_aux_6; label="", c=1, title="Flow split ratio f", size=(400,400))
end

# ╔═╡ 4b62ec6b-fdb8-4e80-830f-a9787a3891e7
begin
	p_pins_aux = Plots.plot(cumsum(TP.time_steps)./60, pfs_aux_in; label="p_in", c=1)
	Plots.plot!(p_pins_aux, cumsum(TP.time_steps)./60, pfs_aux1; label="p_aux1", c=2, linewidth=2, xlabel="time in min", ylabel="pressure in Pa", title="Inlet pressures")
	Plots.plot!(p_pins_aux, cumsum(TP.time_steps)./60, pfs_aux2; label="p_aux2", c=3, size=(600,400), legend=:topleft)
end

# ╔═╡ cf4f3581-556a-40fa-b698-2f96f60349a7
begin
	p_Fd_aux = Plots.plot(cumsum(TP.time_steps)./60, Ffs_aux_o*60e6; label="F_O", c=1)
	Plots.plot!(p_Fd_aux, cumsum(TP.time_steps)./60, Ffs_aux_ms*60e6; label="F_MS", c=2, xlabel="time in min", ylabel="flow in mL/min", title="Flows to detectors")
	p_Faux_aux = Plots.plot(cumsum(TP.time_steps)./60, Faux_1.(cumsum(TP.time_steps)); label="F_aux1", c=3)
	Plots.plot!(p_Faux_aux, cumsum(TP.time_steps)./60, Faux_2.(cumsum(TP.time_steps)); label="F_aux2", c=4, xlabel="time in min", ylabel="flow in mL/min", title="Auxillary flows")
	Plots.plot(p_Fd_aux, p_Faux_aux, size=(800,400), legend=:top)
end

# ╔═╡ 1344b595-f674-484a-9340-a4918f2e7fc5
gr_Faux(t) = GasChromatographySystems.plot_graph_with_flow(sys_est_Faux(t), p2fun_flow_eval, t; lay=stressl, node_size=40, arrow_size=20, arrow_shift=0.5, elabels_fontsize=12, nlabels_fontsize=12)

# ╔═╡ c6c321a4-fd1e-4096-9367-084eae1b9049
gr_Faux(select_t_aux)

# ╔═╡ b601f315-44a2-44a6-ba6a-e8db806e9b9a
gr_est(t) = GasChromatographySystems.plot_graph_with_flow(sys_est, p2fun_flow_eval, t; lay=stressl, node_size=40, arrow_size=20, arrow_shift=0.5, elabels_fontsize=12, nlabels_fontsize=12)

# ╔═╡ cad4d0b4-3479-4d53-92b1-a704316e8edd
gr_est(select_t_est)

# ╔═╡ 2d213e3d-2ae1-4549-8b06-66ca820a53a0
gr(t) = GasChromatographySystems.plot_graph_with_flow(sys_flow, p2fun_flow_eval, t; lay=stressl, node_size=40, arrow_size=20, arrow_shift=0.5, elabels_fontsize=12, nlabels_fontsize=12)

# ╔═╡ b473ea3a-87ed-46b3-afdd-2315d668a99d
gr(select_t)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
GasChromatographySimulator = "dd82b6e2-56ef-419d-b271-0be268cb65f5"
GasChromatographySystems = "2a17fa6e-c440-4af9-acac-c1be9e3a006f"
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
CSV = "~0.10.14"
DataFrames = "~1.6.1"
GasChromatographySimulator = "~0.4.6"
GasChromatographySystems = "~0.2.1"
Graphs = "~1.11.1"
HypertextLiteral = "~0.9.5"
NLsolve = "~4.5.1"
Plots = "~1.40.4"
PlutoUI = "~0.7.59"
Symbolics = "~5.14.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.0"
manifest_format = "2.0"
project_hash = "3ef806f3a09706c0d6e40f6a5bb0015936a9ffc7"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Preferences", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "d7832de8cf7af26abac741f10372080ac6cb73df"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.34.7"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cde29ddf7e5726c9fb511f340244ea3481267608"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.7.2"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["ArrayInterfaceCore", "Compat", "IfElse", "LinearAlgebra", "SnoopPrecompile", "Static"]
git-tree-sha1 = "dedc16cbdd1d32bead4617d27572f582216ccf23"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "6.0.25"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e5f08b5689b1aad068e01751889f2f615c7db36d"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.29"

[[deps.ArrayInterfaceGPUArrays]]
deps = ["Adapt", "ArrayInterfaceCore", "GPUArraysCore", "LinearAlgebra"]
git-tree-sha1 = "fc114f550b93d4c79632c2ada2924635aabfa5ed"
uuid = "6ba088a2-8465-4c0a-af30-387133b534db"
version = "0.2.2"

[[deps.ArrayInterfaceOffsetArrays]]
deps = ["ArrayInterface", "OffsetArrays", "Static"]
git-tree-sha1 = "3d1a9a01976971063b3930d1aed1d9c4af0817f8"
uuid = "015c0d05-e682-4f19-8f0a-679ce4c54826"
version = "0.1.7"

[[deps.ArrayInterfaceStaticArrays]]
deps = ["Adapt", "ArrayInterface", "ArrayInterfaceCore", "ArrayInterfaceStaticArraysCore", "LinearAlgebra", "Static", "StaticArrays"]
git-tree-sha1 = "f12dc65aef03d0a49650b20b2fdaf184928fd886"
uuid = "b0d46f97-bff5-4637-a19a-dd75974142cd"
version = "0.1.5"

[[deps.ArrayInterfaceStaticArraysCore]]
deps = ["ArrayInterfaceCore", "LinearAlgebra", "StaticArraysCore"]
git-tree-sha1 = "01a9f8e6cfc2bfdd01d333f70b8014a04893103c"
uuid = "dd5226c6-a4d4-4bc7-8575-46859f9c95b9"
version = "0.1.4"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "734d86dfee527f71f61370f410fb3f6c26740d33"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "0.8.19"

[[deps.Arrow]]
deps = ["ArrowTypes", "BitIntegers", "CodecLz4", "CodecZstd", "ConcurrentUtilities", "DataAPI", "Dates", "EnumX", "LoggingExtras", "Mmap", "PooledArrays", "SentinelArrays", "Tables", "TimeZones", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "f8d411d1b45459368567dc51f683ed78a919d795"
uuid = "69666777-d1a9-59fb-9406-91d4454c9d45"
version = "2.7.2"

[[deps.ArrowTypes]]
deps = ["Sockets", "UUIDs"]
git-tree-sha1 = "404265cd8128a2515a81d5eae16de90fdef05101"
uuid = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
version = "2.3.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Automa]]
deps = ["PrecompileTools", "TranscodingStreams"]
git-tree-sha1 = "588e0d680ad1d7201d4c6a804dcb1cd9cba79fbb"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.0.3"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "SnoopPrecompile", "SparseArrays"]
git-tree-sha1 = "6ef8fc1d77b60f41041d59ce61ef9eb41ed97a83"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "0.17.18"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bijections]]
git-tree-sha1 = "c9b163bd832e023571e86d0b90d9de92a9879088"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.6"

[[deps.BitFlags]]
git-tree-sha1 = "2dc09997850d68179b69dafb58ae806167a32b1b"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.8"

[[deps.BitIntegers]]
deps = ["Random"]
git-tree-sha1 = "a55462dfddabc34bc97d3a7403a2ca2802179ae6"
uuid = "c3b6d118-76ef-56ca-8cc7-ebb389d030a1"
version = "0.3.1"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "0c5f81f47bbbcf4aea7b2959135713459170798b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.5"

[[deps.BoundaryValueDiffEq]]
deps = ["BandedMatrices", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "NLsolve", "Reexport", "SciMLBase", "SparseArrays"]
git-tree-sha1 = "ed8e837bfb3d1e3157022c9636ec1c722b637318"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "2.11.0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "585a387a490f1c4bd88be67eea15b93da5e85db7"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.5"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "6c834533dc1fabd820c1db03c839bf97e45a3fab"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.14"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["CRC32c", "Cairo", "Colors", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools"]
git-tree-sha1 = "d69c7593fe9d7d617973adcbe4762028c6899b2c"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.11.11"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a2f1c8c668c8e3cb4cca4e57a8efdb09067bb3fd"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+2"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "71acdbf594aab5bbb2cec89b208c41b4c411e49f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.24.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.ChemicalIdentifiers]]
deps = ["Arrow", "Downloads", "Preferences", "Scratch", "UUIDs", "Unicode"]
git-tree-sha1 = "5f86727335b9896bc6f9fa81964a03174be91dbe"
uuid = "fa4ea961-1416-484e-bda2-883ee1634ba5"
version = "0.1.9"

[[deps.CloseOpenIntervals]]
deps = ["ArrayInterface", "Static"]
git-tree-sha1 = "d61300b9895f129f4bd684b2aff97cf319b6c493"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.11"

[[deps.CodecLz4]]
deps = ["Lz4_jll", "TranscodingStreams"]
git-tree-sha1 = "b8aecef9f90530cf322a8386630ec18485c17991"
uuid = "5ba52731-8f18-5e0d-9241-30f10d1ec561"
version = "0.4.3"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "59939d8a997469ee05c4b4944560a820f9ba0d73"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.4"

[[deps.CodecZstd]]
deps = ["TranscodingStreams", "Zstd_jll"]
git-tree-sha1 = "0d0612d8646ed6157adaceff420b3bacbc2510a9"
uuid = "6b39b394-51ab-5f42-8807-6242bab2b4c2"
version = "0.8.3"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "4b270d6465eb21ae89b732182c20dc165f8bf9f2"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.25.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "b1c55339b7c6c350ee89f2c1604299660525b248"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.15.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+1"

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "6cbbd4d241d7e6579ab354737f4dd95ca43946e1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.1"

[[deps.Configurations]]
deps = ["ExproniconLite", "OrderedCollections", "TOML"]
git-tree-sha1 = "4358750bb58a3caefd5f37a4a0c5bfdbbf075252"
uuid = "5218b696-f38b-4ac9-8b61-a12ec717816d"
version = "0.17.6"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "260fd2400ed2dab602a7c15cf10c1933c59930a2"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.5"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.Cubature]]
deps = ["Cubature_jll"]
git-tree-sha1 = "c3f4b3b38abd7b5c3ccf59adab2568212e7530d3"
uuid = "667455a9-e2ce-5579-9412-b964f529a492"
version = "1.5.1"

[[deps.Cubature_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0fe9efb84e3eb7b14f885a95aaa0ed50c7e839c8"
uuid = "7bc98958-0e37-5d67-a6ac-a3a19030071a"
version = "1.0.5+0"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "04c738083f29f86e62c8afc341f0967d8717bdb8"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.6.1"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelaunayTriangulation]]
deps = ["EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "1755070db557ec2c37df2664c75600298b0c1cfc"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.0.3"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SimpleUnPack"]
git-tree-sha1 = "6b0a1066e81289463c19410fa02b585edb371f20"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.41.1"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Deno_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cd6756e833c377e0ce9cd63fb97689a255f12323"
uuid = "04572ae6-984a-583e-9378-9577a1c2574d"
version = "1.33.4+0"

[[deps.DiffEqBase]]
deps = ["ArrayInterfaceCore", "ChainRulesCore", "DataStructures", "DocStringExtensions", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "MuladdMacro", "Parameters", "PreallocationTools", "Printf", "RecursiveArrayTools", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "Static", "StaticArrays", "Statistics", "Tricks", "ZygoteRules"]
git-tree-sha1 = "a31156664a970730bf85438c3f1c4f338a52dcf0"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.115.0"

    [deps.DiffEqBase.extensions]
    DiffEqBaseDistributionsExt = "Distributions"
    DiffEqBaseGeneralizedGeneratedExt = "GeneralizedGenerated"
    DiffEqBaseMPIExt = "MPI"
    DiffEqBaseMeasurementsExt = "Measurements"
    DiffEqBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    DiffEqBaseReverseDiffExt = "ReverseDiff"
    DiffEqBaseTrackerExt = "Tracker"
    DiffEqBaseUnitfulExt = "Unitful"
    DiffEqBaseZygoteExt = "Zygote"

    [deps.DiffEqBase.weakdeps]
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    GeneralizedGenerated = "6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "LinearAlgebra", "Markdown", "NLsolve", "Parameters", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "63b6be7b396ad395825f3cc48c56b53bfaf7e69d"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.26.1"

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "Requires", "ResettableStacks", "SciMLBase", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "65cbbe1450ced323b4b17228ccd96349d96795a7"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.21.0"

    [deps.DiffEqNoiseProcess.extensions]
    DiffEqNoiseProcessReverseDiffExt = "ReverseDiff"

    [deps.DiffEqNoiseProcess.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DifferentialEquations]]
deps = ["BoundaryValueDiffEq", "DelayDiffEq", "DiffEqBase", "DiffEqCallbacks", "DiffEqNoiseProcess", "JumpProcesses", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEq", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "SteadyStateDiffEq", "StochasticDiffEq", "Sundials"]
git-tree-sha1 = "418dad2cfd7377a474326bc86a23cb645fac6527"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.6.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "66c4c81f259586e8f002eacebc177e1fb06363b0"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.11"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "9c405847cc7ecda2dc921ccf18b47ca150d7317e"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.109"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "51b4b84d33ec5e0955b55ff4b748b99ce2c3faa9"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.6.7"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "30a1848c4f4fc35d1d4bbbd125650f6a11b5bc6c"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.5.7"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "b3f2ff58735b5f024c392fde763f29b057e4b025"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.8"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceGPUArrays", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "Printf", "SnoopPrecompile", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "eb58c1e1417a6580b983069f1491ff82c37def2c"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.23.0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.ExpressionExplorer]]
git-tree-sha1 = "0da78bef32ca71276337442389a3d1962a1ee0da"
uuid = "21656369-7473-754a-2065-74616d696c43"
version = "1.0.2"

[[deps.ExproniconLite]]
git-tree-sha1 = "6091a6fc0f16639f43d7f78fee225ba365712612"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.8"

[[deps.Extents]]
git-tree-sha1 = "2140cd04483da90b2da7f99b2add0750504fc39c"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.2"

[[deps.EzXML]]
deps = ["Printf", "XML2_jll"]
git-tree-sha1 = "380053d61bb9064d6aa4a9777413b40429c79901"
uuid = "8f5d6c58-4d21-5cfd-889c-e3ad7ee6a615"
version = "1.2.0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "LinearAlgebra", "Polyester", "Static", "StrideArraysCore"]
git-tree-sha1 = "4bef892787c972913d4d84e7255400759bb650e5"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.2.4"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastLapackInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c1293a93193f0ae94be7cf338d33e162c39d8788"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "1.2.9"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "82d8afa92ecf4b52d78d869f038ebfb881267322"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.3"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport", "Requires"]
git-tree-sha1 = "919d9412dbf53a2e6fe74af62a73ceed0bce0629"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.8.3"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "9f00e42f8d99fdde64d40c8ea5d14269a2e2c1aa"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.21"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "7072f1e3e5a8be51d525d64f63d3ec1287ff2790"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.11"

[[deps.FiniteDiff]]
deps = ["ArrayInterfaceCore", "LinearAlgebra", "Requires", "Setfield", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "04ed1f0029b6b3af88343e439b995141cb0d0b8d"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.17.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "2493cdfd0740015955a8e46de4ef28f49460d8bc"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.3"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.FuzzyCompletions]]
deps = ["REPL"]
git-tree-sha1 = "40ec72c57559a4473961bbcd12c96bcd4c2aaab4"
uuid = "fb4132e2-a121-4a70-b8a1-d5b831dcdcc2"
version = "0.5.4"

[[deps.GLFW]]
deps = ["GLFW_jll"]
git-tree-sha1 = "35dbc482f0967d8dceaa7ce007d16f9064072166"
uuid = "f7f18e0c-5ee9-5ccd-a5bf-e8befd85ed98"
version = "3.4.1"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "ff38ba61beff76b8f4acad8ab0c97ef73bb670cb"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.9+0"

[[deps.GLMakie]]
deps = ["ColorTypes", "Colors", "FileIO", "FixedPointNumbers", "FreeTypeAbstraction", "GLFW", "GeometryBasics", "LinearAlgebra", "Makie", "Markdown", "MeshIO", "ModernGL", "Observables", "PrecompileTools", "Printf", "ShaderAbstractions", "StaticArrays"]
git-tree-sha1 = "3ef8015aefacb449d183201714a9c32d86019acc"
uuid = "e9467ef8-e4e7-5192-8a1a-b1aee30e663a"
version = "0.9.11"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "2d6ca471a6c7b536127afccfa7564b5b39227fe0"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.5"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "ddda044ca260ee324c5fc07edb6d7cf3f0b9c350"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.5"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "278e5e0f820178e8a26df3184fcb2280717c79b1"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.5+0"

[[deps.GasChromatographySimulator]]
deps = ["CSV", "ChemicalIdentifiers", "DataFrames", "ForwardDiff", "HypertextLiteral", "Integrals", "Interpolations", "OrdinaryDiffEq", "Plots", "PlutoUI", "Reexport", "UrlDownload"]
git-tree-sha1 = "6be668427d09146a156f1881c5a7bb22682a9546"
uuid = "dd82b6e2-56ef-419d-b271-0be268cb65f5"
version = "0.4.6"

[[deps.GasChromatographySystems]]
deps = ["CSV", "CairoMakie", "DataFrames", "DifferentialEquations", "EzXML", "ForwardDiff", "GLMakie", "GasChromatographySimulator", "GraphIO", "GraphMakie", "GraphRecipes", "Graphs", "HypertextLiteral", "Integrals", "IntegralsCubature", "Interpolations", "Intervals", "JSServe", "LsqFit", "NetworkLayout", "OrdinaryDiffEq", "Plots", "Pluto", "PlutoUI", "QuadGK", "Reexport", "SpecialFunctions", "Symbolics", "UrlDownload", "Waveforms"]
git-tree-sha1 = "234d834dcb225b1c620c10f4c0c309d128a3d591"
uuid = "2a17fa6e-c440-4af9-acac-c1be9e3a006f"
version = "0.2.1"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "af49a0851f8113fcfae2ef5027c6d49d0acec39b"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.4"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "801aef8228f7f04972e596b09d4dba481807c913"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.4"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "b62f2b2d76cee0d61a2ef2b3118cd2a3215d3134"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.11"

[[deps.GeometryTypes]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "d796f7be0383b5416cd403420ce0af083b0f9b28"
uuid = "4d00f742-c7ba-57c2-abde-4428a4b178cb"
version = "0.8.5"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "7c82e6a6cd34e9d935e9aa4051b66c6ff3af59ba"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.2+0"

[[deps.GraphIO]]
deps = ["DelimitedFiles", "Graphs", "Requires", "SimpleTraits"]
git-tree-sha1 = "bc5b7609e9f4583f303a0ab2a7016ea318464da0"
uuid = "aa1b3936-2fda-51b9-ab35-c553d3a640a2"
version = "0.7.0"

    [deps.GraphIO.extensions]
    GraphIODOTExt = "ParserCombinator"
    GraphIOGEXFExt = "EzXML"
    GraphIOGMLExt = "ParserCombinator"
    GraphIOGraphMLExt = "EzXML"
    GraphIOLGCompressedExt = "CodecZlib"

    [deps.GraphIO.weakdeps]
    CodecZlib = "944b1d66-785c-5afd-91f1-9de20f533193"
    EzXML = "8f5d6c58-4d21-5cfd-889c-e3ad7ee6a615"
    ParserCombinator = "fae87a5f-d1ad-5cf0-8f61-c941e1580b46"

[[deps.GraphMakie]]
deps = ["DataStructures", "GeometryBasics", "Graphs", "LinearAlgebra", "Makie", "NetworkLayout", "PolynomialRoots", "SimpleTraits", "StaticArrays"]
git-tree-sha1 = "bdddc4afd944ccc67afbd81791d88d944c36f410"
uuid = "1ecd5474-83a3-4783-bb4f-06765db800d2"
version = "0.5.10"

[[deps.GraphRecipes]]
deps = ["AbstractTrees", "GeometryTypes", "Graphs", "InteractiveUtils", "Interpolations", "LinearAlgebra", "NaNMath", "NetworkLayout", "PlotUtils", "RecipesBase", "SparseArrays", "Statistics"]
git-tree-sha1 = "10920601dc51d2231bb3d2111122045efed8def0"
uuid = "bd48cda9-67a9-57be-86fa-5b3c104eda73"
version = "0.5.13"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "334d300809ae0a68ceee3444c6e99ded412bf0b3"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.11.1"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "6f93a83ca11346771a93bbde2bdad2f65b61498f"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.10.2"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "ExprTools", "Logging", "MultivariatePolynomials", "PrecompileTools", "PrettyTables", "Primes", "Printf", "Random", "SIMD", "TimerOutputs"]
git-tree-sha1 = "6b505ef15e55bdc5bb3ddbcfebdff1c9e67081e8"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.5.1"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "6df9cd6ee79fc59feab33f63a1b3c9e95e2461d5"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.2"

[[deps.HCubature]]
deps = ["Combinatorics", "DataStructures", "LinearAlgebra", "QuadGK", "StaticArrays"]
git-tree-sha1 = "10f37537bbd83e52c63abf6393f209dbd641fedc"
uuid = "19dc6840-f33b-545b-b366-655c7e3ffd49"
version = "1.6.0"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "d1d712be3164d61d1fb98e7ce9bcbc6cc06b45ed"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.8"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "eb8fed28f4994600e29beef49744639d985a04b2"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.16"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "2e4520d67b0cef90865b3ef727594d2a58e0e1f8"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.11"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "b2a7eaa169c13f5bcae8131a83bc30eff8f71be0"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.2"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "437abb322a41d527c197fa800455f79d414f0a3c"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.8"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "355e2b974f2e3212a75dfb60519de21361ad3cb7"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.9"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.Integrals]]
deps = ["CommonSolve", "HCubature", "LinearAlgebra", "MonteCarloIntegration", "QuadGK", "Reexport", "Requires", "SciMLBase"]
git-tree-sha1 = "44dd71130586e1776fde9041bd81f2e0db15db10"
uuid = "de52edbc-65ea-441a-8357-d3a637375a31"
version = "3.9.0"

    [deps.Integrals.extensions]
    IntegralsFastGaussQuadratureExt = "FastGaussQuadrature"
    IntegralsForwardDiffExt = "ForwardDiff"
    IntegralsZygoteExt = ["Zygote", "ChainRulesCore"]

    [deps.Integrals.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FastGaussQuadrature = "442a2c76-b920-505d-bb47-c5924d526838"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.IntegralsCubature]]
deps = ["Cubature", "Integrals"]
git-tree-sha1 = "2124372650ff7443d35c475a9cb83d283906a030"
uuid = "c31f79ba-6e32-46d4-a52f-182a8ac42a54"
version = "0.2.3"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be50fe8df3acbffa0274a744f1a99d29c45a57f4"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.1.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"
weakdeps = ["Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalArithmetic]]
deps = ["CRlibm_jll", "MacroTools", "RoundingEmulator"]
git-tree-sha1 = "433b0bb201cd76cb087b017e49244f10394ebe9c"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.22.14"
weakdeps = ["DiffRules", "ForwardDiff", "RecipesBase"]

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "ac0aaa807ed5eaf13f67afe188ebc07e828ff640"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.10.0"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "59545b0a2b27208b0650df0a46b8e3019f85055b"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.4"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "a53ebe394b71470c7f97c2e7e170d51df21b17af"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.7"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JSServe]]
deps = ["Base64", "CodecZlib", "Colors", "Dates", "Deno_jll", "HTTP", "Hyperscript", "LinearAlgebra", "Markdown", "MsgPack", "Observables", "RelocatableFolders", "SHA", "Sockets", "Tables", "ThreadPools", "URIs", "UUIDs", "WidgetsBase"]
git-tree-sha1 = "4bcf2a78f7c80c6f3d594267bb4e7ec03ac9c172"
uuid = "824d6782-a2ef-11e9-3a09-e5662e0c26f9"
version = "2.3.1"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "fa6d0bcff8583bac20f1ffa708c3913ca605c611"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.5"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c84a835e1a09b289ffcd2271bf2a337bbdda6637"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.3+0"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays", "TreeViews", "UnPack"]
git-tree-sha1 = "7af8d30e281ce558807917b69ba16575d05f412b"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.5.1"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "884c2968c2e8e7e6bf5956af88cb46aa745c854b"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.4.1"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "7d703202e65efa1369de1279c162b915e245eed1"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.9"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "267dad6b4b7b5d529c76d40ff48d33f7e94cb834"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.6"

[[deps.KrylovKit]]
deps = ["ChainRulesCore", "GPUArraysCore", "LinearAlgebra", "Printf"]
git-tree-sha1 = "5cebb47f472f086f7dd31fb8e738a8db728f1f84"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.6.1"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70c5da094887fd2cae843b8db33920bac4b6f07d"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.LabelledArrays]]
deps = ["ArrayInterfaceCore", "ArrayInterfaceStaticArrays", "ArrayInterfaceStaticArraysCore", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "0a92979c14dfa71adbf892f0cd073e34b7189197"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.13.0"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "e0b5cd21dc1b44ec6e64f351976f961e6f31d6c4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.3"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "ArrayInterfaceOffsetArrays", "ArrayInterfaceStaticArrays", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static"]
git-tree-sha1 = "0ad6f0c51ce004dadc24a28a0dfecfb568e52242"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.13"

[[deps.LazilyInitializedFields]]
git-tree-sha1 = "8f7f3cabab0fd1800699663533b6d5cb3fc0e612"
uuid = "0e77f7df-68c5-4e49-93ce-4cd80f5598bf"
version = "1.2.2"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LevyArea]]
deps = ["LinearAlgebra", "Random", "SpecialFunctions"]
git-tree-sha1 = "56513a09b8e0ae6485f34401ea9e2f31357958ec"
uuid = "2d8b4e74-eb68-11e8-0fb9-d5eb67b50637"
version = "1.0.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearSolve]]
deps = ["ArrayInterfaceCore", "DocStringExtensions", "FastLapackInterface", "GPUArraysCore", "IterativeSolvers", "KLU", "Krylov", "KrylovKit", "LinearAlgebra", "Preferences", "RecursiveFactorization", "Reexport", "SciMLBase", "Setfield", "SnoopPrecompile", "SparseArrays", "Sparspak", "SuiteSparse", "UnPack"]
git-tree-sha1 = "960da8a80f9882fb52a5a199e944d3b86f0d2b94"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "1.34.1"

    [deps.LinearSolve.extensions]
    LinearSolveHYPRE = "HYPRE"

    [deps.LinearSolve.weakdeps]
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "ArrayInterfaceOffsetArrays", "ArrayInterfaceStaticArrays", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDTypes", "SLEEFPirates", "SnoopPrecompile", "Static", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "9696a80c21a56b937e3fd89e972f8db5db3186e2"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.150"
weakdeps = ["ChainRulesCore", "ForwardDiff", "SpecialFunctions"]

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.LsqFit]]
deps = ["Distributions", "ForwardDiff", "LinearAlgebra", "NLSolversBase", "Printf", "StatsAPI"]
git-tree-sha1 = "40acc20cfb253cf061c1a2a2ea28de85235eeee1"
uuid = "2fda8390-95c7-5789-9bda-21331edee243"
version = "0.15.0"

[[deps.Lz4_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6c26c5e8a4203d43b5497be3ec5d4e0c3cde240a"
uuid = "5ced341a-0733-55b8-9ab6-a4889d929147"
version = "1.9.4+0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "80b2833b56d466b3858d565adcd16a4a05f2089b"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.1.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "InteractiveUtils", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "MakieCore", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun"]
git-tree-sha1 = "4d49c9ee830eec99d3e8de2425ff433ece7cc1bc"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.20.10"

[[deps.MakieCore]]
deps = ["Observables", "REPL"]
git-tree-sha1 = "248b7a4be0f92b497f7a331aed02c1e9a878f46b"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.7.3"

[[deps.Malt]]
deps = ["Distributed", "Logging", "RelocatableFolders", "Serialization", "Sockets"]
git-tree-sha1 = "18cf4151e390fce29ca846b92b06baf9bc6e002e"
uuid = "36869731-bdee-424d-aa32-cab38c994e3b"
version = "1.1.1"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "96ca8a313eb6437db5ffe946c457a401bbb8ce1d"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.5.7"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.MeshIO]]
deps = ["ColorTypes", "FileIO", "GeometryBasics", "Printf"]
git-tree-sha1 = "8c26ab950860dfca6767f2bbd90fdf1e8ddc678b"
uuid = "7269a6da-0436-5bbc-96c2-40638cbb6118"
version = "0.4.11"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "bf17d9cb4f0d2882351dfad030598f64286e5936"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.8"

[[deps.ModernGL]]
deps = ["Libdl"]
git-tree-sha1 = "b76ea40b5c0f45790ae09492712dd326208c28b2"
uuid = "66fc600b-dfda-50eb-8b99-91cfa97b1301"
version = "1.1.7"

[[deps.MonteCarloIntegration]]
deps = ["Distributions", "Random"]
git-tree-sha1 = "3f78ebce296c927d5c854e83cccdb5dcb1845629"
uuid = "4886b29c-78c9-11e9-0a6e-41e1f4161f7b"
version = "0.0.3"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.MsgPack]]
deps = ["Serialization"]
git-tree-sha1 = "f5db02ae992c260e4826fe78c942954b48e1d9c2"
uuid = "99f44e22-a591-53d1-9472-aa23ef4bd671"
version = "1.2.1"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "dad7be0c92b688bf8f24af170825ccedc104b116"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.5"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "a3589efe0005fc4718775d8641b2de9060d23f73"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.4.4"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkLayout]]
deps = ["GeometryBasics", "LinearAlgebra", "Random", "Requires", "StaticArrays"]
git-tree-sha1 = "91bb2fedff8e43793650e7a677ccda6e6e6e166b"
uuid = "46757867-2c16-5918-afeb-47bfcb05e46a"
version = "0.4.6"
weakdeps = ["Graphs"]

    [deps.NetworkLayout.extensions]
    NetworkLayoutGraphsExt = "Graphs"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ArrayInterface", "DiffEqBase", "EnumX", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SnoopPrecompile", "SparseArrays", "SparseDiffTools", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "a6000c813371cd3cd9cbbdf8a356fc3a97138d92"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "1.6.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "e64b4f5ea6b7389f6f046d13d4896a8f9c1ba71e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+2"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a028ee3cb5641cccc4c24e90c36b0a4f7707bdf5"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.14+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "d9b79c4eed437421ac4285148fcadf42e0700e89"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.9.4"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.OrdinaryDiffEq]]
deps = ["Adapt", "ArrayInterface", "ArrayInterfaceCore", "ArrayInterfaceGPUArrays", "ArrayInterfaceStaticArrays", "ArrayInterfaceStaticArraysCore", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "LinearAlgebra", "LinearSolve", "Logging", "LoopVectorization", "MacroTools", "MuladdMacro", "NLsolve", "NonlinearSolve", "Polyester", "PreallocationTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLNLSolve", "SimpleNonlinearSolve", "SnoopPrecompile", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "0e9df1f4a5fee8c9a6e1349f363eae68b35bbe00"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.40.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "67186a2bc9a90f9f85ff3cc8277868961fb57cbd"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.3"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "ec3edfe723df33528e085e632414499f26650501"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.0"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cb5a2ab6763464ae0f19c86c56c63d4a2b0f5bda"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.52.2+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "6e55c6841ce3411ccb3457ee52fc48cb698d6fb0"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.2.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "442e1e7ac27dd5ff8825c3fa62fbd1e86397974b"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.4"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Pluto]]
deps = ["Base64", "Configurations", "Dates", "Downloads", "ExpressionExplorer", "FileWatching", "FuzzyCompletions", "HTTP", "HypertextLiteral", "InteractiveUtils", "Logging", "LoggingExtras", "MIMEs", "Malt", "Markdown", "MsgPack", "Pkg", "PlutoDependencyExplorer", "PrecompileSignatures", "PrecompileTools", "REPL", "RegistryInstances", "RelocatableFolders", "Scratch", "Sockets", "TOML", "Tables", "URIs", "UUIDs"]
git-tree-sha1 = "7074b3a8339fadaf8524a9252ae7565b85f648f1"
uuid = "c3e4b0f8-55cb-11ea-2926-15256bba5781"
version = "0.19.42"

[[deps.PlutoDependencyExplorer]]
deps = ["ExpressionExplorer", "InteractiveUtils", "Markdown"]
git-tree-sha1 = "4bc5284f77d731196d3e97f23abb732ad6f2a6e4"
uuid = "72656b73-756c-7461-726b-72656b6b696b"
version = "1.0.4"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

[[deps.PoissonRandom]]
deps = ["Random"]
git-tree-sha1 = "a0f1159c33f846aa77c3f30ebbc69795e5327152"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.4"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "e8e0fabcff4df8686c4267503887202a783d498e"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.2"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "240d7170f5ffdb285f9427b92333c3463bf65bf6"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.1"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PolynomialRoots]]
git-tree-sha1 = "5f807b5345093487f733e520a1b7395ee9324825"
uuid = "3a141323-8675-5d76-9d11-e1df1406c778"
version = "1.0.0"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ForwardDiff", "Requires"]
git-tree-sha1 = "2c7658dd593e3adc118b00429e1048829f1abb8c"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.11"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.PrecompileSignatures]]
git-tree-sha1 = "18ef344185f25ee9d51d80e179f8dad33dc48eb1"
uuid = "91cefc8d-f054-46dc-8f8c-26e11d7c5411"
version = "3.0.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "66b20dd35966a748321d3b2537c4584cf40387c7"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.3.2"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "cb420f77dc474d23ee47ca8d14c90810cafe69e7"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.6"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "763a8ceb07833dd51bb9e3bbca372de32c0605ad"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.0"

[[deps.PtrArrays]]
git-tree-sha1 = "f011fbb92c4d401059b2212c05c0601b70f8b759"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "37b7bb7aabf9a085e0044307e1717436117f2b3b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.3+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9b23c31e76e333e6fb4c1595ae6afa74966a729e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.4"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "4743b43e5a9c4a2ede372de7061eed81795b12e7"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.7.0"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "b8a399e95663485820000f26b6a43c794e166a49"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.4"

[[deps.RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceStaticArraysCore", "ChainRulesCore", "DocStringExtensions", "FillArrays", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "Tables", "ZygoteRules"]
git-tree-sha1 = "a5ce741acddc02f0d4fc6505463ca89697d7fb23"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.32.3"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "PrecompileTools", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "c04dacfc546591d43c39dc529c922d6a06a5a694"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.22"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegistryInstances]]
deps = ["LazilyInitializedFields", "Pkg", "TOML", "Tar"]
git-tree-sha1 = "ffd19052caf598b8653b99404058fce14828be51"
uuid = "2792f1a3-b283-48e8-9a74-f99dce5104f3"
version = "0.1.0"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "256eeeec186fa7f26f2801732774ccf277f05db9"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.1.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d483cd324ce5cf5d61b77930f0bbd6cb61927d21"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.2+0"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "04c968137612c4a5629fa531334bb81ad5680f00"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.13"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "2803cab51702db743f3fda07dd1745aadfbf43bd"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.5.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "3aac6d68c5e57449f5b9b865c9ba50ac2970c4cf"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.42"

[[deps.SciMLBase]]
deps = ["ArrayInterfaceCore", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Preferences", "RecipesBase", "RecursiveArrayTools", "RuntimeGeneratedFunctions", "StaticArraysCore", "Statistics", "Tables"]
git-tree-sha1 = "fe89a8113ea445bcff9ee570077830674babb534"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.81.0"

[[deps.SciMLNLSolve]]
deps = ["DiffEqBase", "LineSearches", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "a8eb97c56cac50c21096582afb2a0110784dc36e"
uuid = "e9a6253c-8580-4d32-9898-8661bb511710"
version = "0.1.6"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "90b4f68892337554d31cdcdbe19e48989f26c7e6"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.3"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "79123bc60c5507f035e6d1d9e563bb2971954ec8"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.4.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ArrayInterface", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "Reexport", "Requires", "SciMLBase", "SnoopPrecompile", "StaticArraysCore"]
git-tree-sha1 = "248cf0cd5648c0c4bcefd7a5c33aa4738154043d"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "0.1.12"

    [deps.SimpleNonlinearSolve.extensions]
    SimpleBatchedNonlinearSolveExt = "NNlib"

    [deps.SimpleNonlinearSolve.weakdeps]
    NNlib = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SparseDiffTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceStaticArrays", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "4245283bee733122a9cb4545748d64e0c63337c0"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.30.0"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "342cf4b449c299d8d1ceaf00b7a49f4fbc7940e7"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.9"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "d2fdac9ff3906e27f7a618d47b676941baa6c80c"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.10"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "6e00379a24597be4ae1ee6b2d882e15392040132"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.5"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "cef0472124fab0695b58ca35a77c6fb942fdab8a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.1"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.SteadyStateDiffEq]]
deps = ["DiffEqBase", "DiffEqCallbacks", "EnumX", "LinearAlgebra", "Markdown", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "a0cdb38cc0207213a95aa6caa40b3c94ab323f05"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "1.12.0"

[[deps.StochasticDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqNoiseProcess", "DocStringExtensions", "FillArrays", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "c6b4b802d4d830e0e958f5f2098d8dea0a935f4b"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.58.0"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "SIMDTypes", "Static", "ThreadingUtilities"]
git-tree-sha1 = "8114ba9c3694827838d45ea3c9f6b9ccb4182cf2"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.4.7"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "f4dc295e983502292c4c3f951dbb4e985e35b3be"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.18"
weakdeps = ["Adapt", "GPUArraysCore", "SparseArrays", "StaticArrays"]

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = "GPUArraysCore"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "Reexport", "SciMLBase", "SnoopPrecompile", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "c033830e3c6fb4260243fc907b1e7e93421e7ae8"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.15.1"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "SuiteSparse_jll", "libblastrampoline_jll"]
git-tree-sha1 = "ba4d38faeb62de7ef47155ed321dce40a549c305"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.2+0"

[[deps.SymbolicIndexingInterface]]
deps = ["MacroTools", "RuntimeGeneratedFunctions"]
git-tree-sha1 = "f7b1fc9fc2bc938436b7684c243be7d317919056"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.11"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TermInterface", "TimerOutputs", "Unityper"]
git-tree-sha1 = "20339c0dd70abdb73494955df4fcd9e9ccaff861"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.6.0"

[[deps.Symbolics]]
deps = ["ArrayInterface", "Bijections", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "PrecompileTools", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "SymbolicUtils"]
git-tree-sha1 = "8d28ebc206dec9e250e21b9502a2662265897650"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.14.1"

    [deps.Symbolics.extensions]
    SymbolicsSymPyExt = "SymPy"

    [deps.Symbolics.weakdeps]
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TZJData]]
deps = ["Artifacts"]
git-tree-sha1 = "1607ad46cf8d642aa779a1d45af1c8620dbf6915"
uuid = "dc5dba14-91b3-4cab-a142-028a31da12f7"
version = "1.2.0+2024a"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "cb76cf677714c095e535e3501ac7954732aeea2d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TermInterface]]
git-tree-sha1 = "6f0cee95e74d1f6891ba6b35b8b219fd3d11b567"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "0.4.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadPools]]
deps = ["Printf", "RecipesBase", "Statistics"]
git-tree-sha1 = "50cb5f85d5646bc1422aa0238aa5bfca99ca9ae7"
uuid = "b189fb0b-2eb5-4ed4-bc0c-d34c51242431"
version = "2.1.1"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "eda08f7e9818eb53661b3deb74e3159460dfbc27"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.2"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "bc7fd5c91041f44636b2c134041f7e5263ce58ae"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.10.0"

[[deps.TimeZones]]
deps = ["Dates", "Downloads", "InlineStrings", "Mocking", "Printf", "Scratch", "TZJData", "Unicode", "p7zip_jll"]
git-tree-sha1 = "a6ae8d7a27940c33624f8c7bde5528de21ba730d"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.17.0"
weakdeps = ["RecipesBase"]

    [deps.TimeZones.extensions]
    TimeZonesRecipesBaseExt = "RecipesBase"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "5a13ae8a41237cff5ecf34f73eb1b8f42fff6531"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.24"

[[deps.TranscodingStreams]]
git-tree-sha1 = "a947ea21087caba0a798c5e494d0bb78e3a1a3a0"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.9"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "7ee8ed8904e7dd5d31bb46294ef5644d9e2e44e4"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.21"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "dd260903fdabea27d9b6021689b3cd5401a57748"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.20.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "25008b734a03736c41e2a7dc314ecb95bd6bbdb0"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.6"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.UrlDownload]]
deps = ["HTTP", "ProgressMeter"]
git-tree-sha1 = "758aeefcf0bfdd04cd88e6e3bc9feb9e8c7d6d70"
uuid = "856ac37a-3032-4c1c-9122-f86d88358c8b"
version = "1.0.1"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static"]
git-tree-sha1 = "4c59c2df8d2676c4691a39fa70495a6db0c5d290"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.58"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Waveforms]]
git-tree-sha1 = "28ed144e2362cf22694fc0edc1a6309eb4a5e903"
uuid = "cb13b1c6-351e-5134-b3ad-d6a530956a82"
version = "0.1.2"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WidgetsBase]]
deps = ["Observables"]
git-tree-sha1 = "30a1d631eb06e8c868c559599f915a62d55c2601"
uuid = "eead4739-05f7-45a1-878c-cee36b57321c"
version = "0.1.4"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "52ff2af32e591541550bd753c0da8b9bc92bb9d9"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.7+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e678132f07ddb5bfa46857f0d7620fb9be675d3b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+0"

[[deps.ZygoteRules]]
deps = ["ChainRulesCore", "MacroTools"]
git-tree-sha1 = "27798139afc0a2afa7b1824c206d5e87ea587a00"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.5"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─f27920a2-00d5-11ef-0209-79647fcea22b
# ╟─64c1c01a-8a09-40b1-bfbd-a7668cde9b0b
# ╟─13569378-fba2-4551-8435-aaf2f13bfc66
# ╟─a8bbc6cc-3bcb-43c8-9b00-f2b3e814c136
# ╠═a6a03024-d77d-499a-9617-5c142dad4db6
# ╟─ec81a4e5-8916-46f5-bf1d-6a4781963bbb
# ╟─7345afa9-27c7-472f-b8dc-a55b2fb34c78
# ╟─38afe7e3-84cc-4304-bcc2-5aabdd71e27d
# ╟─dc2e6dfb-7f5a-462f-8fe8-912afa76478f
# ╠═5bcfa4d6-6c4b-44ac-9473-f12d410b0ac7
# ╟─e2df094e-ab54-43a6-8a9d-778031a576c0
# ╠═71902736-dea3-4c85-b067-e5d8d7c9ccb2
# ╟─b0451d2d-5a1b-4027-affb-f6531d11d1a5
# ╠═6c701a8d-d6eb-4303-b000-ce864f222658
# ╟─b0750420-37d2-4e7f-9c1e-584a4890db00
# ╟─bbeb3406-ab4d-45c4-9a87-714c6cdee3bd
# ╟─be454bf6-e220-410f-b7d1-d7e892b7b641
# ╟─39ea740d-a6de-4e68-aaa1-6517935c3a93
# ╟─fd6973a3-68df-49ea-9988-8fc2f26a25e7
# ╟─64146877-f6b4-4d1a-aeb9-147809b756ca
# ╟─f0a02e52-a0c1-40e9-b479-b8162c63f445
# ╟─2a413dc2-9c4f-4608-8f40-034fba017163
# ╟─2fdf479e-1c3b-46a6-9994-6c5ea465ff10
# ╟─dc33f080-4ccf-42ad-b8a8-b192da84f1ec
# ╟─25687b22-dcc3-4571-af9e-e9f260133c3f
# ╟─e41477bc-247e-454e-864d-68a223ad59a2
# ╟─b473ea3a-87ed-46b3-afdd-2315d668a99d
# ╠═e99c055a-0cac-40a1-8969-3fb5e8ba5825
# ╟─795fee94-00d6-442b-aa32-19a296e9c182
# ╟─a72bdfbf-4bc4-4d25-991b-ee55c29ae57f
# ╟─3f570097-406c-4d84-947f-419678aa279b
# ╟─dfb62248-29ba-4cf3-8c28-0f4ff47e5d46
# ╟─aa9810b9-3d2c-4c35-8267-c93998414b20
# ╠═32cb7766-20b2-45c4-9909-cf97153aa050
# ╠═64055873-b815-43fe-adcc-fea47a82768d
# ╟─deceb6b7-08f7-4c45-9663-1062b58fefb5
# ╟─5d4d3407-0136-44ba-8130-774cd09aa58e
# ╟─e72d0302-b19a-4792-8aa0-3a2df88061a2
# ╟─171f9c3f-d986-479b-ad21-46f890187cd8
# ╟─835ea802-4256-479f-9267-0aaefd33d1b6
# ╟─7b42338e-85ae-4418-aa39-34ddcdd80aaa
# ╟─0c670989-8354-4fc6-8578-8032c79ac085
# ╟─0e476729-f41a-4d70-aeed-04632354dd31
# ╟─cad4d0b4-3479-4d53-92b1-a704316e8edd
# ╟─8674b5b0-37b0-4447-b1bc-4d7bf2cd8c33
# ╟─05b6291f-d6f6-40bf-a090-97c531c10c49
# ╟─fafa9077-92ee-4d02-80e6-9c4db34388b8
# ╟─3d49a31c-ecd9-4eb0-8790-f2543706ead5
# ╟─92b0763f-24cd-4bed-a582-76cea406079e
# ╟─3b92560e-11d1-42a2-8220-5dac046025d6
# ╟─5e07ad07-f3de-45ac-b98c-91ed1ddce44c
# ╠═30e2fc65-642b-4210-bcac-2fdea6d1946f
# ╠═00c438cc-58c5-4e30-bd07-41f928d88634
# ╟─afb33442-9cd1-4378-ac25-f67533da8206
# ╟─faaa3135-0135-47f2-ab1b-9f92cdcb309b
# ╟─3c0680df-79b7-4619-8625-d5cfd76e5779
# ╟─5abb07b4-55d3-4632-8a7b-ca0233454dde
# ╟─c6c321a4-fd1e-4096-9367-084eae1b9049
# ╟─cf4f3581-556a-40fa-b698-2f96f60349a7
# ╟─014b9528-a883-4599-a5dc-f51db29938ed
# ╟─ba9f11c7-7857-4008-bb03-614e488e731f
# ╟─4b62ec6b-fdb8-4e80-830f-a9787a3891e7
# ╟─d3759b64-ebb1-4daa-ba6a-302630d87f85
# ╟─fbd5ee2f-092f-429b-8797-660d1107308f
# ╟─b9a2552d-d8e9-4d94-b0a5-097e4b03a843
# ╟─ed7ea63d-66f7-4fe2-be14-6f916e548d10
# ╟─9eaf6a93-a5d4-4cda-a281-d9e95b8df697
# ╟─385d3632-d091-4fde-8045-545733ce6fee
# ╟─cf3c98f1-3cb7-41ef-b7fc-eabded35f1e1
# ╟─fbd3ff22-c286-4c9e-b60b-3f6373ca40e5
# ╟─2138443d-5a25-4e4d-92ee-9f19f931ea18
# ╟─ea46e236-dbf9-41f8-a661-ffe3c76affec
# ╟─89499fa9-ca50-4fcf-9f06-7f3eca7c5194
# ╟─4d9d9e40-0e9b-4044-b764-3fabc18b0970
# ╟─26064dc0-fe8b-4cbb-9f95-8a8ef3d69b7a
# ╟─ed57c1f5-76de-48a5-bd1e-3e7172b2ad48
# ╟─2dc31d04-abdf-478d-a7cc-eea66b00f2ed
# ╟─7cc3428c-6b3e-4acd-9c72-2a0b19270c4c
# ╟─e82d7276-2601-4ce9-b0f4-d15bbe8214b3
# ╟─1344b595-f674-484a-9340-a4918f2e7fc5
# ╟─b601f315-44a2-44a6-ba6a-e8db806e9b9a
# ╟─2d213e3d-2ae1-4549-8b06-66ca820a53a0
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
