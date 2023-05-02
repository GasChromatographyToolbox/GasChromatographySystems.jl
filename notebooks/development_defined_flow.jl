### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ 7722a5a0-e41d-11ed-3f24-b54f45f9a09b
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
	include(joinpath(dirname(pwd()), "src", "GasChromatographySystems.jl"))
	TableOfContents()
end

# ╔═╡ 1acc69fb-ea06-470a-a6c4-b336c64e7c62
md"""
# Development of the simulation of GC systems for defined flows along column sections
"""

# ╔═╡ 7195fdd4-5924-4033-b794-2a491336cce0
md"""
## Modified `ModuleColumn`
"""

# ╔═╡ 724ac332-c561-416e-83c7-bb5a339f26e9
begin
	struct ModuleColumn<:GasChromatographySystems.AbstractModule
		# Module
		# GC column, gradients are possible
		name::String
		length::Float64
		diameter#::Fd # Function
		a_diameter::Array{Float64,1} # Parameters of diameter function, just for information
		film_thickness#::Fdf # Function
		a_film_thickness::Array{Float64,1} # Parameters of film_thickness function, just for information
		stationary_phase::String
		temperature # a number (constant temperature) or a TemperatureProgram structure
		flow # an number (constant flow) or a Function
	end

	function ModuleColumn(name, L, d, df, sp, tp)
		# function to construct the Column structure
		# for the case of constant diameter and constant film thickness
		# and undefined flow
		col = ModuleColumn(name, L, d, [d], df, [df], sp, tp, NaN)
		return col
	end

	function ModuleColumn(name, L, d, df, sp, tp, flow)
		# function to construct the Column structure
		# for the case of constant diameter and constant film thickness
		col = ModuleColumn(name, L, d, [d], df, [df], sp, tp, flow)
		return col
	end
end

# ╔═╡ 14a9675d-7cb0-42aa-9f8e-352770dd7b23
function ModuleColumn_(name, L, d, df, sp, tp)
		# function to construct the Column structure
		# for the case of constant diameter and constant film thickness
		# and undefined flow
		col = ModuleColumn(name, L, d, [d], df, [df], sp, tp, NaN)
		return col
	end

# ╔═╡ be577666-dc2a-40d1-9a85-643928778949
function ModuleColumn_(name, L, d, df, sp, tp, flow)
		# function to construct the Column structure
		# for the case of constant diameter and constant film thickness
		col = ModuleColumn(name, L, d, [d], df, [df], sp, tp, flow)
		return col
	end

# ╔═╡ 8d101c8e-a7df-46b7-81e0-aaa1de6e3b23
md"""
## Definition of the system
"""

# ╔═╡ c5bca10a-9f4d-4f6b-a162-9e5c58d96a9d
default_TP = GasChromatographySystems.TemperatureProgram([0.0, 1800.0], [40.0, 240.0])

# ╔═╡ 34723a63-e7fa-4d26-a047-843ccf124804
md"""
## Examples
"""

# ╔═╡ fcb30907-6fad-4708-a8f4-f8fab3303b5d
md"""
### Split
"""

# ╔═╡ ec5fd728-8ae9-4193-9013-3e3d94a1fc71
function update_system(sys)
	new_timesteps, new_pressuresteps, new_temperaturesteps, index_module_tempprog = GasChromatographySystems.match_programs(sys)
	new_pp = Array{GasChromatographySystems.PressurePoint}(undef, nv(sys.g))
	for i=1:nv(sys.g)
		new_pp[i] = GasChromatographySystems.PressurePoint(sys.pressurepoints[i].name, new_timesteps, new_pressuresteps[i])
	end
	new_modules = Array{GasChromatographySystems.AbstractModule}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		if typeof(sys.modules[i].temperature) <: Number
			new_modules[i] = sys.modules[i]
		elseif typeof(sys.modules[i].temperature) <: GasChromatographySystems.TemperatureProgram
			# add/modify for gradient
			ii = findfirst(index_module_tempprog.==i)
			new_tp = GasChromatographySystems.TemperatureProgram(new_timesteps, new_temperaturesteps[ii])
			new_modules[i] = ModuleColumn(sys.modules[i].name, sys.modules[i].length, sys.modules[i].diameter, sys.modules[i].film_thickness, sys.modules[i].stationary_phase, new_tp, sys.modules[i].flow)
		end
	end
	new_sys = GasChromatographySystems.System(sys.g, new_pp, new_modules, sys.options)
    return new_sys
end

# ╔═╡ 9736b84a-16a4-43df-a735-a5bc4c98c670
begin
	g_split = SimpleDiGraph(4)
	add_edge!(g_split, 1, 2) # Inj -> GC column -> Split point
	add_edge!(g_split, 2, 3) # Split point -> TL column -> Det 1
	add_edge!(g_split, 2, 4) # Split point -> TL column -> Det 2
	# pressure points
	pp_split = Array{GasChromatographySystems.PressurePoint}(undef, nv(g_split))
	pp_split[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 1800.0], [NaN, NaN]) # inlet 
	pp_split[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 1800.0], [NaN, NaN]) # 
	pp_split[3] = GasChromatographySystems.PressurePoint("p₃", [0.0, 1800.0], [eps(Float64), eps(Float64)]) # outlet 1 
	pp_split[4] = GasChromatographySystems.PressurePoint("p₄", [0.0, 1800.0], [101300.0, 101300.0]) # outlet 2
	# modules
	modules_split = Array{GasChromatographySystems.AbstractModule}(undef, ne(g_split))
	modules_split[1] = ModuleColumn_("GC column", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP, 1.0/60e6)
	modules_split[2] = ModuleColumn_("TL column 1", 0.5, 0.1e-3, 0.0e-6, "", 300.0, NaN)
	modules_split[3] = ModuleColumn_("TL column 2", 1.0, 0.15e-3, 0.15e-6, "SLB5ms", 300.0, NaN)
	# system
	sys_split = update_system(GasChromatographySystems.System(g_split, pp_split, modules_split, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ d7e5df4f-d23d-4ab0-a6d2-830e5e24251e
GasChromatographySystems.plot_graph(sys_split; color=:green)

# ╔═╡ 126ab6c1-632c-4d27-a457-38b8759ae4a2
md"""
### Series
"""

# ╔═╡ f9acd758-c86a-4ecf-bf1f-f8a26a4b95de
begin
	g_series = SimpleDiGraph(4)
	add_edge!(g_series, 1, 2) # Inj -> GC column -> Split point
	add_edge!(g_series, 2, 3) # Split point -> TL column -> Det 1
	add_edge!(g_series, 3, 4) # Split point -> TL column -> Det 2
	# pressure points
	pp_series = Array{GasChromatographySystems.PressurePoint}(undef, nv(g_series))
	pp_series[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 1800.0], [NaN, NaN]) # inlet 
	pp_series[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 1800.0], [NaN, NaN]) # 
	pp_series[3] = GasChromatographySystems.PressurePoint("p₃", [0.0, 1800.0], [NaN, NaN]) # outlet 1 
	pp_series[4] = GasChromatographySystems.PressurePoint("p₄", [0.0, 1800.0], [101300.0, 101300.0]) # outlet 2
	# modules
	modules_series = Array{GasChromatographySystems.AbstractModule}(undef, ne(g_series))
	modules_series[1] = ModuleColumn_("GC column", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP, 1.0/60e6)
	modules_series[2] = ModuleColumn_("TL column 1", 0.5, 0.1e-3, 0.0e-6, "", 300.0, NaN)
	modules_series[3] = ModuleColumn_("TL column 2", 1.0, 0.15e-3, 0.15e-6, "SLB5ms", 300.0, NaN)
	# system
	sys_series = update_system(GasChromatographySystems.System(g_series, pp_series, modules_series, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ 963c47a9-e6e9-4f2a-817d-e29c55b974e4
GasChromatographySystems.plot_graph(sys_series; color=:green)

# ╔═╡ 9c337183-adad-463b-af88-9d0844dee46e
md"""
### Complex System
"""

# ╔═╡ 6fda756c-3a32-4a99-b9df-1e40c750f8ec
begin
	g_cs = SimpleDiGraph(8)
	add_edge!(g_cs, 1, 2) # TL1
	add_edge!(g_cs, 2, 3) # GC1
	add_edge!(g_cs, 2, 5) # TL3
	add_edge!(g_cs, 3, 4) # TL2
	add_edge!(g_cs, 5, 6) # GC2
	add_edge!(g_cs, 5, 7) # GC3
	add_edge!(g_cs, 6, 8) # TL4
	add_edge!(g_cs, 7, 8) # TL5
	# pressure points
	pp_cs = Array{GasChromatographySystems.PressurePoint}(undef, nv(g_cs))
	pp_cs[1] = GasChromatographySystems.PressurePoint("p₁", [0.0, 1800.0], [NaN, NaN]) # inlet 
	pp_cs[2] = GasChromatographySystems.PressurePoint("p₂", [0.0, 1800.0], [NaN, NaN]) # 
	pp_cs[3] = GasChromatographySystems.PressurePoint("p₃", [0.0, 1800.0], [NaN, NaN]) #  
	pp_cs[4] = GasChromatographySystems.PressurePoint("p₄", [0.0, 1800.0], [101300.0, 101300.0]) # outlet 1
	pp_cs[5] = GasChromatographySystems.PressurePoint("p5", [0.0, 1800.0], [NaN, NaN]) # 
	pp_cs[6] = GasChromatographySystems.PressurePoint("p6", [0.0, 1800.0], [NaN, NaN]) # 
	pp_cs[7] = GasChromatographySystems.PressurePoint("p7", [0.0, 1800.0], [NaN, NaN]) # 
	pp_cs[8] = GasChromatographySystems.PressurePoint("p8", [0.0, 1800.0], [eps(Float64), eps(Float64)]) # outlet 2
	# modules
	modules_cs = Array{GasChromatographySystems.AbstractModule}(undef, ne(g_cs))
	modules_cs[1] = ModuleColumn_("TL1", 0.30, 0.25e-3, 0.25e-6, "", 300, NaN)
	modules_cs[2] = ModuleColumn_("GC1", 30.0, 0.25e-3, 0.25e-6, "SLB5ms", default_TP, NaN)
	modules_cs[3] = ModuleColumn_("TL3", 1.0, 0.15e-3, 0.15e-6, "", 300.0, NaN)#1.5/60e6)
	modules_cs[4] = ModuleColumn_("TL2", 1.0, 0.15e-3, 0.15e-6, "SLB5ms", 300.0, 1.0/60e6)
	modules_cs[5] = ModuleColumn_("GC2", 10.0, 0.1e-3, 0.1e-6, "SLB5ms", default_TP, NaN)
	modules_cs[6] = ModuleColumn_("GC3", 10.0, 0.15e-3, 0.15e-6, "Wax", default_TP, NaN)
	modules_cs[7] = ModuleColumn_("TL4", 1.0, 0.1e-3, 0.1e-6, "SLB5ms", 300.0, NaN)
	modules_cs[8] = ModuleColumn_("TL5", 1.0, 0.15e-3, 0.15e-6, "Wax", 300.0, NaN)
	# system
	sys_cs = update_system(GasChromatographySystems.System(g_cs, pp_cs, modules_cs, GasChromatographySystems.Options(ng=true)))
end

# ╔═╡ 9dc540a5-dd43-4b75-b174-4fb2453d3694
GasChromatographySystems.plot_graph(sys_cs; color=:green)

# ╔═╡ bf770d0a-7347-47a1-95f8-8669efde4716
md"""
## Pressure and flow calculation
"""

# ╔═╡ 9a3a7286-95eb-424b-8617-8a1f6058d7c8
md"""
- first construct the flow balance equations only using the flows over the edges
- replace the unkown flows with the equation using pressures ``p_{in}``, ``p_{out}`` and the restrictions ``κ``
"""

# ╔═╡ 6f30dc72-4928-475e-86b4-261256f77ab5
md"""
- ``P_1`` is not explicitly included in the flow balance equations
- but it can be calculated once the other unknown pressures ``P_2`` and ``P_3`` are determinated, from the known ``F_1``
- -> add this as additional equation!
"""

# ╔═╡ b7b5069c-ad71-4dcc-b460-2cdd057c8da7
GasChromatographySystems.flow_balance(sys_split)

# ╔═╡ 2a32dfdb-d805-4437-9658-6fba5365883c
g_split

# ╔═╡ eaa0b056-4a6b-4f6c-9ff8-a0e4d0b8ddb8
E = collect(edges(g_split))

# ╔═╡ 5f4eccb4-2a86-4705-b8cf-245dbbf78c8c
srcE = src.(E)

# ╔═╡ 272eb44b-9103-49d8-b1b4-811a0122ae6e
dstE = dst.(E)

# ╔═╡ 9d60fb4d-d013-45af-890a-3474a521f283
inner_V = GasChromatographySystems.inner_vertices(g_split)

# ╔═╡ ed44d3b8-b8ab-49a7-9fe3-2a54ac6de92b
i_src = findall(srcE.==inner_V)

# ╔═╡ be251bb9-a4a8-4ca7-a7ca-9a4cd3c3409d
i_dst = findall(dstE.==inner_V)

# ╔═╡ 870b5176-94f8-4915-83a1-77c4dbcdbb0a
@variables A, P²[1:nv(sys_split.g)], κ[1:ne(sys_split.g)], F[1:ne(sys_split.g)]

# ╔═╡ 814655cb-ff2b-4042-a1ee-09cb4dfc173d
F[1]

# ╔═╡ b3054fae-a2ea-47dc-9bb6-ce526e9e9b39
begin
	balance = 0
	for j=1:length(i_dst) # ingoing flows
		balance = balance + F[i_dst[j]]
	end
	for j=1:length(i_src) # outgoing flows
		balance = balance - F[i_src[j]]
	end
end

# ╔═╡ 1942d1cc-c573-44f2-8055-fc7140762a90
balance

# ╔═╡ 38be367b-3d58-46a7-8062-49771bd52952
function flow_balance(g, i_n, F)
	#@variables P²[1:nv(g)], κ[1:ne(g)]
	E = collect(edges(g)) # all edges
	srcE = src.(E) # source nodes of the edges
	dstE = dst.(E) # destination nodes of the edges
	# find edges, where node `i_n` is the source
	i_src = findall(srcE.==i_n)
	# find edges, where node `i_n` is the destination
	i_dst = findall(dstE.==i_n)
	balance = 0
	for j=1:length(i_dst) # ingoing flows
		balance = balance + F[i_dst[j]]
	end
	for j=1:length(i_src) # outgoing flows
		balance = balance - F[i_src[j]]
	end
	return balance
end

# ╔═╡ 7be3d2e8-31cb-4a8f-80a0-18d75f9fb6b6
# first construct the flow balance equations only using the flows over the edges
function flow_balance(sys)
	@variables F[1:ne(sys.g)]
	inner_V = GasChromatographySystems.inner_vertices(sys.g) # one flow balance equation for every inner node
	bal_eq = Array{Symbolics.Equation}(undef, length(inner_V))
	for i=1:length(inner_V)
		bal_eq[i] = flow_balance(sys.g, inner_V[i], F) ~ 0
	end
	return bal_eq
end

# ╔═╡ 1fef6d42-b51f-4a4c-b09e-1750351b66eb
flow_balance(g_split, inner_V, F)

# ╔═╡ 82fbba32-886a-4f6f-a855-7d0947c58afe
fbal = flow_balance(sys_split)

# ╔═╡ 3a886889-6655-4b6d-b2f5-a50b2053da36
GasChromatographySystems.flow_balance(sys_split)

# ╔═╡ bf8c00b0-10c5-4de1-a4ab-27105e2b06b1
function unknown_F(sys)
	i_unknown_flows = Int[]
	for i=1:ne(sys.g)
		if isnan(sys.modules[i].flow)
			push!(i_unknown_flows, i)
		end
	end
	return i_unknown_flows
end

# ╔═╡ fa9eca1c-682b-495c-afc4-16bee0d58ce9
unknown_F(sys_split)

# ╔═╡ b28571c6-10c4-430a-8e62-c39571b0a38b
GasChromatographySystems.unknown_p(sys_split)

# ╔═╡ 0b2a67f3-bcdb-4c51-8322-42b51d43b71c
function substitute_unknown_flows(sys)
	@variables A, P²[1:nv(sys.g)], κ[1:ne(sys.g)], F[1:ne(sys.g)]
	E = collect(edges(sys.g)) # all edges
	srcE = src.(E) # source nodes of the edges
	dstE = dst.(E) # destination nodes of the edges
	i_unknown_F = unknown_F(sys) # indices of the modules with an unknown flow
	# create dictionary for the substitution of the unknown flows
	sub_dict = Dict()
	for i=1:length(i_unknown_F)
		j = i_unknown_F[i]
		sub_dict = merge(sub_dict, Dict(F[j] => A*(P²[srcE[j]]-P²[dstE[j]])/κ[j]))
	end
	# index of the known flows
	i_known_F = collect(1:length(edges(sys.g)))[Not(unknown_F(sys))]
	# substitute the unknown flows in all balance equations
	bal_eq = flow_balance(sys)
	sub_bal_eq = Array{Equation}(undef, length(bal_eq)+length(i_known_F))
	for i=1:length(bal_eq)
		sub_bal_eq[i] = substitute(bal_eq[i], sub_dict)
	end
	for i=1:length(i_known_F)
		j = i_known_F[i]
		sub_bal_eq[length(bal_eq)+i] = F[j] - A*(P²[srcE[j]]-P²[dstE[j]])/κ[j] ~ 0
	end
	return sub_bal_eq
end

# ╔═╡ 4702ad20-72e9-4e7d-b5ef-4dab809183fc
collect(1:length(edges(sys_split.g)))[Not(unknown_F(sys_split))]

# ╔═╡ b4e1ae34-2ec0-433d-9664-21e7e1fa564b
substitute_unknown_flows(sys_split)

# ╔═╡ 71c1e3e2-cdcf-426b-8409-2a6cbd3f966c
flow_balance(sys_series)

# ╔═╡ dd729291-42a7-40dd-b383-1728a1f13ee8
unknown_F(sys_series)

# ╔═╡ 289628a2-9ac4-449a-8890-163b4a75c804
GasChromatographySystems.unknown_p(sys_series)

# ╔═╡ dd010e56-1a71-4bc4-9a39-15233588158f
substitute_unknown_flows(sys_series)

# ╔═╡ e7d3833f-7cac-425f-b1a7-468afafe5492
substitute_unknown_flows(sys_cs)

# ╔═╡ 89f63619-d588-4358-8259-dfb0eaf0aad5
length(unknown_F(sys_series))

# ╔═╡ ddcd2f49-b3ce-4dd3-bc65-a08a60f46756
length(GasChromatographySystems.unknown_p(sys_series))

# ╔═╡ e1a387b7-4e54-4f1d-bd20-6dcdba4d9d4a
length(substitute_unknown_flows(sys_series))

# ╔═╡ bf15d0c9-c328-4f91-9cd1-2c0b7ff9c93a
function solve_balance(sys)
	@variables A, P²[1:nv(sys.g)], κ[1:ne(sys.g)], F[1:ne(sys_split.g)]
	i_unknown_p = GasChromatographySystems.unknown_p(sys) # indices of the nodes with unknown pressure
	i_unknown_F = unknown_F(sys) # indices of the edges with an unknown flow
	bal_eq = substitute_unknown_flows(sys)
	#num_use_eq = GasChromatographySystems.unknowns_in_flow_balances(sys)[2]
	if length(i_unknown_p) == length(bal_eq)
		sol = Symbolics.solve_for(bal_eq, [P²[i_unknown_p[i]] for i=1:length(i_unknown_p)])
	#elseif length(i_unknown_p) == 1
	#	sol = Symbolics.solve_for(F_balance[i_unknown_p[1]-1], [P²[i_unknown_p[1]]])
	#elseif length(i_unknown_p) == (length(F_balance) - 1)
		#if length(findall(minimum(num_use_eq).==num_use_eq)) == 1
	#		leave_out_eq = findlast(num_use_eq.==minimum(num_use_eq))
		#else
		#	leave_out_eq = 
		#end
	#	sol = Symbolics.solve_for(F_balance[Not(leave_out_eq)], [P²[i_unknown_p[i]] for i=1:length(i_unknown_p)])
	#elseif length(i_unknown_p) < (length(F_balance) - 1)
	#	error("More flow balance equations than unknown pressures. ToDo: leave equations out.")
	elseif length(i_unknown_p) > length(bal_eq)
		error("More unknown pressures than flow balance equations.")
	else # loop of the i_unknown_p
		# identifie inner nodes which are unknown pressures, there index in the inner_vertices() is the index of the flow balance equations to use 
		# this works only for unknown_p which are inner nodes, unknown_p at outer nodes (e.g. inlet pressure), lead to error 
		inner_V = GasChromatographySystems.inner_vertices(sys.g) 
		bal_eq_i = Array{Int}(undef, length(i_unknown_p))
		for i=1:length(i_unknown_p)
			# if the unknown pressure is not an inner node, than add a equation from the end of the list of balance equation, which should be a flow definition.
			if isnothing(findfirst(i_unknown_p[i].==inner_V))
				bal_eq_i[i] = length(bal_eq)-(i+0)
			else
				bal_eq_i[i] = findfirst(i_unknown_p[i].==inner_V)
			end
		end
		sol = Symbolics.solve_for([bal_eq[x] for x in bal_eq_i], [P²[i_unknown_p[i]] for i=1:length(i_unknown_p)])
	end
	return sol
end

# ╔═╡ b31e2dec-8712-41c6-b6b1-6a0fb00d7035
solve_balance(sys_series)

# ╔═╡ 59f3a04e-df80-4bad-a898-ed21ce06f953
solve_balance(sys_split)

# ╔═╡ 7a7a2347-4852-46f6-ace3-7ebafeb03564
solve_balance(sys_cs)

# ╔═╡ f5a56163-f78e-4b4c-920c-34ec2a7f50aa
function solve_pressure(sys)
	@variables A, P²[1:nv(sys.g)], κ[1:ne(sys.g)], F[1:ne(sys.g)]
	#balance = flow_balance(sys)
	i_unknown_p = GasChromatographySystems.unknown_p(sys)
	i_known_F = collect(1:length(edges(sys.g)))[Not(unknown_F(sys))]
	solutions = solve_balance(sys)
	a = π/256 * GasChromatographySystems.Tn/GasChromatographySystems.pn
	κs = GasChromatographySystems.flow_restrictions(sys)
	p²s = GasChromatographySystems.pressures_squared(sys)
	p_solution = Array{Function}(undef, length(i_unknown_p))
	for i=1:length(i_unknown_p)
		sub_dict(t) = merge(Dict((P²[j] => p²s[j](t) for j=setdiff(1:nv(sys.g), i_unknown_p))), Dict((κ[j] => κs[j](t) for j=1:ne(sys.g))), Dict(A => a), Dict(F[j] => sys.modules[j].flow for j in i_known_F))
		f(t) = sqrt(substitute(solutions[i], sub_dict(t)))
		p_solution[i] = f
	end
	return p_solution, i_unknown_p
end

# ╔═╡ d5f3d0ec-6d9c-458a-97ed-bf2b3003dd9c
sol = solve_pressure(sys_series)

# ╔═╡ d39a02a8-3b58-45a7-a6d2-f8d9374047a3
sol[1][1](0)

# ╔═╡ 11a0cb05-1f14-4d51-bdd7-a6b6577e2a6f
function pressure_functions(sys)
	pres, unk = solve_pressure(sys)
	#p²s = pressures_squared(sys)
	p_func = Array{Any}(undef, nv(sys.g))
	for i=1:nv(sys.g)
		f = if i in unk
			pres[findfirst(unk.==i)]
		else
			GasChromatographySimulator.steps_interpolation(sys.pressurepoints[i].timesteps, identity.(sys.pressurepoints[i].pressure_steps))
		end
	p_func[i] = f
	end
	return p_func
end

# ╔═╡ 9f16c03e-8d0c-4d88-a5de-1ee23a5d15b7
function flow_functions(sys)
	p_func = pressure_functions(sys)
	F_func = Array{Function}(undef, ne(sys.g))
	E = collect(edges(sys.g))
	srcE = src.(E)
	dstE = dst.(E)
	for i=1:ne(sys.g)
		pin(t) = p_func[srcE[i]](t)
		pout(t) = p_func[dstE[i]](t)
		if typeof(sys.modules[i].temperature) <: GasChromatographySystems.TemperatureProgram
			T_itp = GasChromatographySimulator.temperature_interpolation(sys.modules[i].temperature.timesteps, sys.modules[i].temperature.temperaturesteps, sys.modules[i].temperature.gradient_function, sys.modules[i].length)
		elseif typeof(sys.modules[i].temperature) <: Number
			gf(x) = [zero(x), zero(x)]
			T_itp = GasChromatographySimulator.temperature_interpolation([0.0, 3600.0], [sys.modules[i].temperature, sys.modules[i].temperature], gf, sys.modules[i].length)
		end
		f(t) = GasChromatographySimulator.flow(t, T_itp, pin, pout, sys.modules[i].length, sys.modules[i].diameter, sys.options.mobile_phase; ng=sys.options.ng, vis=sys.options.vis, control=sys.options.control)
		F_func[i] = f
	end
	return F_func
end

# ╔═╡ 0452b0ef-b677-420f-93af-0acc7c77935d
function plot_graph_with_flow(sys, t; lay = Spring(), color=:lightblue, node_size=80, arrow_size=20, arrow_shift=0.8, dataaspect=false, nlabels_fontsize=14, elabels_fontsize=14, elabels_distance = 20)
	p_func = pressure_functions(sys)
	F_func = flow_functions(sys)
	fig, ax, p = GraphMakie.graphplot(sys.g, 
						layout=lay,
						nlabels = [sys.pressurepoints[i].name*"\n$(trunc(Int, p_func[i](t)/1000))kPa" for i in 1:nv(sys.g)],
						nlabels_align=(:center,:center),
						nlabels_fontsize = nlabels_fontsize,
						node_size = [node_size for i in 1:nv(sys.g)],
						node_color = [color for i in 1:nv(sys.g)],
						elabels = [sys.modules[i].name*"\n $(round(F_func[i](t)*60*1e6; sigdigits=2)) mL/min" for i in 1:ne(sys.g)],
						elabels_distance = elabels_distance,
						elabels_fontsize = elabels_fontsize,
						arrow_size = arrow_size,
						arrow_shift = arrow_shift
					)
	hidedecorations!(ax)
	hidespines!(ax)
	if dataaspect == true
		ax.aspect = DataAspect()
	end
	return fig
end

# ╔═╡ dedf7afb-a36d-472b-82d0-754e703bc8a9
md"""
### Split
"""

# ╔═╡ 27c84842-a39a-4b47-bb4b-c79232fe98c5
plot_graph_with_flow(sys_split,0)

# ╔═╡ d9df5e5a-4501-4b6f-8eb5-b7f82dfd5bcf
function plot_flow_over_time(sys; dt=60.0)
	#plotly()
	F_func = flow_functions(sys)
	com_timesteps = GasChromatographySystems.common_timesteps(sys)
	trange = 0:dt:sum(com_timesteps)
	p_flow = Plots.plot(xlabel="time in s", ylabel="flow in mL/min")
	for i=1:ne(sys.g)
		Plots.plot!(p_flow, trange, F_func[i].(trange).*60e6, label="F_$(sys.modules[i].name)")
	end
	return p_flow
end

# ╔═╡ fc04fdec-8706-493d-b1a3-d4129d92295f
plot_flow_over_time(sys_split)

# ╔═╡ ef0dafce-81b0-4bc5-88ed-1acd351876c4
function plot_pressure_over_time(sys; dt=60.0)
	#plotly()
	p_func = pressure_functions(sys)
	com_timesteps = GasChromatographySystems.common_timesteps(sys)
	trange = 0:dt:sum(com_timesteps)
	p_pres = Plots.plot(xlabel="time in s", ylabel="pressure in Pa(a)", legend=:topleft)
	for i=1:nv(sys.g)
		Plots.plot!(p_pres, trange, p_func[i].(trange), label="$(sys.pressurepoints[i].name)")
	end
	return p_pres
end

# ╔═╡ b9b26d01-f614-47bf-b7cc-f55410e87cd8
plot_pressure_over_time(sys_split)

# ╔═╡ cd059b8c-cf58-4ec8-93e6-55f45f89e037
md"""
### Series
"""

# ╔═╡ 57f97467-1b96-4ff1-83c4-4c11c083dadc
plot_graph_with_flow(sys_series,0)

# ╔═╡ 61f01373-00ff-4c62-a690-54a475d50853
plot_flow_over_time(sys_series)

# ╔═╡ 3ede7f5d-8423-4f1b-9eac-24d55885841a
plot_pressure_over_time(sys_series)

# ╔═╡ c536018b-0ae9-461d-b14c-06e6dd0b8945
md"""
### Complex System
"""

# ╔═╡ 4098cf91-bee6-4aa7-ac8d-028a09a9b5af
substitute_unknown_flows(sys_cs)

# ╔═╡ 83236f28-33ca-489b-991b-9ef4316d64fb
solve_balance(sys_cs)

# ╔═╡ 0130dc96-8bf8-4ac2-823c-c7956c5c8e23
plot_graph_with_flow(sys_cs,1500)

# ╔═╡ f695a4a4-6a20-49b5-bd52-9154526ec75c
plot_flow_over_time(sys_cs)

# ╔═╡ 25d11148-4ce7-4bb0-91b6-ac13bbd5d3eb
plot_pressure_over_time(sys_cs)

# ╔═╡ d5a39d20-d3ee-4ec1-b12d-277031e6ff55
sys_cs

# ╔═╡ Cell order:
# ╠═7722a5a0-e41d-11ed-3f24-b54f45f9a09b
# ╠═1acc69fb-ea06-470a-a6c4-b336c64e7c62
# ╠═7195fdd4-5924-4033-b794-2a491336cce0
# ╠═724ac332-c561-416e-83c7-bb5a339f26e9
# ╠═14a9675d-7cb0-42aa-9f8e-352770dd7b23
# ╠═be577666-dc2a-40d1-9a85-643928778949
# ╠═8d101c8e-a7df-46b7-81e0-aaa1de6e3b23
# ╠═c5bca10a-9f4d-4f6b-a162-9e5c58d96a9d
# ╠═34723a63-e7fa-4d26-a047-843ccf124804
# ╠═fcb30907-6fad-4708-a8f4-f8fab3303b5d
# ╠═9736b84a-16a4-43df-a735-a5bc4c98c670
# ╠═ec5fd728-8ae9-4193-9013-3e3d94a1fc71
# ╠═d7e5df4f-d23d-4ab0-a6d2-830e5e24251e
# ╠═126ab6c1-632c-4d27-a457-38b8759ae4a2
# ╠═f9acd758-c86a-4ecf-bf1f-f8a26a4b95de
# ╠═963c47a9-e6e9-4f2a-817d-e29c55b974e4
# ╠═9c337183-adad-463b-af88-9d0844dee46e
# ╠═6fda756c-3a32-4a99-b9df-1e40c750f8ec
# ╠═9dc540a5-dd43-4b75-b174-4fb2453d3694
# ╠═bf770d0a-7347-47a1-95f8-8669efde4716
# ╠═9a3a7286-95eb-424b-8617-8a1f6058d7c8
# ╠═6f30dc72-4928-475e-86b4-261256f77ab5
# ╠═b7b5069c-ad71-4dcc-b460-2cdd057c8da7
# ╠═2a32dfdb-d805-4437-9658-6fba5365883c
# ╠═eaa0b056-4a6b-4f6c-9ff8-a0e4d0b8ddb8
# ╠═5f4eccb4-2a86-4705-b8cf-245dbbf78c8c
# ╠═272eb44b-9103-49d8-b1b4-811a0122ae6e
# ╠═9d60fb4d-d013-45af-890a-3474a521f283
# ╠═ed44d3b8-b8ab-49a7-9fe3-2a54ac6de92b
# ╠═be251bb9-a4a8-4ca7-a7ca-9a4cd3c3409d
# ╠═870b5176-94f8-4915-83a1-77c4dbcdbb0a
# ╠═814655cb-ff2b-4042-a1ee-09cb4dfc173d
# ╠═b3054fae-a2ea-47dc-9bb6-ce526e9e9b39
# ╠═1942d1cc-c573-44f2-8055-fc7140762a90
# ╠═38be367b-3d58-46a7-8062-49771bd52952
# ╠═1fef6d42-b51f-4a4c-b09e-1750351b66eb
# ╠═7be3d2e8-31cb-4a8f-80a0-18d75f9fb6b6
# ╠═82fbba32-886a-4f6f-a855-7d0947c58afe
# ╠═3a886889-6655-4b6d-b2f5-a50b2053da36
# ╠═bf8c00b0-10c5-4de1-a4ab-27105e2b06b1
# ╠═fa9eca1c-682b-495c-afc4-16bee0d58ce9
# ╠═b28571c6-10c4-430a-8e62-c39571b0a38b
# ╠═0b2a67f3-bcdb-4c51-8322-42b51d43b71c
# ╠═4702ad20-72e9-4e7d-b5ef-4dab809183fc
# ╠═b4e1ae34-2ec0-433d-9664-21e7e1fa564b
# ╠═71c1e3e2-cdcf-426b-8409-2a6cbd3f966c
# ╠═dd729291-42a7-40dd-b383-1728a1f13ee8
# ╠═289628a2-9ac4-449a-8890-163b4a75c804
# ╠═dd010e56-1a71-4bc4-9a39-15233588158f
# ╠═e7d3833f-7cac-425f-b1a7-468afafe5492
# ╠═89f63619-d588-4358-8259-dfb0eaf0aad5
# ╠═ddcd2f49-b3ce-4dd3-bc65-a08a60f46756
# ╠═e1a387b7-4e54-4f1d-bd20-6dcdba4d9d4a
# ╠═bf15d0c9-c328-4f91-9cd1-2c0b7ff9c93a
# ╠═b31e2dec-8712-41c6-b6b1-6a0fb00d7035
# ╠═59f3a04e-df80-4bad-a898-ed21ce06f953
# ╠═7a7a2347-4852-46f6-ace3-7ebafeb03564
# ╠═f5a56163-f78e-4b4c-920c-34ec2a7f50aa
# ╠═d5f3d0ec-6d9c-458a-97ed-bf2b3003dd9c
# ╠═d39a02a8-3b58-45a7-a6d2-f8d9374047a3
# ╠═11a0cb05-1f14-4d51-bdd7-a6b6577e2a6f
# ╠═9f16c03e-8d0c-4d88-a5de-1ee23a5d15b7
# ╠═0452b0ef-b677-420f-93af-0acc7c77935d
# ╠═dedf7afb-a36d-472b-82d0-754e703bc8a9
# ╠═27c84842-a39a-4b47-bb4b-c79232fe98c5
# ╠═d9df5e5a-4501-4b6f-8eb5-b7f82dfd5bcf
# ╠═fc04fdec-8706-493d-b1a3-d4129d92295f
# ╠═ef0dafce-81b0-4bc5-88ed-1acd351876c4
# ╠═b9b26d01-f614-47bf-b7cc-f55410e87cd8
# ╠═cd059b8c-cf58-4ec8-93e6-55f45f89e037
# ╠═57f97467-1b96-4ff1-83c4-4c11c083dadc
# ╠═61f01373-00ff-4c62-a690-54a475d50853
# ╠═3ede7f5d-8423-4f1b-9eac-24d55885841a
# ╠═c536018b-0ae9-461d-b14c-06e6dd0b8945
# ╠═4098cf91-bee6-4aa7-ac8d-028a09a9b5af
# ╠═83236f28-33ca-489b-991b-9ef4316d64fb
# ╠═0130dc96-8bf8-4ac2-823c-c7956c5c8e23
# ╠═f695a4a4-6a20-49b5-bd52-9154526ec75c
# ╠═25d11148-4ce7-4bb0-91b6-ac13bbd5d3eb
# ╠═d5a39d20-d3ee-4ec1-b12d-277031e6ff55
