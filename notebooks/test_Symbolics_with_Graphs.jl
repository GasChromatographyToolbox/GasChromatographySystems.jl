### A Pluto.jl notebook ###
# v0.19.17

using Markdown
using InteractiveUtils

# ╔═╡ 00a70e44-28ff-11ed-3a74-b939bdc0aa6a
begin
	##import Pkg
    # activate the shared project environment
    ##Pkg.activate(Base.current_project())
	#using DataFrames, CSV, Interpolations, Plots, QuadGK, DifferentialEquations, ForwardDiff, Intervals
	using Plots
	##using GasChromatographySystems
	using Graphs
	using GraphRecipes
	using PlutoUI
	# using GLMakie
	using CairoMakie
	using GraphMakie
	using NetworkLayout
	using Symbolics
	TableOfContents()
end

# ╔═╡ 62270453-ed6c-4e07-afb8-bb2a72e1fb54
md"""
# Test of Symbolics.jl with Graphs
"""

# ╔═╡ 1e6aaf53-831d-4717-8296-a3db8221020e
md"""
## Defining a Graph

Using the example from the abstract for Liege 2023 as graph.
"""

# ╔═╡ c486ef03-a6df-4e41-b5c2-f63f59934e06
begin
	# setup of the Graph
	g = SimpleDiGraph(9)
	add_edge!(g, 1, 2) # Inj->GC1->pressure regulator 
	add_edge!(g, 2, 3) # pressure regulator -> GC2 -> split
	add_edge!(g, 3, 4) # Split -> TL1 -> Det1
	add_edge!(g, 3, 5) # Split -> TL2 -> Modulator inlet
	add_edge!(g, 5, 6) # Modulator inlet -> Modulator -> GC3 inlet
	add_edge!(g, 6, 7) # GC3 inlet -> GC3 -> Split
	add_edge!(g, 7, 8) # Split -> TL3 -> Det2
	add_edge!(g, 7, 9) # Split -> TL4 -> Det3
	edge_labels = ["GC₁", "GC₂", "TL₁", "TL₂", "Modulator", "GC₃", "TL₃", "TL₄"]
	node_labels = ["p₁", "p₂", "p₃", "p₄", "p₅", "p₆", "p₇", "p₈", "p₉"]
	lay = Stress()
	fig, ax, p = GraphMakie.graphplot(g, 
						layout=lay,
						nlabels=node_labels, 
						nlabels_align=(:center,:center),
						node_size = [40 for i in 1:nv(g)],
						node_color = [:lightblue for i in 1:nv(g)],
						elabels = edge_labels
					)
	hidedecorations!(ax)
	hidespines!(ax)
	ax.aspect = DataAspect()
	fig
end

# ╔═╡ d57b71a0-d105-4ce5-af27-11f660c1d02a
E = collect(edges(g))

# ╔═╡ 90318630-d3db-4c34-aff2-f64451160c5d
V = collect(vertices(g))

# ╔═╡ 223cb31c-bce7-412d-af77-f184de61de29
deg = Graphs.degree(g)

# ╔═╡ dbcd4234-6d2f-4664-beb3-7be07f6c64ba
md"""
## Create variables
"""

# ╔═╡ e640ba26-2cd2-4668-a2ed-03db48805c02
md"""
Define the square of the pressure as variable -> balance equations are linear for these.
"""

# ╔═╡ 4b35d9b6-12ed-49d8-8e86-481b60bfc7e5
@variables P²[1:nv(g)], κ[1:ne(g)]

# ╔═╡ 3f554c27-8658-4414-bd92-2cbd1a94e018
md"""
## Create the equations of flow balance from the Graph
"""

# ╔═╡ 0fd7f836-ac86-4d93-aedb-39cac6ff333e
ii = findall(Graphs.degree(g).>1) # indices of `inner nodes`

# ╔═╡ 65d48d0d-1f47-44a1-beb1-20a5f9aac925
ii[2]

# ╔═╡ 3121f8e9-ba3a-4761-9e1b-827299a52392
srcE = src.(E) # the source nodes of the edges 

# ╔═╡ 4cd7955f-afb9-4196-b8ac-12b574bb8417
i_src = findall(srcE.==ii[2])  # indices (of `srcE`) of the nodes, which `ii[2]=3`, -> 3th and 4th edge have node 3 as source

# ╔═╡ 2354dbf3-b90e-49a8-b2db-277b201afffa
dstE = dst.(E) # the destination nodes of the edges

# ╔═╡ 40074400-5281-4f41-a288-83bebbbefe67
i_dst = findall(dstE.==ii[2]) # indices (of `dstE`) of the nodes, which `ii[2]=3`, -> 2nd edge has node 3 as source

# ╔═╡ 26739daf-fc72-48e4-afd4-814bdde40335
(P²[srcE[i_dst[1]]]-P²[dstE[i_dst[1]]])/κ[findall(dstE.==ii[2])[1]] # construction of a symbolic expression for the flow into node 3

# ╔═╡ 0239c209-bb5a-4fa7-a4b9-71b4f85f1b03
(P²[srcE[i_src[2]]]-P²[dstE[i_src[2]]])/κ[findall(srcE.==ii[2])[2]] # construction of a symbolic expression for the flow from node 3 to node 5 (the second destination of node 3)

# ╔═╡ 5ecbbaf1-e8a9-4c97-9054-c613d7e87306
begin
	i_inner = findall(Graphs.degree(g).>1)
	balance = Array{Num}(undef, length(i_inner)) # correlation of the balance equations to the array does not work, only the last element is correct
	for i=length(i_inner)
		ii = i_inner[i]
		# find edges, where node `i` is the source
		i_src = findall(src.(E).==ii)
		# find edges, where node `i` is the destination
		i_dst = findall(dst.(E).==ii)
		balance[i] = 0
		for j=1:length(i_dst)
			balance[i] = balance[i] + (P²[srcE[i_dst[j]]]-P²[dstE[i_dst[j]]])/κ[findfirst(dstE.==ii)[j]]
		end
		for j=1:length(i_src)
			balance[i] = balance[i] - (P²[srcE[i_src[j]]]-P²[dstE[i_src[j]]])/κ[findall(srcE.==ii)[j]]
		end
		#push!(balance, balance_)
		#balance[i] = balance_
	end
end

# ╔═╡ 087bfc21-de03-46ee-99f1-be247da061f5
balance

# ╔═╡ 5b159783-5a44-4c31-b333-2512581c6608
balance[5] ~ 0

# ╔═╡ 570c4a24-f31f-455d-bc9f-42416e2c4099
inner_V = findall(Graphs.degree(g).>1)

# ╔═╡ 2f8406f6-ada7-483a-abd8-c31fe363b8a3
bal_eq = Array{Symbolics.Equation}(undef, length(inner_V))

# ╔═╡ 89662190-13e1-486e-920f-e047e69cf46b
md"""
System of equations for the pressures (the squares of the pressures)
"""

# ╔═╡ adea4142-52df-4229-819b-6a44f5b18973
bal_eq

# ╔═╡ a686c1d5-bd07-42ba-828b-024f757c6cb6
bal_eq_rhs = Array{Symbolics.Num}(undef, length(inner_V))

# ╔═╡ 66cc1647-e274-4871-adec-772178e7e93f
bal_eq_rhs

# ╔═╡ 6d0e1d5f-9f8a-4e96-989c-0b5131b04400
md"""
With ` ~ 0` the symbolic expresion of type `Symbolics.Num` becomes an equation of the type `Symbolics.Equation`.

Some manipulations are possible with `Symbolics.Num` but not with `Symbolics.Equation`, e.g. adding to expressions/equations:
"""

# ╔═╡ c35ee1dd-2514-4851-b229-71b60052d0e2
bal_eq[3] + bal_eq[4]

# ╔═╡ b0b708c7-b24a-42f2-a9d2-6ead62d7558b
bal_eq_rhs[3] + bal_eq_rhs[4]

# ╔═╡ 58836887-701b-49b0-aa17-16205175738d
md"""
## Solve the linear balance equations
"""

# ╔═╡ db0ffead-08ec-477f-b54b-e629b464ac4e
# unknown pressures:
unknownP² = P²[3], P²[5], P²[6], P²[7]

# ╔═╡ 19464015-f3e8-4158-9290-944c2f957e8d
md"""
First equation has only on unknown: ``P_3^2``.
"""

# ╔═╡ 970d759a-4f2a-4351-beac-3e9e067668d1
Symbolics.solve_for(bal_eq[1], P²[3])

# ╔═╡ 47556faf-29b6-42df-8ca7-c2307d09e328
Symbolics.solve_for([bal_eq[1], bal_eq[2]], [P²[3], P²[5]])

# ╔═╡ 5797032b-94c0-4c5a-9137-52460bfe8849
Symbolics.solve_for([bal_eq[1], bal_eq[2], bal_eq[3]], [P²[3], P²[5], P²[6]])

# ╔═╡ 4789d19a-2933-41c8-9703-fa60d3f0e859
Symbolics.solve_for([bal_eq[1], bal_eq[2], bal_eq[3], bal_eq[5]], [P²[3], P²[5], P²[6], P²[7]])

# ╔═╡ f73a20a3-2e23-46ac-a4b6-837913a57a86
md"""
Balance equations (3) and (4) can be combined into one, by adding them.

**Balance equations in following nodes without splitting (the degree of the neighboring nodes is 2) can be added to reduce the number of equations.**
"""

# ╔═╡ bb52f00e-fa4e-457c-b982-0a0c20cb8b44
md"""
## Another example
"""

# ╔═╡ c3a5dabf-0d45-4ab8-b850-60e0d6fd983f
begin
	g2 = SimpleDiGraph(14)
	add_edge!(g2, 1, 2)
	add_edge!(g2, 2, 3)
	add_edge!(g2, 2, 5)
	add_edge!(g2, 2, 6)
	add_edge!(g2, 3, 4)
	add_edge!(g2, 4, 10)
	add_edge!(g2, 5, 8)
	add_edge!(g2, 6, 9)
	add_edge!(g2, 6, 10)
	add_edge!(g2, 7, 11)
	add_edge!(g2, 9, 12)
	add_edge!(g2, 10, 13)
	add_edge!(g2, 11, 13)
	add_edge!(g2, 13, 14)

	fig2, ax2, p2 = GraphMakie.graphplot(g2, 
						layout=Stress(),
						nlabels=["$(i)" for i in 1:nv(g2)], 
						nlabels_align=(:center,:center),
						node_size = [40 for i in 1:nv(g2)],
						node_color = [:lightblue for i in 1:nv(g2)]
					)
	hidedecorations!(ax2)
	hidespines!(ax2)
	ax2.aspect = DataAspect()
	fig2
end

# ╔═╡ 5a05bf2b-cd76-4f80-8d0c-a175ec23e352
V2 = collect(vertices(g2))

# ╔═╡ df28b32f-b8f1-45c7-b370-dd901aa737a4
inner_V2 = findall(Graphs.degree(g2).>1)

# ╔═╡ 4127c3c2-f78b-4dcf-bc90-cf08c2b7db9e
bal_eq2 = Array{Symbolics.Equation}(undef, length(inner_V2))

# ╔═╡ e1a4d66d-8a49-45bf-814a-7f7b44ff2b93
bal_eq2

# ╔═╡ 15afefc5-dfb9-4544-a9ab-c90a134ca475
@variables P²__[1:nv(g2)], κ__[1:ne(g2)]

# ╔═╡ d1d87ffb-5253-40f4-a522-f8f7e76be1fb
unknown_P²__ = [P²__[2], P²__[3], P²__[4], P²__[5], P²__[6], P²__[9], P²__[10], P²__[11], P²__[13]]

# ╔═╡ aee2c695-3b86-4cf3-b5c4-6168e2378b27
#sol = Symbolics.solve_for(bal_eq2, unknown_P²__)

# ╔═╡ 2024ed58-bd26-49cb-892c-9dbad661dce8
md"""
# End
"""

# ╔═╡ 80874883-e4b0-4e94-9e2c-574f3be348da
function flow_balance(g, i_n, P², κ)
	#@variables P²[1:nv(g)], κ[1:ne(g)]
	E = collect(edges(g))
	srcE = src.(E)
	dstE = dst.(E)
	# find edges, where node `i_n` is the source
	i_src = findall(srcE.==i_n)
	# find edges, where node `i_n` is the destination
	i_dst = findall(dstE.==i_n)
	balance = 0
	for j=1:length(i_dst) # ingoing flows
		balance = balance + (P²[srcE[i_dst[j]]]-P²[dstE[i_dst[j]]])/κ[i_dst[j]]
	end
	for j=1:length(i_src) # outgoing flows
		balance = balance - (P²[srcE[i_src[j]]]-P²[dstE[i_src[j]]])/κ[i_src[j]]
	end
	return balance
end

# ╔═╡ 8944e501-a6ca-4233-a58f-b49ca86ae47d
flow_balance(g, 2, P², κ)

# ╔═╡ 7539ec79-a944-4783-b045-21d68c68fb83
for i=1:length(inner_V)
	bal_eq[i] = flow_balance(g, inner_V[i], P², κ) ~ 0
end

# ╔═╡ 7c24d387-8664-4998-9476-fa35f0fc0755
for i=1:length(inner_V)
	bal_eq_rhs[i] = flow_balance(g, inner_V[i], P², κ)
end

# ╔═╡ 23c25c5e-cf47-42fa-a833-d1b4e801def1
for i=1:length(inner_V2)
	bal_eq2[i] = flow_balance(g2, inner_V2[i], P²__, κ__) ~ 0
end

# ╔═╡ Cell order:
# ╠═00a70e44-28ff-11ed-3a74-b939bdc0aa6a
# ╠═62270453-ed6c-4e07-afb8-bb2a72e1fb54
# ╠═1e6aaf53-831d-4717-8296-a3db8221020e
# ╠═c486ef03-a6df-4e41-b5c2-f63f59934e06
# ╠═d57b71a0-d105-4ce5-af27-11f660c1d02a
# ╠═90318630-d3db-4c34-aff2-f64451160c5d
# ╠═223cb31c-bce7-412d-af77-f184de61de29
# ╠═dbcd4234-6d2f-4664-beb3-7be07f6c64ba
# ╠═e640ba26-2cd2-4668-a2ed-03db48805c02
# ╠═4b35d9b6-12ed-49d8-8e86-481b60bfc7e5
# ╠═3f554c27-8658-4414-bd92-2cbd1a94e018
# ╠═0fd7f836-ac86-4d93-aedb-39cac6ff333e
# ╠═65d48d0d-1f47-44a1-beb1-20a5f9aac925
# ╠═3121f8e9-ba3a-4761-9e1b-827299a52392
# ╠═4cd7955f-afb9-4196-b8ac-12b574bb8417
# ╠═2354dbf3-b90e-49a8-b2db-277b201afffa
# ╠═40074400-5281-4f41-a288-83bebbbefe67
# ╠═26739daf-fc72-48e4-afd4-814bdde40335
# ╠═0239c209-bb5a-4fa7-a4b9-71b4f85f1b03
# ╠═5ecbbaf1-e8a9-4c97-9054-c613d7e87306
# ╠═087bfc21-de03-46ee-99f1-be247da061f5
# ╠═5b159783-5a44-4c31-b333-2512581c6608
# ╠═8944e501-a6ca-4233-a58f-b49ca86ae47d
# ╠═570c4a24-f31f-455d-bc9f-42416e2c4099
# ╠═2f8406f6-ada7-483a-abd8-c31fe363b8a3
# ╠═89662190-13e1-486e-920f-e047e69cf46b
# ╠═7539ec79-a944-4783-b045-21d68c68fb83
# ╠═adea4142-52df-4229-819b-6a44f5b18973
# ╠═a686c1d5-bd07-42ba-828b-024f757c6cb6
# ╠═7c24d387-8664-4998-9476-fa35f0fc0755
# ╠═66cc1647-e274-4871-adec-772178e7e93f
# ╠═6d0e1d5f-9f8a-4e96-989c-0b5131b04400
# ╠═c35ee1dd-2514-4851-b229-71b60052d0e2
# ╠═b0b708c7-b24a-42f2-a9d2-6ead62d7558b
# ╠═58836887-701b-49b0-aa17-16205175738d
# ╠═db0ffead-08ec-477f-b54b-e629b464ac4e
# ╠═19464015-f3e8-4158-9290-944c2f957e8d
# ╠═970d759a-4f2a-4351-beac-3e9e067668d1
# ╠═47556faf-29b6-42df-8ca7-c2307d09e328
# ╠═5797032b-94c0-4c5a-9137-52460bfe8849
# ╠═4789d19a-2933-41c8-9703-fa60d3f0e859
# ╠═f73a20a3-2e23-46ac-a4b6-837913a57a86
# ╠═bb52f00e-fa4e-457c-b982-0a0c20cb8b44
# ╠═c3a5dabf-0d45-4ab8-b850-60e0d6fd983f
# ╠═5a05bf2b-cd76-4f80-8d0c-a175ec23e352
# ╠═df28b32f-b8f1-45c7-b370-dd901aa737a4
# ╠═4127c3c2-f78b-4dcf-bc90-cf08c2b7db9e
# ╠═23c25c5e-cf47-42fa-a833-d1b4e801def1
# ╠═e1a4d66d-8a49-45bf-814a-7f7b44ff2b93
# ╠═15afefc5-dfb9-4544-a9ab-c90a134ca475
# ╠═d1d87ffb-5253-40f4-a522-f8f7e76be1fb
# ╠═aee2c695-3b86-4cf3-b5c4-6168e2378b27
# ╟─2024ed58-bd26-49cb-892c-9dbad661dce8
# ╠═80874883-e4b0-4e94-9e2c-574f3be348da
