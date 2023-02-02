### A Pluto.jl notebook ###
# v0.19.17

using Markdown
using InteractiveUtils

# ╔═╡ e7111e94-1fd0-11ed-3423-c1975e0240ea
begin
	import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
	#using DataFrames, CSV, Interpolations, Plots, QuadGK, DifferentialEquations, ForwardDiff, Intervals
	using Plots
	using GasChromatographySystems
	using Graphs
	using GraphRecipes
	using PlutoUI
	# using GLMakie
	using CairoMakie
	using GraphMakie
	using NetworkLayout
	TableOfContents()
end

# ╔═╡ 00cf5a16-8a1e-431d-a40e-ef620b1c1508
md"""
# Test directed Graphs

Using the flow calculator for the FF-TG-GC as complex example.
"""

# ╔═╡ d53a1637-9705-4d32-9738-2521fc4119f8
g = SimpleDiGraph(14) # start a graph with 14 vertices

# ╔═╡ 1da8f8d8-2c51-41b1-b14a-18deb6e0f50a
# adding edges
begin
	add_edge!(g, 1, 2)
	add_edge!(g, 2, 3)
	add_edge!(g, 2, 5)
	add_edge!(g, 2, 6)
	add_edge!(g, 3, 4)
	add_edge!(g, 4, 10)
	add_edge!(g, 5, 8)
	add_edge!(g, 6, 9)
	add_edge!(g, 6, 10)
	add_edge!(g, 7, 11)
	add_edge!(g, 9, 12)
	add_edge!(g, 10, 13)
	add_edge!(g, 11, 13)
	add_edge!(g, 13, 14)
end

# ╔═╡ 06d4b0b5-9280-43d5-a51c-253266e1d989
g

# ╔═╡ 96f6b288-ad6a-44b1-bb42-a63e9560e0cb
md"""
## Plotting with GraphRecipes.jl
"""

# ╔═╡ 1d03023b-176c-4057-b92d-94257b052a8d
GraphRecipes.graphplot(g; curves=false)

# ╔═╡ ea9453a4-3171-43a4-a297-3f5378ad4370
md"""
## Plotting with GraphMakie.jl
"""

# ╔═╡ a9d942c5-05d9-439e-b651-59a9a27b30fc
GraphMakie.graphplot(g)

# ╔═╡ 13fff928-4655-4e6f-9d66-ad21cf66261e
md"""
## Interaction Examples GraphMakie.jl

Using examples from the documentation of GraphMakie.jl.

Seems not to work in the notebook.
"""

# ╔═╡ dcbdbc38-e338-48bf-a1ca-78a3669b0ed0
begin
	f, ax, p = GraphMakie.graphplot(g,
									edge_width = [2.0 for i in 1:ne(g)],
									edge_color = [colorant"gray" for i in 1:ne(g)],
									node_size = [10 for i in 1:nv(g)],
									node_color = [colorant"red" for i in 1:nv(g)])
	hidedecorations!(ax)
	hidespines!(ax)
	ax.aspect = DataAspect()
	# later want to enable drag interactions -> disable default :rectanglezoom interaction
	deregister_interaction!(ax, :rectanglezoom)
	f
end

# ╔═╡ 509e3f76-5afa-4c13-8449-58aab1d3c3a2
begin # hover interaction
	function node_hover_action(state, idx, event, axis)
		p.node_size[][idx] = state ? 20 : 10
		p.node_size[] = p.node_size[] # trigger observable
	end

	nhover = NodeHoverHandler(node_hover_action)
	register_interaction!(ax, :nhover, nhover)
	f
end

# ╔═╡ bf84d2b0-4f67-49fd-badd-6ade14c2ec70
f

# ╔═╡ e45e2c05-1497-4ac7-80ce-60fa824808e5
md"""
## Idea to build up a graph

### 1. Define the graph by the number of nodes and label for the nodes (here just the numbers)
"""

# ╔═╡ a4851f51-c596-4ef0-9fbb-f818918eccb4
new_g = SimpleDiGraph(9)

# ╔═╡ 5353737d-725f-47d4-bef6-7f56e346ce07
md"""
### 2. Plot the graph
"""

# ╔═╡ a1e3cb5c-9442-4321-962d-a277c29d7f67
GraphMakie.graphplot(new_g, 
						nlabels=repr.(1:nv(new_g)), 
						nlabels_align=(:center,:center),
						node_size = [25 for i in 1:nv(new_g)],
						node_color = [:lightblue for i in 1:nv(new_g)]
					)

# ╔═╡ e70c78e6-9599-4187-aa11-6832ce72b208
md"""
It seems, I cannot plot a graph with 0 edges with GraphRecipes.jl.
"""

# ╔═╡ 77da05e2-c23b-4a22-bc8d-2afebd83e4bf
md"""
### 3. Add the edges one after the other, starting from the injector
"""

# ╔═╡ 83ac58c9-1772-4e2b-8792-d0c03bcb68a5
add_edge!(new_g, 1, 2) # Inj->GC1->pressure regulator 

# ╔═╡ 5f8f507c-e408-4d9e-8d26-435b21993ffe
GraphMakie.graphplot(new_g, 
						nlabels=repr.(1:nv(new_g)), 
						nlabels_align=(:center,:center),
						node_size = [25 for i in 1:nv(new_g)],
						node_color = [:lightblue for i in 1:nv(new_g)]
					)

# ╔═╡ 05892523-78bd-416c-abb6-190d3cee75ed
add_edge!(new_g, 2, 3) # pressure regulator -> GC2 -> split

# ╔═╡ 946066e9-992f-4dc5-b8df-d500843e3c2f
GraphMakie.graphplot(new_g, 
						nlabels=repr.(1:nv(new_g)), 
						nlabels_align=(:center,:center),
						node_size = [25 for i in 1:nv(new_g)],
						node_color = [:lightblue for i in 1:nv(new_g)]
					)

# ╔═╡ 7483bc49-c071-4691-a390-48010001ac5f
add_edge!(new_g, 3, 4) # Split -> TL1 -> Det1

# ╔═╡ f544661e-0c3a-48c3-8045-93cf729d2db1
GraphMakie.graphplot(new_g, 
						nlabels=repr.(1:nv(new_g)), 
						nlabels_align=(:center,:center),
						node_size = [25 for i in 1:nv(new_g)],
						node_color = [:lightblue for i in 1:nv(new_g)]
					)

# ╔═╡ af5862ec-a8f8-4e70-85de-803d03a59c8a
add_edge!(new_g, 3, 5) # Split -> TL2 -> Modulator inlet

# ╔═╡ 5f8253d9-dc39-46cd-bd0b-768fd769249f
GraphMakie.graphplot(new_g, 
						nlabels=repr.(1:nv(new_g)), 
						nlabels_align=(:center,:center),
						node_size = [25 for i in 1:nv(new_g)],
						node_color = [:lightblue for i in 1:nv(new_g)]
					)

# ╔═╡ 3be14798-252a-4bdf-b275-a95fb01efc02
add_edge!(new_g, 5, 6) # Modulator inlet -> Modulator -> GC3 inlet

# ╔═╡ 003468a6-a16d-4c04-9fb0-e656a489df9c
GraphMakie.graphplot(new_g, 
						nlabels=repr.(1:nv(new_g)), 
						nlabels_align=(:center,:center),
						node_size = [25 for i in 1:nv(new_g)],
						node_color = [:lightblue for i in 1:nv(new_g)]
					)

# ╔═╡ 281807fd-1053-4d63-94b2-9748705396c8
add_edge!(new_g, 6, 7) # GC3 inlet -> GC3 -> Split

# ╔═╡ dcda3695-f93a-4a91-a14d-0ed2de88e2e3
GraphMakie.graphplot(new_g, 
						nlabels=repr.(1:nv(new_g)), 
						nlabels_align=(:center,:center),
						node_size = [25 for i in 1:nv(new_g)],
						node_color = [:lightblue for i in 1:nv(new_g)]
					)

# ╔═╡ f3cb6ba0-7ded-4912-9163-2984b3860832
add_edge!(new_g, 7, 8) # Split -> TL3 -> Det2

# ╔═╡ 3a4b7809-29a5-4128-8819-dd16ebeb627f
GraphMakie.graphplot(new_g, 
						nlabels=repr.(1:nv(new_g)), 
						nlabels_align=(:center,:center),
						node_size = [25 for i in 1:nv(new_g)],
						node_color = [:lightblue for i in 1:nv(new_g)]
					)

# ╔═╡ 7c9c3609-a10c-460a-b345-869f7651fa62
add_edge!(new_g, 7, 9) # Split -> TL4 -> Det3

# ╔═╡ 943df79a-10a2-44c7-9c2b-fdbebfc467d5
GraphMakie.graphplot(new_g, 
						nlabels=repr.(1:nv(new_g)), 
						nlabels_align=(:center,:center),
						node_size = [25 for i in 1:nv(new_g)],
						node_color = [:lightblue for i in 1:nv(new_g)]
					)

# ╔═╡ 4edb80ac-068c-4b5d-ac4d-2c365d1a36cc
md"""
### 4. Add Edge labels

respectivly add the meta-data for the edges (column lengths, diameters, stationary phase, temperatures/temperature programs)
"""

# ╔═╡ 2099f8e5-d37c-4e35-9de0-5ca7a13acd8f
edge_labels = ["GC₁", "GC₂", "TL₁", "TL₂", "Modulator", "GC₃", "TL₃", "TL₄"]

# ╔═╡ 456220ed-7cbb-43ba-933c-688662bf21bb
GraphMakie.graphplot(new_g, 
						nlabels=repr.(1:nv(new_g)), 
						nlabels_align=(:center,:center),
						node_size = [25 for i in 1:nv(new_g)],
						node_color = [:lightblue for i in 1:nv(new_g)],
						elabels = edge_labels
					)

# ╔═╡ 6b052008-d4ad-4eba-9686-fcae9fff582b
md"""
### 5. Change Node labels to the pressure points
"""

# ╔═╡ c9f3413d-5802-4f89-ae8d-3bb74c23c762
node_labels = ["p₁", "p₂", "p₃", "p₄", "p₅", "p₆", "p₇", "p₈", "p₉"]

# ╔═╡ 595763b4-7518-42fd-a6f1-6d90bf36533e
md"""
### 6. Adjust the layout of the graph plot
"""

# ╔═╡ 56d0ec78-1a7d-4197-a74c-a02fa62922fc
lay = Stress()

# ╔═╡ 5ae20243-f3ff-4256-a8d2-56bdbf673deb
begin
new_f, new_ax, new_p = GraphMakie.graphplot(new_g, 
						layout=lay,
						nlabels=node_labels, 
						nlabels_align=(:center,:center),
						node_size = [40 for i in 1:nv(new_g)],
						node_color = [:lightblue for i in 1:nv(new_g)],
						elabels = edge_labels
					)
	hidedecorations!(new_ax)
	hidespines!(new_ax)
	new_ax.aspect = DataAspect()
	new_f
end

# ╔═╡ 63a0e7aa-7885-44e1-a907-48e7a1acd9c0
#Makie.save("graph.svg", new_f)

# ╔═╡ 55cfb28a-9b37-4e49-b580-317bda8d3c3d
GraphRecipes.graphplot(new_g; curves=false)

# ╔═╡ c35a470c-95eb-4577-9a20-795b13018134
md"""
## Do something with the constructed graph
"""

# ╔═╡ 80884498-94a6-4e46-bfb8-fb38085aa30d
md"""
### Adjacency Matrix
"""

# ╔═╡ 6ab8d15c-5689-4ad8-bd66-95746971f4af
A = adjacency_matrix(new_g) 

# ╔═╡ 8fb59d0a-95c1-4d7b-8e86-973d6ba65dc4
begin
	rowsum = Array{Int}(undef, size(A)[1])
	for i=1:size(A)[1]
		rowsum[i] = sum(A[i,:])
	end
	rowsum
end

# ╔═╡ 833112a9-a09e-4a02-bfb4-d3f0261b8431
begin
	colsum = Array{Int}(undef, size(A)[2])
	for i=1:size(A)[2]
		colsum[i] = sum(A[:,i])
	end
	colsum
end

# ╔═╡ e27f3b2d-373d-4846-812f-9ddf3a93d32a
md"""
#### Observations
- `rowsum[j]` -> number of outflow from the node `j`
- `colsum[j]` -> number of inflow into node `j`


- `rowsum[j] = 0` for outlets, e.g. detectors
- `rowsum[j] > 1` split point
- `colsum[j] = 1` inlet (primary)
"""

# ╔═╡ 2c5538c6-fbc4-43db-ac28-63a51ee84810
md"""
### Rules

- the pressures of outer nodes (nodes with only one in-/outflow, `degree=1`) must be known
- number of inner nodes (nodes with at least one inflow and one outflow `degree>1`) is the number of flow balance equations

to setup the flow balance equations:
- take the edges connected to an inner node `j`
- flow of the edges where the node `j` is the destination are positive
- flow of the edges where the node `j` is the source are negative
- the sum of the flows in a node `j` is equal to Zero
- e.g. `E=[(1,2), (2,3)]` => `F₁₂ - F₂₃ = 0`
- the flow `Fᵢⱼ ∝ (pᵢ²-pⱼ²)/κᵢⱼ`, where `κᵢⱼ` the flow resistance of the edge `(i,j)` (or the module, e.g. GC column)
"""

# ╔═╡ d3a4189a-4134-44ec-baba-8029e1656025
Graphs.degree(new_g) # number of neighbors

# ╔═╡ 3fec9da2-0975-43d6-aaf3-feb13c7632f7
neighbors(new_g, 3)

# ╔═╡ 1878b77e-892c-4e41-85e1-2069cde40210
all_neighbors(new_g, 3)

# ╔═╡ fedac331-6e41-447e-a374-9fec674963d7
E = collect(edges(new_g))

# ╔═╡ ad8ee901-8daf-4111-8aee-3609d4acda1a
src(E[1]), dst(E[1])

# ╔═╡ 02f967c4-c2d3-43bd-8192-3e1f00c6ec1e
V = collect(vertices(new_g))

# ╔═╡ 590cefc6-3088-40b7-ab58-94f7f81672a9
md"""
# End
"""

# ╔═╡ Cell order:
# ╠═e7111e94-1fd0-11ed-3423-c1975e0240ea
# ╠═00cf5a16-8a1e-431d-a40e-ef620b1c1508
# ╠═d53a1637-9705-4d32-9738-2521fc4119f8
# ╠═1da8f8d8-2c51-41b1-b14a-18deb6e0f50a
# ╠═06d4b0b5-9280-43d5-a51c-253266e1d989
# ╠═96f6b288-ad6a-44b1-bb42-a63e9560e0cb
# ╠═1d03023b-176c-4057-b92d-94257b052a8d
# ╠═ea9453a4-3171-43a4-a297-3f5378ad4370
# ╠═a9d942c5-05d9-439e-b651-59a9a27b30fc
# ╠═13fff928-4655-4e6f-9d66-ad21cf66261e
# ╠═dcbdbc38-e338-48bf-a1ca-78a3669b0ed0
# ╠═509e3f76-5afa-4c13-8449-58aab1d3c3a2
# ╠═bf84d2b0-4f67-49fd-badd-6ade14c2ec70
# ╠═e45e2c05-1497-4ac7-80ce-60fa824808e5
# ╠═a4851f51-c596-4ef0-9fbb-f818918eccb4
# ╠═5353737d-725f-47d4-bef6-7f56e346ce07
# ╠═a1e3cb5c-9442-4321-962d-a277c29d7f67
# ╠═e70c78e6-9599-4187-aa11-6832ce72b208
# ╠═77da05e2-c23b-4a22-bc8d-2afebd83e4bf
# ╠═83ac58c9-1772-4e2b-8792-d0c03bcb68a5
# ╠═5f8f507c-e408-4d9e-8d26-435b21993ffe
# ╠═05892523-78bd-416c-abb6-190d3cee75ed
# ╠═946066e9-992f-4dc5-b8df-d500843e3c2f
# ╠═7483bc49-c071-4691-a390-48010001ac5f
# ╠═f544661e-0c3a-48c3-8045-93cf729d2db1
# ╠═af5862ec-a8f8-4e70-85de-803d03a59c8a
# ╠═5f8253d9-dc39-46cd-bd0b-768fd769249f
# ╠═3be14798-252a-4bdf-b275-a95fb01efc02
# ╠═003468a6-a16d-4c04-9fb0-e656a489df9c
# ╠═281807fd-1053-4d63-94b2-9748705396c8
# ╠═dcda3695-f93a-4a91-a14d-0ed2de88e2e3
# ╠═f3cb6ba0-7ded-4912-9163-2984b3860832
# ╠═3a4b7809-29a5-4128-8819-dd16ebeb627f
# ╠═7c9c3609-a10c-460a-b345-869f7651fa62
# ╠═943df79a-10a2-44c7-9c2b-fdbebfc467d5
# ╠═4edb80ac-068c-4b5d-ac4d-2c365d1a36cc
# ╠═2099f8e5-d37c-4e35-9de0-5ca7a13acd8f
# ╠═456220ed-7cbb-43ba-933c-688662bf21bb
# ╠═6b052008-d4ad-4eba-9686-fcae9fff582b
# ╠═c9f3413d-5802-4f89-ae8d-3bb74c23c762
# ╠═595763b4-7518-42fd-a6f1-6d90bf36533e
# ╠═56d0ec78-1a7d-4197-a74c-a02fa62922fc
# ╠═5ae20243-f3ff-4256-a8d2-56bdbf673deb
# ╠═63a0e7aa-7885-44e1-a907-48e7a1acd9c0
# ╠═55cfb28a-9b37-4e49-b580-317bda8d3c3d
# ╠═c35a470c-95eb-4577-9a20-795b13018134
# ╠═80884498-94a6-4e46-bfb8-fb38085aa30d
# ╠═6ab8d15c-5689-4ad8-bd66-95746971f4af
# ╠═8fb59d0a-95c1-4d7b-8e86-973d6ba65dc4
# ╠═833112a9-a09e-4a02-bfb4-d3f0261b8431
# ╟─e27f3b2d-373d-4846-812f-9ddf3a93d32a
# ╠═2c5538c6-fbc4-43db-ac28-63a51ee84810
# ╠═d3a4189a-4134-44ec-baba-8029e1656025
# ╠═3fec9da2-0975-43d6-aaf3-feb13c7632f7
# ╠═1878b77e-892c-4e41-85e1-2069cde40210
# ╠═fedac331-6e41-447e-a374-9fec674963d7
# ╠═ad8ee901-8daf-4111-8aee-3609d4acda1a
# ╠═02f967c4-c2d3-43bd-8192-3e1f00c6ec1e
# ╠═590cefc6-3088-40b7-ab58-94f7f81672a9
