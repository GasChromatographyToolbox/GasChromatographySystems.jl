### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 03885f52-0e61-11ed-1db1-87cf0916dcdd
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
	TableOfContents()
end

# ╔═╡ cc1fbc1e-a684-4e1b-b3dc-2e326b09a6f8
md"""
# Test Graphs.jl
"""

# ╔═╡ ffd2f954-0263-48a0-a172-622657b89422
gc = path_graph(6) # Defines a graph with 11 vertices

# ╔═╡ cd3223d2-0b9c-467e-8975-1a05e84d0ca9
nv(gc) # number of vertices

# ╔═╡ 27272425-9ee1-4330-a94e-06f7f787a3b4
ne(gc) # number of edges

# ╔═╡ ba7eabd8-45a4-4cef-a99a-1b220c390b25
graphplot(gc)

# ╔═╡ a23080ae-372a-4204-9571-d55e16690aa8
md"""
## Creating graphs
"""

# ╔═╡ 9c384647-1281-4ad7-bb80-cbf19fde5c62
g = SimpleGraph(2, 1) # start with 2 vertices and 1 edge

# ╔═╡ ba1e0b92-afe8-46b9-9de0-11c78c95c01c
graphplot(g)

# ╔═╡ 66c1485a-8b6e-4b6a-b2fe-fbf889d0a925
add_vertex!(g) # add one vertex

# ╔═╡ ffd75d2f-ae5b-49c6-ad50-66057b7752a6
add_edge!(g, 2, 3) # add edge from vertex 2 to vertex 3 

# ╔═╡ 65f7706e-cb30-409d-92c4-715c7e0dcf17
graphplot!(g)

# ╔═╡ 892fe528-797f-4802-ac52-60860e102b0f
add_vertices!(g, 6) # add 6 vertices

# ╔═╡ e68cb8bc-88e7-4e3d-b2b1-7528025d7339
begin
	# connect the new vertices
	add_edge!(g, 2, 4)
	add_edge!(g, 3, 5)
	add_edge!(g, 3, 6)
	add_edge!(g, 5, 8)
	add_edge!(g, 6, 9)
	add_edge!(g, 4, 7)
end

# ╔═╡ a8bcb649-14e2-4060-bc2a-3016ac88e75e
md"""
## Plotting with GraphRecipes.jl
"""

# ╔═╡ ac1d2730-0e47-4e89-88d9-fd73361ef6ec
graphplot(g; curves=false)

# ╔═╡ 0c43193b-fe7a-41fb-ab0e-90ea1b754073
names = ["PP", "TL", "GC", "GC", "TL", "TL", "PP", "PP", "PP"]

# ╔═╡ b50e5d2f-dbc3-4c37-b729-4ca225e5a515
graphplot(g; names=names, curves=false)

# ╔═╡ f37538e2-6c74-4955-bcea-e594f274f18e
x = [0, 1, 2, 2, 3, 3, 3, 4, 4]

# ╔═╡ aa29fb3c-fe64-46a4-9f2a-00a537bbe48d
y = [0, 0, 0.5, -0.5, 1, 0, -0.5, 1, 0]

# ╔═╡ 7998deab-115a-422e-a46d-28b21df9885d
graphplot(g; names=names, curves=false, x=x, y=y, nodesize=0.5, nodeshape=:circle)

# ╔═╡ 00c70d8a-dbee-4a93-a3d9-628648ec8457
md"""
## Access Graphs
"""

# ╔═╡ 04e4d769-2f93-4c2a-802c-fd10a4cd4f8e
nv(g) # number of vertices

# ╔═╡ 5080c54b-f93a-4716-8c57-f840d3ecc5a9
ne(g) # number of edges

# ╔═╡ 4ca130d3-fb2d-49e3-a427-2ad2f3087500
v = vertices(g) # returns an iterable object containing all the vertices in g

# ╔═╡ 1581841c-2845-4cc2-82ea-046fa3ad7930
collect(v)

# ╔═╡ 6dc5c6ab-45e5-4e11-9cc3-ae2a928330a4
e = edges(g) # returns an iterable object containing all the edges in g

# ╔═╡ 07e9aa3c-4078-45c0-944f-aac1ccbdc74e
collect(e)

# ╔═╡ 7834c6db-a171-4a1b-b6bc-d5fa646c68bb


# ╔═╡ 920ceeab-5be7-421c-9b83-9d7a31ecb0e9
src.(e), dst.(e)

# ╔═╡ f73c4b58-8d09-4c39-8a82-5f2a771cfa57
a_star(g, 1, 9) # shortest path from source (1) to target (9)
# could be used to define the different paths from inlet (1) to all the outlets

# ╔═╡ e60d44f5-168d-4f4f-94c3-e8bc9ec2b289
Ag = adjacency_matrix(g)

# ╔═╡ 58a12cfa-8307-4732-b391-6168a17280dd
sum(Ag[1,:]) # -> start/end point

# ╔═╡ d3efdcd1-748b-4453-9b4f-0293c5be50e9
sum(Ag[2,:]) # -> spilt point

# ╔═╡ e8bbe03a-06a2-46c8-84d1-efff78c35964
sum(Ag[4,:]) # -> straight connection

# ╔═╡ f1fae2c2-2c64-4c11-9127-ffefa4e1b447
incidence_matrix(g) # ???

# ╔═╡ a2d9563d-2abe-4842-829b-d28f7e91f5b9
Lg = laplacian_matrix(g) # ??? -> diagonal has the sum(Ag[i,:]) -> property of the vertices

# ╔═╡ 6c7bae1c-4b16-4cec-b0ed-4dc19196672e
md"""
## Properties of the GC system
"""

# ╔═╡ d628a3de-d7c6-4df8-97b2-5a4b30147e84
Ag[1,:]

# ╔═╡ 323ee5d9-be5a-4661-81b2-45bc8dd62b2a
Ag[9,:]

# ╔═╡ b41cb994-2337-4811-869a-4814b0a25951
Ag[:,9]

# ╔═╡ 825b5d6f-e395-455c-921a-aaf7a6d00cfb
a = findfirst(dst.(e).==9)

# ╔═╡ 71c5d2ba-6460-4bbc-93c0-9c9736538491
typeof(a)

# ╔═╡ f6ad8a60-33bd-4d62-b26e-6563582c2114
length(a)

# ╔═╡ 8ac80ce8-9e19-40d8-9279-16b171755231
isa(a,Int)

# ╔═╡ 8b69fcc9-5feb-48b8-9ca5-92b33bda17d3
isnothing(a)

# ╔═╡ c2462e99-a45a-4e7d-8b2c-995168fb058c
begin
	diag = Array{Int}(undef, nv(g))
	feat = Array{String}(undef, nv(g))
	Lg_ = laplacian_matrix(g)
	src_e = src.(e)
	dst_e = dst.(e)
	start_v = Int[]
	end_v = Int[]	
	split_v = Int[]
	connection_v = Int[]
	for i=1:nv(g)
		diag[i] = Lg_[i,i]
		if diag[i] == 1 # end points
			if isnothing(findfirst(src.(e).==i)) && isa(findfirst(dst.(e).==i), Int)# it is a end point if the correlared edge has this vertex not as source (but as destination )
				feat[i] = "end"
				push!(end_v, i)
			else
				feat[i] = "start"
				push!(start_v, i)
			end
		elseif diag[i] == 2 # straight connection
			feat[i] = "connection"
			push!(connection_v, i)
		elseif diag[i] == 3 # splitter
			feat[i] = "split"
			push!(split_v, i)
		end
	end
	diag, feat, start_v, end_v, connection_v, split_v
end

# ╔═╡ 994512a7-36f2-47bc-9796-b8a91b087a84
inneighbors(g, 1)

# ╔═╡ 30d835b3-2d29-4253-ac1c-c764d3ee0ab7
md"""
### System of Balances

At an end point only modules of type `Pressure_Point` are allowd, which have a defined pressure value ``p_i``.

At straight connections the flow from on module ``i`` (vertex/node) is the same as in the connected module ``j``:

``
F_i = F_j
``

At a split point three modules ``i``, ``j`` and ``k`` are connected and for the flow it is:

``
F_i - F_j -F_k = 0
``

The normalized volumetric flow is calculated as:

``
F = A \frac{p_{in}^2-p{out}^2}{κ_L}
``

with ``A = π/256 T_n/p_n`` a constant, the inlet pressure ``p_{in}``, the outlet pressure ``p_{out}`` and the flow restriction ``κ_L = \int_0^L ηT/d^4 dy = ηTL/d^4`` (conventional GC). Therefore every module (besides `Pressure_Point`) is defined by the three quantities ``p_{in}``, ``p_{out}`` and ``κ_L``. 

``κ_L`` is defined by the geometry of the module (``L``, ``d``) the mobile phase gas (``η``) and the temperture program (``T``, ``η``)).
"""

# ╔═╡ 9459f5f8-f20d-4e5d-a50a-1e4de8e3cbb8
md"""
The pressures are depending on the pressures applied on the end of the GC-system (the graph), which are defined by the pressures at the `Pressure Points` and are calculated by solving the system of equations defined by the balance equations of the graph.  

In the example above we should get this system of equations:

``
F_2 - F_3 - F_4 = 0
``

``
F_3 - F_5 - F_6 = 0
``

The pressures of the end points 1, 7, 8, 9 define pressures for the connected modules:

``
p_1 = p_{in,2}
``

``
p_7 = p_{out,4}
``

``
p_8 = p_{out,6}
``

``
p_9 = p_{out,5}
``

The desicion about inlet or oultlet pressure is based on the edge `(s,d)`. If the number of the module is the source `s`, than it is the inlet pressure. If the number of the module is the destination `d`, than it is the outlet pressure.

In the same way, the `Pressure_Points` can be defined as global inlet pressure if the number of the `Pressure_Point` is the source `s` of the edge and as global outlet pressure if the number is the destination `d` of the edge.

"""

# ╔═╡ 83fe322c-3252-4522-85ce-9fcdd0d5b5ed
md"""
Also, the outlet pressure of one module (if it is the source `s` in the edge), is the inlet pressure of the connected module (the destination `d` of the same edge). Therfore we have:

``
p_{out,2} = p_{in,3} = p_{in,4}
``

``
p_{out,3} = p_{in,5} = p_{in,6}
``

If we can calculate the pressures ``p_{out,2}`` and ``p_{out,3}``, which are the pressures at the spilt points, all inlet and outlet pressures of all modules are known and the flows through the whole system can be calculated.

For the system of equations, using the definition of the flow ``F``, and the relations of the pressures we get two equations for two unknowns:

``
\frac{p_1^2 - p_{out,2}^2}{κ_2} - \frac{p_{out,2}^2 - p_{out,3}^2}{κ_3} - \frac{p_{out,2}^2 - p_7^2}{κ_4} = 0
``

``
\frac{p_{out,2}^2 - p_{out,3}^2}{κ_3} - \frac{p_{out,3}^2 - p_9^2}{κ_5} - \frac{p_{out,3}^2 - p_8^2}{κ_6} = 0
``
"""

# ╔═╡ d824ad62-d8d4-4cea-bf80-f43e1eb82e8d
md"""
It should be possible, to extract thes systems of equations from the graphs representing the GC system.

With the right package, the 'automatic' solving of the system of equations should be possible.
"""

# ╔═╡ d7294fff-a6f9-4ff5-985f-513444b95973
f(pin,pout,κ) = (pin^2-pout^2)/κ

# ╔═╡ 7ff721fc-0aa5-46e8-9ed2-cc8bcad5a357
begin
	# go throug the graph and find the vertices neigboring the end points

# ╔═╡ 7e53f22e-d924-4c31-b41a-e79add9108f1
md"""
# End
"""

# ╔═╡ Cell order:
# ╠═03885f52-0e61-11ed-1db1-87cf0916dcdd
# ╠═cc1fbc1e-a684-4e1b-b3dc-2e326b09a6f8
# ╠═ffd2f954-0263-48a0-a172-622657b89422
# ╠═cd3223d2-0b9c-467e-8975-1a05e84d0ca9
# ╠═27272425-9ee1-4330-a94e-06f7f787a3b4
# ╠═ba7eabd8-45a4-4cef-a99a-1b220c390b25
# ╠═a23080ae-372a-4204-9571-d55e16690aa8
# ╠═9c384647-1281-4ad7-bb80-cbf19fde5c62
# ╠═ba1e0b92-afe8-46b9-9de0-11c78c95c01c
# ╠═66c1485a-8b6e-4b6a-b2fe-fbf889d0a925
# ╠═ffd75d2f-ae5b-49c6-ad50-66057b7752a6
# ╠═65f7706e-cb30-409d-92c4-715c7e0dcf17
# ╠═892fe528-797f-4802-ac52-60860e102b0f
# ╠═e68cb8bc-88e7-4e3d-b2b1-7528025d7339
# ╠═a8bcb649-14e2-4060-bc2a-3016ac88e75e
# ╠═ac1d2730-0e47-4e89-88d9-fd73361ef6ec
# ╠═0c43193b-fe7a-41fb-ab0e-90ea1b754073
# ╠═b50e5d2f-dbc3-4c37-b729-4ca225e5a515
# ╠═f37538e2-6c74-4955-bcea-e594f274f18e
# ╠═aa29fb3c-fe64-46a4-9f2a-00a537bbe48d
# ╠═7998deab-115a-422e-a46d-28b21df9885d
# ╠═00c70d8a-dbee-4a93-a3d9-628648ec8457
# ╠═04e4d769-2f93-4c2a-802c-fd10a4cd4f8e
# ╠═5080c54b-f93a-4716-8c57-f840d3ecc5a9
# ╠═4ca130d3-fb2d-49e3-a427-2ad2f3087500
# ╠═1581841c-2845-4cc2-82ea-046fa3ad7930
# ╠═6dc5c6ab-45e5-4e11-9cc3-ae2a928330a4
# ╠═07e9aa3c-4078-45c0-944f-aac1ccbdc74e
# ╠═7834c6db-a171-4a1b-b6bc-d5fa646c68bb
# ╠═920ceeab-5be7-421c-9b83-9d7a31ecb0e9
# ╠═f73c4b58-8d09-4c39-8a82-5f2a771cfa57
# ╠═e60d44f5-168d-4f4f-94c3-e8bc9ec2b289
# ╠═58a12cfa-8307-4732-b391-6168a17280dd
# ╠═d3efdcd1-748b-4453-9b4f-0293c5be50e9
# ╠═e8bbe03a-06a2-46c8-84d1-efff78c35964
# ╠═f1fae2c2-2c64-4c11-9127-ffefa4e1b447
# ╠═a2d9563d-2abe-4842-829b-d28f7e91f5b9
# ╠═6c7bae1c-4b16-4cec-b0ed-4dc19196672e
# ╠═d628a3de-d7c6-4df8-97b2-5a4b30147e84
# ╠═323ee5d9-be5a-4661-81b2-45bc8dd62b2a
# ╠═b41cb994-2337-4811-869a-4814b0a25951
# ╠═825b5d6f-e395-455c-921a-aaf7a6d00cfb
# ╠═71c5d2ba-6460-4bbc-93c0-9c9736538491
# ╠═f6ad8a60-33bd-4d62-b26e-6563582c2114
# ╠═8ac80ce8-9e19-40d8-9279-16b171755231
# ╠═8b69fcc9-5feb-48b8-9ca5-92b33bda17d3
# ╠═c2462e99-a45a-4e7d-8b2c-995168fb058c
# ╠═994512a7-36f2-47bc-9796-b8a91b087a84
# ╟─30d835b3-2d29-4253-ac1c-c764d3ee0ab7
# ╟─9459f5f8-f20d-4e5d-a50a-1e4de8e3cbb8
# ╟─83fe322c-3252-4522-85ce-9fcdd0d5b5ed
# ╟─d824ad62-d8d4-4cea-bf80-f43e1eb82e8d
# ╠═d7294fff-a6f9-4ff5-985f-513444b95973
# ╠═7ff721fc-0aa5-46e8-9ed2-cc8bcad5a357
# ╠═7e53f22e-d924-4c31-b41a-e79add9108f1
