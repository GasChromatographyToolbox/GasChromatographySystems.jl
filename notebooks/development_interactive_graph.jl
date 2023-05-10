### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ 6fa4b912-e4d3-11ed-31fb-c1dd2964b917
begin
	import Pkg
    # activate the shared project environment
	Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
	#Pkg.upgrade_manifest()
	#Pkg.precompile()
    Pkg.instantiate()
	using JSServe
	Page()
end

# ╔═╡ d47b602d-091d-40f0-a1b2-87c4388014d3
begin
	using CSV, DataFrames
	using Plots, GLMakie, GraphMakie
	using Graphs, NetworkLayout, Symbolics
	using GasChromatographySimulator
	using PlutoUI
	include(joinpath(dirname(pwd()), "src", "GasChromatographySystems.jl"))
	TableOfContents()
end

# ╔═╡ 57aa64ff-5de1-4110-8a66-faa38031a467
md"""
The following code is from: https://discourse.julialang.org/t/interactive-network-visualization/49054/3
"""

# ╔═╡ 97e9f95d-25d7-442e-a9a4-d5916e07e83e
function networkplot(g)
    fig, ax, p = graphplot(g)
    
    hidedecorations!(ax)
    hidespines!(ax)

    function node_drag_action(state, idx, event, axis)
        p[:node_pos][][idx] = event.data
        p[:node_pos][] = p[:node_pos][]
    end

    ndrag = NodeDragHandler(node_drag_action)

    deregister_interaction!(ax, :rectanglezoom)
    register_interaction!(ax, :ndrag, ndrag)

    fig
end

# ╔═╡ d4cbd3cd-9acc-4412-a939-fbff03c76efb
g = erdos_renyi(10, 20)

# ╔═╡ 2752a07d-278a-45ac-a868-462d7b71fb77
networkplot(g)

# ╔═╡ Cell order:
# ╠═6fa4b912-e4d3-11ed-31fb-c1dd2964b917
# ╠═d47b602d-091d-40f0-a1b2-87c4388014d3
# ╠═57aa64ff-5de1-4110-8a66-faa38031a467
# ╠═97e9f95d-25d7-442e-a9a4-d5916e07e83e
# ╠═d4cbd3cd-9acc-4412-a939-fbff03c76efb
# ╠═2752a07d-278a-45ac-a868-462d7b71fb77
