# plotting of graphs and system functions

# begin - plotting of the graphs and functions
function plot_graph(g, edge_labels, node_labels; lay = Spring(), color=:lightblue, node_size=40, arrow_size=20, arrow_shift=0.8, dataaspect=false, nlabels_fontsize=10, elabels_fontsize=10)
	fig, ax, p = GraphMakie.graphplot(g, 
						layout=lay,
						nlabels=node_labels,
						nlabels_fontsize = nlabels_fontsize, 
						nlabels_align=(:center,:center),
						node_size = [node_size for i in 1:Graphs.nv(g)],
						node_color = [color for i in 1:Graphs.nv(g)],
						elabels = edge_labels,
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

function plot_graph(sys; lay = Spring(), color=:lightblue, node_size=40, arrow_size=20, arrow_shift=0.8, dataaspect=false, nlabels_fontsize=20, elabels_fontsize=20)
	fig, ax, p = GraphMakie.graphplot(sys.g, 
						layout=lay,
						nlabels = [sys.pressurepoints[i].name for i in 1:Graphs.nv(sys.g)],
						nlabels_align=(:center,:center),
						nlabels_fontsize = nlabels_fontsize,
						node_size = [node_size for i in 1:Graphs.nv(sys.g)],
						node_color = [color for i in 1:Graphs.nv(sys.g)],
						elabels = [sys.modules[i].name for i in 1:Graphs.ne(sys.g)],
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

#=
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
=#

function plot_graph_with_flow(sys, p2fun, t; lay = GasChromatographySystems.Spring(), color=:lightblue, node_size=80, arrow_size=20, arrow_shift=0.8, dataaspect=false, nlabels_fontsize=14, elabels_fontsize=14, elabels_distance = 20)
	p_func = GasChromatographySystems.pressure_functions(sys, p2fun)
	F_func = GasChromatographySystems.flow_functions(sys, p2fun)
	fig, ax, p = GasChromatographySystems.GraphMakie.graphplot(sys.g, 
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
	GasChromatographySystems.hidedecorations!(ax)
	GasChromatographySystems.hidespines!(ax)
	if dataaspect == true
		ax.aspect = DataAspect()
	end
	return fig
end

#=
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
=#

function plot_flow_over_time(sys, p2fun; dt=60.0)
	#plotly()
	F_func = GasChromatographySystems.flow_functions(sys, p2fun)
	com_timesteps = GasChromatographySystems.common_timesteps(sys)
	trange = 0:dt:sum(com_timesteps)
	p_flow = Plots.plot(xlabel="time in s", ylabel="flow in mL/min")
	for i=1:ne(sys.g)
		Plots.plot!(p_flow, trange, F_func[i].(trange).*60e6, label="F_$(sys.modules[i].name)")
	end
	return p_flow
end

#=
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
=#

function plot_pressure_over_time(sys, p2fun; dt=60.0)
	#plotly()
	p_func = GasChromatographySystems.pressure_functions(sys, p2fun)
	com_timesteps = GasChromatographySystems.common_timesteps(sys)
	trange = 0:dt:sum(com_timesteps)
	p_pres = Plots.plot(xlabel="time in s", ylabel="pressure in Pa(a)", legend=:topleft)
	for i=1:nv(sys.g)
		Plots.plot!(p_pres, trange, p_func[i].(trange), label="$(sys.pressurepoints[i].name)")
	end
	return p_pres
end

#=
function plot_holdup_time_path_over_time(sys, num_paths; dt=60.0)
	#plotly()
	tMp_func = holdup_time_path(sys, num_paths)
	com_timesteps = GasChromatographySystems.common_timesteps(sys)
	trange = 0:dt:sum(com_timesteps)
	p_tM = Plots.plot(xlabel="time in s", ylabel="hold-up time in s")
	for i=1:num_paths
		Plots.plot!(p_tM, trange, tMp_func[i].(trange), label="path: $(i)")
	end
	return p_tM
end
=#

function plot_holdup_time_path_over_time(sys, p2fun, num_paths; dt=60.0)
	#plotly()
	tMp_func = holdup_time_path(sys, p2fun, num_paths)
	com_timesteps = GasChromatographySystems.common_timesteps(sys)
	trange = 0:dt:sum(com_timesteps)
	p_tM = Plots.plot(xlabel="time in s", ylabel="hold-up time in s")
	for i=1:num_paths
		Plots.plot!(p_tM, trange, tMp_func[i].(trange), label="path: $(i)")
	end
	return p_tM
end
# end - plotting of the graphs and functions

# old? 
function plot_GCxGC(pl_GCxGC, sys; categories = String[])
	t1end = sum(GasChromatographySystems.common_timesteps(sys))
	t2end = ceil(maximum(pl_GCxGC.tR2.-pl_GCxGC.tR1))
	p_gcxgc = Plots.scatter(pl_GCxGC.tR1, pl_GCxGC.tR2.-pl_GCxGC.tR1, xlims=(0.0,t1end), ylims=(0.0, t2end), label="", xlabel=string(sys.modules[1].sp, ", t¹ in s"), ylabel=string(sys.modules[2].sp, ", t² in s"))
	if !isempty(categories)
		for i=1:length(categories) 
			pl_f = filter([:Cat] => x -> categories[i] in x, pl_GCxGC)
			Plots.scatter!(p_gcxgc, pl_f.tR1, pl_f.tR2.-pl_f.tR1, label=categories[i], legend=:bottomright)
		end
	end
	return p_gcxgc
end

