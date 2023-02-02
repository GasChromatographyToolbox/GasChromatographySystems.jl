using Graphs
using GLMakie
using GraphMakie
using GLMakie.Colors

# example from the GraphMakie documentation
g = wheel_graph(10)
f, ax, p = graphplot(g,
                     edge_width = [2.0 for i in 1:ne(g)],
                     edge_color = [:gray for i in 1:ne(g)],
                     node_size = [10 for i in 1:nv(g)],
                     node_color = [:red for i in 1:nv(g)])
hidedecorations!(ax); hidespines!(ax)
ax.aspect = DataAspect()

# to use drag interactions, disable the default rectangle zoom interaction
deregister_interaction!(ax, :rectanglezoom)

# hover interactions
function node_hover_action(state, idx, event, axis)
    p.node_size[][idx] = state ? 20 : 10
    p.node_size[] = p.node_size[] # trigger observable
end
nhover = NodeHoverHandler(node_hover_action)
register_interaction!(ax, :nhover, nhover)

function edge_hover_action(state, idx, event, axis)
    p.edge_width[][idx]= state ? 5.0 : 2.0
    p.edge_width[] = p.edge_width[] # trigger observable
end
ehover = EdgeHoverHandler(edge_hover_action)
register_interaction!(ax, :ehover, ehover)

# click interactions, click changes the color randomly
#function node_click_action(idx, args...)
#    p.node_color[][idx] = rand(RGB)
#    p.node_color[] = p.node_color[]
#end
#nclick = NodeClickHandler(node_click_action)
#register_interaction!(ax, :nclick, nclick)

#function edge_click_action(idx, args...)
#    p.edge_color[][idx] = rand(RGB)
#    p.edge_color[] = p.edge_color[]
#end
#eclick = EdgeClickHandler(edge_click_action)
#register_interaction!(ax, :eclick, eclick)

# drag interaction
function node_drag_action(state, idx, event, axis)
    p[:node_pos][][idx] = event.data
    p[:node_pos][] = p[:node_pos][]
end
ndrag = NodeDragHandler(node_drag_action)
register_interaction!(ax, :ndrag, ndrag)

#---------
using GLMakie
using GraphMakie
using Graphs
g = wheel_graph(10)
f, ax, p = graphplot(g, edge_width=[3 for i in 1:ne(g)],
                     node_size=[10 for i in 1:nv(g)])

deregister_interaction!(ax, :rectanglezoom)
register_interaction!(ax, :nhover, NodeHoverHighlight(p))
register_interaction!(ax, :ehover, EdgeHoverHighlight(p))
register_interaction!(ax, :ndrag, NodeDrag(p))
register_interaction!(ax, :edrag, EdgeDrag(p))


#------------
# from Makie documentation
# draw a line from the point where left mouse button is clicked to the point where the left mouse button is released
using GLMakie

points = Observable(Point2f[])

scene = Scene(camera = campixel!)
linesegments!(scene, points, color = :black)
scatter!(scene, points, color = :gray)

on(events(scene).mousebutton) do event
    if event.button == Mouse.left
        if event.action == Mouse.press || event.action == Mouse.release
            mp = events(scene).mouseposition[]
            push!(points[], mp)
            notify(points)
        end
    end
end

scene
# this could perhaps be used to create an edge from the node next to the first point to the node next to second point
# deletion of an edge could be done by clicking on an edge (possibly in combination with a certain key)