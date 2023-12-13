# Definition of a GC system

In this segment the structures and function, which are used to define a complex GC system, wil be presented.

The structure of a complex GC system consists of four parts:

1. Representation of the GC system as a simple directional graph.
2. Pressure points defining the pressures at the vertices of the graph (connections of the capillaries).
3. Modules defining the features of the capillaries, e.g. dimensions and temperature.
4. Additional options, e.g. typ of the mobile phase gas.

```@docs
GasChromatographySystems.System
```

## GC system as a graph

A complex GC system consisting of multiple connected capillaries can be abstracted and represented by a [Graph](https://en.wikipedia.org/wiki/Graph_(discrete_mathematics)). A Graph ``G(V,E)`` is mathematical construct consisting of a set of vertices ``V`` and a set of edges ``E``, used to represent connections between multiple objects. In the case of the complex GC system the vertices, also called nodes, represent the connections of the capillaries and are related to the inlet and outlet pressures of the capillaries. The edges are the capillaries. A simple directed graph is used. In this context, "simple" means, that each possible connection between vertices exists only once and "directed" means, that edges have direction going from one vertex ``i`` to vertex ``j`` is different from an edge connecting vertices ``j`` and ``i``. This is used to represent the preferred flow direction through the capillaries.

The Julia package [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) is used for the implementation of the simple directed graphs.

In the structure `GasChromatographySystems.System(g, pressurepoints, modules, options)` the graph is defined in `g`. For example the graph representing 3 capillaries in series is constructed by

First, a simple directed graph with 4 vertices is defined:

```@example ex_3series
    a = 1 + 1
    #g = Graphs.SimpleDiGraph(4)
    #nothing # hide
```

In following steps the vertices are connected by the edges:

```@example ex_3series
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 3, 4)
    #nothing # hide
```

This graph can be visualized:

```@example ex_3series
    plot_graph(g, ["$(i)=>$(i+1)" for i=1:3], ["p$(i)" for i=1:4])
    #nothing # hide
```
## Pressure points

## Modules

## Additional options

