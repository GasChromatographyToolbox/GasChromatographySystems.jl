# Definition of a GC system

In this segment the structures and function, which are used to define a complex GC system, will be presented.

The structure of a complex GC system consists of four parts:

1. Representation of the GC system as a simple directional graph.
2. Pressure points defining the pressures at the vertices of the graph (connections of the capillaries).
3. Modules defining the features of the capillaries, e.g. dimensions and temperature.
4. Additional options, e.g. type of the mobile phase gas.

These parts are collected in `GasChromatographySystems.System(g, pressurepoints, modules, options)`.

## GC system as a graph

A complex GC system consisting of multiple connected capillaries can be abstracted and represented by a [Graph](https://en.wikipedia.org/wiki/Graph_(discrete_mathematics)). A Graph ``G(V,E)`` is mathematical construct consisting of a set of vertices ``V`` and a set of edges ``E``, used to represent connections between multiple objects. In the case of the complex GC system the vertices, also called nodes, represent the connections of the capillaries and are related to the inlet and outlet pressures of the capillaries. The edges are the capillaries. A simple directed graph is used. In this context, "simple" means, that each possible connection between vertices exists only once and "directed" means, that edges have direction going from one vertex ``i`` to vertex ``j`` is different from an edge connecting vertices ``j`` and ``i``. This is used to represent the preferred flow direction through the capillaries.

The Julia package [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) is used for the implementation of the simple directed graphs.

In the structure `GasChromatographySystems.System(g, pressurepoints, modules, options)` the graph is defined in `g`. For example the graph representing 3 capillaries in series is constructed by

First, a simple directed graph with 4 vertices is defined:

```@example ex_3series
    using GasChromatographySystems
    g = SimpleDiGraph(4)
    nothing # hide
```

In following steps the vertices are connected by the edges:

```@example ex_3series
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 3, 4)
    nothing # hide
```

This graph can be visualized:

```@example ex_3series
    using CairoMakie # hide
    f = GasChromatographySystems.plot_graph(g, ["$(i)=>$(i+1)" for i=1:3], ["p$(i)" for i=1:4])
    CairoMakie.save("plot_graph_3series.svg", f) #hide
    nothing # hide
```

![](plot_graph_3series.svg)

## Pressure points

The vertices ``V`` represent the connections of the capillaries of the GC system. To fully define the flow of the mobile phase through the GC system the pressure at these connections must be known or calculated. To store this information about the pressures, two structures are defined. One to define the vertex as pressure point and one to store information about a pressure program.

The first structure is `GasChromatographySystems.PressurePoint(name, PP)`. It is composed of the name of the pressure point resp. vertex and the value of the pressure, either as a number or as a pressure program. The pressure program is defined by the second structure of `GasChromatographySystems.PressureProgram(time_steps, pressure_steps)`, to store the `pressure_steps` at the times of the corresponding `time_steps`. In between the `time_steps` the pressure values are linearly interpolated.

For our example of three capillaries in series we define a pressure program for the inlet pressure, with a pressure of 150000 Pa held for the first two minutes and than increasing the pressure over a time span of ten minutes to 230000 Pa and hold it there for a further minute:

```@example ex_3series
    PPin = GasChromatographySystems.PressureProgram([0.0, 120.0, 600.0, 60.0], [150000.0, 150000.0, 230000.0, 230000.0])
nothing # hide
```

An array of type `GasChromatographySystems.PressurePoint` with the length of number of vertices of the graph for our example is created with the defined pressure program for the inlet, a constant pressure of 101300 Pa for the outlet. The pressures in between are unknown and `NaN` (not a number) will be assigned to these pressure points.

```@example ex_3series
    pp = Array{GasChromatographySystems.PressurePoint}(undef, nv(g))
	pp[1] = GasChromatographySystems.PressurePoint("p₁", PPin) # inlet 
	pp[2] = GasChromatographySystems.PressurePoint("p₂", NaN)
	pp[3] = GasChromatographySystems.PressurePoint("p₃", NaN) 
	pp[4] = GasChromatographySystems.PressurePoint("p₄", 101300.0) #outlet 
    nothing # hide
```

## Modules

Edges of the graph are **modules** — typically a chromatographic capillary (`ModuleColumn`) or a thermal modulator segment (`ModuleTM`). Arrays of modules must match `ne(g)`.

### `ModuleColumn`

A column edge stores length ``L``, diameter ``d``, film thickness ``df``, stationary phase name ``sp``, temperature ``T`` (constant or `TemperatureProgram`), optional flow ``F``, and `ModuleColumnOptions`.

**Uniform vs profiled ``d`` / ``df``.** Pass a `Number` for a constant value, or a function of position ``x`` (m) for a non-uniform profile. Convenience constructors set `a_d` / `a_df` to `[value]` for numbers and to empty vectors for callables. Evaluate either form with [`GasChromatographySimulator.profile`](https://GasChromatographyToolbox.github.io/GasChromatographySimulator.jl/stable/). Variable diameter needs `ng=false` so flow restriction integrates along ``x``.

```@example ex_3series
    # Uniform film
    col_u = GasChromatographySystems.ModuleColumn(
        "1 → 2", 10.0, 0.25e-3, 0.25e-6, "Rxi5ms",
        GasChromatographySystems.default_TP(),
    )
    # Non-uniform film (example linear taper)
    df_fun(x) = 0.25e-6 * (1 + 0.5 * x / 10.0)
    col_df = GasChromatographySystems.ModuleColumn(
        "trap body", 10.0, 0.25e-3, df_fun, "Rxi5ms",
        GasChromatographySystems.default_TP();
        ng=false,
    )
    nothing # hide
```

**Migration domain.** `ModuleColumnOptions(domain=:spatial)` (default) integrates ``t(z)``. Set `domain=:temporal` to integrate ``x(t)`` via GasChromatographySimulator (reversal / stop / trap-style edges). Vacuum-like outlets (`pout/pin < 0.1`) stay spatial. See [Simulation](@ref) and [Segment handoff](@ref).

```@example ex_3series
    opt_td = GasChromatographySystems.ModuleColumnOptions(domain=:temporal)
    col_td = GasChromatographySystems.ModuleColumn(
        "trap", 0.5, 0.25e-3, 0.25e-6, "Rxi5ms",
        GasChromatographySystems.TemperatureProgram([0.0, 60.0], [40.0, 40.0]),
        opt_td,
    )
    nothing # hide
```

For the three-capillary series graph, one typical assignment is:

```@example ex_3series
    modules = Array{GasChromatographySystems.AbstractModule}(undef, ne(g))
    modules[1] = GasChromatographySystems.ModuleColumn("1 → 2", 10.0, 0.25e-3, 0.25e-6, "Rxi5ms", GasChromatographySystems.default_TP())
    modules[2] = GasChromatographySystems.ModuleColumn("2 → 3", 5.0, 0.25e-3, 0.25e-6, "Rxi5ms", GasChromatographySystems.default_TP())
    modules[3] = GasChromatographySystems.ModuleColumn("3 → 4", 2.0, 0.25e-3, 0.1e-6, "Wax", GasChromatographySystems.default_TP())
    nothing # hide
```

Helper constructors such as `SeriesSystem`, `SplitSystem`, and `GCxGC_TM` build common graphs for you (see [Examples](@ref)).

### `ModuleTM`

Thermal modulator edges add modulation period ``PM``, cold/hot duty, and jet temperatures on top of a short column segment. Options live in `ModuleTMOptions`. Path simulation uses the modulator path in `simulate_along_paths` (slicing / simplified TM), not `domain=:temporal`.

## Additional options

System-wide defaults are stored in `GasChromatographySystems.Options` (mobile-phase gas, viscosity model, pressure vs flow control, ODE system flag, retention-factor threshold, …). Module-level solver tolerances and `domain` stay on `ModuleColumnOptions` / `ModuleTMOptions`.

```@example ex_3series
    opt = GasChromatographySystems.Options(gas="He", control="Pressure")
    sys = GasChromatographySystems.System("series3", g, pp, modules, opt)
    nothing # hide
```

After construction, `update_system(sys)` synchronizes temperature and pressure programs onto common time steps when any module or pressure point uses a program.

