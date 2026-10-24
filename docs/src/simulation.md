# Simulation

Path simulation walks the capillaries of a GC `System` (graph) in topological order. Each chromatographic segment is solved with [GasChromatographySimulator.jl](https://github.com/GasChromatographyToolbox/GasChromatographySimulator.jl). After a segment finishes, its peak list supplies the initials for the next segment on that path.

## Handoff between segments

For ordinary column → column connections:

- ``t_0`` of the next segment = ``t_R`` of the previous
- ``\tau_0`` of the next segment = ``\tau_R`` of the previous

Programs share one absolute clock. Details, exceptions (same capillary / trap phases, open-end release, modulator slicing), and the map into time- vs position-domain solvers are collected under [Segment handoff](@ref).

## Flow and pressures

Before migration, inlet and outlet pressures on every edge come from the flow calculator (graph balance). See [Flow calculation](@ref) and [Definition of a GC system](@ref).
