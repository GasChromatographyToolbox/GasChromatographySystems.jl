# Examples

## Built-in system helpers

Convenience constructors build common graphs:

- `SeriesSystem` — capillaries in series
- `SplitSystem` — one inlet split into two outlets
- `GCxGC_TM` / related helpers — comprehensive 2D GC with thermal modulator segments

```julia
using GasChromatographySystems
sys = SeriesSystem(sps = ["SLB5ms", "Wax"])
```

## Temporal `ModuleColumn` (reverse-capable)

Set `ModuleColumnOptions(domain=:temporal)` (or pass `domain=:temporal` into helpers that forward module options). Forward elution matches the spatial solver within a few percent; reverse flush uses `pin < pout` and continues from a parked ``(x,\sigma)`` with `change_initial_spatial`.

See the Pluto/notebook demo:

- [`notebooks/temporal_ModuleColumn.jl`](https://github.com/GasChromatographyToolbox/GasChromatographySystems.jl/blob/main/notebooks/temporal_ModuleColumn.jl)

Tests covering the same path live in `test/runtests.jl` (`temporal ModuleColumn` and `ModuleColumn functional d and df`).

## Flow-calculator demos

Notebooks under `notebooks/FlowCalcPaper/` illustrate multi-capillary pressure/flow balance (research demos; not Documenter `@example` blocks):

- [`4-Way-Splitter_Demo.jl`](https://github.com/GasChromatographyToolbox/GasChromatographySystems.jl/blob/main/notebooks/FlowCalcPaper/4-Way-Splitter_Demo.jl)
- [`DeansSwitch_Demo.jl`](https://github.com/GasChromatographyToolbox/GasChromatographySystems.jl/blob/main/notebooks/FlowCalcPaper/DeansSwitch_Demo.jl)

## GCxGC notebooks

Thermal-modulation / comprehensive 2D GC demos live under `notebooks/GCxGC-TM_Paper/` (uses `GCxGC_TM` and related helpers):

- [`GCxGC_TM_Demo.jl`](https://github.com/GasChromatographyToolbox/GasChromatographySystems.jl/blob/main/notebooks/GCxGC-TM_Paper/GCxGC_TM_Demo.jl)
