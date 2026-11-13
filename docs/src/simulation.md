# Simulation

Path simulation walks the capillaries of a GC `System` (graph) along chosen paths. Each chromatographic segment is solved with [GasChromatographySimulator.jl](https://github.com/GasChromatographyToolbox/GasChromatographySimulator.jl). After a segment finishes, its peak list supplies the initials for the next segment on that path.

Requires **GasChromatographySimulator ≥ 0.6.2** (temporal `Options.domain`, `peaklist_temporal`, `profile` for ``d``/``df``).

## Typical workflow

1. **Define the system** — graph, `PressurePoint`s, `ModuleColumn` / `ModuleTM`, options (see [Definition of a GC system](@ref)).
2. **Resolve flows** — solve the pressure/flow balance to obtain pressures (and thus programs) on every edge ([Flow calculation](@ref)).
3. **Build segment parameters** — `graph_to_parameters(sys, p2fun, db_dataframe, selected_solutes; …)` produces a `Parameters` object per edge (column, program, substances, options). Module `domain` is copied into GasChromatographySimulator `Options.domain`.
4. **Simulate along paths** — `simulate_along_paths(sys, p2fun, paths, …)` walks each path:
   - first segment: `simulate_ModuleColumn(segment_par, t₀, τ₀)` (or a `ModuleTM` equivalent);
   - later segments: hand off the previous peak list, e.g. `simulate_ModuleColumn(segment_par, prev_peaklist)`.

The handoff helper is `change_initial(par, prev_pl)`, which sets ``t_0 = t_R`` and ``\tau_0 = \tau_R`` from the upstream peak list. For focussed reinjection after a thermal modulator, use `change_initial_focussed` / `refocus` with an explicit ``\tau_0`` (often zeros).

## Spatial vs temporal column edges

| `ModuleColumnOptions.domain` | ODE | Typical use |
|---|---|---|
| `:spatial` (default) | ``t(z)``, ``\tau^2(z)`` | Ordinary forward elution |
| `:temporal` | ``x(t)``, ``\sigma^2(t)`` | Stop / reverse / trap-style load–flush on one edge |

Both go through GasChromatographySimulator `simulate`. Temporal exits use `peaklist_temporal` (``t_R``, ``\tau_R``, ``\sigma_R``, ``exit_end``). Reverse flow is allowed on `:temporal` edges when `pin` and `pout` cross (`path_possible` / `allows_reverse_flow`); `:spatial` edges keep the forward-flow check. Vacuum-like `pout/pin < 0.1` falls back to spatial.

Same-capillary phase continuation (band still inside the tube) uses `change_initial_spatial` with ``(t_0, x_0, \sigma_0)``. Demo: `notebooks/temporal_ModuleColumn.jl`.

## Handoff between segments

For ordinary column → column connections:

- ``t_0`` of the next segment = ``t_R`` of the previous
- ``\tau_0`` of the next segment = ``\tau_R`` of the previous

Programs share one absolute clock. Details, exceptions (same capillary / trap phases, open-end release, modulator refocus), and the map into time- vs position-domain solvers are under [Segment handoff](@ref).

## Flow and pressures

Before migration, inlet and outlet pressures on every edge come from the flow calculator (graph balance). See [Flow calculation](@ref) and [Definition of a GC system](@ref).
