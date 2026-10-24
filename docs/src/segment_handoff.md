# Segment handoff

GasChromatographySystems.jl chains capillaries on a graph: each edge is simulated with [GasChromatographySimulator.jl](https://github.com/GasChromatographyToolbox/GasChromatographySimulator.jl), and the result of one segment becomes the initial condition of the next. The rule depends on the **connection type**, not on whether a given segment is integrated in position or in time.

All modules on a path must share the **same absolute program time** so ``T(t)`` and pressures are consistent at the handoff instant.

## Default: column → column (junction)

This is what `change_initial` / path simulation do for ordinary `ModuleColumn` (and analogous) segments.

From the upstream peak list:

- ``t_0^\mathrm{next} = t_R^\mathrm{prev}``
- ``\tau_0^\mathrm{next} = \tau_R^\mathrm{prev}``

The next segment starts at **its own inlet**. The next column only sees a temporal injection profile; spatial width is recovered locally if needed via ``\sigma_0 = \tau_0 |u_s|`` at that inlet (time-domain solver). Do **not** copy upstream ``\sigma_R`` across a junction when ``u_s``, diameter, or film differ.

Modulator / valve junctions that **slice** a peak into reinjection windows still follow the temporal picture: each slice gets a ``t_0`` on the open window and a ``\tau_0`` from `valve_initial_width` (`:inherit` from the upstream ``\tau_R``, or a fixed focus width). That is a junction handoff with optional refocus, not a continuous spatial continuation.

## Same capillary, changing phase

If one physical tube is solved in several program phases (e.g. trap load → hold → reverse flush) **without** leaving the capillary, carry position and spatial band width ``(x, \sigma)`` between phases. That pattern lives mainly in GasChromatographySimulator time-domain workflows; Systems path simulation assumes each graph edge is entered at the inlet with ``(t_0, \tau_0)``.

## Trap or open-end → analytical column

Treat as a **junction**: estimate the temporal release window (e.g. ``t_{50}`` and ``\tau_0`` from ``t_{05}``–``t_{95}`` of mass past the open end), then set the first analytical segment’s ``t_0`` and ``\tau_0``. Outlet spatial ``\sigma`` while the band is still exiting is not the column ``\sigma_0``.

## Solver independence

| Next segment solver | After a junction, set |
|---|---|
| Position-domain (`simulate`, …) | ``t_0``, ``\tau_0`` |
| Time-domain (`solving_odesystem_time`, …) | ``t_0``, ``x_0 = 0``, ``\sigma_0 = \tau_0\|u_s(0,t_0)\|`` |

The choice of what to transfer is fixed by the connection; only the representation of that transfer changes with the solver.

Canonical write-up with equations: [GasChromatographySimulator — Segment handoff](https://GasChromatographyToolbox.github.io/GasChromatographySimulator.jl/dev/segment_handoff/) (once published on `dev`/`stable`; source: `docs/src/segment_handoff.md` in that repository).
