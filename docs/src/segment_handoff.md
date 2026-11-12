# Segment handoff

GasChromatographySystems.jl chains capillaries on a graph: each edge is simulated with [GasChromatographySimulator.jl](https://github.com/GasChromatographyToolbox/GasChromatographySimulator.jl), and the result of one segment becomes the initial condition of the next. The rule depends on the **connection type**, not on whether a given segment is integrated in position or in time.

All modules on a path must share the **same absolute program time** so ``T(t)`` and pressures are consistent at the handoff instant.

## Default: column → column (junction)

This is what `change_initial` / `simulate_ModuleColumn` / `simulate_along_paths` do for ordinary `ModuleColumn` segments.

From the upstream peak list:

- ``t_0^\mathrm{next} = t_R^\mathrm{prev}``
- ``\tau_0^\mathrm{next} = \tau_R^\mathrm{prev}``

The next segment starts at **its own inlet**. The next column only sees a temporal injection profile; spatial width is recovered locally if needed via ``\sigma_0 = \tau_0 |u_s|`` at that inlet (time-domain solver). Do **not** copy upstream ``\sigma_R`` across a junction when ``u_s``, diameter, or film differ.

Thermal modulators (`ModuleTM`) still follow the temporal picture. With `refocus=false`, the downstream segment inherits ``\tau_R`` (via `change_initial`). With `refocus=true` (or `change_initial_focussed`), each reinjected band can start with a chosen ``\tau_0`` (often near zero) at the slice time — optional focusing, not a continuous spatial continuation of the upstream ``\sigma``.

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

## Future: trap module & bidirectional flow

Capillary-trap workflows (load → hold → reverse flush → GC column → restrictor) are only
partly covered by today’s graph path simulation. Tracked needs before a first-class
**`ModuleTrap`** (and before Experiments can drop its stand-in series balance):

1. **GasChromatographySimulator compatibility** — stay aligned with current GCSim APIs
   (time-domain solvers, `Program` / `Substance` validation, segment-handoff conventions).
   Flow calculator and `graph_to_parameters` must keep working on the GCSim line used by
   the toolbox shared environment.
2. **Bidirectional / reverse flows** — load and reverse flush reverse the preferred flow on
   the same physical capillary. Directed graph edges and ``F = A\lambda(p_i^2-p_j^2)`` need an
   explicit reverse-flow story (negative ``F``, flipped edge, or phase-dependent orientation)
   without forcing a spatial remapping of ``T(x)`` / ``d(x)`` for every phase.
3. **Entry point beyond the inlet** — Systems currently assumes each edge is entered at
   ``x=0`` with ``(t_0,\tau_0)``. Trap phases that never leave the tube must carry
   ``(x_0,\sigma_0)`` (and absolute program time) between phases. Define a handoff mode for
   “same capillary, new phase” in the Systems API, not only junction reinjection.
4. **`ModuleTrap`** — package heated zones / body programs, open-end release metrics
   (``t_{05}``–``t_{95}``, mass left), and multi-phase scheduling so trap → column is a
   supported path, with pneumatics from the flow calculator staying coupled as ``κ(t)``
   changes (series balance, not a one-shot set-point).

Until that lands, [GasChromatographyExperiments](https://github.com/JanLeppert/GasChromatographyExperiments)
rebalances trap→column→restrictor pressures over the release clock in
`HotInletTrapSystem` (`release_series=:coupled`); see that repo’s
`docs/coupled_release_pneumatics.md`.
