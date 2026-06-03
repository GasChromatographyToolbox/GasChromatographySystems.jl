# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Valve module core structures in `Structures.jl`: `ValveProgram`, `ModuleValveOptions`, and `ModuleValve` with docstrings and convenience constructors.
- `AbstractValveProgram` supertype; `PeriodicValveProgram` for compact periodic modulation (O(1) `valve_state`); `expand_valve_program` and `ValveProgram(mp, t_closed, t_end; ...)` expand to explicit segments when needed.
- `default_periodic_ValveProgram()` returns `PeriodicValveProgram(10.0, 2.0, 1800.0)`.
- Piecewise-constant valve state helper `valve_state(vp, t)` for open/closed switching without linear interpolation of boolean states.
- Test coverage for `ValveProgram`, `PeriodicValveProgram`, `ModuleValveOptions`, and `ModuleValve` constructors/defaults, including periodic/inverted schedules and agreement with expanded programs.
- Regression tests for `common_timesteps`, `match_programs`, and `update_system` with mismatched column, valve, and pressure program grids (tee fixture; includes post-balance `flow_functions` on valve edge).
- `index_modules_with_valve_program(sys)` — edge indices for `ModuleValve` with `AbstractValveProgram` `state` (not used by `match_programs`; helper for future tooling).
- Phase 5.3 (partial): `is_simulation_segment`, `path_is_chromatographic`, chromatographic `all_paths(g, modules)` / `all_paths(sys)` (exclude paths through `ModuleValve`); tests in `Chromatographic paths (exclude ModuleValve)`.
- `graph_to_parameters` placeholder `Parameters` on valve edges (`sp = ""`, `d_open`, default solver options); `all_stationary_phases` skips modules without `sp`.
- Docstrings and typed signatures for `index_parameter`, `common_edges`, `positive_flow`, and `path_possible` in `SolvingSystems.jl`.
- Valve junction transport (`ValveJunction.jl`): `incident_valve_modules`, `slice_peaks_by_valve`, `simulate_valve_junction`, `apply_valve_junctions_at_vertex`, `select_valve_initial_width`; `simulate_along_paths` splits peaks at path vertices with time-varying incident valves (before downstream column/TM).
- `ValveInitialWidth` type alias and `ModuleValveOptions.valve_initial_width` — `:inherit` (default) or `(:fixed, width_s)` for per-slice initial peak width after the junction (simulator `τ₀`; `0.0` = sharp band, TM `refocus` analogue).
- `slicing` optional `ann_prefix` keyword (default `"s"`; valve slices use `"v"`).
- Tests: `Valve junction slicing` (area conservation, slice ordering, `valve_initial_width` / `select_valve_initial_width`, sharp-band `change_initial`); `change_initial finite-row guard`.
- `GCxGC_DPM` builder: default `Vern9()`, `abstol=1e-10`, `reltol=1e-8`, `opt_valve` with `valve_initial_width=(:fixed, 0.0)`; docstring notes on `graph_to_parameters` `dt` and verifying `p₂(t)`.

### Changed
- `graph_to_parameters` now enforces `GasChromatographySimulator.Options(control="Pressure")` (with warning when `sys.options.control != "Pressure"`), to match pressure-balanced network simulations (`solve_balance`/`build_pressure_squared_functions`) and avoid accidental over-driving with flow-control.
- Near-zero-flow diagnostics in `simulate_along_paths`: logs path/segment/module context for non-finite peak rows and downstream handoff failures (`change_initial`).
- Valve junction API: `valve_slicing_schedule` returns `(mp, t_closed, phase_shift)` instead of TM-style `(PM, ratio, shift)`; `slice_peaks_by_valve` / `simplified_valve_junction` use valve phase names; added `t_start_next_open_window`.
- CI: upgraded Codecov upload to `codecov/codecov-action@v5` with `files: lcov.info` and `CODECOV_TOKEN` (replaces deprecated v1 uploader and `CODECOV_SECRET`).
- README: fixed CI and Codecov badge links (`GasChromatographyToolbox` org name).
- `ValveProgram` documentation now states the intended stepwise semantics (segment durations in `time_steps`) and references `valve_state` for evaluation.
- `module_temperature` now accepts `ModuleValve` using the same constant/program temperature handling used for `ModuleColumn`.
- Flow balance helpers now use `edge_restriction(...)` dispatch for `ModuleColumn`, `ModuleTM`, and `ModuleValve` instead of direct `.d` access.
- `common_timesteps`, `match_programs`, and `update_system` synchronize pressure and temperature only; `ModuleValve` `state` stays compact (`valve_state` at runtime).
- `flow_functions(sys, p2fun)` and `holdup_time_functions(sys, p2fun)` dispatch on `ModuleValve` using `valve_state` and instantaneous `d_open`/`d_closed` with GCSim `flow` / `holdup_time` (same pattern on both edges).
- Docstrings for `edge_restriction`, `flow_restrictions`, `flow_permeabilities`, `flow_functions`, and `holdup_time_functions`.
- `all_paths` requires `modules` (or `sys`) so valve edges can be filtered; `holdup_time_path` uses `all_paths(sys, num_paths)`.
- Removed valve program synchronization: `ValveProgram` / `PeriodicValveProgram` are no longer merged or resampled in `common_timesteps` / `update_system` (fixes large grids and slow `pressure_functions` on periodic DPM systems).

### Fixed
- `simulate_along_paths`: downstream column after a valve junction uses sliced `Parameters` from `apply_valve_junctions_at_vertex` (`par_in`), not the path template `par_sys`.
- `change_initial`: when peak-list annotations match valve slice rows (`"v…"`), downstream `τ₀` is taken from `par.sub` (honours `valve_initial_width=(:fixed, 0.0)`); otherwise keeps upstream `τR` (TM → column and column → column handoffs unchanged).
- `ValveProgram(time_steps, state_steps)` now validates matching vector lengths in the inner constructor (prevents inconsistent instances from the default typed constructor path).
- Implemented valve hydraulics in permeability/restriction evaluation using `σ(t) = valve_state(...)` with open/closed restrictions (`d_open`/`d_closed`), enabling `ModuleValve` edges in flow solves.
- `update_system` no longer mis-handles `ModuleValve` inside the column/TM temperature branch (constant-`T` valves unchanged; valves with `TemperatureProgram` `T` resample `T` only).
- `edge_restriction` for `ModuleValve`: correct `flow_restriction` arguments (`d_open` / `d_closed`), `module_.state` instead of `mod`, and scalar κ blend at `t`.

### Removed
- Unused duplicate GitHub Actions workflows under `data/.github/workflows/` (only `.github/workflows/` at the repo root is used).

## [0.2.7] - 2026-05-25

### Added
- Regression tests for finite `Program` pressure steps and zero injection times (`t₀`, `τ₀`) after `graph_to_parameters` on a series system with `NaN` junction pressures.

### Changed
- Updated `GasChromatographySimulator` compatibility to the `0.6` line in `Project.toml`.
- Updated ODE solver documentation in `Structures.jl` to match the solver set exposed by `GasChromatographySimulator` (`OwrenZen3`, `OwrenZen4`, `OwrenZen5`, `Tsit5`, `Vern9`, `BS5`, `DP5`).
- Aligned module option defaults with `GasChromatographySimulator` by keeping `OwrenZen5()` as the default algorithm in `ModuleTMOptions` and related `ModuleTM` constructors (previously `Vern9()` for thermal modulators).

### Fixed
- `graph_to_parameters` now builds inlet/outlet pressure step vectors by evaluating the resolved flow-balance pressure functions at each module's `time_steps`, instead of copying `NaN` placeholders from internal `PressurePoint`s. Restores compatibility with `GasChromatographySimulator` 0.6 `Program` validation while keeping `NaN` markers in system definitions for unknown pressures.
- `graph_to_parameters` passes finite default initial conditions (`t₀ = 0`, `τ₀ = 0`) to `load_solute_database` instead of `NaN`, matching `simulate_along_paths` and `GasChromatographySimulator` 0.6 `Substance` validation. Later graph segments still receive updated values via `change_initial` in `SolvingSystems.jl`.
