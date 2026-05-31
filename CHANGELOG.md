# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Valve module core structures in `Structures.jl`: `ValveProgram`, `ModuleValveOptions`, and `ModuleValve` with docstrings and convenience constructors.
- Periodic valve program builder `ValveProgram(mp, t_closed, t_end; inverted=false, t_start=0.0)` and `default_periodic_ValveProgram()`.
- Piecewise-constant valve state helper `valve_state(vp, t)` for open/closed switching without linear interpolation of boolean states.
- Test coverage for `ValveProgram`, `ModuleValveOptions`, and `ModuleValve` constructors/defaults, including periodic/inverted schedules.
- `index_modules_with_valve_program(sys)` for edges carrying a `ModuleValve` with `ValveProgram` state.
- Program synchronization for valves: `common_timesteps` merges `ValveProgram.time_steps`; `match_programs` resamples open/closed `state_steps` onto the common segment grid; `update_system` rebuilds `ModuleValve` modules (constant or programmed `T`, synchronized `state`).

### Changed
- CI: upgraded Codecov upload to `codecov/codecov-action@v5` with `files: lcov.info` and `CODECOV_TOKEN` (replaces deprecated v1 uploader and `CODECOV_SECRET`).
- README: fixed CI and Codecov badge links (`GasChromatographyToolbox` org name).
- `ValveProgram` documentation now states the intended stepwise semantics (segment durations in `time_steps`) and references `valve_state` for evaluation.
- `module_temperature` now accepts `ModuleValve` using the same constant/program temperature handling used for `ModuleColumn`.
- Flow balance helpers now use `edge_restriction(...)` dispatch for `ModuleColumn`, `ModuleTM`, and `ModuleValve` instead of direct `.d` access.
- `common_timesteps`, `match_programs`, and `update_system` docstrings now describe valve programs alongside pressure and temperature programs.
- `flow_functions(sys, p2fun)` and `holdup_time_functions(sys, p2fun)` dispatch on `ModuleValve` using `valve_state` and instantaneous `d_open`/`d_closed` with GCSim `flow` / `holdup_time` (same pattern on both edges).
- Docstrings for `edge_restriction`, `flow_restrictions`, `flow_permeabilities`, `flow_functions`, and `holdup_time_functions`.

### Fixed
- `ValveProgram(time_steps, state_steps)` now validates matching vector lengths in the inner constructor (prevents inconsistent instances from the default typed constructor path).
- Implemented valve hydraulics in permeability/restriction evaluation using `σ(t) = valve_state(...)` with open/closed restrictions (`d_open`/`d_closed`), enabling `ModuleValve` edges in flow solves.
- `update_system` no longer mis-handles `ModuleValve` inside the column/TM temperature branch (constant-`T` valves now receive a synchronized `ValveProgram`; valves with `TemperatureProgram` `T` update both programs).
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
