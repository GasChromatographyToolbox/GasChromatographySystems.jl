# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.7] - 2026-05-25

### Changed
- Updated `GasChromatographySimulator` compatibility to the `0.6` line in `Project.toml`.
- Updated ODE solver documentation in `Structures.jl` to match the solver set exposed by `GasChromatographySimulator` (`OwrenZen3`, `OwrenZen4`, `OwrenZen5`, `Tsit5`, `Vern9`, `BS5`, `DP5`).
- Aligned module option defaults with `GasChromatographySimulator` by keeping `OwrenZen5()` as the default algorithm in `ModuleTMOptions` and related `ModuleTM` constructors (previously `Vern9()` for thermal modulators).

### Fixed
- `graph_to_parameters` now builds inlet/outlet pressure step vectors by evaluating the resolved flow-balance pressure functions at each module's `time_steps`, instead of copying `NaN` placeholders from internal `PressurePoint`s. Restores compatibility with `GasChromatographySimulator` 0.6 `Program` validation while keeping `NaN` markers in system definitions for unknown pressures.
- `graph_to_parameters` passes finite default initial conditions (`t₀ = 0`, `τ₀ = 0`) to `load_solute_database` instead of `NaN`, matching `simulate_along_paths` and `GasChromatographySimulator` 0.6 `Substance` validation. Later graph segments still receive updated values via `change_initial` in `SolvingSystems.jl`.
