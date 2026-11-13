# GasChromatographySystems.jl

Documentation for [GasChromatographySystems.jl](https://github.com/GasChromatographyToolbox/GasChromatographySystems)

## Introduction

The package GasChromatographySystems.jl simulates the separation of substances in a [gas chromatographic (GC) system](https://en.wikipedia.org/wiki/Gas_chromatography). The core of the simulation is realized with the package [GasChromatographySimulator.jl](https://github.com/GasChromatographyToolbox/GasChromatographySimulator.jl), which uses ordinary differential equations (ODE) to model solute migration and peak broadening — by default as ``t(z)`` and ``\tau^2(z)``, or optionally as ``x(t)`` and ``\sigma^2(t)`` on selected column edges (`ModuleColumnOptions(domain=:temporal)`). GasChromatographySystems.jl extends this to complex GC systems of multiple capillaries. Each segment is solved with GasChromatographySimulator and the result feeds the next segment (typically ``t_0 = t_R`` and ``\tau_0 = \tau_R`` at a junction — see [Segment handoff](@ref)).

To model the separation in the complex GC system, it is necessary to know the inlet and outlet pressures of every capillary of this system, which define the flows ``F`` and hold-up times ``t_M`` along these capillaries. Therefore, the GasChromatographySystems.jl package also includes a flow calculator for complex GC systems. This flow calculator uses a  [Graph](https://en.wikipedia.org/wiki/Graph_(discrete_mathematics)) representation of the system of capillaries to automatically calculate the pressures at the connection points between the capillaries.

## Content

```@contents
Pages = [
    "installation.md",
    "definition_sys.md",
    "flowcalc.md",
    "simulation.md",
    "segment_handoff.md",
    "examples.md",
    "docstrings.md"
]
Depth = 2
```

## Contribution

Please open an issue if you:
- want to report a bug
- have problems using the package (please first look at the documentation)
- have ideas for new features or ways to improve the usage of this package

You can contribute (e.g. fix bugs, add new features, add to the documentation) to this package by Pull Request:
- first discuss your contributions in a new issue
- ensure that all tests pass locally before starting the pull request
- new features should be included in `runtests.jl`
- add description to the pull request, link to corresponding issues by `#` and issue number
- the pull request will be reviewed
