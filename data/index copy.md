# GasChromatographySystems.jl

Documentation for GasChromatographySystems.jl



## Introduction

The package GasChromatographySystems.jl simulates the separation of substances in a [gas chromatographic (GC) system](https://en.wikipedia.org/wiki/Gas_chromatography). The core of the simulation is realized with the package [GasChromatographySimulator.jl](https://github.com/JanLeppert/GasChromatographySimulator.jl), which uses ordinary differential equations (ODE) to model the migration ``t(x)`` of a substance through the GC system and the development of the peak variance during this migration ``τ²(x)``. In difference to this package, GasChromatographySystems.jl simulates this separation not only for one GC column/capillary, but for a complex GC system consisting of multiple capillaries. Each separation on the single capillaries is modeled using the GasChromatographySimulator.jl package and uses the results of one capillary as initial values for the next capillary.

To model the separation in the complex GC system, it is necessary to know the inlet and outlet pressures of every capillary of this system, which define the flows ``F`` and hold-up times ``t_M`` along these capillaries. Therefore, the GasChromatographySystems.jl package also includes a flow calculator for complex GC systems. This flow calculator uses a [Graph](https://en.wikipedia.org/wiki/Graph_(discrete_mathematics)) representation of the system of capillaries to automatically calculate the pressures at the connection points between the capillaries.  

## Content

The manual is structured as followed:

```@contents
Pages = [
    "installation.md",
    "definition_sys.md",
    "flowcalc.md",
    "simulation.md,
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




