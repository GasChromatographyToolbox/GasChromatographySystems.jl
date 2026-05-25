# GasChromatographySystems.jl

[![DOI](https://zenodo.org/badge/423777500.svg)](https://zenodo.org/doi/10.5281/zenodo.10982222)
[![CI](https://github.com/GasChromatographyToolbox/GasChromatographySystems.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/GasChromatographyToolbox/GasChromatographySystems.jl/actions/workflows/ci.yml)
[![codecov.io](http://codecov.io/github/GasChromatographyToolbox/GasChromatographySystems.jl/coverage.svg?branch=main)](http://codecov.io/github/GasChromatographyToolbox/GasChromatographySystems.jl?branch=main)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://gaschromatographytoolbox.github.io/GasChromatographySystems.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://gaschromatographytoolbox.github.io/GasChromatographySystems.jl/dev)

A package for the simulation of complex gas chromatography (GC) systems with multiple columns, flow splitters, and thermal modulators. The package extends [GasChromatographySimulator.jl](https://github.com/GasChromatographyToolbox/GasChromatographySimulator.jl) to handle:

- **Multi-column systems** with series, parallel, and split configurations
- **Flow calculations** for complex capillary networks using graph theory
- **Thermal modulation** for comprehensive two-dimensional GC (GC×GC)
- **Pressure and temperature programs** synchronized across all system components
- **Deans switching** and multi-way splitting for advanced GC configurations

## Installation

To install the package type:

```julia
julia> ] add GasChromatographySystems
```

To use the package type:

```julia
julia> using GasChromatographySystems
```

## Documentation

Please read the [documentation page](https://gaschromatographytoolbox.github.io/GasChromatographySystems.jl/dev/) for more information.

## Notebooks

In the folder [notebooks](https://github.com/GasChromatographyToolbox/GasChromatographySystems.jl/tree/main/notebooks) several notebooks, using [Pluto.jl](https://github.com/fonsp/Pluto.jl), for the simulation of complex GC systems are available.

To use these notebooks [Julia, v1.6 or above,](https://julialang.org/downloads/#current_stable_release) must be installed and **Pluto** must be added:

```julia
julia> ]
(v1.7) pkg> add Pluto
```

To run Pluto, use the following commands:

```julia
julia> using Pluto
julia> Pluto.run()
```

Pluto will open your browser. In the field `Open from file` the URL of a notebook or the path to a locally downloaded notebook can be insert and the notebook will open and load the necessary packages.

### Overview of notebooks

- `FlowCalcPaper/4-Way-Splitter_Demo.jl` - Flow calculation of a 4-way splitter system demonstrating complex pressure balancing in multi-outlet configurations
- `FlowCalcPaper/DeansSwitch_Demo.jl` - Flow calculation of a Deans switch system showing how to model switching between different column configurations
- `GCxGC-TM_Paper/GCxGC_TM_Demo.jl` - Comprehensive two-dimensional GC simulation with thermal modulation, demonstrating GC×GC separations

## Key Features

### System Modeling
- **Graph-based representation** of GC systems using the Graphs.jl package
- **Automatic pressure calculation** at connection points between capillaries
- **Flow balance equations** solved using symbolic computation (Symbolics.jl)
- **Synchronized programs** for temperature and pressure across all modules

### Module Types
- **ModuleColumn**: Standard GC columns with temperature programs and gradients
- **ModuleTM**: Thermal modulators for GC×GC with periodic temperature modulation
- **PressurePoint**: System nodes with constant or programmed pressure values

### Advanced Configurations
- **Series systems**: Multiple columns in sequence
- **Split systems**: Flow splitting to multiple detectors
- **Deans switching**: Dynamic column switching during analysis
- **Multi-way splitters**: Complex flow distribution networks

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

## Citation

```
@software{leppert2024flowcalc,
  title = {Generalized flow calculation of the gas flow in a network of capillaries used in gas chromatography},
  author = {Leppert, Jan and Brehmer, Tillman and Boeker, Peter and Wüst, Matthias},
  volume = {47},
  number = {16},
  journal = {Journal of Separation Science},
  issn = {1615-9306, 1615-9314},
  year = {2024},
  pages = {2400419},
  doi = {10.1002/jssc.202400419},
  
}
```
