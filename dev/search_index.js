var documenterSearchIndex = {"docs":
[{"location":"simulation/#Simulation","page":"Simulation","title":"Simulation","text":"","category":"section"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"flowcalc/#Flow-calculation","page":"Flow calculation","title":"Flow calculation","text":"","category":"section"},{"location":"flowcalc/#Paths","page":"Flow calculation","title":"Paths","text":"","category":"section"},{"location":"flowcalc/#Flow-balance","page":"Flow calculation","title":"Flow balance","text":"","category":"section"},{"location":"flowcalc/#Hold-up-times","page":"Flow calculation","title":"Hold-up times","text":"","category":"section"},{"location":"installation/#Installation","page":"Instalation","title":"Installation","text":"","category":"section"},{"location":"installation/","page":"Instalation","title":"Instalation","text":"To use GasChromatographicystems.jl, you need to install Julia 1.9 or greater first (official Julia website) and than add the package:","category":"page"},{"location":"installation/","page":"Instalation","title":"Instalation","text":"julia> using Pkg; Pkg.add(GasChromatographySystems)","category":"page"},{"location":"installation/","page":"Instalation","title":"Instalation","text":"To use the package type:","category":"page"},{"location":"installation/","page":"Instalation","title":"Instalation","text":"using GasChromatographySystems","category":"page"},{"location":"definition_sys/#Definition-of-a-GC-system","page":"Definition of a GC system","title":"Definition of a GC system","text":"","category":"section"},{"location":"definition_sys/","page":"Definition of a GC system","title":"Definition of a GC system","text":"In this segment the structures and function, which are used to define a complex GC system, will be presented.","category":"page"},{"location":"definition_sys/","page":"Definition of a GC system","title":"Definition of a GC system","text":"The structure of a complex GC system consists of four parts:","category":"page"},{"location":"definition_sys/","page":"Definition of a GC system","title":"Definition of a GC system","text":"Representation of the GC system as a simple directional graph.\nPressure points defining the pressures at the vertices of the graph (connections of the capillaries).\nModules defining the features of the capillaries, e.g. dimensions and temperature.\nAdditional options, e.g. typ of the mobile phase gas.","category":"page"},{"location":"definition_sys/","page":"Definition of a GC system","title":"Definition of a GC system","text":"(ref to GasChromatographySystems.System)","category":"page"},{"location":"definition_sys/#GC-system-as-a-graph","page":"Definition of a GC system","title":"GC system as a graph","text":"","category":"section"},{"location":"definition_sys/","page":"Definition of a GC system","title":"Definition of a GC system","text":"A complex GC system consisting of multiple connected capillaries can be abstracted and represented by a Graph. A Graph G(VE) is mathematical construct consisting of a set of vertices V and a set of edges E, used to represent connections between multiple objects. In the case of the complex GC system the vertices, also called nodes, represent the connections of the capillaries and are related to the inlet and outlet pressures of the capillaries. The edges are the capillaries. A simple directed graph is used. In this context, \"simple\" means, that each possible connection between vertices exists only once and \"directed\" means, that edges have direction going from one vertex i to vertex j is different from an edge connecting vertices j and i. This is used to represent the preferred flow direction through the capillaries.","category":"page"},{"location":"definition_sys/","page":"Definition of a GC system","title":"Definition of a GC system","text":"The Julia package Graphs.jl is used for the implementation of the simple directed graphs.","category":"page"},{"location":"definition_sys/","page":"Definition of a GC system","title":"Definition of a GC system","text":"In the structure GasChromatographySystems.System(g, pressurepoints, modules, options) the graph is defined in g. For example the graph representing 3 capillaries in series is constructed by","category":"page"},{"location":"definition_sys/","page":"Definition of a GC system","title":"Definition of a GC system","text":"First, a simple directed graph with 4 vertices is defined:","category":"page"},{"location":"definition_sys/","page":"Definition of a GC system","title":"Definition of a GC system","text":"    using GasChromatographySystems\n    g = SimpleDiGraph(4)\n    nothing # hide","category":"page"},{"location":"definition_sys/","page":"Definition of a GC system","title":"Definition of a GC system","text":"In following steps the vertices are connected by the edges:","category":"page"},{"location":"definition_sys/","page":"Definition of a GC system","title":"Definition of a GC system","text":"    add_edge!(g, 1, 2)\n    add_edge!(g, 2, 3)\n    add_edge!(g, 3, 4)\n    nothing # hide","category":"page"},{"location":"definition_sys/","page":"Definition of a GC system","title":"Definition of a GC system","text":"This graph can be visualized:","category":"page"},{"location":"definition_sys/","page":"Definition of a GC system","title":"Definition of a GC system","text":"    using CairoMakie # hide\n    f = GasChromatographySystems.plot_graph(g, [\"$(i)=>$(i+1)\" for i=1:3], [\"p$(i)\" for i=1:4])\n    CairoMakie.save(\"plot_graph_3series.svg\", f) #hide\n    nothing # hide","category":"page"},{"location":"definition_sys/","page":"Definition of a GC system","title":"Definition of a GC system","text":"(Image: )","category":"page"},{"location":"definition_sys/#Pressure-points","page":"Definition of a GC system","title":"Pressure points","text":"","category":"section"},{"location":"definition_sys/","page":"Definition of a GC system","title":"Definition of a GC system","text":"The vertices V represent the connections of the capillaries of the GC system. To fully define the flow of the mobile phase through the GC system the pressure at these connections must be known or calculated. To store this information about the pressures, two structures are defined. One to define the vertex as pressure point and one to store information about a pressure program.","category":"page"},{"location":"definition_sys/","page":"Definition of a GC system","title":"Definition of a GC system","text":"The first structure is GasChromatographySystems.PressurePoint(name, PP). It is composed of the name of the pressure point resp. vertex and the value of the pressure, either as a number or as a pressure program. The pressure program is defined by the second structure of GasChromatographySystems.PressureProgram(time_steps, pressure_steps), to store the pressure_steps at the times of the corresponding time_steps. In between the time_steps the pressure values are linearly interpolated.","category":"page"},{"location":"definition_sys/","page":"Definition of a GC system","title":"Definition of a GC system","text":"For our example of three capillaries in series we define a pressure program for the inlet pressure, with a pressure of 150000 Pa held for the first two minutes and than increasing the pressure over a time span of ten minutes to 230000 Pa and hold it there for a further minute:","category":"page"},{"location":"definition_sys/","page":"Definition of a GC system","title":"Definition of a GC system","text":"    PPin = GasChromatographySystems.PressureProgram([0.0, 120.0, 600.0, 60.0], [150000.0, 150000.0, 230000.0, 230000.0])\nnothing # hide","category":"page"},{"location":"definition_sys/","page":"Definition of a GC system","title":"Definition of a GC system","text":"An array of type GasChromatographySystems.PressurePoint with the length of number of vertices of the graph for our example is created with the defined pressure program for the inlet, a constant pressure of 101300 Pa for the outlet. The pressures in between are unknown and NaN (not a number) will be assigned to these pressure points.","category":"page"},{"location":"definition_sys/","page":"Definition of a GC system","title":"Definition of a GC system","text":"    pp = Array{GasChromatographySystems.PressurePoint}(undef, nv(g))\n\tpp[1] = GasChromatographySystems.PressurePoint(\"p₁\", PPin) # inlet \n\tpp[2] = GasChromatographySystems.PressurePoint(\"p₂\", NaN)\n\tpp[3] = GasChromatographySystems.PressurePoint(\"p₃\", NaN) \n\tpp[4] = GasChromatographySystems.PressurePoint(\"p₄\", 101300.0) #outlet \n    nothing # hide","category":"page"},{"location":"definition_sys/#Modules","page":"Definition of a GC system","title":"Modules","text":"","category":"section"},{"location":"definition_sys/#Additional-options","page":"Definition of a GC system","title":"Additional options","text":"","category":"section"},{"location":"#GasChromatographySystems.jl","page":"Home","title":"GasChromatographySystems.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for GasChromatographySystems.jl","category":"page"},{"location":"#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package GasChromatographySystems.jl simulates the separation of substances in a gas chromatographic (GC) system. The core of the simulation is realized with the package GasChromatographySimulator.jl, which uses ordinary differential equations (ODE) to model the migration t(x) of a substance through the GC system and the development of the peak variance during this migration τ²(x). In difference to this package, GasChromatographySystems.jl simulates this separation not only for one GC column/capillary, but for a complex GC system consisting of multiple capillaries. Each separation on the single capillaries is modeled using the GasChromatographySimulator.jl package and uses the results of one capillary as initial values for the next capillary.","category":"page"},{"location":"","page":"Home","title":"Home","text":"To model the separation in the complex GC system, it is necessary to know the inlet and outlet pressures of every capillary of this system, which define the flows F and hold-up times t_M along these capillaries. Therefore, the GasChromatographySystems.jl package also includes a flow calculator for complex GC systems. This flow calculator uses a  Graph representation of the system of capillaries to automatically calculate the pressures at the connection points between the capillaries.  ","category":"page"},{"location":"#Content","page":"Home","title":"Content","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The manual is structured as followed:","category":"page"},{"location":"#Contribution","page":"Home","title":"Contribution","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Please open an issue if you:","category":"page"},{"location":"","page":"Home","title":"Home","text":"want to report a bug \nhave problems using the package (please first look at the documentation)\nhave ideas for new features or ways to improve the usage of this package ","category":"page"},{"location":"","page":"Home","title":"Home","text":"You can contribute (e.g. fix bugs, add new features, add to the documentation) to this package by Pull Request: ","category":"page"},{"location":"","page":"Home","title":"Home","text":"first discuss your contributions in a new issue\nensure that all tests pass locally before starting the pull request\nnew features should be included in runtests.jl\nadd description to the pull request, link to corresponding issues by # and issue number\nthe pull request will be reviewed","category":"page"},{"location":"docstrings/#Module-Index","page":"Docstrings","title":"Module Index","text":"","category":"section"},{"location":"docstrings/#Detailed-API","page":"Docstrings","title":"Detailed API","text":"","category":"section"},{"location":"docstrings/","page":"Docstrings","title":"Docstrings","text":"Modules = [GasChromatographySystems]\nOrder   = [:constant, :type, :function, :macro]","category":"page"},{"location":"docstrings/#GasChromatographySystems.ModuleColumn","page":"Docstrings","title":"GasChromatographySystems.ModuleColumn","text":"ModuleColumn(; Tcold_abs=true, sflank=40, tflank=20, alg=Vern9(), abstol=1e-10, reltol=1e-8, dtinit=1e-6, ng=true, Tcontrol=\"inlet\")\n\nStructure describing the module of a gas chromatographic column, the basic buliding block of a gas chromatographic system. These are edges of the graph representation of a GC system. \n\nArguments\n\nname: Name of the module.\nL: Length of the column in m.\nd: Internal diameter of the column in m. This can be a function of column position, if a non-uniform diameter is used (not tested yet).\na_d: Parameters of the diameter function d. Empty, if d is not a function. \ndf: Film thickness of the column in m. This can be a function of column position, if a non-uniform film thickness is used (not tested yet).\na_df: Parameters of the film thickness function df. Empty, if df is not a function. \nsp: The name of the stationary phase.\nT: Temperature of the module. This can be a value (in °C) or a temperature program as a TemperatureProgram structure. \nF: Flow of the mobile phase through the column in mL/min. This can be an number, a fuction or a NaN (if it is not defined yet an should be calculated based on values from other modules or pressure points in the system).\nopt: Options for this module. This is a structure of the form ModuleColumnOptions.\n\nFour methods to construct this structure exist:\n\nModuleColumn(name, L, d, df, sp, T, opt::ModuleColumnOptions): The flow is not defined, uniform d and df.\nModuleColumn(name, L, d, df, sp, T, F, opt::ModuleColumnOptions): The flow is defined by F and uniform d and df.\nModuleColumn(name, L, d, df, sp, T; alg=OwrenZen5(), abstol=1e-8, reltol=1e-6, ng=false, Tcontrol=\"inlet\"): The flow is not defined, uniform d and df. Options according to ModuleColumnOptions.\nModuleColumn(name, L, d, df, sp, tp, flow; alg=OwrenZen5(), abstol=1e-8, reltol=1e-6, ng=false, Tcontrol=\"inlet\"): The flow is defined by F and uniform d and df. Options according to ModuleColumnOptions.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#GasChromatographySystems.ModuleColumnOptions","page":"Docstrings","title":"GasChromatographySystems.ModuleColumnOptions","text":"ModuleColumnOptions(; alg=OwrenZen5(), abstol=1e-8, reltol=1e-5, ng=false, Tcontrol=\"inlet\")\n\nStructure describing options for column modules. \n\nArguments\n\nalg: The algorithm used for the ODE solver. The algorithms OwrenZen3(), OwrenZen4() and OwrenZen5() are recommended.\nabstol: The absolute tolerance for the ODE solver. Recommended value 1e-6 to 1e-8.\nreltol: The relative tolerance for the ODE solver. Recommended value 1e-3 to 1e-5. \nng: Option to calculate the simulation without a gradient (ng = true, default) or with a gradient (ng = false).\nTcontrol: Option defining at which point of the column the temperature program is calculated. The options are inlet (x=0) and outlet (x=L).\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#GasChromatographySystems.ModuleTM","page":"Docstrings","title":"GasChromatographySystems.ModuleTM","text":"ModuleTM(; Tcold_abs=true, sflank=40, tflank=20, alg=Vern9(), abstol=1e-10, reltol=1e-8, dtinit=1e-6, ng=true, Tcontrol=\"inlet\")\n\nStructure describing the module of a thermal modulator, a short section of a column with a periodic repeating temperature modulation between a cold and hot jet. These are edges of the graph representation of a GC system.\n\nArguments\n\nname: Name of the module.\nL: Length of the column in m.\nd: Internal diameter of the column in m. This can be a function of column position, if a non-uniform diameter is used (not tested yet).\na_d: Parameters of the diameter function d. Empty, if d is not a function. \ndf: Film thickness of the column in m. This can be a function of column position, if a non-uniform film thickness is used (not tested yet).\na_df: Parameters of the film thickness function df. Empty, if df is not a function. \nsp: The name of the stationary phase.\nT: Basic temperature of the module. On top of this temperature the modulation occurs. This can be a value (in °C) or a temperature program as a TemperatureProgram structure. \nshift: Time shift in s of the periodic modulation relative to the start of the chromatogram.\nPM: Modulation period in s.\nratio: The ratio of the duration of cold and hot jet.\nThot: Temperature in °C, by which T is increased while the hot jet is active. \nTcold: Temperature in °C, by which T is decrease while the cold jet is active (if Tcold_abs = false in ModuleColumnOptions) or to which the column is cooled down (if Tcold_abs = true in ModuleTMOptions)\nF: Flow of the mobile phase through the column in mL/min. This can be an number, a fuction or a NaN (if it is not defined yet an should be calculated based on values from other modules or pressure points in the system).\nopt: Options for this module. This is a structure of the form ModuleTMOptions.\n\nFour methods to construct this structure exist:\n\nModuleTM(name, L, d, df, sp, T, shift, PM, ratio, Thot, Tcold, opt::ModuleTMOptions): The flow is not defined, uniform d and df.\nModuleTM(name, L, d, df, sp, T, shift, PM, ratio, Thot, Tcold, F, opt::ModuleTMOptions): The flow is defined by F and uniform d and df.\nModuleTM(name, L, d, df, sp, T, shift, PM, ratio, Thot, Tcold; Tcold_abs=true, sflank=40, tflank=20, alg=Vern9(), abstol=1e-10, reltol=1e-8, dtinit=1e-6, ng=false, Tcontrol=\"inlet\"): The flow is not defined, uniform d and df. Options according to ModuleTMOptions.\nModuleTM(name, L, d, df, sp, tp, shift, pm, ratio, Thot, Tcold, F; Tcold_abs=true, sflank=40, tflank=20, alg=Vern9(), abstol=1e-10, reltol=1e-8, dtinit=1e-6, ng=false, Tcontrol=\"inlet\"): The flow is defined by F and uniform d and df. Options according to ModuleTMOptions.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#GasChromatographySystems.ModuleTMOptions","page":"Docstrings","title":"GasChromatographySystems.ModuleTMOptions","text":"ModuleTMOptions(; Tcold_abs=true, sflank=40, tflank=20, alg=Vern9(), abstol=1e-10, reltol=1e-8, dtinit=1e-6, ng=true, Tcontrol=\"inlet\")\n\nStructure describing options for thermal modulator modules. \n\nArguments\n\nTcold_abs: Calculate the low temperature at the modulator as the absolute value of the defined parameter Tcold (Tcold_abs = true) or as relative temperature difference from the defined oven temperature (Tcold_abs = false).\nsflank: Flank factor of the smoothed rectangle temperature function in space. A higher factor results in a steeper slope at the edges of the modulator point. Values between 12 and 100, or Inf if ng = true. Recommend sflank = 40. Only relevante, if option ng = false. \ntflank: Flank factor of the smoothed rectangle temperature function in time. A higher factor results in a steeper slope at the begining and end of the hot jet. Values between 12 and 100. Recommend tflank = 20.\nalg: The algorithm used for the ODE solver. For the thermal modulator module the algorithm Vern9() is recommend.\nabstol: The absolute tolerance for the ODE solver. Recommended value 1e-6 to 1e-12.\nreltol: The relative tolerance for the ODE solver. Recommended value 1e-3 to 1e-10.\ndtinit: The initial step width for the ODE solver. A value of L*1e-6, with L the length of the modulator module, is recommend. \nng: Option to calculate the simulation without a gradient (ng = true, default) or with a gradient (ng = false).\nTcontrol: Option defining at which point of the column the temperature program is calculated. The options are inlet (x=0) and outlet (x=L).\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#GasChromatographySystems.Options","page":"Docstrings","title":"GasChromatographySystems.Options","text":"Options(; gas=\"He\", odesys=true, vis=\"Blumberg\", control=\"Pressure\", k_th=1e12)\n\nStructure describing general options for the simulation of the system. \n\nArguments\n\ngas: The type of gas used as mobile phase in the gas chromatographic system. Allowed values: He, H2 or N2.\nodesys: Combine the ODEs for migration and peak-width into a system of ODEs (odesys = true) or solve the two ODEs separately (odesys = false).\nvis: Used model of viscosity. HP is a second-order polynomial taken from the HP flow calculator. Blumberg is an emperical formula according to the book   Temperature-programmed Gas Chromatography by Leonid M. Blumberg (2010, Wiley-VCH).\ncontrol: Control of the \"Flow\" or of the \"Pressure\" (at column inlet) during the program.\nk_th: Threshold for the maximum of the retention factor. If the calculated retention factor is bigger than k_th than the retention factor is set to the value k_th.   This is done to avoid to small step widths in the solver for highly retained soultes at the beginning of a GC program. \n\n\n\n\n\n","category":"type"},{"location":"docstrings/#GasChromatographySystems.PressurePoint","page":"Docstrings","title":"GasChromatographySystems.PressurePoint","text":"PressurePoint(name, PP)\n\nStructure describing the pressure program at the connection points of the modules. These are the vertices of the graph representation of a GC system.\n\nArguments\n\nname: Name of the PressurePoint.\nP: Pressure program as structure GasChromatographySystems.PressureProgram or as a number for constant pressure.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#GasChromatographySystems.PressureProgram","page":"Docstrings","title":"GasChromatographySystems.PressureProgram","text":"PressureProgram(time_steps, pressure_steps)\n\nStructure describing the pressure program. \n\nArguments\n\ntime_steps: Time steps in s, after which the corresponding pressure in pressure_steps is reached. \npressure_steps: Pressure steps in Pa.\n\nA default temperature program is avaliable:\n\ndefault_PP(): Pressure increase from 100.000 Pa to 200.000 in 1800 s (30 min). \n\n\n\n\n\n","category":"type"},{"location":"docstrings/#GasChromatographySystems.System","page":"Docstrings","title":"GasChromatographySystems.System","text":"System(g, pressurepoints, modules, options)\n\nStructure describing a GC system. \n\nArguments\n\nname: Name of the GC system.\ng: The graph representation of the GC system, using SimpleDiGraph from Graphs.jl. \npressurepoints: The vertices of graph g and their pressure programs. These are structures of type PressurePoint. \nmodules: The edges of graph g. These are column segments of type ModuleColumn or thermal modulator points of type ModuleTM.\n\nExamples\n\nSimple GC:\n\nDefinition of the graph:\n\n    g = SimpleDiGraph(1)\n    add_edge!(g, 1, 2)\n\nDefinition of the two pressure points, here column inlet has the default pressure program and the outlet pressure is constant:\n\n    pp = Array{GasChromatographySystems.PressurePoint}(undef, nv(g))\n    pp[1] = GasChromatographySystems.PressurePoint(\"p_in\", GasChromatographySystems.default_PP())\n    pp[2] = GasChromatographySystems.PressurePoint(\"p_out\", pout)\n\nDefinition of the column module with default temperature program default_TP(), flow is unknown and is calculated from the pressures, default module options:\n\n    modules = Array{GasChromatographySystems.AbstractModule}(undef, ne(g))  \n    modules[1] = GasChromatographySystems.ModuleColumn(\"column\", 30.0, 0.25*1e-3, 0.25*1e-6, \"Rxi5ms\", GasChromatographySystems.default_TP())\n\nCombination to construct the System using the function update_system() and default system options:\n\n    sys = GasChromatographySystems.update_system(GasChromatographySystems.System(g, pp, modules, GasChromatographySystems.Options()))\n\nSeries of n columns:\n\nDefinition of the graph:\n\n    g = SimpleDiGraph(n+1)\n\tfor i=1:n\n\t\tadd_edge!(g, i, i+1) \n\tend\n\nSplit of Columns:\n\nGCxGC with thermal modulation and loop:\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#GasChromatographySystems.TemperatureProgram","page":"Docstrings","title":"GasChromatographySystems.TemperatureProgram","text":"TemperatureProgram(time_steps, temp_steps, gf, a_gf)\n\nStructure describing the temperature program. \n\nArguments\n\ntime_steps: Time steps in s, after which the corresponding temperature in temp_steps is reached. \ntemp_steps: Temperature steps in °C.\ngf: Gradient function gf(x, a_gf), describes the thermal gradient.\na_gf:Parameters of the gradient function.\n\nTwo method to construct this structure exist for a uniform temperatur program (same for all column positions x):\n\nTemperatureProgram(time_steps, temp_steps)\nTemperatureProgram(CP), with CP beeing a conventional notation of a temperature program in the form of an array with the pattern [T1, t1, r1, T2, t2, r2, ..., Tn, tn], with Ti temperature niveaus, ti holding times for the corresponding tempertures, ri the heating ramp between temperatures Ti and Ti+1.\n\nA default temperature program is avaliable:\n\ndefault_TP(): Heating from 40°C to 340°C in 1800 s (30 min). \n\n\n\n\n\n","category":"type"},{"location":"docstrings/#GasChromatographySystems.build_pressure_squared_functions-Tuple{Any}","page":"Docstrings","title":"GasChromatographySystems.build_pressure_squared_functions","text":"build_pressure_squared_functions(sys; mode=\"λ\")\n\nConstruct array of functions of the solutions for the unkown squared pressures of the flow balance equations of the system of capillaries sys.\n\nThe arguments for the build functions are arrays of the ordered known squared pressures p^2, the ordered known flow permabilities λ resp. flow restrictions κ, the ordered known flows F and constant A = π256 p_nT_n.  \n\nArguments\n\nsys: System structure of the capillary system for which the flow balance is set up.\nmode: Mode for flow equations to use flow permeabilities λ (mode = λ; default) or flow restrictions κ (mode = κ)\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#GasChromatographySystems.check_expressions_λ_κ-Tuple{Any}","page":"Docstrings","title":"GasChromatographySystems.check_expressions_λ_κ","text":"check_expressions_λ_κ(sol; mode=\"λ\", n=100)\n\nChecks the expressions of the array sol (solutions to the flow balance equations) if they use the flow permeabilities λ or the flow restrictions κ and substitutes them if needed. \n\nArguments\n\nsol: Symbolic expressions (of the solutions for the flow balance equations) \nmode: Mode for flow equations to use flow permeabilities λ (mode = λ; default) or flow restrictions κ (mode = κ)\nn: Maximum of the number of expected symbols λ resp. κ. Could be replaced by the number of edges of the used system n = ne(sys.g). \n\n\n\n\n\n","category":"method"},{"location":"docstrings/#GasChromatographySystems.flow_balance-Tuple{Any}","page":"Docstrings","title":"GasChromatographySystems.flow_balance","text":"flow_balance(sys)\n\nConstructing the flow balance equations of the capillary system sys in the form of an array of symbolic equations.\n\nFor every inner vertice, the sum of ingoing flows (positive) and of outgoing flows (negative) are equated to zero.\n\nSum  F_in + Sum F_out = 0\n\nArguments\n\nsys: System structure of the capillary system for which the flow balance is set up.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#GasChromatographySystems.flow_functions-Tuple{Any, Any}","page":"Docstrings","title":"GasChromatographySystems.flow_functions","text":"flow_functions(sys, p2fun; mode=\"λ\")\n\nCollects the flow functions as functions of time t for all edges of the system of capillaries sys.\n\nThe flow over edge i => j is calculated as\n\nF_ij = fracAκ_ij left(p_i^2-p_j^2right)\n\nwith flow restriction κ_ij = int_0^L_ij η(T_ij)T_ijd_ij^4 dx, pressures p_i resp. p_j at the vertices i resp. j, temperature T_ij, capillary length L_ij and diameter d_ij of the edge i => j.\n\nArguments\n\nsys: System structure of the capillary system for which the flow balance is set up.\np2fun: Julia function of the solutions of the flow balance equations from build_pressure_squared_functions(sys; mode=\"λ\")\nmode: Mode for flow equations to use flow permeabilities λ (mode = λ; default) or flow restrictions κ (mode = κ)\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#GasChromatographySystems.flow_permeabilities-Tuple{Any}","page":"Docstrings","title":"GasChromatographySystems.flow_permeabilities","text":"flow_permeabilities(sys)\n\nCalculates the flow permeabilities λ of all edges (capliaries) of a system of capillaries.\n\nArguments\n\nsys: System structure of the capillary system for which the flow balance is set up.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#GasChromatographySystems.flow_restrictions-Tuple{Any}","page":"Docstrings","title":"GasChromatographySystems.flow_restrictions","text":"flow_restrictions(sys)\n\nCalculates the flow restrictions κ of all edges (capliaries) of a system of capillaries.\n\nArguments\n\nsys: System structure of the capillary system for which the flow balance is set up.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#GasChromatographySystems.holdup_time_functions-Tuple{Any, Any}","page":"Docstrings","title":"GasChromatographySystems.holdup_time_functions","text":"holdup_time_functions(sys, p2fun; mode=\"λ\")\n\nCollects the hold-up time functions as functions of time t for all edges of the system of capillaries sys.\n\nThe hold-up time over edge i => j is calculated as\n\nt_M_ij = frac1283 η(T_ij) fracL_ij^2d_ij^2 fracp_i^3-p_j^3left(p_i^2-p_j^2right)^2\n\nwith flow restriction, pressures p_i resp. p_j at the vertices i resp. j, temperature T_ij, capillary length L_ij and diameter d_ij of the edge i => j.\n\nArguments\n\nsys: System structure of the capillary system for which the flow balance is set up.\np2fun: Julia function of the solutions of the flow balance equations from build_pressure_squared_functions(sys; mode=\"λ\")\nmode: Mode for flow equations to use flow permeabilities λ (mode = λ; default) or flow restrictions κ (mode = κ)\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#GasChromatographySystems.holdup_time_path-Tuple{Any, Any, Any}","page":"Docstrings","title":"GasChromatographySystems.holdup_time_path","text":"holdup_time_path(sys, p2fun, numpaths; mode=\"λ\")\n\nCalculates the hold-up times of the numpaths paths as functions of time t of the system of capillaries sys.\n\nThe hold-up time over edge i => j is calculated as\n\nArguments\n\nsys: System structure of the capillary system for which the flow balance is set up.\np2fun: Julia function of the solutions of the flow balance equations from build_pressure_squared_functions(sys; mode=\"λ\")\nnumpaths: Number of the different paths between inlet and outlets of system sys.\nmode: Mode for flow equations to use flow permeabilities λ (mode = λ; default) or flow restrictions κ (mode = κ)\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#GasChromatographySystems.interpolate_pressure_functions-Tuple{Any, Any}","page":"Docstrings","title":"GasChromatographySystems.interpolate_pressure_functions","text":"interpolate_pressure_functions(sys, p2fun; dt=1, mode=\"λ\")\n\nInterpolates (linearly) all pressure funtions at the vertices of the system of capillaries sys between the time steps dt. For the speed of the simulation these interpolated functions are faster than the pure solution functions of the flow balance equations.\n\nArguments\n\nsys: System structure of the capillary system for which the flow balance is set up.\np2fun: Julia function of the solutions of the flow balance equations from build_pressure_squared_functions(sys; mode=\"λ\")\ndt: time steps, where the original pressure function is evaluated. Inbetween these time steps the pressure function is linearly interpolated. \nmode: Mode for flow equations to use flow permeabilities λ (mode = λ; default) or flow restrictions κ (mode = κ)\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#GasChromatographySystems.pressure_functions-Tuple{Any, Any}","page":"Docstrings","title":"GasChromatographySystems.pressure_functions","text":"pressure_functions(sys, p2fun; mode=\"λ\")\n\nCollect all pressure functions as functions of time t at the vertices of the capillary system sys, either from defined input values or from the solutions of the flow balance equations. \n\nArguments\n\nsys: System structure of the capillary system for which the flow balance is set up.\np2fun: Julia function of the solutions of the flow balance equations from build_pressure_squared_functions(sys; mode=\"λ\")\nmode: Mode for flow equations to use flow permeabilities λ (mode = λ; default) or flow restrictions κ (mode = κ)\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#GasChromatographySystems.save_build_pressure_squared_functions-Tuple{Any, Any}","page":"Docstrings","title":"GasChromatographySystems.save_build_pressure_squared_functions","text":"save_build_pressure_squared_functions(sys, solution; filename=pwd()*\"/p2fun_\"*sys.name, mode=\"λ\")\n\nConstructs and saves the array of functions of the solutions solution for the unkown squared pressures of the flow balance equations of the system of capillaries sys.\n\nArguments\n\nsys: System structure of the capillary system for which the flow balance is set up.\nfilename: Filename, where the solution functions are saved. Default pwd()*\"/p2fun_\"*sys.name attached with \"_λ.jl\" or \"_κ.jl\", depending on mode.\nmode: Mode for flow equations to use flow permeabilities λ (mode = λ; default) or flow restrictions κ (mode = κ)\n\nOutput\n\nA dictionary with the folowing keys is saved in a .jl file:\n\n\"p2fun\": the build function of the solutions\n\"i_known_p\": the indices of the known pressures.\n\"i_known_λ\": the indices of the known flow permeabilities.\n\"i_known_F\": the indices of the known flows.\n\"mode\": mode of the functions using flow permeabilities λ or flow restrictions κ.\n\nLoading\n\nThe saved dictionary can easily be loaded into Julia by\n\np2fun_load = include(\"p2fun_saved.jl\")\n\nThe build functions have to be evaluated by eval.(p2fun_load[\"p2fun\"]) before usage. The arguments for the squared pressure functions are the ordered known squared pressures p^2, the ordered known flow permabilities λ resp. flow restrictions κ, the ordered known flows F and constant A = π256 p_nT_n. \n\n\n\n\n\n","category":"method"},{"location":"docstrings/#GasChromatographySystems.solve_balance-Tuple{Any}","page":"Docstrings","title":"GasChromatographySystems.solve_balance","text":"solve_balance(sys; mode=\"λ\", bal_eq = flow_balance(sys))\n\nSolves the substituted flow balance equations of the capillary system sys for the squared pressures of vertices with undefined pressures as an array of symbolic expressions.\n\nArguments\n\nsys: System structure of the capillary system for which the flow balance is set up.\nmode: Mode for flow equations to use flow permeabilities λ (mode = λ; default) or flow restrictions κ (mode = κ)\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#GasChromatographySystems.substitute_pressure_squared_functions-Tuple{Any, Any}","page":"Docstrings","title":"GasChromatographySystems.substitute_pressure_squared_functions","text":"substitute_pressure_squared_functions(p2fun, sys; mode=\"λ\")\n\nSubstitutes the the pressure functions (solutions to the flow balance equations) with the known quantities of pressures, flows, flow restictions/permabilities.\n\nThis results in an array of pressures p at vertices without defined pressure as function of time t. \n\nArguments\n\np2fun: Julia function of the solutions of the flow balance equations from build_pressure_squared_functions(sys; mode=\"λ\")\nsys: System structure of the capillary system for which the flow balance is set up.\nmode: Mode for flow equations to use flow permeabilities λ (mode = λ; default) or flow restrictions κ (mode = κ)\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#GasChromatographySystems.substitute_unknown_flows-Tuple{Any}","page":"Docstrings","title":"GasChromatographySystems.substitute_unknown_flows","text":"substitute_unknown_flows(sys; mode=\"λ\", bal_eq = flow_balance(sys))\n\nSubstitutes flow equations for undefined flows over edges in the capillary system sys in the flow balance equation system. \n\nFlow equation over edge ji between vertice ìandjwith flow permabilityλ_{i,j}`:\n\nF_ij = A λ_ij left(p_i^2-p_j^2right)\n\nFor flow resistance κ_ij:\n\nF_ij = fracAκ_ij left(p_i^2-p_j^2right)\n\nArguments\n\nsys: System structure of the capillary system for which the flow balance is set up.\nmode: Mode for flow equations to use flow permeabilities λ (mode = λ; default) or flow restrictions κ (mode = κ)\nbal_eq: Array of the symbolic flow balance equations; defaults to flow_balance(sys)\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#GasChromatographySystems.unknown_F-Tuple{Any}","page":"Docstrings","title":"GasChromatographySystems.unknown_F","text":"unknown_F(sys)\n\nExtract the index of the edges of the graph of system sys for which the flows F are not defined. \n\nArguments\n\nsys: System structure of the capillary system for which the flow balance is set up.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#GasChromatographySystems.unknown_p-Tuple{Any}","page":"Docstrings","title":"GasChromatographySystems.unknown_p","text":"unknown_p(sys)\n\nExtract the index of the vertices of the graph of system sys for which the pressure p are not defined. \n\nArguments\n\nsys: System structure of the capillary system for which the flow balance is set up.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#GasChromatographySystems.unknown_λ-Tuple{Any}","page":"Docstrings","title":"GasChromatographySystems.unknown_λ","text":"unknown_λ(sys)\n\nExtract the index of the edges of the graph of system sys for which the flow permeabilities λ are not defined. \n\nArguments\n\nsys: System structure of the capillary system for which the flow balance is set up.\n\n\n\n\n\n","category":"method"}]
}
