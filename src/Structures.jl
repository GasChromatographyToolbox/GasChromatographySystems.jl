# definition of structures and related methods

# structures containing options
"""
    Options(; gas="He", odesys=true, vis="Blumberg", control="Pressure", k_th=1e12)

Structure describing general options for the simulation of the system. 

# Arguments
* `gas`: The type of gas used as mobile phase in the gas chromatographic system. Allowed values: He, H2 or N2.
* `odesys`: Combine the ODEs for migration and peak-width into a system of ODEs (`odesys = true`) or solve the two ODEs separately (`odesys = false`).
* `vis`: Used model of viscosity. `HP` is a second-order polynomial taken from the HP flow calculator. `Blumberg` is an emperical formula according to the book
    `Temperature-programmed Gas Chromatography` by Leonid M. Blumberg (2010, Wiley-VCH).
* `control`: Control of the "Flow" or of the "Pressure" (at column inlet) during the program.
* `k_th`: Threshold for the maximum of the retention factor. If the calculated retention factor is bigger than `k_th` than the retention factor is set to the value `k_th`.
    This is done to avoid to small step widths in the solver for highly retained soultes at the beginning of a GC program. 
"""
struct Options
    gas::String# gas of the mobile phase
    odesys::Bool  		# calculate the two ODEs (migration and peak-width) separately (false) or 
                        # combined as a system of ODEs (true)                        
    vis::String         # viscosity model 'HP' or 'Blumberg'
    control::String     # control of the 'Flow' or of the inlet 'Pressure' during the program
    k_th                # threshold for the max. possible retention factor
end

function Options(;gas="He", odesys=true, vis="Blumberg", control="Pressure", k_th=1e12)
    opt = Options(gas, odesys, vis, control, k_th)
    return opt
end

"""
    ModuleColumnOptions(; alg=OwrenZen5(), abstol=1e-8, reltol=1e-5, ng=false, Tcontrol="inlet")

Structure describing options for column modules. 

# Arguments
* `alg`: The algorithm used for the ODE solver. The algorithms `OwrenZen3()`, `OwrenZen4()` and `OwrenZen5()` are recommended.
* `abstol`: The absolute tolerance for the ODE solver. Recommended value 1e-6 to 1e-8.
* `reltol`: The relative tolerance for the ODE solver. Recommended value 1e-3 to 1e-5. 
* `ng`: Option to calculate the simulation without a gradient (`ng = true`) or with a gradient (`ng = false`).
* `Tcontrol`: Option defining at which point of the column the temperature program is calculated. The options are `inlet` (x=0) and `outlet` (x=L).
"""
struct ModuleColumnOptions
	alg                 # algorithmen for the ODE solver
    abstol              # absolute tolerance for ODE solver
    reltol              # relative tolerance for ODE solver 
	ng::Bool            # non-gradient calculation, ignores a defined spatial change of d, df or T
	Tcontrol::String    # temperature control at 'inlet' (top) or 'outlet' (bottom) of the column
end

function ModuleColumnOptions(; alg=OwrenZen5(), abstol=1e-8, reltol=1e-6, ng=false, Tcontrol="inlet")
	ModuleColumnOptions(alg, abstol, reltol, ng, Tcontrol)
end

"""
    ModuleTMOptions(; Tcold_abs=true, sflank=40, tflank=20, alg=Vern9(), abstol=1e-10, reltol=1e-8, dtinit=1e-6, ng=true, Tcontrol="inlet")

Structure describing options for thermal modulator modules. 

# Arguments
* `Tcold_abs`: Calculate the low temperature at the modulator as the absolute value of the defined parameter `Tcold` (`Tcold_abs = true`) or as relative temperature difference from the defined oven temperature (`Tcold_abs = false`).
* `sflank`: Flank factor of the smoothed rectangle temperature function in space. A higher factor results in a steeper slope at the edges of the modulator point. Values between 12 and 100, or Inf if `ng = true`. Recommend `sflank = 40`. Only relevante, if option `ng = false`. 
* `tflank`: Flank factor of the smoothed rectangle temperature function in time. A higher factor results in a steeper slope at the begining and end of the hot jet. Values between 12 and 100. Recommend `tflank = 20`.
* `alg`: The algorithm used for the ODE solver. For the thermal modulator module the algorithm `Vern9()` is recommend.
* `abstol`: The absolute tolerance for the ODE solver. Recommended value 1e-6 to 1e-12.
* `reltol`: The relative tolerance for the ODE solver. Recommended value 1e-3 to 1e-10.
* `dtinit`: The initial step width for the ODE solver. A value of `L*1e-6`, with `L` the length of the modulator module, is recommend. 
* `ng`: Option to calculate the simulation without a gradient (`ng = true`) or with a gradient (`ng = false`).
* `Tcontrol`: Option defining at which point of the column the temperature program is calculated. The options are `inlet` (x=0) and `outlet` (x=L).
"""
struct ModuleTMOptions
	Tcold_abs::Bool 	# Tcold as absolute value (`true`) or as relative difference to the program temperature (`false`) 
	sflank::Float64 	# flank factor for the spatial smoothed rectangle, values between 12 and 100, or Inf if `ng = true`. 
	tflank::Float64 	# flank factor for the temporal smoothed rectangle, values between 12 and 100
	alg 				# solver algorithm for the simulation, typical `Vern9()`, `OwrenZen5()`. With the option "simplifiedTM" no solving of ODEs is used for the modulator spot but an approximation is used, assuming a rectangle modulation.
	abstol::Float64 	# absolute tolerance for the solving of the ODEs, typical value 1e-10
	reltol::Float64 	# relativ tolerance for solving of the ODEs, typical value 1e-8
	dtinit::Float64 	# initial step width for the solving of the ODEs, typical value `module length * 1e-6
	ng::Bool 			# model modulation spot also as smoothed rectangle over the length (`false`) or as uniform (`true`)
	Tcontrol::String    # temperature control at 'inlet' (top) or 'outlet' (bottom) of the column
end

function ModuleTMOptions(; Tcold_abs=true, sflank=40, tflank=20, alg=Vern9(), abstol=1e-10, reltol=1e-8, dtinit=1e-6, ng=true, Tcontrol="inlet")
	ModuleTMOptions(Tcold_abs, sflank, tflank, alg, abstol, reltol, dtinit, ng, Tcontrol)
end



# structures of the system modules 
abstract type AbstractModule end

"""
    ModuleColumn(; Tcold_abs=true, sflank=40, tflank=20, alg=Vern9(), abstol=1e-10, reltol=1e-8, dtinit=1e-6, ng=true, Tcontrol="inlet")

Structure describing the module of a gas chromatographic column, the basic buliding block of a gas chromatographic system. These are edges of the graph representation of a GC system. 

# Arguments
* `name`: Name of the module.
* `L`: Length of the column in m.
* `d`: Internal diameter of the column in m. This can be a function of column position, if a non-uniform diameter is used (not tested yet).
* `a_d`: Parameters of the diameter function `d`. Empty, if `d` is not a function. 
* `df`: Film thickness of the column in m. This can be a function of column position, if a non-uniform film thickness is used (not tested yet).
* `a_df`: Parameters of the film thickness function `df`. Empty, if `df` is not a function. 
* `sp`: The name of the stationary phase.
* `T`: Temperature of the module. This can be a value (in °C) or a temperature program as a `TemperatureProgram` structure. 
* `F`: Flow of the mobile phase through the column in mL/min. This can be an number, a fuction or a NaN (if it is not defined yet an should be calculated based on values from other modules or pressure points in the system).
* `opt`: Options for this module. This is a structure of the form `ModuleColumnOptions`.

Four methods to construct this structure exist:
* `ModuleColumn(name, L, d, df, sp, T, opt::ModuleColumnOptions)`: The flow is not defined, uniform `d` and `df`.
* `ModuleColumn(name, L, d, df, sp, T, F, opt::ModuleColumnOptions)`: The flow is defined by `F` and uniform `d` and `df`.
* `ModuleColumn(name, L, d, df, sp, T; alg=OwrenZen5(), abstol=1e-8, reltol=1e-6, ng=false, Tcontrol="inlet")`: The flow is not defined, uniform `d` and `df`. Options according to `ModuleColumnOptions`.
* `ModuleColumn(name, L, d, df, sp, tp, flow; alg=OwrenZen5(), abstol=1e-8, reltol=1e-6, ng=false, Tcontrol="inlet")`: The flow is defined by `F` and uniform `d` and `df`. Options according to `ModuleColumnOptions`.
"""
struct ModuleColumn<:GasChromatographySystems.AbstractModule
	# Module
	# GC column, gradients are possible
	name::String
	L::Float64
	d#::Fd # Function
	a_d::Array{Float64,1} # Parameters of diameter function, just for information
	df#::Fdf # Function
	a_df::Array{Float64,1} # Parameters of film_thickness function, just for information
	sp::String
	T # a number (constant temperature) or a TemperatureProgram structure
	F # an number (constant flow) or a Function
	opt::ModuleColumnOptions
end

function ModuleColumn(name, L, d, df, sp, T, opt::ModuleColumnOptions)
	# function to construct the Column structure
	# for the case of constant diameter and constant film thickness
	# and undefined flow
	col = ModuleColumn(name, L, d, [d], df, [df], sp, T, NaN, opt)
	return col
end

function ModuleColumn(name, L, d, df, sp, T, F, opt::ModuleColumnOptions)
	# function to construct the Column structure
	# for the case of constant diameter and constant film thickness
	col = ModuleColumn(name, L, d, [d], df, [df], sp, T, F, opt)
	return col
end

function ModuleColumn(name, L, d, df, sp, tp; alg=OwrenZen5(), abstol=1e-8, reltol=1e-6, ng=false, Tcontrol="inlet")
	opt = ModuleColumnOptions(; alg=alg, abstol=abstol, reltol=reltol, ng=ng, Tcontrol=Tcontrol)
	TM = ModuleColumn(name, L, d, [d], df, [df], sp, tp, NaN, opt)
	return TM
end

function ModuleColumn(name, L, d, df, sp, tp, flow; alg=OwrenZen5(), abstol=1e-8, reltol=1e-6, ng=false, Tcontrol="inlet")
	opt = ModuleColumnOptions(; alg=alg, abstol=abstol, reltol=reltol, ng=ng, Tcontrol=Tcontrol)
	TM = ModuleColumn(name, L, d, [d], df, [df], sp, tp, flow, opt)
	return TM
end

"""
    ModuleTM(; Tcold_abs=true, sflank=40, tflank=20, alg=Vern9(), abstol=1e-10, reltol=1e-8, dtinit=1e-6, ng=true, Tcontrol="inlet")

Structure describing the module of a thermal modulator, a short section of a column with a periodic repeating temperature modulation between a cold and hot jet. These are edges of the graph representation of a GC system.

# Arguments
* `name`: Name of the module.
* `L`: Length of the column in m.
* `d`: Internal diameter of the column in m. This can be a function of column position, if a non-uniform diameter is used (not tested yet).
* `a_d`: Parameters of the diameter function `d`. Empty, if `d` is not a function. 
* `df`: Film thickness of the column in m. This can be a function of column position, if a non-uniform film thickness is used (not tested yet).
* `a_df`: Parameters of the film thickness function `df`. Empty, if `df` is not a function. 
* `sp`: The name of the stationary phase.
* `T`: Basic temperature of the module. On top of this temperature the modulation occurs. This can be a value (in °C) or a temperature program as a `TemperatureProgram` structure. 
* `shift`: Time shift in s of the periodic modulation relative to the start of the chromatogram.
* `PM`: Modulation period in s.
* `ratio`: The ratio of the duration of cold and hot jet.
* `Thot`: Temperature in °C, by which `T` is increased while the hot jet is active. 
* `Tcold`: Temperature in °C, by which `T` is decrease while the cold jet is active (if `Tcold_abs = false` in `ModuleColumnOptions`) or to which the column is cooled down (if `Tcold_abs = true` in `ModuleTMOptions`)
* `F`: Flow of the mobile phase through the column in mL/min. This can be an number, a fuction or a NaN (if it is not defined yet an should be calculated based on values from other modules or pressure points in the system).
* `opt`: Options for this module. This is a structure of the form `ModuleTMOptions`.

Four methods to construct this structure exist:
* `ModuleTM(name, L, d, df, sp, T, shift, PM, ratio, Thot, Tcold, opt::ModuleTMOptions)`: The flow is not defined, uniform `d` and `df`.
* `ModuleTM(name, L, d, df, sp, T, shift, PM, ratio, Thot, Tcold, F, opt::ModuleTMOptions)`: The flow is defined by `F` and uniform `d` and `df`.
* `ModuleTM(name, L, d, df, sp, T, shift, PM, ratio, Thot, Tcold; Tcold_abs=true, sflank=40, tflank=20, alg=Vern9(), abstol=1e-10, reltol=1e-8, dtinit=1e-6, ng=false, Tcontrol="inlet")`: The flow is not defined, uniform `d` and `df`. Options according to `ModuleTMOptions`.
* `ModuleTM(name, L, d, df, sp, tp, shift, pm, ratio, Thot, Tcold, F; Tcold_abs=true, sflank=40, tflank=20, alg=Vern9(), abstol=1e-10, reltol=1e-8, dtinit=1e-6, ng=false, Tcontrol="inlet")`: The flow is defined by `F` and uniform `d` and `df`. Options according to `ModuleTMOptions`.
"""
struct ModuleTM<:GasChromatographySystems.AbstractModule
	# Module
	# thermal modulator
	name::String
	L::Float64
	d#::Fd # Function
	a_d::Array{Float64,1} # Parameters of diameter function, just for information
	df#::Fdf # Function
	a_df::Array{Float64,1} # Parameters of film_thickness function, just for information
	sp::String
	T # a number (constant temperature) or a TemperatureProgram structure
	shift::Float64
	PM::Float64 # a number, modulation periode 
	ratio::Float64 # a number, ratio of the duration between cold and hot jet, approx. as rectangular function
	Thot::Float64 # heating with hot jet
	Tcold::Float64 # cooling with cold jet
	F # an number (constant flow) or a Function
	opt::ModuleTMOptions
end

function ModuleTM(name, L, d, df, sp, T, shift, PM, ratio, Thot, Tcold, opt::ModuleTMOptions)
	TM = ModuleTM(name, L, d, [d], df, [df], sp, T, shift, PM, ratio, Thot, Tcold, NaN, opt)
	return TM
end

function ModuleTM(name, L, d, df, sp, T, shift, PM, ratio, Thot, Tcold, F, opt::ModuleTMOptions)
	TM = ModuleTM(name, L, d, [d], df, [df], sp, T, shift, PM, ratio, Thot, Tcold, F, opt)
	return TM
end

function ModuleTM(name, L, d, df, sp, T, shift, PM, ratio, Thot, Tcold; Tcold_abs=true, sflank=40, tflank=20, alg=Vern9(), abstol=1e-10, reltol=1e-8, dtinit=1e-6, ng=false, Tcontrol="inlet")
	opt = ModuleTMOptions(; Tcold_abs=Tcold_abs, sflank=sflank, tflank=tflank, alg=alg, abstol=abstol, reltol=reltol, dtinit=dtinit, ng=ng, Tcontrol=Tcontrol)
	TM = ModuleTM(name, L, d, df, sp, T, shift, PM, ratio, Thot, Tcold, NaN, opt)
	return TM
end

function ModuleTM(name, L, d, df, sp, tp, shift, pm, ratio, Thot, Tcold, F; Tcold_abs=true, sflank=40, tflank=20, alg=Vern9(), abstol=1e-10, reltol=1e-8, dtinit=1e-6, ng=false, Tcontrol="inlet")
	opt = ModuleTMOptions(; Tcold_abs=Tcold_abs, sflank=sflank, tflank=tflank, alg=alg, abstol=abstol, reltol=reltol, dtinit=dtinit, ng=ng, Tcontrol=Tcontrol)
	TM = ModuleTM(name, L, d, df, sp, tp, shift, pm, ratio, Thot, Tcold, F, opt)
	return TM
end

# add a flow modulator module

# temperature program structure
"""
    TemperatureProgram(time_steps, temp_steps, gf, a_gf)

Structure describing the temperature program. 

# Arguments
* `time_steps`: Time steps in s, after which the corresponding temperature in `temp_steps` is reached. 
* `temp_steps`: Temperature steps in °C.
* `gf`: Gradient function `gf(x, a_gf)`, describes the thermal gradient.
* `a_gf`:Parameters of the gradient function.

One method to construct this structure exist for a uniform temperatur program (same for all column positions `x`):
* `TemperatureProgram(time_steps, temp_steps)`

A default temperature program is avaliable:
* `default_TP()`: Heating from 40°C to 340°C in 1800 s (30 min). 
"""
struct TemperatureProgram{F<:Function}
    time_steps::Array{<:Real,1}
    temp_steps::Array{<:Real,1}
    gf::F
    a_gf::Array{<:Real,2} # Parameters of the gradient_function, just for information
    TemperatureProgram(ts,Ts,gf,a) = (length(ts)!=length(Ts) || length(ts)!=length(gf(0.0))) ? error("Mismatch between length(time_steps) = $(length(ts)), length(temp_steps) = $(length(Ts)) and length(gf(0.0)) = $(length(gf(0.0)))") : new{typeof(gf)}(ts,Ts,gf,a)
end

function TemperatureProgram(time_steps, temp_steps)
    # function to construct the TEmperature Program structure
    # without a thermal gradient
    # using as gradient function the exponential model 'gf_exp(x,a_gf,Tcontrol)'
    a_gf = [zeros(length(time_steps)) zeros(length(time_steps)) ones(length(time_steps)) zeros(length(time_steps))]
    gf(x) = GasChromatographySimulator.gradient(x, a_gf)
    prog = TemperatureProgram(time_steps, temp_steps, gf, a_gf)
    return prog
end

default_TP() = GasChromatographySystems.TemperatureProgram([0.0, 1800.0], [40.0, 340.0])

# pressure point structure
"""
    PressurePoint(name, time_steps, pressure_steps)

Structure describing the pressure program at the connection points of the modules. These are the vertices of the graph representation of a GC system.

# Arguments
* `name`: Name of the `PressurePoint`.
* `time_steps`: Time steps in , after which the corresponding pressure in `pressure_steps` is reached. 
* `pressure_steps`: Pressure steps in kPa(a)?.
"""
struct PressurePoint
	# Pressure program, same structure for inlet and outlet
	name::String
	time_steps::Array{Float64,1}
	pressure_steps::Array{Float64,1}
	PressurePoint(n,ts,ps) = length(ts)!=length(ps) ? error("Mismatch between length(time_steps) = $(length(ts)), length(pressure_steps) = $(length(ps))") : new(n,ts,ps)
end

# system structure
"""
    System(g, pressurepoints, modules, options)

Structure describing a GC system. 

# Arguments
* `g`: The graph representation of the GC system, using `SimpleDiGraph` from `Graphs.jl`. 
* `pressurepoints`: The vertices of graph `g` and their pressure programs. These are structures of type `PressurePoint`. 
* `modules`: The edges of graph `g`. These are column segments of type `ModuleColumn` or thermal modulator points of type `ModuleTM`.

# Examples
## Simple GC:

Definition of the graph:
```julia
    g = SimpleDiGraph(1)
    add_edge!(g, 1, 2)
```

Definition of the two pressure points, here column inlet and outlet with constant pressure:
```julia
    pp = Array{GasChromatographySystems.PressurePoint}(undef, nv(g))
    pp[1] = GasChromatographySystems.PressurePoint("p_in", [0.0, 3600.0], [pin, pin])
    pp[2] = GasChromatographySystems.PressurePoint("p_out", [0.0, 3600.0], [pout, pout])
```

Definition of the column module with default temperature program `default_TP()`, flow is unknown and is calculated from the pressures, default module options:
```julia
    modules = Array{GasChromatographySystems.AbstractModule}(undef, ne(g))  
    modules[1] = GasChromatographySystems.ModuleColumn("column", 30.0, 0.25*1e-3, 0.25*1e-6, "Rxi5ms", GasChromatographySystems.default_TP())
```

Combination to construct the System using the function `update_system()` and default system options:
```julia
    sys = GasChromatographySystems.update_system(GasChromatographySystems.System(g, pp, modules, GasChromatographySystems.Options()))
```

## Series of `n` columns:

Definition of the graph:
```julia
    g = SimpleDiGraph(n+1)
	for i=1:n
		add_edge!(g, i, i+1) 
	end
```

## Split of Columns:

## GCxGC with thermal modulation and loop:

"""
struct System
	g::Graphs.SimpleDiGraph{Int}
	pressurepoints::Array{PressurePoint}
	modules::Array{AbstractModule}
	options::Options
	System(g_,pressurepoints_,modules_,options_) = nv(g_)!=length(pressurepoints_) || ne(g_)!=length(modules_) ? error("Mismatch between number of nodes ($(nv(g_))) and number of pressure points ($(length(pressurepoints_))) and/or mismatch between number of edges ($(ne(g_))) and number of modules ($(length(modules_))).") : new(g_,pressurepoints_,modules_,options_)
end