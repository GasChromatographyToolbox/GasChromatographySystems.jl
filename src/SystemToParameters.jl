# begin - system to parameters
# transform system to GasChromatographySimulator.Parameters
function is_simulation_segment(module_::AbstractModule)
	if module_ isa ModuleColumn || module_ isa ModuleTM
		return true
	else
		return false
	end
end

"""
	all_stationary_phases(sys)

Returns the stationary phases of all segments which have a stationary phase (entry 'sp') in the system 'sys'.
# Arguments
- `sys`: The GC system structure containing the network of modules

# Returns
- Array of stationary phases
"""
function all_stationary_phases(sys)
	stat_phases = String[]
	for i=1:ne(sys.g)
		if :sp in fieldnames(typeof(sys.modules[i]))
			push!(stat_phases, sys.modules[i].sp)
		end
	end
	return stat_phases
end

"""
	common_solutes(db, sys)

Returns the solutes with common stationary phases between the database 'db' and the GC-System 'sys'.
# Arguments
- `db`: The database structure containing the solute properties
- `sys`: The GC system structure containing the network of modules

# Returns
- DataFrame of solutes with common stationary phases
"""
function common_solutes(db, sys)
	# gives the solutes with common stationary phases between the database 'db' and the
	# GC-System 'GCsys'
	usp = setdiff(unique(GasChromatographySystems.all_stationary_phases(sys)), [""])
	if length(usp)==0 # no stationary phase
		common_solutes = DataFrame(Name=db.Name, CAS=db.CAS)
	else
		filter_db = Array{DataFrame}(undef, length(usp))
		for i=1:length(usp)
			filter_db[i] = filter([:Phase] => x -> x==usp[i], db)
		end
		if length(usp)==1 # only one stationary phase
			common_db = filter_db[1]
		else
			common_db = innerjoin(filter_db[1], filter_db[2], on=:CAS, makeunique=true)
			if length(usp)>2 # more than two stationary phases
				for i=3:length(usp)
				common_db = innerjoin(common_db, filter_db[i], on=:CAS, makeunique=true)
				end
			end
		end
		CAS = unique(common_db.CAS)
		Name = Array{String}(undef, length(CAS))
		for i=1:length(CAS)
			ii = findfirst(common_db.CAS.==CAS[i])
			Name[i] = common_db.Name[ii]
		end
		common_solutes = DataFrame(Name=Name, CAS=CAS)
	end	
	return common_solutes
end

"""
    graph_to_parameters(sys, p2fun, db_dataframe, selected_solutes; interp=true, dt=1, mode="λ")

Convert a gas chromatography system graph into simulation parameters for each module.

This function processes a GC system graph and generates the necessary parameters for simulating
solute transport through each module in the system. It handles column modules (ModuleColumn) and thermal
modulators (ModuleTM), setting up temperature programs, pressure functions, and substance parameters.
 For valve modules (ModuleValve) placeholder parameters are used, as ModuleValve is not used for simulation.

# Arguments
- `sys`: The GC system structure containing the network of modules
- `p2fun`: Pressure-squared solutions from `build_pressure_squared_functions(sys, solve_balance(sys))` (or an equivalent saved function)
- `db_dataframe`: Database containing solute properties
- `selected_solutes`: List of solutes to include in the simulation

# Keyword Arguments
- `interp`: Whether to use interpolated pressure functions (default: true). On tee graphs with [`ModuleValve`](@ref), uses valve-phase [`steps_interpolation`](@ref) (see [`interpolate_pressure_functions`](@ref)).
- `dt`: Spacing (s) for the auxiliary uniform pressure grid when `interp=true` (default: 1). [`PeriodicValveProgram`](@ref) phase boundaries are always included; `dt ≪ mp` still recommended for resolving slow column/T programs.
- `mode`: Mode for flow calculations ("λ" for permeability or "κ" for restriction) (default: "λ")

# Returns
- Array of `GasChromatographySimulator.Parameters` objects, one for each module in the system

# Notes
- For each module edge, inlet/outlet pressure step vectors (`Fpin_steps`, `pout_steps`) are sampled from the resolved pressure functions at that module's `time_steps`, matching the `Fpin_itp` and `pout_itp` passed to `Program`.
- When this function is used after network pressure-balance (`solve_balance` / `build_pressure_squared_functions`),
  simulation options are forced to `control="Pressure"` because inlet/outlet pressures are already defined by the
  solved pressure field. Using `control="Flow"` here would re-impose flow control on top of pressure-defined programs.
- Handles both constant and programmed temperature/pressure conditions
- Sets up column parameters including length, diameter, and stationary phase
- Configures temperature programs with interpolation functions
- Loads solute properties from the database for the specified stationary phase
- Applies module-specific options including numerical solver settings
- Supports both ModuleColumn and ModuleTM (thermal modulator) types
- ModuleValve is not used for simulation (only placeholder), use 'd_open' as diameter and 0.0 as film thickness and "" as stationary phase. Simulation specific options are set to default values.
"""
function graph_to_parameters(sys, p2fun, db_dataframe, selected_solutes; interp=true, dt=1, mode="λ")
	E = collect(edges(sys.g))
	srcE = src.(E) # source indices
	dstE = dst.(E) # destination indices
	control_mode = "Pressure"
	if sys.options.control != "Pressure"
		@warn "graph_to_parameters: overriding sys.options.control='$(sys.options.control)' to 'Pressure' (pressure-balanced network uses solved pin/pout programs)."
	end
	if interp == true # linear interpolation of pressure functions with step width dt
		p_func = interpolate_pressure_functions(sys, p2fun; dt=dt, mode=mode)
	else
		p_func = pressure_functions(sys, p2fun; mode=mode)
	end
	parameters = Array{GasChromatographySimulator.Parameters}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		# column parameters
		if is_simulation_segment(sys.modules[i])
			col = GasChromatographySimulator.Column(sys.modules[i].L, sys.modules[i].d, [sys.modules[i].d], sys.modules[i].df, [sys.modules[i].df], sys.modules[i].sp, sys.options.gas)
		else
			# ModuleValve is not used for simulation (only placeholder), use 'd_open' as diameter and 0.0 as film thickness and "" as stationary phase
			col = GasChromatographySimulator.Column(sys.modules[i].L, sys.modules[i].d_open, [sys.modules[i].d_open], 0.0, [0.0], "", sys.options.gas)
		end

		# program parameters
		time_steps, temp_steps, gf, a_gf, T_itp = module_temperature(sys.modules[i], sys)
		pin_itp = p_func[srcE[i]]
		pout_itp = p_func[dstE[i]]	
		pin_steps = Float64[pin_itp(t) for t in time_steps]
		pout_steps = Float64[pout_itp(t) for t in time_steps]
		
		prog = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)

		# substance parameters
		n_sub = length(selected_solutes)
		sub = GasChromatographySimulator.load_solute_database(db_dataframe, col.sp, sys.options.gas, selected_solutes, zeros(n_sub), zeros(n_sub))

		# option parameters
		if is_simulation_segment(sys.modules[i])
			opt = GasChromatographySimulator.Options(alg=sys.modules[i].opt.alg, abstol=sys.modules[i].opt.abstol, reltol=sys.modules[i].opt.reltol, Tcontrol=sys.modules[i].opt.Tcontrol, odesys=sys.options.odesys, ng=sys.modules[i].opt.ng, vis=sys.options.vis, control=control_mode, k_th=sys.options.k_th)
		else
			# ModuleValveOptions has only 'ng' entry, use default values for other options, as ModuleValve is not used for simulation (only placeholder)
			opt = GasChromatographySimulator.Options(odesys=sys.options.odesys, ng=sys.modules[i].opt.ng, vis=sys.options.vis, control=control_mode, k_th=sys.options.k_th)
		end

		parameters[i] = GasChromatographySimulator.Parameters(col, prog, sub, opt)
	end
	return parameters
end

# end - system to parameters
