# begin - system to parameters
# transform system to GasChromatographySimulator.Parameters
function all_stationary_phases(sys)
	stat_phases = Array{String}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		stat_phases[i] = sys.modules[i].sp
	end
	return stat_phases
end

function common_solutes(db, sys)
	# gives the soultes with common stationary phases between the database 'db' and the
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
solute transport through each module in the system. It handles both column modules and thermal
modulators, setting up temperature programs, pressure functions, and substance parameters.

# Arguments
- `sys`: The GC system structure containing the network of modules
- `p2fun`: Pressure functions for the system
- `db_dataframe`: Database containing solute properties
- `selected_solutes`: List of solutes to include in the simulation

# Keyword Arguments
- `interp`: Whether to use interpolated pressure functions (default: true)
- `dt`: Time step for pressure function interpolation (default: 1)
- `mode`: Mode for flow calculations ("λ" for permeability or "κ" for restriction) (default: "λ")

# Returns
- Array of `GasChromatographySimulator.Parameters` objects, one for each module in the system

# Notes
- Handles both constant and programmed temperature/pressure conditions
- Sets up column parameters including length, diameter, and stationary phase
- Configures temperature programs with interpolation functions
- Loads solute properties from the database for the specified stationary phase
- Applies module-specific options including numerical solver settings
- Supports both ModuleColumn and ModuleTM (thermal modulator) types
"""
function graph_to_parameters(sys, p2fun, db_dataframe, selected_solutes; interp=true, dt=1, mode="λ")
	# ng should be taken from the separat module options 
	E = collect(edges(sys.g))
	srcE = src.(E) # source indices
	dstE = dst.(E) # destination indices
	if interp == true # linear interpolation of pressure functions with step width dt
		p_func = interpolate_pressure_functions(sys, p2fun; dt=dt, mode=mode)
	else
		p_func = pressure_functions(sys, p2fun; mode=mode)
	end
	parameters = Array{GasChromatographySimulator.Parameters}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		# column parameters
		col = GasChromatographySimulator.Column(sys.modules[i].L, sys.modules[i].d, [sys.modules[i].d], sys.modules[i].df, [sys.modules[i].df], sys.modules[i].sp, sys.options.gas)

		# program parameters
		time_steps, temp_steps, gf, a_gf, T_itp = module_temperature(sys.modules[i], sys)
		pin_steps = if typeof(sys.pressurepoints[srcE[i]].P) <: PressureProgram
			sys.pressurepoints[srcE[i]].P.pressure_steps
		else
			fill(sys.pressurepoints[srcE[i]].P, length(time_steps))
		end
		pout_steps = if typeof(sys.pressurepoints[dstE[i]].P) <: PressureProgram
			sys.pressurepoints[dstE[i]].P.pressure_steps
		else
			fill(sys.pressurepoints[dstE[i]].P, length(time_steps))
		end
		pin_itp = p_func[srcE[i]]
		pout_itp = p_func[dstE[i]]	
		
		prog = GasChromatographySimulator.Program(time_steps, temp_steps, pin_steps, pout_steps, gf, a_gf, T_itp, pin_itp, pout_itp)

		# substance parameters
		sub = GasChromatographySimulator.load_solute_database(db_dataframe, sys.modules[i].sp, sys.options.gas, selected_solutes, NaN.*ones(length(selected_solutes)), NaN.*ones(length(selected_solutes)))

		# option parameters
		#opt = if typeof(sys.modules[i]) == GasChromatographySystems.ModuleTM 
		opt = GasChromatographySimulator.Options(alg=sys.modules[i].opt.alg, abstol=sys.modules[i].opt.abstol, reltol=sys.modules[i].opt.reltol, Tcontrol=sys.modules[i].opt.Tcontrol, odesys=sys.options.odesys, ng=sys.modules[i].opt.ng, vis=sys.options.vis, control=sys.options.control, k_th=sys.options.k_th)
		#else
			# put options for ModuleColumn in separat options structure?
		#	GasChromatographySimulator.Options(alg=sys.modules[i].opt.alg, abstol=sys.modules[i].opt.abstol, reltol=sys.modules[i].opt.reltol, Tcontrol=sys.modules[i].opt.Tcontrol, odesys=sys.options.odesys, ng=sys.modules[i].opt.ng, vis=sys.options.vis, control=sys.options.control, k_th=sys.options.k_th)
		#end

		parameters[i] = GasChromatographySimulator.Parameters(col, prog, sub, opt)
	end
	return parameters
end

# end - system to parameters
