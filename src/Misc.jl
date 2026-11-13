# misc helpers: program sync, update_system, traces

# common programs
function common_timesteps(sys; default=[0.0, 36000.0])
	com_timesteps = []
	for i=1:nv(sys.g)
		if typeof(sys.pressurepoints[i].P) <: PressureProgram
			com_timesteps = GasChromatographySimulator.common_time_steps(com_timesteps, sys.pressurepoints[i].P.time_steps)
		end
	end
	for i=1:ne(sys.g)
		if typeof(sys.modules[i].T) <: TemperatureProgram
			com_timesteps = GasChromatographySimulator.common_time_steps(com_timesteps, sys.modules[i].T.time_steps)
		end
	end
	if isempty(com_timesteps)
		com_timesteps = default
	end
	return com_timesteps
end

function index_modules_with_temperature_program(sys)
	i_tempprog = Int[]
	for i=1:ne(sys.g)
		if typeof(sys.modules[i].T) <: TemperatureProgram
			push!(i_tempprog, i)
		end
	end
	return i_tempprog
end

function index_pressurepoints_with_pressure_program(sys)
	i_pressprog = Int[]
	for i=1:nv(sys.g)
		if typeof(sys.pressurepoints[i].P) <: PressureProgram
			push!(i_pressprog, i)
		end
	end
	return i_pressprog
end

"""
    match_programs(sys)

Synchronize and match temperature and pressure programs across all modules in a gas chromatography system.

This function ensures that all temperature and pressure programs in the system use the same time steps
by interpolating values to common time points. It handles both pressure programs at pressure points
and temperature programs in modules, including thermal gradients.

# Arguments
- `sys`: A `System` structure containing the gas chromatography system configuration

# Returns
- `com_times`: Array of common time steps for all programs
- `new_press_steps`: Array of interpolated pressure values for each pressure point with a pressure program
- `new_temp_steps`: Array of interpolated temperature values for each module with a temperature program
- `new_a_gf`: Array of interpolated gradient parameters for each module with a temperature program
- `i_pressprog`: Indices of pressure points that have pressure programs
- `i_tempprog`: Indices of modules that have temperature programs

# Notes
- For modules with thermal gradients (`ng=false`), the function interpolates all gradient parameters (ΔT, x0, L0, α)
- For modules without thermal gradients (`ng=true`), the gradient parameters are set to default values
- All programs are synchronized to common time steps to ensure consistent simulation
"""
function match_programs(sys)
	com_times = GasChromatographySystems.common_timesteps(sys)
	i_pressprog = GasChromatographySystems.index_pressurepoints_with_pressure_program(sys)
	new_press_steps = Array{Array{Float64,1}}(undef, length(i_pressprog))
	for i=1:length(i_pressprog)
		new_press_steps[i] = GasChromatographySimulator.new_value_steps(sys.pressurepoints[i_pressprog[i]].P.pressure_steps, sys.pressurepoints[i_pressprog[i]].P.time_steps, com_times)
	end
	i_tempprog = GasChromatographySystems.index_modules_with_temperature_program(sys)
	new_temp_steps = Array{Array{Float64,1}}(undef, length(i_tempprog))
	new_a_gf = Array{Array{Float64, 2}}(undef, length(i_tempprog))
	for i=1:length(i_tempprog)
		new_temp_steps[i] = GasChromatographySimulator.new_value_steps(sys.modules[i_tempprog[i]].T.temp_steps, sys.modules[i_tempprog[i]].T.time_steps, com_times)
		if sys.modules[i_tempprog[i]].opt.ng == false
			# with gradient
			ΔT = GasChromatographySimulator.new_value_steps(sys.modules[i_tempprog[i]].T.a_gf[:,1], sys.modules[i_tempprog[i]].T.time_steps, com_times)
			x0 = GasChromatographySimulator.new_value_steps(sys.modules[i_tempprog[i]].T.a_gf[:,2], sys.modules[i_tempprog[i]].T.time_steps, com_times)
			L0 = GasChromatographySimulator.new_value_steps(sys.modules[i_tempprog[i]].T.a_gf[:,3], sys.modules[i_tempprog[i]].T.time_steps, com_times)
			alpha = GasChromatographySimulator.new_value_steps(sys.modules[i_tempprog[i]].T.a_gf[:,4], sys.modules[i_tempprog[i]].T.time_steps, com_times)
			new_a_gf[i] = [ΔT x0 L0 alpha]
		else
			# without gradient
			new_a_gf[i] = [zeros(length(com_times)) zeros(length(com_times)) ones(length(com_times)) zeros(length(com_times))]
		end
	end
	return com_times, new_press_steps, new_temp_steps, new_a_gf, i_pressprog, i_tempprog
end

"""
    update_system(sys)

Update and synchronize all temperature and pressure programs in a gas chromatography system to use common time steps.

This function ensures that all modules and pressure points in the system use synchronized time steps by:
1. Finding common time steps across all programs
2. Interpolating pressure and temperature values to these common time points
3. Creating new pressure points and modules with synchronized programs
4. Preserving constant temperature/pressure values where applicable

# Arguments
- `sys`: A `System` structure containing the gas chromatography system configuration

# Returns
- A new `System` structure with synchronized programs

# Notes
- For pressure points:
  - Constant pressure values are preserved unchanged
  - Pressure programs are interpolated to common time steps
- For modules:
  - Constant temperature values are preserved unchanged
  - Temperature programs are interpolated to common time steps
  - Thermal gradients are handled differently based on the `ng` option:
    - With gradient (`ng=false`): All gradient parameters are interpolated
    - Without gradient (`ng=true`): Uses default gradient parameters
- The function maintains all other module properties (length, diameter, etc.)
- The graph structure and system options remain unchanged
"""
function update_system(sys)
	new_timesteps, new_pressuresteps, new_temperaturesteps, new_a_gf, index_pp_pressprog, index_module_tempprog = match_programs(sys)
	new_pp = Array{GasChromatographySystems.PressurePoint}(undef, nv(sys.g))
	for i=1:nv(sys.g)
		if typeof(sys.pressurepoints[i].P) <: Number
			new_pp[i] = sys.pressurepoints[i]
		elseif typeof(sys.pressurepoints[i].P) <: GasChromatographySystems.PressureProgram
			ii = findfirst(index_pp_pressprog.==i)
			new_presprog = GasChromatographySystems.PressureProgram(new_timesteps, new_pressuresteps[ii])
			new_pp[i] = GasChromatographySystems.PressurePoint(sys.pressurepoints[i].name, new_presprog)
		end
	end
	new_modules = Array{GasChromatographySystems.AbstractModule}(undef, ne(sys.g))
	for i=1:ne(sys.g)
		if typeof(sys.modules[i].T) <: Number
			new_modules[i] = sys.modules[i]
		elseif typeof(sys.modules[i].T) <: GasChromatographySystems.TemperatureProgram
			ii = findfirst(index_module_tempprog.==i)
			if sys.modules[i].opt.ng == false 
				# with gradient
				gf(x) = GasChromatographySimulator.gradient(x, new_a_gf[ii]) 
				new_tp = GasChromatographySystems.TemperatureProgram(new_timesteps, new_temperaturesteps[ii], gf, new_a_gf[ii])
			else 
				# without gradient
				new_tp = GasChromatographySystems.TemperatureProgram(new_timesteps, new_temperaturesteps[ii])
			end
			if typeof(sys.modules[i]) == GasChromatographySystems.ModuleColumn
				new_modules[i] = GasChromatographySystems.ModuleColumn(sys.modules[i].name, sys.modules[i].L, sys.modules[i].d, sys.modules[i].df, sys.modules[i].sp, new_tp, sys.modules[i].F, sys.modules[i].opt)
			elseif typeof(sys.modules[i]) == GasChromatographySystems.ModuleTM
				new_modules[i] = GasChromatographySystems.ModuleTM(sys.modules[i].name, sys.modules[i].L, sys.modules[i].d, sys.modules[i].df, sys.modules[i].sp, new_tp, sys.modules[i].shift, sys.modules[i].PM, sys.modules[i].ratio, sys.modules[i].Thot, sys.modules[i].Tcold, sys.modules[i].F, sys.modules[i].opt)
			end
		end
	end
	new_sys = GasChromatographySystems.System(sys.name, sys.g, new_pp, new_modules, sys.options)
    return new_sys
end

"""
    check_temperature_gradient(T)

Check if a temperature program includes a thermal gradient.

This function determines whether a given temperature specification includes a thermal gradient
by examining the gradient parameters. It handles both constant temperature values and
temperature programs.

# Arguments
- `T`: Either a constant temperature value (Number) or a TemperatureProgram structure

# Returns
- `gradient`: Boolean indicating whether a thermal gradient is present
  - `false` for constant temperature values
  - `false` for temperature programs with no gradient (all gradient parameters zero)
  - `true` for temperature programs with non-zero gradient parameters

# Notes
- For temperature programs, the presence of a gradient is determined by checking if
  the first column of gradient parameters (a_gf[:, 1]) contains any non-zero values
- This assumes that the first parameter of a_gf is non-zero when a temperature gradient exists
- usage for the option ng as: `ng=!check_temperature_gradient(TPs[i])`
"""
function check_temperature_gradient(T)
	if typeof(T) <: Number
		gradient = false
	elseif typeof(T) <: GasChromatographySystems.TemperatureProgram
		if maximum(abs.(T.a_gf[:, 1])) > 0.0
			# assuming, that the first parameter of a_gf is not zero for temperature gradients
			gradient = true
		else
			gradient = false
		end
	end
	return gradient
end

function traces(sol, par, i_select)
	z = sol[i_select].t
	
		tt = Array{Float64}(undef, length(sol[i_select].t))
		ττ = Array{Float64}(undef, length(sol[i_select].t))
		TT = Array{Float64}(undef, length(sol[i_select].t))
		kk = Array{Float64}(undef, length(sol[i_select].t))
		for j=1:length(sol[i_select].t)
			tt[j] = sol[i_select].u[j][1]
			ττ[j] = sol[i_select].u[j][2]
			TT[j] = par.prog.T_itp(sol[i_select].t[j], sol[i_select].u[j][1])
			kk[j] = GasChromatographySimulator.retention_factor(sol[i_select].t[j], sol[i_select].u[j][1], par.col, par.prog, par.sub[i_select], par.opt)
		end
	return trace = DataFrame(z=z, t=tt, τ²=ττ, T=TT, k=kk)
end
