module GasChromatographySystems

using Reexport
@reexport using GasChromatographySimulator
using GasChromatographyTools

# some constants
Tst = 273.15            # K
R = 8.31446261815324    # J mol⁻¹ K⁻¹
Tn = 25.0 + Tst         # K
pn = 101300             # Pa

# structures for the GC-system moduls 
Base.@kwdef struct Options
    # Options for the whole GC system
	mobile_phase::String
	alg
	abstol::Float64
	reltol::Float64
	Tcontrol::String
	odesys::Bool
    #nongrad::Bool = false
end

struct Transferline
    # Module
	# transferline, no gradient
	length::Float64
	diameter::Float64
	film_thickness::Float64
	stationary_phase::String
	temperature::Float64
end

struct Temperature_Program{F<:Function}
	timesteps::Array{Float64,1}
	temperaturesteps::Array{Float64,1}
	gradient_function::F
    a_gradient_function::Array{Float64,2} # Parameters of the gradient_function, just for information
	Temperature_Program(ts,Ts,gf,a) = (length(ts)!=length(Ts) || length(ts)!=length(gf(0.0))) ? error("Mismatch between length(timesteps) = $(length(ts)), length(temperaturesteps) = $(length(Ts)) and length(gradient_function(0.0)) = $(length(gf(0.0)))") : new{typeof(gf)}(ts,Ts,gf,a)
end

struct Column{Fd<:Function, Fdf<:Function}
    # Module
	# GC column, gradients are possible
	length::Float64
	diameter::Fd # Function
    a_diameter::Array{Float64,1} # Parameters of diameter function, just for information
	film_thickness::Fdf # Function
    a_film_thickness::Array{Float64,1} # Parameters of film_thickness function, just for information
	stationary_phase::String
	temperature_program::Temperature_Program
end

struct Pressure_Point
	# Pressure program, same structure for inlet and outlet
	timesteps::Array{Float64,1}
	pressure_steps::Array{Float64,1}
	Pressure_Point(ts,ps) = length(ts)!=length(ps) ? error("Mismatch between length(timesteps) = $(length(ts)), length(pressure_steps) = $(length(ps))") : new(ts,ps)
end
#----------end-structur-------------------

#---several test functions----------------
function test_of_GCsys(GCsys)
	# test function for the GCsys-tupel
	press_i = pressure_points_index(GCsys)
	err = String[]
	if typeof(GCsys[1])==Pressure_Point && typeof(GCsys[end])==Pressure_Point && diff(press_i)>ones(length(press_i)-1)
		test = true
	else
		test = false # -> sollte einen Fehler schmeißen 
		if typeof(GCsys[1])!=Pressure_Point
			push!(err, "First element of GC-system must be a Pressure_Point.")
		end
		if typeof(GCsys[end])!=Pressure_Point
			push!(err, "Last element of GC-system must be a Pressure_Point.")
		end
		if diff(press_i)>ones(length(press_i)-1)
			push!(err, "Between two Pressure_Point elements must be at least one module.")
		end
		error(join(err, " "))
	end
	return test
end

function modules_with_timesteps(GCsys)
	# check which modules have timesteps or temperature_program (which has timesteps):
	index = Vector{Int}()
	for i=1:length(GCsys)
		if (hasproperty(GCsys[i], :timesteps) || hasproperty(GCsys[i], :temperature_program))
			push!(index, i)
		end
	end
	return index
end

function modules_with_diameter(GCsys)
	# check which modules have diameter (Transferline or Column):
	index = Vector{Int}()
	for i=1:length(GCsys)
		if hasproperty(GCsys[i], :diameter)
			push!(index, i)
		end
	end
	return index
end

function modules_with_film_thickness(GCsys)
	# check which modules have film_thickness (Transferline or Column):
	index = Vector{Int}()
	for i=1:length(GCsys)
		if hasproperty(GCsys[i], :film_thickness)
			push!(index, i)
		end
	end
	return index
end

function modules_index(GCsys)
	index = Vector{Int}()
	for i=1:length(GCsys)
		if !isa(GCsys[i], Pressure_Point)
			push!(index, i)
		end
	end
	return index
end

function pressure_points_index(GCsys)
	index = Vector{Int}()
	for i=1:length(GCsys)
		if isa(GCsys[i], Pressure_Point)
			push!(index, i)
		end
	end
	return index
end
#---end-several test functions------------

#---translate-GCsys-into-GasChromatographySimulator.Parameters-----
function temperature_matrix(Module::Transferline)
	nx = [0.0, Module.length]
	nt = [0.0, 1.0e6]
	Tmat = Array{Float64}(undef, length(nx), length(nt))
	for j=1:length(nt)
		for i=1:length(nx)
			Tmat[i,j] = Module.temperature + 273.15
		end
	end
	return Tmat, nx, nt
end

function temperature_matrix(Module::Column)
	T(x) = Module.temperature_program.temperaturesteps .+ Module.temperature_program.gradient_function(x)
	nx = 0.0:Module.length/10000:Module.length # mm exact -> stepwidth 1e-3
	nt = cumsum(Module.temperature_program.timesteps)
	Tmat = Array{Float64}(undef, length(nx), length(nt))
	for j=1:length(nt)
		for i=1:length(nx)
			Tmat[i,j] = T(nx[i])[j] + 273.15
		end
	end
	return Tmat, collect(nx), nt
end

function temperature_interpolation(GCsys, Option)
	L = Float64[]
	Tmat = Array{Float64, 2}[]
	nx = Array{Float64,1}[]
	nt = Array{Float64,1}[]
	xshift = Float64[]
	T_itp = []
	for i=1:length(GCsys)
		if typeof(GCsys[i])!=Pressure_Point
			push!(L, GCsys[i].length)
			Tmat_i, nx_i, nt_i = temperature_matrix(GCsys[i])
			push!(Tmat, Tmat_i)
			push!(nx, nx_i)
			push!(nt, nt_i)
			if length(L)==1
				xstart = 0.0
			else
				xstart = sum(L[1:end-1])
			end
			push!(xshift, xstart)
			push!(T_itp, LinearInterpolation((nx_i, nt_i), Tmat_i, extrapolation_bc=Flat()))
		end
	end
	return T_itp, xshift, L
end

function pressure_interpolation(Modul::Pressure_Point)
	if length(Modul.timesteps)==1
		nt = [0.0, 1.0e6]
		p = Modul.pressure_steps.*ones(2)
	else
		nt = cumsum(Modul.timesteps)
		p = Modul.pressure_steps
	end
	p_itp = LinearInterpolation((nt, ), p, extrapolation_bc=Flat())
	return p_itp
end

function new_time_steps(GCsys)
	# constructs the new timesteps common for all moduls
	index_mod = modules_with_timesteps(GCsys)
	timesteps = Array{Array{Float64},1}(undef, length(index_mod))
	ctselements = Float64[]
	for i=1:length(index_mod)
		try
			timesteps[i] = GCsys[index_mod[i]].timesteps
			for j=1:length(timesteps[i])
				push!(ctselements, cumsum(timesteps[i])[j])
			end
		catch
			timesteps[i] = GCsys[index_mod[i]].temperature_program.timesteps
			for j=1:length(timesteps[i])
				push!(ctselements, cumsum(timesteps[i])[j])
			end
		end
	end
	new_timesteps = [0.0; diff(sort(unique(ctselements)))]
	return new_timesteps
end

function new_temperature_steps(GCsys, Option)
	T_itp, xshift, L = temperature_interpolation(GCsys, Option)
	new_times = new_time_steps(GCsys)
	new_temperatures = Array{Float64}(undef, length(new_times), length(L))
	for j=1:length(L)
		for i=1:length(new_times)
			if Option.Tcontrol == "inlet"
				new_temperatures[i,j] = T_itp[j](0.0, cumsum(new_times)[i])-Tst
			elseif Option.Tcontrol == "outlet"
				new_temperatures[i,j] = T_itp[j](L[j], cumsum(new_times)[i])-Tst
			end
		end
	end
	return new_temperatures
end

function new_pressure_steps(GCsys)
	new_times = new_time_steps(GCsys)
	pp_index = pressure_points_index(GCsys)
	new_pressures = Array{Float64}(undef, length(new_times), length(pp_index))
	for j=1:length(pp_index)
		for i=1:length(new_times)
			ii = findfirst(cumsum(new_times)[i].==cumsum(GCsys[pp_index[j]].timesteps))
			if isa(ii, Int)
				new_pressures[i,j] = GCsys[pp_index[j]].pressure_steps[ii]
			else
				pp_itp = pressure_interpolation(GCsys[pp_index[j]])
				new_pressures[i,j] = pp_itp(cumsum(new_times)[i])
			end
		end
	end
	return new_pressures
end

function new_gradient_parameter_steps(Module::Column, new_times)
	old_times = Module.temperature_program.timesteps
	a_gf = Module.temperature_program.a_gradient_function
	new_a_gf = Array{Float64}(undef, length(new_times), size(a_gf)[2])
	for i=1:size(a_gf)[2]
		itp = LinearInterpolation((cumsum(old_times), ), a_gf[:,i], extrapolation_bc=Flat())
		for j=1:length(new_times)
			new_a_gf[j,i] = itp(cumsum(new_times)[j])
		end
	end
	return new_a_gf
end

function pressure_at_moduls(GCsys, Option)
	# new version for multiple Pressure_Points, which can be switched on and off
	T_itp, xshift, L = temperature_interpolation(GCsys, Option)
	T_sys(x,t) = GasChromatographyTools.general_step(x, L, TT(x,t,T_itp,xshift))
	modul_index = modules_index(GCsys)
	pp_index = pressure_points_index(GCsys)
	new_times = new_time_steps(GCsys)
	new_pressures = new_pressure_steps(GCsys)
	pp_itp = Array{Any}(undef, length(pp_index))
	for i=1:length(pp_index)
		#f_p(t) = pressure_time(t, Pressure_interpolation(GCsys[pp_index[i]]),
		#new_times, new_pressures[:,i])
		# -> pressure_time() was used for cases of pressure-values on 'NaN'
		# 		if Pressure_point is switched off
		# -> for now the case of 'NaN' is not used 
		pp_itp[i] = pressure_interpolation(GCsys[pp_index[i]])
	end
	d_sys(x) = GasChromatographyTools.general_step(x, L, dd(x,GCsys,xshift))
	gas = Option.mobile_phase
	xL = [0.0; cumsum(L)]

	n = length(L) # number of non-PressurePoints moduls
	m = length(new_times) # number of the time steps

	subsys_index = Array{Int64}(undef, n)
	press = Array{Float64}(undef, n, m, 2)
	
	for j=1:m
		# estimate here the number of subsystems
		# n_subsys = number_of_subsystems(cumsum(new_times)[j], new_times,
		# new_pressures, length(pp_index))
		# -> number_of_subsystems() was used for cases of pressure-values on 'NaN'
		# 		if Pressure_point is switched off
		# -> for now the case of 'NaN' is not used 
		n_subsys = length(pp_index)-1
		L_subsys = zeros(n_subsys)
		pi = Array{Any}(undef, n_subsys)
		po = Array{Any}(undef, n_subsys)
		# determine the length of the subsystems, inlet and outlet pressure functions and corresponding subsystem number
		if n_subsys==1
			L_subsys[1] = sum(L)
			pi[1] = pp_itp[1]
			po[1] = pp_itp[end]
			subsys_index = ones(Int, n)
		else
			# indices of pressure_points which are switched on
			pp_on_index = Int[]
			for jj=1:length(pp_index)
				if !isnan(new_pressures[:,jj][j])
					push!(pp_on_index, pp_index[jj])
				end
			end
			# loop over the subsystems
			for jj=1:n_subsys
				for i=1:n # determine the subsystem to which a modul belongs
					if modul_index[i] in pp_on_index[jj]..pp_on_index[jj+1]
						L_subsys[jj] += L[i]
						pi[jj] = pp_itp[jj]
						po[jj] = pp_itp[jj+1]
						subsys_index[i] = jj
					end
				end
			end
		end
		for i=1:n
			for k=1:2
				#println("i=$(i), j=$(j), k=$(k)")
				jj = subsys_index[i]
				if jj>1
					xx = xL[i+k-1]-cumsum(L_subsys)[jj-1]
					# d_sys and T_sys must be modified to account for the shift of the positions in the different subsystems
					T_subsys(x,t) = GasChromatographyTools.general_step(x, L[findall(subsys_index.==jj)], TT(x,t,T_itp,xshift)[findall(subsys_index.==jj)])
					#T_subsys = fT_subsys
					d_subsys = x -> GasChromatographyTools.general_step(x, L[findall(subsys_index.==jj)], dd(x,GCsys,xshift)[findall(subsys_index.==jj)])
				else
					xx = xL[i+k-1]
					d_subsys = d_sys
					T_subsys = T_sys
				end
				#println("xx=$(xx), subsys_index=$(subsys_index[i])") 
				press[i,j,k] = GasChromatographySimulator.pressure(xx,cumsum(new_times)[j], T_subsys, pi[jj], po[jj], L_subsys[jj], d_subsys, gas)	
				
			end
		end
	end
	# k = 1 ... inlet pressures
	# k = 2 ... outlet pressures
	return press
end

function pressure_at_moduls_itp(L, new_times, press) # -> GasChromatographySystems
	pi_itp = Array{Any}(undef, length(L))
	po_itp = Array{Any}(undef, length(L))
	for i=1:length(L)
		pi_itp[i] = LinearInterpolation((cumsum(new_times), ), press[i,:,1], extrapolation_bc=Flat())
		po_itp[i] = LinearInterpolation((cumsum(new_times), ), press[i,:,2], extrapolation_bc=Flat())
	end
	return pi_itp, po_itp
end

function TT(x, t, T_itp, xshift)
	# helper function to collect an array of interpolated temperatures with variables `x`and `t`, `x` is shifted to account the correct x-position of the module in the GC-system
	n = length(xshift) # number of modules
	TTarray = Array{Float64}(undef, n)
	for i=1:n
		TTarray[i] = T_itp[i](x-xshift[i],t)
	end
	return TTarray
end

function dd(x, GCsys, xshift)
	# helper function to collect an array of diameter-functions
	index_d = modules_with_diameter(GCsys)
	n = length(index_d)
	ddarray = Array{Float64}(undef, n)
	for i=1:n
		if isa(GCsys[index_d[i]].diameter, Function)
			ddarray[i] = GCsys[index_d[i]].diameter(x-xshift[i])
		else
			ddarray[i] = GCsys[index_d[i]].diameter
		end
	end
	return ddarray
end

function ddff(x, GCsys, xshift)
	# helper function to collect an array of fim_thickness-functions
	index_d = modules_with_film_thickness(GCsys)
	n = length(index_d)
	ddffarray = Array{Float64}(undef, n)
	for i=1:n
		if isa(GCsys[index_d[i]].film_thickness, Function)
			ddffarray[i] = GCsys[index_d[i]].film_thickness(x-xshift[i])
		else
			ddffarray[i] = GCsys[index_d[i]].film_thickness
		end
	end
	return ddffarray
end

function initilize_parameters(GCsys, Option, solutes, db)
	T_itp, xshift, L = temperature_interpolation(GCsys, Option)#
	#xL = [0.0; cumsum(L)]
	#Lsys = sum(L)
	new_times = new_time_steps(GCsys)#
	new_temperatures = new_temperature_steps(GCsys, Option)#
	#new_pressures = new_pressure_steps(GCsys, Option)
	modul_i = modules_index(GCsys)#
	#press_i = moduls_are_pressure_point(GCsys)
	
	TTsys(x,t) = GasChromatographyTools.general_step(x, L, TT(x,t,T_itp,xshift))#
	#pi = Pressure_interpolation(GCsys[press_i[1]])
	#po = Pressure_interpolation(GCsys[press_i[end]])
	d_sys(x) = GasChromatographyTools.general_step(x, L, dd(x,GCsys,xshift))#
	df_sys(x) = GasChromatographyTools.general_step(x, L, ddff(x,GCsys,xshift))#
	# need a gf_sys function?
	
	press = pressure_at_moduls(GCsys, Option)#
	piitp, poitp = pressure_at_moduls_itp(L, new_times, press)#
	
	parameters = Array{GasChromatographySimulator.Parameters}(undef, length(L))
	for i=1:length(L)
		ii = modul_i[i]
		if isa(GCsys[ii], Transferline)
			a_d = [GCsys[ii].diameter]
			d(x) = GasChromatographySimulator.gradient(x,a_d)
			a_df = [GCsys[ii].film_thickness]
			df(x) = GasChromatographySimulator.gradient(x,a_df)
			system = GasChromatographySimulator.System(L[i], d, a_d, df, a_df, GCsys[ii].stationary_phase, Option.mobile_phase)
		elseif isa(GCsys[ii], Column)
			system = GasChromatographySimulator.System(L[i], 
														GCsys[ii].diameter, 
														GCsys[ii].a_diameter, 
														GCsys[ii].film_thickness, 
														GCsys[ii].a_film_thickness, 
														GCsys[ii].stationary_phase, 
														Option.mobile_phase)
		end
		# entrys of timesteps, temperaturesteps, pressuresteps and gradientfunction in ProgramGC2 have only a placeholder function, all the information about the program is in T_itp, piitp, poitp -> check these entrys
		if hasproperty(GCsys[ii], :temperature_program)
			## TODO -> gradient function and its parameters array have to be adapted for new_times_steps!!!
			## perhaps, I have to change the different gf functions into one, where the structure of the parameters
			## decides the actual function
			new_a_gf = new_gradient_parameter_steps(GCsys[ii], new_times)
			gf(x) = GasChromatographySimulator.gradient(x, new_a_gf; Tcontrol=Option.Tcontrol)
			program = GasChromatographySimulator.Program(new_times, 
															new_temperatures[:,i], 
															press[i,:,1], 
															press[i,:,2], 
															gf, 
															new_a_gf,
															T_itp[i], 
															piitp[i], 
															poitp[i])
		elseif typeof(GCsys[ii])==Transferline
			gfzero(x) = GasChromatographySimulator.gradient(x, zeros(length(new_times),1))
			program = GasChromatographySimulator.Program(new_times, 
															GCsys[ii].temperature.*ones(length(new_times)), 
															press[i,:,1], 
															press[i,:,2], 
															gfzero, 
															zeros(length(new_times),1),
															T_itp[i], 
															piitp[i], 
															poitp[i])
		end

		substance = GasChromatographySimulator.load_solute_database(db, GCsys[ii].stationary_phase, Option.mobile_phase, solutes, zeros(length(solutes)), zeros(length(solutes)))
		options = GasChromatographySimulator.Options(alg = Option.alg, 
                                						abstol = Option.abstol, 
                                						reltol = Option.reltol, 
                                						Tcontrol = Option.Tcontrol, 
                                						odesys = Option.odesys)
		parameters[i] = GasChromatographySimulator.Parameters(system, program, substance, options)
	end
	return parameters
end
#----end-translate-GCsys-into-GasChromatographySimulator.Parameters----- 

#----run-the-simulation-for-a-GC-system---
function linear_GC_system_simulation(GCSys, Option, solutes, db; ng=false)
	parameters = initilize_parameters(GCSys, Option, solutes, db)
	if Option.odesys==true
		solution = Array{Array{Any}}(undef, length(parameters))
		solution[1] = GasChromatographySimulator.solve_system_multithreads(parameters[1]; ng=ng)
		for i=2:length(parameters)
			t₀ = Array{Float64}(undef, length(solutes))
			τ₀ = Array{Float64}(undef, length(solutes))
			for j=1:length(solutes)
				t₀[j] = solution[i-1][j].u[end][1]
				τ₀[j] = sqrt(solution[i-1][j].u[end][2])
			end
			newpar = GasChromatographyTools.change_initial(parameters[i], τ₀, t₀)
			solution[i] = GasChromatographySimulator.solve_system_multithreads(newpar; ng=ng)
		end
	end
	return parameters, solution
end
#----end-run-the-simulation-for-a-GC-system---
end # module
