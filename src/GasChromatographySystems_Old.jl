module GasChromatographySystems

using Reexport
@reexport using GasChromatographySimulator
using Plots
using Intervals

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

struct Column#{Fd<:Function, Fdf<:Function}
    # Module
	# GC column, gradients are possible
	length::Float64
	diameter#::Fd # Function
    a_diameter::Array{Float64,1} # Parameters of diameter function, just for information
	film_thickness#::Fdf # Function
    a_film_thickness::Array{Float64,1} # Parameters of film_thickness function, just for information
	stationary_phase::String
	temperature_program::Temperature_Program
end

# add Cold-Spot as Module

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

function unique_stationary_phases(GCsys)
    # gives the stationary phases of the GC system as an array
	n = length(GCsys)
	sp = String[] 
	for i=1:n
		# test if the module has an entry 'stationary_phase' (Moduls 'Transferline' and 'Column')
		if (typeof(GCsys[i])==Transferline || typeof(GCsys[i])==Column)
			push!(sp, GCsys[i].stationary_phase)
		end
	end
	sp
	return unique(sp)
end

function common_solutes(db, GCsys)
	# gives the soultes with common stationary phases between the database 'db' and the
	# GC-System 'GCsys'
	usp = setdiff(unique_stationary_phases(GCsys), [""])
	if length(usp)==0 # no stationary phase
		common_solutes = db.Name
	else
		filter_db = Array{DataFrame}(undef, length(usp))
		for i=1:length(usp)
			filter_db[i] = filter([:Phase] => x -> x==usp[i], db)
		end
		if length(usp)==1 # only one stationary phase
			common_db = filter_db[1]
		else
			common_db = innerjoin(filter_db[1], filter_db[2], on=:Name, makeunique=true)
			if length(usp)>2 # more than two stationary phases
				for i=3:length(usp)
				common_db = innerjoin(common_db, filter_db[i], on=:Name, makeunique=true)
				end
			end
		end
		common_solutes = common_db.Name
	end	
	return common_solutes
end
#---end-several test functions------------

#---translate-GCsys-into-GasChromatographySimulator.Parameters-----
function length_vector(GCsys)
	L = Float64[]
	for i=1:length(GCsys)
		if typeof(GCsys[i])!=Pressure_Point
			push!(L, GCsys[i].length)
		end
	end
	return L
end

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
			push!(T_itp, linear_interpolation((nx_i, nt_i), Tmat_i, extrapolation_bc=Flat()))
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
	p_itp = linear_interpolation((nt, ), p, extrapolation_bc=Flat())
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
	new_time_steps = [0.0; diff(sort(unique(ctselements)))]
	return new_time_steps
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
		itp = linear_interpolation((cumsum(old_times), ), a_gf[:,i], extrapolation_bc=Flat())
		for j=1:length(new_times)
			new_a_gf[j,i] = itp(cumsum(new_times)[j])
		end
	end
	return new_a_gf
end

"""
	general_step(x, L, a)

A generalized step function

# Arguments
* `x::Float64`: variable (e.g. x-position)
* `L::Array{Float64,1}`: Array with different length of the constant segments
* `a::Array{<:Any,1}`: Array with the values for the segments, can also be
  functions.
"""
function general_step(x::Float64, L::Array{Float64,1}, a::Array{<:Any,1})
    if length(L)!=length(a)
        error("Parameters `L` and `a` must have the same length.")
        return
    end
	intervals = Array{Intervals.Interval}(undef, length(L))
	cumL = round.(cumsum([0; L]), digits=4)
	for i=1:length(intervals)
		if i==length(intervals)
			intervals[i] = Interval{Closed, Closed}(cumL[i],cumL[i+1])
		else
			intervals[i]=Interval{Closed, Open}(cumL[i],cumL[i+1])
		end
		if x in intervals[i]
			return a[i]
		end
	end
end

function pressure_at_moduls(GCsys, Option)
	# new version for multiple Pressure_Points, which can be switched on and off
	T_itp, xshift, L = temperature_interpolation(GCsys, Option)
	T_sys(x,t) = general_step(x, L, TT(x,t,T_itp,xshift))
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
	d_sys(x) = general_step(x, L, dd(x,GCsys,xshift))
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
					T_subsys(x,t) = general_step(x, L[findall(subsys_index.==jj)], TT(x,t,T_itp,xshift)[findall(subsys_index.==jj)])
					#T_subsys = fT_subsys
					d_subsys = x -> general_step(x, L[findall(subsys_index.==jj)], dd(x,GCsys,xshift)[findall(subsys_index.==jj)])
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
		pi_itp[i] = linear_interpolation((cumsum(new_times), ), press[i,:,1], extrapolation_bc=Flat())
		po_itp[i] = linear_interpolation((cumsum(new_times), ), press[i,:,2], extrapolation_bc=Flat())
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

function pp(x, t, T_itp, xshift, piitp, poitp, L, gas, GCsys)
	n = length(L)
	pparray = Array{Float64}(undef, n)
	for i=1:n
		d_func(x) = dd(x, GCsys, xshift)[i] # diameter function and T_itp has to be modified for the x-shift
		if x-xshift[i]<0 || x-xshift[i]>cumsum(L)[i]
			pparray[i] = NaN
		else
			pparray[i] = GasChromatographySimulator.pressure(x-xshift[i], t, T_itp[i], piitp[i], poitp[i], L[i], d_func, gas)
		end
	end
	return pparray
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
	
	TTsys(x,t) = general_step(x, L, TT(x,t,T_itp,xshift))#
	#pi = Pressure_interpolation(GCsys[press_i[1]])
	#po = Pressure_interpolation(GCsys[press_i[end]])
	d_sys(x) = general_step(x, L, dd(x,GCsys,xshift))#
	df_sys(x) = general_step(x, L, ddff(x,GCsys,xshift))#
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
			system = GasChromatographySimulator.Column(L[i], d, a_d, df, a_df, GCsys[ii].stationary_phase, Option.mobile_phase)
		elseif isa(GCsys[ii], Column)
			system = GasChromatographySimulator.Column(L[i], 
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
"""
	change_initial(par, init_t, init_τ)

Change the inital time and peak widths for the substances in a defined GC-system
`par` to the values `init_t` and `, init_τ`.
""" 
function change_initial(par::GasChromatographySimulator.Parameters, init_t, init_τ)
	# copys the parameters `par` and changes the values of par.sub[i].τ₀ and par.sub[i].t₀ to init_τ[i] resp. init_t[i]
	newsub = Array{GasChromatographySimulator.Substance}(undef, length(par.sub))
	for i=1:length(par.sub)
		newsub[i] = GasChromatographySimulator.Substance(par.sub[i].name, par.sub[i].CAS, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀, par.sub[i].ann, par.sub[i].Cag, init_t[i], init_τ[i])
	end
	newpar = GasChromatographySimulator.Parameters(par.col, par.prog, newsub, par.opt)
	return newpar
end

function linear_GC_system_simulation(GCSys, Option, solutes, db)
	parameters = initilize_parameters(GCSys, Option, solutes, db)
	if Option.odesys==true
		solution = Array{Array{Any}}(undef, length(parameters))
		peaklist = Array{DataFrame}(undef, length(parameters))
		#solution[1] = GasChromatographySimulator.solve_system_multithreads(parameters[1])
		peaklist[1], solution[1] = GasChromatographySimulator.simulate(parameters[1])
		for i=2:length(parameters)
			t₀ = Array{Float64}(undef, length(solutes))
			τ₀ = Array{Float64}(undef, length(solutes))
			for j=1:length(solutes)
				t₀[j] = solution[i-1][j].u[end][1]
				τ₀[j] = sqrt(solution[i-1][j].u[end][2])
			end
			newpar = change_initial(parameters[i], t₀, τ₀)
			peaklist[i], solution[i] = GasChromatographySimulator.simulate(newpar)
		end
	end
	return parameters, peaklist, solution
end
#----end-run-the-simulation-for-a-GC-system---

# a non-multithreads-version (should be in GasChromatographySimulator)
function solve_system(par)
	n = length(par.sub)
	sol = Array{Any}(undef, n)
	for i=1:n
		sol[i] = GasChromatographySimulator.solving_odesystem_r(par.col, par.prog, par.sub[i], par.opt)
	end
	return sol
end

function linear_GC_system_simulation_nmt(GCSys, Option, solutes, db)
	parameters = initilize_parameters(GCSys, Option, solutes, db)
	if Option.odesys==true
		solution = Array{Array{Any}}(undef, length(parameters))
		solution[1] = solve_system(parameters[1])
		for i=2:length(parameters)
			t₀ = Array{Float64}(undef, length(solutes))
			τ₀ = Array{Float64}(undef, length(solutes))
			for j=1:length(solutes)
				t₀[j] = solution[i-1][j].u[end][1]
				τ₀[j] = sqrt(solution[i-1][j].u[end][2])
			end
			newpar = change_initial(parameters[i], τ₀, t₀)
			solution[i] = solve_system(newpar)
		end
	end
	return parameters, solution
end

#----begin-plot-functions---------------------
function pressure_plot(GCsys, Option, plot_selector; x₀=0.0, t₀=0.0)
	new_times = new_time_steps(GCsys)
	T_itp, xshift, L = temperature_interpolation(GCsys, Option)
	press = pressure_at_moduls(GCsys, Option)
	piitp, poitp = pressure_at_moduls_itp(L, new_times, press)
	p_func(x,t) = general_step(x, L, pp(x, t, T_itp, xshift, piitp, poitp, L, Option.mobile_phase, GCsys))
	if plot_selector=="p(x)"
		nx = 0.0:sum(L)/1000:sum(L)
		p = p_func.(nx, t₀)
		pplot = plot(nx, p, xlabel="x in m", ylabel="pressure in Pa(a)", legend=:top, label="t=$(t₀)s", size=(800,500))
	elseif plot_selector=="p(t)"
		nt = 0.0:sum(new_times)/1000:sum(new_times)
		p = p_func.(x₀, nt)
		pplot = plot(nt, p, xlabel="t in s", ylabel="pressure in Pa(a)", legend=:top, label="x=$(x₀)m", size=(800,500))
	elseif plot_selector=="p(x,t)"
		nx = 0.0:sum(L)/100:sum(L)
		nt = 0.0:sum(new_times)/100:sum(new_times)
		p = Array{Float64}(undef, length(nx), length(nt))
		for j=1:length(nt)
			for i=1:length(nx)
				p[i,j] = p_func(nx[i], nt[j])
			end
		end
		pplot = plot(nx, nt, p', st=:surface, xlabel="x in m", ylabel="t in s", zlabel="pressure in Pa(a)", size=(800,500))
	end
	return pplot, p
end

function temperature_plot(GCsys, Option, plot_selector; x₀=0.0, t₀=0.0)
	new_times = new_time_steps(GCsys)
	T_itp, xshift, L = temperature_interpolation(GCsys, Option)
	TTsys(x,t) = general_step(x, L, TT(x,t,T_itp,xshift))
	if plot_selector=="T(x)"
		nx = 0.0:sum(L)/1000:sum(L)
		T = TTsys.(nx, t₀).-273.15
		Tplot = plot(nx,  T, xlabel="x in m", ylabel="T in °C", legend=:top, label="t=$(t₀)s", size=(800,500))
	elseif plot_selector=="T(t)"
		nt = 0.0:sum(new_times)/1000:sum(new_times)
		T = TTsys.(x₀, nt).-273.15
		Tplot = plot(nt, T, xlabel="t in s", ylabel="T in °C", legend=:top, label="x=$(x₀)m", size=(800,500))
	elseif plot_selector=="T(x,t)"
		nx = 0.0:sum(L)/100:sum(L)
		nt = 0.0:sum(new_times)/100:sum(new_times)
		T = Array{Float64}(undef, length(nx), length(nt))
		for j=1:length(nt)
			for i=1:length(nx)
				T[i,j] = TTsys(nx[i], nt[j])-273.15
			end
		end
		Tplot = plot(nx, nt, T', st=:surface, xlabel="x in m", ylabel="t in s", zlabel="T in °C", size=(800,500))
	end
	return Tplot, T
end

function flow_plot(GCsys, Option)
	time, flow, κL = flow_system(GCsys, Option)
	pflow = plot(time, flow[:,1], xlabel="t in s", ylabel="flow in mL/min", label="subsys 1", legend=:bottomright)
	if size(flow)[2]>1
		for i=2:size(flow)[2]
			plot!(pflow, time, flow[:,i], label="subsys $(i)")
		end
	end
	return pflow, flow, time
end

function flow_system(GCsys, Option)
	new_times = new_time_steps(GCsys)
	new_pressures = new_pressure_steps(GCsys)
	pp_index = pressure_points_index(GCsys)
	modul_index = modules_index(GCsys)
	tt = 0.0:sum(new_times)/1000:sum(new_times)
	n_subsys_at_new_times = Array{Int}(undef, length(new_times))
	for i=1:length(new_times)
		n_subsys_at_new_times[i] = length(pp_index)-1 # number_of_subsystems(cumsum(new_times)[i], new_times, new_pressures, length(pp_index))
	end
	max_n_subsys = maximum(n_subsys_at_new_times)
	T_itp, xshift, L = temperature_interpolation(GCsys, Option)
	press = pressure_at_moduls(GCsys, Option)
	pi_itp, po_itp = pressure_at_moduls_itp(L, new_times, press)
	F_subsys = fill(NaN, length(tt), max_n_subsys)
	κL_subsys = fill(NaN, length(tt), max_n_subsys)
	for i=1:length(tt)
		F, κL = flow(tt[i], new_times, new_pressures, pp_index, modul_index, pi_itp, po_itp, T_itp, xshift, L, Option.mobile_phase, GCsys)
		for j=1:length(F)
			F_subsys[i,j] = F[j]
			κL_subsys[i,j] = κL[j]
		end
	end
	return tt, F_subsys, κL_subsys
end

function flow(t, new_times, new_pressures, pp_index, modul_index, pi_itp, po_itp, T_itp, xshift, L, gas, GCsys)
	n = length(L) # number of moduls
	n_subsys = length(pp_index)-1#number_of_subsystems(t, new_times, new_pressures, length(pp_index))
	if n_subsys==1
		pi_itp_sys = pi_itp[1]
		po_itp_sys = po_itp[end]
		T_sys(x,t) = general_step(x, L, TT(x,t,T_itp,xshift))
		d_sys(x) = general_step(x, L, dd(x,GCsys,xshift))
		κL = [GasChromatographySimulator.flow_restriction(sum(L), t, T_sys, d_sys, gas)]
		F = [π/256*Tn/pn*(pi_itp_sys(t)^2-po_itp_sys(t)^2)/κL[1]*60*1e6]
		#pin = [pi_itp_sys(t)]
		#pout = [po_itp_sys(t)]
		subsys_index = ones(length(modul_index))
		#p_index_first = [1]
		#p_index_last = [n]
	else
		L_subsys = zeros(n_subsys)
		subsys_index = Array{Int64}(undef, n)
		κL = Array{Float64}(undef, n_subsys)
		F = Array{Float64}(undef, n_subsys)
		#pin = Array{Float64}(undef, n_subsys)
		#pout = Array{Float64}(undef, n_subsys)
		#p_index_first = Array{Float64}(undef, n_subsys)
		#p_index_last = Array{Float64}(undef, n_subsys)
		# time index
		j = findfirst(t.<=cumsum(new_times))
		# indices of pressure_points which are switched on
		pp_on_index = Int[]
		for jj=1:length(pp_index)
			if !isnan(new_pressures[j,jj])
				push!(pp_on_index, pp_index[jj])
			end
		end
		# loop over the subsystems
		for jj=1:n_subsys
			for i=1:n # determine the subsystem to which a modul belongs
				if modul_index[i] in pp_on_index[jj]..pp_on_index[jj+1]
					L_subsys[jj] += L[i]
					subsys_index[i] = jj
				end
			end
			#p_index_first[jj] = findfirst(subsys_index.==jj)
			#p_index_last[jj] = findlast(subsys_index.==jj)
			pi_itp_subsys = pi_itp[findfirst(subsys_index.==jj)]
			po_itp_subsys = po_itp[findlast(subsys_index.==jj)]
			# new temperature and diameter functions of the subsystems:
			fT_subsys(x,t) = general_step(x, L[findall(subsys_index.==jj)], TT(x,t,T_itp,xshift)[findall(subsys_index.==jj)])
			T_subsys = fT_subsys
			d_subsys = x -> general_step(x, L[findall(subsys_index.==jj)], dd(x,GCsys,xshift)[findall(subsys_index.==jj)])
			# calculate the total flow resistance of the subsystem
			κL[jj] = GasChromatographySimulator.flow_restriction(L_subsys[jj], t, T_subsys, d_subsys, gas)
			# calculate the flow of the subsystem
			F[jj] = π/256*Tn/pn*(pi_itp_subsys(t)^2-po_itp_subsys(t)^2)/κL[jj]*60*1e6
			#pin[jj] = pi_itp_subsys(t)
			#pout[jj] = po_itp_subsys(t)
		end
	end
	return F, κL#, pin, pout, subsys_index, p_index_first, p_index_last
end
#----end-plot-functions-----------------------

end # module
