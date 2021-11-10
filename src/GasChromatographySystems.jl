module GasChromatographySystems

using GasChromatographySimulator

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
	press_i = moduls_are_pressure_point(GCsys)
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
	# check which modules have timesteps (Pressure_Point or Column (which has a temperature_program)):
	index = Vector{Int}()
	for i=1:length(GCsys)
		#if (typeof(GCsys[i])==Pressure_Point || hasproperty(GCsys[i], :temperature_program))
		if (typeof(GCsys[i])==Pressure_Point || typeof(GCsys[i])==Column)
			push!(index, i)
		end
	end
	return index
end

function modules_with_diameter(GCsys)
	# check which modules have diameter (Transferline or Column):
	index = Vector{Int}()
	for i=1:length(GCsys)
		#if (typeof(GCsys[i])==Pressure_Point || hasproperty(GCsys[i], :temperature_program))
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

function moduls_not_pressure_point(GCsys)
	index = Vector{Int}()
	for i=1:length(GCsys)
		if !isa(GCsys[i], Pressure_Point)
			push!(index, i)
		end
	end
	return index
end

function moduls_are_pressure_point(GCsys)
	index = Vector{Int}()
	for i=1:length(GCsys)
		if isa(GCsys[i], Pressure_Point)
			push!(index, i)
		end
	end
	return index
end

end # module
