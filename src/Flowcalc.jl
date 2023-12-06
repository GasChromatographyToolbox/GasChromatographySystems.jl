# definition of all functions related to flow calculations


# first construct the flow balance equations only using the flows over the edges
"""
    flow_balance(sys)

Constructing the flow balance equations of the capillary system `sys` in the form of an array of symbolic equations.

For every inner vertice, the sum of ingoing flows (positive) and of outgoing flows (negative) are equated to zero.

```math
\\Sum  F_{in} + \\Sum F_{out} = 0
``` 

# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
"""
function flow_balance(sys)
	@variables F[1:ne(sys.g)]
	inner_V = GasChromatographySystems.inner_vertices(sys.g) # one flow balance equation for every inner node
	bal_eq = Array{Symbolics.Equation}(undef, length(inner_V))
	for i=1:length(inner_V)
		bal_eq[i] = flow_balance(sys.g, inner_V[i], F) ~ 0
	end
	return bal_eq
end

function flow_balance(g, i_n, F)
	@variables A
	E = collect(edges(g)) # all edges
	srcE = src.(E) # source nodes of the edges
	dstE = dst.(E) # destination nodes of the edges
	# find edges, where node `i_n` is the source
	i_src = findall(srcE.==i_n)
	# find edges, where node `i_n` is the destination
	i_dst = findall(dstE.==i_n)
	balance = 0
	for j=1:length(i_dst) # ingoing flows
		balance = balance + F[i_dst[j]]/A
	end
	for j=1:length(i_src) # outgoing flows
		balance = balance - F[i_src[j]]/A
	end
	return balance
end

# second substitute the unknowns in the flow balance equations and add equations for the known flows
"""
    substitute_unknown_flows(sys; mode="λ", bal_eq = flow_balance(sys))

Substitutes flow equations for undefined flows over edges in the capillary system `sys` in the flow balance equation system. 

Flow equation over edge ``j,i`` between vertice `ì`` and ``j`` with flow permability ``λ_{i,j}``:

```math
F_{i,j} = A λ_{i,j} \\left(p_i^2-p_j^2\\right)
```

For flow resistance ``κ_{i,j}``:

```math
F_{i,j} = \\frac{A}{κ_{i,j}} \\left(p_i^2-p_j^2\\right)
```

# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
* `mode`: Mode for flow equations to use flow permeabilities λ (`mode = λ`; default) or flow restrictions κ (`mode = κ`)
* `bal_eq`: Array of the symbolic flow balance equations; defaults to `flow_balance(sys)`
"""
function substitute_unknown_flows(sys; mode="λ", bal_eq = flow_balance(sys))
    if mode == "λ"
	    return substitute_unknown_flows_λ(sys, bal_eq)
    elseif mode == "κ"
        return substitute_unknown_flows_κ(sys, bal_eq)
    else
        error("Unknown `mode` selection. Use `mode = λ` for flow permeabilities or `mode = κ` for flow restrictions.")
	end
end

function substitute_unknown_flows_λ(sys, bal_eq)
	@variables A, P²[1:nv(sys.g)], λ[1:ne(sys.g)], F[1:ne(sys.g)]
	E = collect(edges(sys.g)) # all edges
	srcE = src.(E) # source nodes of the edges
	dstE = dst.(E) # destination nodes of the edges
	i_unknown_F = GasChromatographySystems.unknown_F(sys) # indices of the modules with an unknown flow
	# create dictionary for the substitution of the unknown flows
	sub_dict = Dict()
	for i=1:length(i_unknown_F)
		j = i_unknown_F[i]
		sub_dict = merge(sub_dict, Dict(F[j] => A*λ[j]*(P²[srcE[j]]-P²[dstE[j]])))
	end
	# index of the known flows
	i_known_F = collect(1:length(edges(sys.g)))[Not(i_unknown_F)]
	# substitute the unknown flows in all balance equations
#	bal_eq = flow_balance(sys)
	sub_bal_eq = Array{Equation}(undef, length(bal_eq)+length(i_known_F))
	for i=1:length(bal_eq)
		sub_bal_eq[i] = substitute(bal_eq[i], sub_dict)
	end
	for i=1:length(i_known_F)
		j = i_known_F[i]
		sub_bal_eq[length(bal_eq)+i] = F[j]/A - λ[j]*(P²[srcE[j]]-P²[dstE[j]]) ~ 0
	end
	return sub_bal_eq
end

function substitute_unknown_flows_κ(sys, bal_eq)
	@variables A, P²[1:nv(sys.g)], κ[1:ne(sys.g)], F[1:ne(sys.g)]
	E = collect(edges(sys.g)) # all edges
	srcE = src.(E) # source nodes of the edges
	dstE = dst.(E) # destination nodes of the edges
	i_unknown_F = GasChromatographySystems.unknown_F(sys) # indices of the modules with an unknown flow
	# create dictionary for the substitution of the unknown flows
	sub_dict = Dict()
	for i=1:length(i_unknown_F)
		j = i_unknown_F[i]
		sub_dict = merge(sub_dict, Dict(F[j] => A/κ[j]*(P²[srcE[j]]-P²[dstE[j]])))
	end
	# index of the known flows
	i_known_F = collect(1:length(edges(sys.g)))[Not(i_unknown_F)]
	# substitute the unknown flows in all balance equations
#	bal_eq = GasChromatographySystems.flow_balance(sys)
	sub_bal_eq = Array{Equation}(undef, length(bal_eq)+length(i_known_F))
	for i=1:length(bal_eq)
		sub_bal_eq[i] = substitute(bal_eq[i], sub_dict)
	end
	for i=1:length(i_known_F)
		j = i_known_F[i]
		sub_bal_eq[length(bal_eq)+i] = F[j]/A - (P²[srcE[j]]-P²[dstE[j]])/κ[j] ~ 0
	end
	return sub_bal_eq
end

"""
   unknown_F(sys)

Extract the index of the edges of the graph of system `sys` for which the flows `F` are not defined. 

# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
"""
function unknown_F(sys)
	i_unknown_flows = Int[]
	for i=1:ne(sys.g)
		if isnan(sys.modules[i].F)
			push!(i_unknown_flows, i)
		end
	end
	return i_unknown_flows
end

"""
   unknown_p(sys)

Extract the index of the vertices of the graph of system `sys` for which the pressure `p` are not defined. 

# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
"""
function unknown_p(sys)
	i_unknown_pressures = Int[]
	for i=1:nv(sys.g)
		if typeof(sys.pressurepoints[i].P) <: GasChromatographySystems.PressureProgram
			if isnan(sys.pressurepoints[i].P.pressure_steps[1])
				push!(i_unknown_pressures, i)
			end
		elseif typeof(sys.pressurepoints[i].P) <: Number
			if isnan(sys.pressurepoints[i].P)
				push!(i_unknown_pressures, i)
			end
		end
	end
	return i_unknown_pressures
end

"""
   unknown_λ(sys)

Extract the index of the edges of the graph of system `sys` for which the flow permeabilities `λ` are not defined. 

# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
"""
function unknown_λ(sys)
	i_unknown_permeability = Int[]
	for i=1:ne(sys.g)
		if isnan(sys.modules[i].L) || isnan(sys.modules[i].d)
			push!(i_unknown_permeability, i)
		end
	end
	return i_unknown_permeability
end

# not needed???
#=
function unknowns_in_flow_balances(sys)
	@variables P²[1:nv(sys.g)], λ[1:ne(sys.g)]
	i_unknown_p = GasChromatographySystems.unknown_p(sys)
	F_balance = GasChromatographySystems.flow_balance(sys)
	in_equation = Array{Array{Int64,1}}(undef, length(i_unknown_p))
	touched = zeros(Int, length(F_balance))
	for i=1:length(i_unknown_p)
		var_unknown = Symbolics.get_variables(P²[i_unknown_p[i]])
		in_equation_ = Int[]
		for j=1:length(F_balance)
			var_equation = Symbolics.get_variables(F_balance[j])
			counter_touched = touched[j]
			for k=1:length(var_equation)
				if isequal(var_unknown[1], var_equation[k])
					push!(in_equation_, j)
					counter_touched += 1
				end
			end
			touched[j] = counter_touched
		end
		in_equation[i] = in_equation_
	end
	# in_equation ... Array of Arrays, lists the idices of the flow balances in which the unknown is in it.
	# touched ... total number of different unknowns in the flow balance equations
	return in_equation, touched
end
=#

"""
    solve_balance(sys; mode="λ", bal_eq = flow_balance(sys))

Solves the substituted flow balance equations of the capillary system `sys` for the squared pressures of vertices with undefined pressures as an array of symbolic expressions.

# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
* `mode`: Mode for flow equations to use flow permeabilities λ (`mode = λ`; default) or flow restrictions κ (`mode = κ`)
"""
function solve_balance(sys; mode="λ", bal_eq = flow_balance(sys))
	if mode == "λ"
		return solve_balance_λ(sys, substitute_unknown_flows_λ(sys, bal_eq))
	elseif mode == "κ"
		return solve_balance_κ(sys, substitute_unknown_flows_κ(sys, bal_eq))
	else
        error("Unknown `mode` selection. Use `mode = λ` for flow permeabilities or `mode = κ` for flow restrictions.")
	end
end

function solve_balance_λ(sys, subst_bal_eq_λ) # should be standard
	# unkown permeabilities λ ???
	@variables P²[1:nv(sys.g)]
	i_unknown_p = GasChromatographySystems.unknown_p(sys) # indices of the nodes with unknown pressure
	#i_unknown_F = unknown_F(sys) # indices of the edges with an unknown flow
#	bal_eq = substitute_unknown_flows_λ(sys)
	#num_use_eq = GasChromatographySystems.unknowns_in_flow_balances(sys)[2]
	if length(i_unknown_p) == length(subst_bal_eq_λ)
		a, b, lin = Symbolics.linear_expansion(subst_bal_eq_λ , [P²[i_unknown_p[i]] for i=1:length(i_unknown_p)])
		if lin == true
			sol = (inv(a)*-b)
		else
			error("System of flow balance equations is not linear for the unknown parameters.")
		end
	elseif length(i_unknown_p) > length(subst_bal_eq_λ)
		error("More unknown pressures than flow balance equations.")
	else # loop of the i_unknown_p -> is this really used???
		# identifie inner nodes which are unknown pressures, there index in the inner_vertices() is the index of the flow balance equations to use 
		# this works only for unknown_p which are inner nodes, unknown_p at outer nodes (e.g. inlet pressure), lead to error 
		inner_V = GasChromatographySystems.inner_vertices(sys.g) 
		bal_eq_i = Array{Int}(undef, length(i_unknown_p))
		for i=1:length(i_unknown_p)
			# if the unknown pressure is not an inner node, than add a equation from the end of the list of balance equation, which should be a flow definition.
			if isnothing(findfirst(i_unknown_p[i].==inner_V))
				bal_eq_i[i] = length(subst_bal_eq_λ)-(i+0)
			else
				bal_eq_i[i] = findfirst(i_unknown_p[i].==inner_V)
			end
		end
		a, b, lin = Symbolics.linear_expansion([subst_bal_eq_λ[x] for x in bal_eq_i], [P²[i_unknown_p[i]] for i=1:length(i_unknown_p)])
		if lin == true
			sol = (inv(a)*-b)
		else
			error("System of flow balance equations is not linear for the unknown parameters.")
		end
	end
	return sol#Symbolics.simplify.(sol) # simplification seems to hang or take forever in some cases
end

function solve_balance_κ(sys, subst_bal_eq_κ)
	# unkown permeabilities λ ???
	@variables P²[1:nv(sys.g)]
	i_unknown_p = GasChromatographySystems.unknown_p(sys) # indices of the nodes with unknown pressure
	#i_unknown_F = unknown_F(sys) # indices of the edges with an unknown flow
#	bal_eq = substitute_unknown_flows_κ(sys)
	#num_use_eq = GasChromatographySystems.unknowns_in_flow_balances(sys)[2]
	if length(i_unknown_p) == length(subst_bal_eq_κ)
		a, b, lin = Symbolics.linear_expansion(subst_bal_eq_κ , [P²[i_unknown_p[i]] for i=1:length(i_unknown_p)])
		if lin == true
			sol = (inv(a)*-b)
		else
			error("System of flow balance equations is not linear for the unknown parameters.")
		end
	elseif length(i_unknown_p) > length(subst_bal_eq_κ)
		error("More unknown pressures than flow balance equations.")
	else # loop of the i_unknown_p -> is this really used???
		# identifie inner nodes which are unknown pressures, there index in the inner_vertices() is the index of the flow balance equations to use 
		# this works only for unknown_p which are inner nodes, unknown_p at outer nodes (e.g. inlet pressure), lead to error 
		inner_V = GasChromatographySystems.inner_vertices(sys.g) 
		bal_eq_i = Array{Int}(undef, length(i_unknown_p))
		for i=1:length(i_unknown_p)
			# if the unknown pressure is not an inner node, than add a equation from the end of the list of balance equation, which should be a flow definition.
			if isnothing(findfirst(i_unknown_p[i].==inner_V))
				bal_eq_i[i] = length(subst_bal_eq_κ)-(i+0)
			else
				bal_eq_i[i] = findfirst(i_unknown_p[i].==inner_V)
			end
		end
		a, b, lin = Symbolics.linear_expansion([subst_bal_eq_κ[x] for x in bal_eq_i], [P²[i_unknown_p[i]] for i=1:length(i_unknown_p)])
		if lin == true
			sol = (inv(a)*-b)
		else
			error("System of flow balance equations is not linear for the unknown parameters.")
		end
	end
	return sol#Symbolics.simplify.(sol)
end

"""
	build_pressure_squared_functions(sys; mode="λ")

Construct array of functions of the solutions for the unkown squared pressures of the flow balance equations of the system of capillaries `sys`.

The arguments for the build functions are arrays of the ordered known squared pressures ``p^2``, the ordered known flow permabilities ``λ`` resp. flow restrictions ``κ``, the ordered known flows ``F`` and constant ``A = π/256 p_n/T_n``.  
	
# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
* `mode`: Mode for flow equations to use flow permeabilities λ (`mode = λ`; default) or flow restrictions κ (`mode = κ`)
"""
function build_pressure_squared_functions(sys; mode="λ")
	if mode == "λ"
		return build_pressure_squared_functions_λ(sys, solve_balance(sys; mode="λ"))
	elseif mode == "κ"
		return build_pressure_squared_functions_κ(sys, solve_balance(sys; mode="κ"))
	else
        error("Unknown `mode` selection. Use `mode = λ` for flow permeabilities or `mode = κ` for flow restrictions.")
	end
end

function build_pressure_squared_functions_λ(sys, solutions)
	#build the function for the squared pressures at the unknown vertices from the symbolic expressions of the solutions of the solved balance equations
	@variables A, P²[1:nv(sys.g)], λ[1:ne(sys.g)], F[1:ne(sys.g)]

	i_unknown_p = GasChromatographySystems.unknown_p(sys)
	i_unknown_F = GasChromatographySystems.unknown_F(sys)
	i_unknown_λ = GasChromatographySystems.unknown_λ(sys)

	i_known_p = collect(1:nv(sys.g))[Not(i_unknown_p)]
	i_known_F = collect(1:ne(sys.g))[Not(i_unknown_F)]
	i_known_λ = collect(1:ne(sys.g))[Not(i_unknown_λ)]
	
	#solutions = GasChromatographySystems.solve_balance_λ(sys)
	p_solution = Array{Function}(undef, length(i_unknown_p))
	for i=1:length(i_unknown_p)
		pfun = build_function(solutions[i], [[P²[j] for j=i_known_p]; [λ[j] for j=i_known_λ]; [F[j] for j=i_known_F]; A];
               expression = Val{false},
               target = Symbolics.JuliaTarget(),
               parallel=nothing)
		p_solution[i] = pfun
	end
	return p_solution
end

function build_pressure_squared_functions_κ(sys, solutions)
	#build the function for the squared pressures at the unknown vertices from the symbolic expressions of the solutions of the solved balance equations
	@variables A, P²[1:nv(sys.g)], κ[1:ne(sys.g)], F[1:ne(sys.g)]

	i_unknown_p = GasChromatographySystems.unknown_p(sys)
	i_unknown_F = GasChromatographySystems.unknown_F(sys)
	i_unknown_λ = GasChromatographySystems.unknown_λ(sys)

	i_known_p = collect(1:nv(sys.g))[Not(i_unknown_p)]
	i_known_F = collect(1:ne(sys.g))[Not(i_unknown_F)]
	i_known_λ = collect(1:ne(sys.g))[Not(i_unknown_λ)]
	
	#solutions = GasChromatographySystems.solve_balance_κ(sys)
	p_solution = Array{Function}(undef, length(i_unknown_p))
	for i=1:length(i_unknown_p)
		pfun = build_function(solutions[i], [[P²[j] for j=i_known_p]; [κ[j] for j=i_known_λ]; [F[j] for j=i_known_F]; A];
               expression = Val{false},
               target = Symbolics.JuliaTarget(),
               parallel=nothing)
		p_solution[i] = pfun
	end
	return p_solution
end

"""
	substitute_pressure_squared_functions(p2fun, sys; mode="λ")

Substitutes the the pressure functions (solutions to the flow balance equations) with the known quantities of pressures, flows, flow restictions/permabilities.

This results in an array of pressures ``p`` at vertices without defined pressure as function of time ``t``. 

# Arguments
* `p2fun`: Julia function of the solutions of the flow balance equations from `build_pressure_squared_functions(sys; mode="λ")`
* `sys`: System structure of the capillary system for which the flow balance is set up.
* `mode`: Mode for flow equations to use flow permeabilities λ (`mode = λ`; default) or flow restrictions κ (`mode = κ`)
"""
function substitute_pressure_squared_functions(p2fun, sys; mode="λ")
	if mode == "λ"
		return substitute_pressure_squared_functions_λ(p2fun, sys)
	elseif mode == "κ"
		return substitute_pressure_squared_functions_κ(p2fun, sys)
	else
        error("Unknown `mode` selection. Use `mode = λ` for flow permeabilities or `mode = κ` for flow restrictions.")
	end
end

function substitute_pressure_squared_functions_λ(p2fun, sys)
	# substitute the actual values of the system `sys` into the p^2 functions to result in an array of p(t) functions (unknown pressures)
	i_unknown_p = GasChromatographySystems.unknown_p(sys)
	i_unknown_F = GasChromatographySystems.unknown_F(sys)
	i_unknown_λ = GasChromatographySystems.unknown_λ(sys)
	i_known_p = collect(1:nv(sys.g))[Not(i_unknown_p)]
	i_known_F = collect(1:ne(sys.g))[Not(i_unknown_F)]
	i_known_λ = collect(1:ne(sys.g))[Not(i_unknown_λ)]
	a = π/256 * GasChromatographySystems.Tn/GasChromatographySystems.pn
	λs = GasChromatographySystems.flow_permeabilities(sys)
	p²s = GasChromatographySystems.pressures_squared(sys)
	p_solution = Array{Function}(undef, length(p2fun))
	for i=1:length(p2fun)
		f(t) = sqrt(p2fun[i]([[p²s[j](t) for j=i_known_p]; [λs[j](t) for j=i_known_λ]; [sys.modules[j].F for j=i_known_F]; a]))
		p_solution[i] = f
	end
	return p_solution
end

function substitute_pressure_squared_functions_κ(p2fun, sys)
	# substitute the actual values of the system `sys` into the p^2 functions to result in an array of p(t) functions (unknown pressures)
	i_unknown_p = GasChromatographySystems.unknown_p(sys)
	i_unknown_F = GasChromatographySystems.unknown_F(sys)
	i_unknown_λ = GasChromatographySystems.unknown_λ(sys)
	i_known_p = collect(1:nv(sys.g))[Not(i_unknown_p)]
	i_known_F = collect(1:ne(sys.g))[Not(i_unknown_F)]
	i_known_λ = collect(1:ne(sys.g))[Not(i_unknown_λ)]
	a = π/256 * GasChromatographySystems.Tn/GasChromatographySystems.pn
	κs = GasChromatographySystems.flow_restrictions(sys)
	p²s = GasChromatographySystems.pressures_squared(sys)
	p_solution = Array{Function}(undef, length(p2fun))
	for i=1:length(p2fun)
		f(t) = sqrt(p2fun[i]([[p²s[j](t) for j=i_known_p]; [κs[j](t) for j=i_known_λ]; [sys.modules[j].F for j=i_known_F]; a]))
		p_solution[i] = f
	end
	return p_solution
end

"""
	solve_pressure(sys; mode="λ")

Solves for the unkown pressures of the system of capillaries as functions of time ``t``. It also returns the indices of the vertices wit undefined pressures.

# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
* `mode`: Mode for flow equations to use flow permeabilities λ (`mode = λ`; default) or flow restrictions κ (`mode = κ`)
"""
function solve_pressure(sys; mode="λ")
	p2fun = build_pressure_squared_functions(sys; mode=mode)
	p_solution = substitute_pressure_squared_functions(p2fun, sys; mode=mode)
	i_unknown_p = GasChromatographySystems.unknown_p(sys)
	return p_solution, i_unknown_p
end

#=
function solve_pressure(sys, solutions)
	@variables A, P²[1:nv(sys.g)], λ[1:ne(sys.g)], F[1:ne(sys.g)]
	#balance = flow_balance(sys)
	i_unknown_p = GasChromatographySystems.unknown_p(sys)
	i_unknown_F = GasChromatographySystems.unknown_F(sys)
	i_unknown_λ = unknown_λ(sys)

	i_known_p = collect(1:nv(sys.g))[Not(i_unknown_p)]
	i_known_F = collect(1:ne(sys.g))[Not(i_unknown_F)]
	i_known_λ = collect(1:ne(sys.g))[Not(i_unknown_λ)]
	
	a = π/256 * GasChromatographySystems.Tn/GasChromatographySystems.pn
	λs = GasChromatographySystems.flow_permeabilities(sys)
	p²s = GasChromatographySystems.pressures_squared(sys)
	p_solution = Array{Function}(undef, length(i_unknown_p))
	for i=1:length(i_unknown_p)
		#sub_dict(t) = merge(Dict((P²[Not(i_unknown_p)] => p²s[j](t) for j=setdiff(1:nv(sys.g), i_unknown_p))), Dict((λ[j] => λs[j](t) for j=1:ne(sys.g))), Dict(A => a), Dict(F[j] => sys.modules[j].F for j in i_known_F))
		pfun = build_function(solutions[i], [[P²[j] for j=i_known_p]; [λ[j] for j=i_known_λ]; [F[j] for j=i_known_F]; A];
               expression = Val{false},
               target = Symbolics.JuliaTarget(),
               parallel=nothing)
		f(t) = sqrt(eval(pfun)([[p²s[j](t) for j=i_known_p]; [λs[j](t) for j=i_known_λ]; [sys.modules[j].F for j=i_known_F]; a]))
		p_solution[i] = f
	end
	return p_solution, i_unknown_p
end
=#

"""
	pressure_functions(sys, p2fun)

Collect all pressure functions as functions of time t at the vertices of the capillary system `sys`, either from defined input values or from the solutions of the flow balance equations. 

# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
* `p2fun`: Julia function of the solutions of the flow balance equations from `build_pressure_squared_functions(sys; mode="λ")`
"""
function pressure_functions(sys, p2fun)
	# collect all pressure functions of all the vertices into an array of functions of time t. If the pressures are known (defined by time_steps and pressure_steps), than a linear interpolation 
	# is used, if the pressures are unkown, the symbolic solution functions `p2fun` of these pressures are used with substituted values of the system.
	# This function should be used, if parameters of the system are to be changes, e.g. column length or diameter, but the structure of the system is the same (same grape, same unknown pressures/flows)
	pres = substitute_pressure_squared_functions(p2fun, sys)
	i_unknown_p = GasChromatographySystems.unknown_p(sys)
	#p²s = pressures_squared(sys)
	p_func = Array{Any}(undef, nv(sys.g))
	for i=1:nv(sys.g)
		f = if i in i_unknown_p
			pres[findfirst(i_unknown_p.==i)]
		elseif typeof(sys.pressurepoints[i].P) <: GasChromatographySystems.PressureProgram
			GasChromatographySimulator.steps_interpolation(sys.pressurepoints[i].P.time_steps, identity.(sys.pressurepoints[i].P.pressure_steps))
		elseif typeof(sys.pressurepoints[i].P) <: Number
			GasChromatographySimulator.steps_interpolation([0.0, 36000.0], identity.(fill(sys.pressurepoints[i].P, 2)))
		end
	p_func[i] = f
	end
	return p_func
end

"""
	pressure_functions(sys)

Collect all pressure functions as functions of time t at the vertices of the capillary system `sys`, either from defined input values or from the solutions of the flow balance equations. 

# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
"""
function pressure_functions(sys)
	pres, unk = solve_pressure(sys)
	#p²s = pressures_squared(sys)
	p_func = Array{Any}(undef, nv(sys.g))
	for i=1:nv(sys.g)
		f = if i in unk
			pres[findfirst(unk.==i)]
		elseif typeof(sys.pressurepoints[i].P) <: GasChromatographySystems.PressureProgram
			GasChromatographySimulator.steps_interpolation(sys.pressurepoints[i].P.time_steps, identity.(sys.pressurepoints[i].P.pressure_steps))
		elseif typeof(sys.pressurepoints[i].P) <: Number
			GasChromatographySimulator.steps_interpolation([0.0, 36000.0], identity.(fill(sys.pressurepoints[i].P, 2)))
		end
	p_func[i] = f
	end
	return p_func
end

"""
	interpolate_pressure_functions(sys; dt=1)

Interpolates (linearly) all pressure funtions at the vertices of the system of capillaries `sys` between the time steps `dt`. For the speed of the simulation these interpolated functions are faster than the pure solution functions of the flow balance equations.

# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
* `dt`: time steps, where the original pressure function is evaluated. Inbetween these time steps the pressure function is linearly interpolated. 
"""
function interpolate_pressure_functions(sys; dt=1)
	tsteps = GasChromatographySystems.common_timesteps(sys)
	tend = sum(tsteps)
	if all(degree(sys.g).<3) # no split/merge nodes, straight graph -> pressure at nodes is linear
		trange = cumsum(tsteps)
	else
		trange = 0:dt:tend
	end
	i_unknown_p = GasChromatographySystems.unknown_p(sys)
	p_func = GasChromatographySystems.pressure_functions(sys)
	p_itp = Array{Any}(undef, length(p_func))
	for i=1:length(p_func)
		#if all(isnan.(sys.pressurepoints[i].pressure_steps))
		if i ∈ i_unkonwn_p
			p_itp[i] = LinearInterpolation((trange, ), p_func[i].(trange), extrapolation_bc=Flat())
		else
			p_itp[i] = p_func[i]
		end # no additional interpolation, if the pressure is allready defined by a linear program
	end
	return p_itp
end

"""
	flow_functions(sys)

Collects the flow functions as functions of time t for all edges of the system of capillaries `sys`.

The flow over edge `i => j` is calculated as

```math
F_{i,j} = \\frac{A}{κ_{i,j}} \\left(p_i^2-p_j^2\\right)
```

with flow restriction ``κ_{i,j} = \\int_0^{L_{i,j}} η(T_{i,j})T_{i,j}/d_{i,j}^4 dx``, pressures ``p_i`` resp. ``p_j`` at the vertices ``i`` resp. ``j``, temperature ``T_{i,j}``, capillary length ``L_{i,j}`` and diameter ``d_{i,j}`` of the edge `i => j`.

# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
"""
function flow_functions(sys)
	p_func = pressure_functions(sys)
	F_func = Array{Function}(undef, ne(sys.g))
	E = collect(edges(sys.g))
	srcE = src.(E)
	dstE = dst.(E)
	for i=1:ne(sys.g)
		pin(t) = p_func[srcE[i]](t)
		pout(t) = p_func[dstE[i]](t)
		T_itp = GasChromatographySystems.module_temperature(sys.modules[i], sys)[5]
		f(t) = GasChromatographySimulator.flow(t, T_itp, pin, pout, sys.modules[i].L, sys.modules[i].d, sys.options.gas; ng=sys.modules[i].opt.ng, vis=sys.options.vis, control=sys.options.control)
		F_func[i] = f
	end
	return F_func
end

"""
	holdup_time_functions(sys)

Collects the hold-up time functions as functions of time t for all edges of the system of capillaries `sys`.

The hold-up time over edge `i => j` is calculated as

```math
t_{M_{i,j}} = \\frac{128}{3} η(T_{i,j}) \\frac{L_{i,j}^2}{d_{i,j}^2} \\frac{p_i^3-p_j^3}{\\left(p_i^2-p_j^2\\right)^2}
```
	
with flow restriction, pressures ``p_i`` resp. ``p_j`` at the vertices ``i`` resp. ``j``, temperature ``T_{i,j}``, capillary length ``L_{i,j}`` and diameter ``d_{i,j}`` of the edge `i => j`.
	
# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
"""
function holdup_time_functions(sys)
	p_func = GasChromatographySystems.pressure_functions(sys)
	tM_func = Array{Function}(undef, GasChromatographySystems.ne(sys.g))
	E = collect(GasChromatographySystems.edges(sys.g))
	srcE = GasChromatographySystems.src.(E)
	dstE = GasChromatographySystems.dst.(E)
	for i=1:GasChromatographySystems.ne(sys.g)
		pin(t) = p_func[srcE[i]](t)
		pout(t) = p_func[dstE[i]](t)
		T_itp = GasChromatographySystems.module_temperature(sys.modules[i], sys)[5]
		f(t) = GasChromatographySimulator.holdup_time(t, T_itp, pin, pout, sys.modules[i].L, sys.modules[i].d, sys.options.gas; ng=sys.modules[i].opt.ng, vis=sys.options.vis, control=sys.options.control)
		tM_func[i] = f
	end
	return tM_func
end

function holdup_time_functions(sys, p2fun)
	# collecting the hold-up time functions of every edge as function of time t for system `sys` and the squared pressure solution functions `p2fun`. This function should be used, if parameters of the system are to be changes, e.g. column length or diameter, but the structure of the system is the same (same grape, same unknown pressures/flows)
	p_func = pressure_functions(sys, p2fun)
	tM_func = Array{Function}(undef, GasChromatographySystems.ne(sys.g))
	E = collect(GasChromatographySystems.edges(sys.g))
	srcE = GasChromatographySystems.src.(E)
	dstE = GasChromatographySystems.dst.(E)
	for i=1:GasChromatographySystems.ne(sys.g)
		pin(t) = p_func[srcE[i]](t)
		pout(t) = p_func[dstE[i]](t)
		T_itp = GasChromatographySystems.module_temperature(sys.modules[i], sys)[5]
		f(t) = GasChromatographySimulator.holdup_time(t, T_itp, pin, pout, sys.modules[i].L, sys.modules[i].d, sys.options.gas; ng=sys.modules[i].opt.ng, vis=sys.options.vis, control=sys.options.control)
		tM_func[i] = f
	end
	return tM_func
end

"""
	holdup_time_path(sys, numpaths)

Calculates the hold-up times of the `numpaths` paths as functions of time t of the system of capillaries `sys`.

The hold-up time over edge `i => j` is calculated as
	
# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
* `numpaths`: Number of the different paths between inlet and outlets of system `sys`.
"""
function holdup_time_path(sys, num_paths)
	tM = holdup_time_functions(sys)
	paths = all_paths(sys.g, num_paths)[2]
	
	tMp = Array{Function}(undef, num_paths)
	for i=1:num_paths
		i_paths = GasChromatographySystems.index_parameter(sys.g, paths[i])
		f(t) = sum([tM[x](t) for x in i_paths])
		tMp[i] = f
	end
	return tMp
end

function holdup_time_path(sys, p2fun, num_paths)
	# collecting the hold-up time functions of every path as function of time t for system `sys` and the squared pressure solution functions `p2fun`. This function should be used, if parameters of the system are to be changes, e.g. column length or diameter, but the structure of the system is the same (same grape, same unknown pressures/flows)
	tM = holdup_time_functions(sys, p2fun)
	paths = all_paths(sys.g, num_paths)[2]
	
	tMp = Array{Function}(undef, num_paths)
	for i=1:num_paths
		i_paths = GasChromatographySystems.index_parameter(sys.g, paths[i])
		f(t) = sum([tM[x](t) for x in i_paths])
		tMp[i] = f
	end
	return tMp
end