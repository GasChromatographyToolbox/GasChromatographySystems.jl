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

function export_balance_equations_for_Mathematica(subst_bal_eq; filename="bal_eq_for_Mathematica.txt")
	# change format for Mathematica
	bal_eq_char = collect(replace(string(subst_bal_eq), "Symbolics.Equation[" => "{", "~" => "=="))
	bal_eq_char[end] = '}'
	# add the unknown pressure for which the balance equations should be solved
	unknown_p_char = collect(replace(string([P²[i] for i ∈ GasChromatographySystems.unknown_p(sys)]),  "Symbolics.Num[" => "{"))
	unknown_p_char[end] = '}'
	write(filename, "{"*join(bal_eq_char)*", "*join(unknown_p_char)*"}")
end

function import_solution_from_Mathmatica(file)
	# do the Symbols have to be defined before?
	sol_str = read(file, String)
	# change the format for Julia
	translate_sol_str = replace(sol_str, r"P²\[([0-9]+)\] -> " => "", "{" => "Symbolics.Num[", "}" => "]")
	# parse the string to convert it to symbolic expressions
	solution = eval(Meta.parse(translate_sol_str))
	return solution
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
		# convert the equations of `subst_bal_eq_λ` into matrix form a*b=0, where the rows of matrix a contain the factors for P²[i_unknown] in the order of therm_mod
		# which here is just the order of the pressure points in the system
		# solution `sol` follows the same order
		# a different order of the flow balance equations `subst_bal_eq_λ` (e.g. because balance equations for unknown outer pressures points are appended at the end of the flow balance equation array)
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
		# convert the equations of `subst_bal_eq_κ` into matrix form a*b=0, where the rows of matrix a contain the factors for P²[i_unknown] in the order of therm_mod
		# which here is just the order of the pressure points in the system
		# solution `sol` follows the same order
		# a different order of the flow balance equations `subst_bal_eq_κ` (e.g. because balance equations for unknown outer pressures points are appended at the end of the flow balance equation array)
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
		# order of the parameters defined here, if the system is changed with definition of pressures or flows, this also has to change => reevaluate this equation
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
    check_expressions_λ_κ(sol; mode="λ", n=100)

Checks the expressions of the array `sol` (solutions to the flow balance equations) if they use the flow permeabilities λ or the flow restrictions κ and substitutes them if needed. 

# Arguments
* `sol`: Symbolic expressions (of the solutions for the flow balance equations) 
* `mode`: Mode for flow equations to use flow permeabilities λ (`mode = λ`; default) or flow restrictions κ (`mode = κ`)
* `n`: Maximum of the number of expected symbols λ resp. κ. Could be replaced by the number of edges of the used system `n = ne(sys.g)`. 
"""
function check_expressions_λ_κ(sol; mode="λ", n=100)
	symbols = unique(vcat(Symbolics.get_variables.(sol)...)) # all unique symbols in the array of expressions `sol`
	@variables λ[1:n], κ[1:n] # how big should they be??? -> count them in `symbols`
	count_λ = count([any(isequal.(λ[i], symbols)) for i=1:n])
	count_κ = count([any(isequal.(κ[i], symbols)) for i=1:n])
	if count_λ == 0 && count_κ == 0
		error("Solutions do not contain symbols λ or κ.")
	elseif count_λ != 0 && count_κ != 0
		error("Solutions contain both symbols λ and κ.")
	elseif count_λ != 0 && count_κ == 0 && mode == "λ"
		# everything as it should be
		return true, sol
	elseif count_λ == 0 && count_κ != 0 && mode == "κ"
		# everything as it should be
		return true, sol
	elseif count_λ == 0 && count_κ != 0 && mode == "λ"
		# substitute κ with λ
		sub_dict = Dict([κ[i] => 1/λ[i] for i=1:count_κ])
		sub_sol = substitute(sol, sub_dict)
		return false, sub_sol
	elseif count_λ != 0 && count_κ == 0 && mode == "κ"
		# substitute λ with κ
		sub_dict = Dict([λ[i] => 1/κ[i] for i=1:count_λ])
		sub_sol = substitute(sol, sub_dict)
		return false, sub_sol
	end
end

"""
	save_build_pressure_squared_functions(sys, solution; filename=pwd()*"/p2fun_"*sys.name, mode="λ")

Constructs and saves the array of functions of the solutions `solution` for the unkown squared pressures of the flow balance equations of the system of capillaries `sys`.

# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
* `filename`: Filename, where the solution functions are saved. Default `pwd()*"/p2fun_"*sys.name` attached with `"_λ.jl"` or `"_κ.jl"`, depending on `mode`.
* `mode`: Mode for flow equations to use flow permeabilities λ (`mode = λ`; default) or flow restrictions κ (`mode = κ`)

# Output
A dictionary with the folowing keys is saved in a `.jl` file:
* `"p2fun"`: the build function of the solutions
* `"i_known_p"`: the indices of the known pressures.
* `"i_known_λ"`: the indices of the known flow permeabilities.
* `"i_known_F"`: the indices of the known flows.
* `"mode"`: mode of the functions using flow permeabilities λ or flow restrictions κ.

# Loading
The saved dictionary can easily be loaded into Julia by
```Julia
p2fun_load = include("p2fun_saved.jl")
```
The build functions have to be evaluated by `eval.(p2fun_load["p2fun"])` before usage. The arguments for the squared pressure functions are the ordered known squared pressures ``p^2``, the ordered known flow permabilities ``λ`` resp. flow restrictions ``κ``, the ordered known flows ``F`` and constant ``A = π/256 p_n/T_n``. 
"""
function save_build_pressure_squared_functions(sys, solution; filename=pwd()*"/p2fun_"*sys.name, mode="λ")
	# make here a check for the symbols used in the solutions expressions and select the mode based on this
	if mode == "λ"
		check_sol = check_expressions_λ_κ(solution; mode="λ", n=100)[2]
		return save_build_pressure_squared_functions_λ(sys, filename, check_sol)
	elseif mode == "κ"
		check_sol = check_expressions_λ_κ(solution; mode="κ", n=100)[2]
		return save_build_pressure_squared_functions_κ(sys, filename, check_sol)
	else
		error("Unknown `mode` selection. Use `mode = λ` for flow permeabilities or `mode = κ` for flow restrictions.")
	end
end

function save_build_pressure_squared_functions_λ(sys, filename, solutions)
	#build the function for the squared pressures at the unknown vertices from the symbolic expressions of the solutions of the solved balance equations
	@variables A, P²[1:nv(sys.g)], λ[1:ne(sys.g)], F[1:ne(sys.g)]

	i_unknown_p = GasChromatographySystems.unknown_p(sys)
	i_unknown_F = GasChromatographySystems.unknown_F(sys)
	i_unknown_λ = GasChromatographySystems.unknown_λ(sys)

	i_known_p = collect(1:nv(sys.g))[Not(i_unknown_p)]
	i_known_F = collect(1:ne(sys.g))[Not(i_unknown_F)]
	i_known_λ = collect(1:ne(sys.g))[Not(i_unknown_λ)]
	
	#solutions = GasChromatographySystems.solve_balance_λ(sys)
	p_solution = Array{Expr}(undef, length(i_unknown_p))
	for i=1:length(i_unknown_p)
		# order of the parameters defined here, if the system is changed with definition of pressures or flows, this also has to change => reevaluate this equation
		pfun = build_function(solutions[i], [[P²[j] for j=i_known_p]; [λ[j] for j=i_known_λ]; [F[j] for j=i_known_F]; A];
			   expression = Val{true}, # must be set to true for saving
			   target = Symbolics.JuliaTarget(),
			   parallel=nothing)
		p_solution[i] = pfun
	end
	file = filename*"_λ.jl"
	p2fun_dict = Dict("p2fun" => p_solution, "i_known_p" => i_known_p, "i_known_λ" => i_known_λ, "i_known_F" => i_known_F, "mode" => "λ")
	write(file, string(p2fun_dict))
	return p2fun_dict, file
end

function save_build_pressure_squared_functions_κ(sys, filename, solutions)
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
			   expression = Val{true}, # must be set to true for saving
			   target = Symbolics.JuliaTarget(),
			   parallel=nothing)
		p_solution[i] = pfun
	end
	file = filename*"_κ.jl"
	p2fun_dict = Dict("p2fun" => p_solution, "i_known_p" => i_known_p, "i_known_λ" => i_known_λ, "i_known_F" => i_known_F, "mode" => "κ")
	write(file, string(p2fun_dict))
	return p2fun_dict, file
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
	pressure_functions(sys, p2fun; mode="λ")

Collect all pressure functions as functions of time t at the vertices of the capillary system `sys`, either from defined input values or from the solutions of the flow balance equations. 

# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
* `p2fun`: Julia function of the solutions of the flow balance equations from `build_pressure_squared_functions(sys; mode="λ")`
* `mode`: Mode for flow equations to use flow permeabilities λ (`mode = λ`; default) or flow restrictions κ (`mode = κ`)
"""
function pressure_functions(sys, p2fun; mode="λ")
	# collect all pressure functions of all the vertices into an array of functions of time t. If the pressures are known (defined by time_steps and pressure_steps), than a linear interpolation 
	# is used, if the pressures are unkown, the symbolic solution functions `p2fun` of these pressures are used with substituted values of the system.
	# This function should be used, if parameters of the system are to be changes, e.g. column length or diameter, but the structure of the system is the same (same grape, same unknown pressures/flows)
	pres = substitute_pressure_squared_functions(p2fun, sys; mode=mode)
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
	pressure_functions(sys; mode="λ")

Collect all pressure functions as functions of time t at the vertices of the capillary system `sys`, either from defined input values or from the solutions of the flow balance equations. 

# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
* `mode`: Mode for flow equations to use flow permeabilities λ (`mode = λ`; default) or flow restrictions κ (`mode = κ`)
"""
function pressure_functions(sys; mode="λ")
	pres, unk = solve_pressure(sys; mode=mode)
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
	interpolate_pressure_functions(sys, p2fun; dt=1, mode="λ")

Interpolates (linearly) all pressure funtions at the vertices of the system of capillaries `sys` between the time steps `dt`. For the speed of the simulation these interpolated functions are faster than the pure solution functions of the flow balance equations.

# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
* `p2fun`: Julia function of the solutions of the flow balance equations from `build_pressure_squared_functions(sys; mode="λ")`
* `dt`: time steps, where the original pressure function is evaluated. Inbetween these time steps the pressure function is linearly interpolated. 
* `mode`: Mode for flow equations to use flow permeabilities λ (`mode = λ`; default) or flow restrictions κ (`mode = κ`)
"""
function interpolate_pressure_functions(sys, p2fun; dt=1, mode="λ")
	tsteps = GasChromatographySystems.common_timesteps(sys)
	tend = sum(tsteps)
	if all(degree(sys.g).<3) # no split/merge nodes, straight graph -> pressure at nodes is linear
		trange = cumsum(tsteps)
	else
		trange = 0:dt:tend
	end
	i_unknown_p = GasChromatographySystems.unknown_p(sys)
	p_func = GasChromatographySystems.pressure_functions(sys, p2fun; mode=mode) #!!! mode "λ" "κ" !!!
	p_itp = Array{Any}(undef, length(p_func))
	for i=1:length(p_func)
		#if all(isnan.(sys.pressurepoints[i].pressure_steps))
		if i ∈ i_unknown_p
			p_itp[i] = LinearInterpolation((trange, ), p_func[i].(trange), extrapolation_bc=Flat())
		else
			p_itp[i] = p_func[i]
		end # no additional interpolation, if the pressure is allready defined by a linear program
	end
	return p_itp
end

"""
	interpolate_pressure_functions(sys; dt=1, mode="λ")

Interpolates (linearly) all pressure funtions at the vertices of the system of capillaries `sys` between the time steps `dt`. For the speed of the simulation these interpolated functions are faster than the pure solution functions of the flow balance equations.

# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
* `dt`: time steps, where the original pressure function is evaluated. Inbetween these time steps the pressure function is linearly interpolated. 
* `mode`: Mode for flow equations to use flow permeabilities λ (`mode = λ`; default) or flow restrictions κ (`mode = κ`)
"""
function interpolate_pressure_functions(sys; dt=1, mode="λ")
	tsteps = GasChromatographySystems.common_timesteps(sys)
	tend = sum(tsteps)
	if all(degree(sys.g).<3) # no split/merge nodes, straight graph -> pressure at nodes is linear
		trange = cumsum(tsteps)
	else
		trange = 0:dt:tend
	end
	i_unknown_p = GasChromatographySystems.unknown_p(sys)
	p_func = GasChromatographySystems.pressure_functions(sys; mode=mode) #!!! mode "λ" "κ" !!!
	p_itp = Array{Any}(undef, length(p_func))
	for i=1:length(p_func)
		#if all(isnan.(sys.pressurepoints[i].pressure_steps))
		if i ∈ i_unknown_p
			p_itp[i] = LinearInterpolation((trange, ), p_func[i].(trange), extrapolation_bc=Flat())
		else
			p_itp[i] = p_func[i]
		end # no additional interpolation, if the pressure is allready defined by a linear program
	end
	return p_itp
end

"""
	flow_functions(sys; mode="λ")

Collects the flow functions as functions of time t for all edges of the system of capillaries `sys`.

The flow over edge `i => j` is calculated as

```math
F_{i,j} = \\frac{A}{κ_{i,j}} \\left(p_i^2-p_j^2\\right)
```

with flow restriction ``κ_{i,j} = \\int_0^{L_{i,j}} η(T_{i,j})T_{i,j}/d_{i,j}^4 dx``, pressures ``p_i`` resp. ``p_j`` at the vertices ``i`` resp. ``j``, temperature ``T_{i,j}``, capillary length ``L_{i,j}`` and diameter ``d_{i,j}`` of the edge `i => j`.

# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
* `mode`: Mode for flow equations to use flow permeabilities λ (`mode = λ`; default) or flow restrictions κ (`mode = κ`)
"""
function flow_functions(sys; mode="λ")
	p_func = pressure_functions(sys; mode=mode) #!!! mode "λ" "κ" !!!
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

# flow_functions version with p2fun
function flow_functions(sys, p2fun; mode="λ")
	p_func = pressure_functions(sys, p2fun, mode=mode)
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
	holdup_time_functions(sys; mode="λ")

Collects the hold-up time functions as functions of time t for all edges of the system of capillaries `sys`.

The hold-up time over edge `i => j` is calculated as

```math
t_{M_{i,j}} = \\frac{128}{3} η(T_{i,j}) \\frac{L_{i,j}^2}{d_{i,j}^2} \\frac{p_i^3-p_j^3}{\\left(p_i^2-p_j^2\\right)^2}
```
	
with flow restriction, pressures ``p_i`` resp. ``p_j`` at the vertices ``i`` resp. ``j``, temperature ``T_{i,j}``, capillary length ``L_{i,j}`` and diameter ``d_{i,j}`` of the edge `i => j`.
	
# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
* `mode`: Mode for flow equations to use flow permeabilities λ (`mode = λ`; default) or flow restrictions κ (`mode = κ`)
"""
function holdup_time_functions(sys; mode="λ")
	p_func = GasChromatographySystems.pressure_functions(sys; mode=mode) #!!! mode "λ" "κ" !!!
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

function holdup_time_functions(sys, p2fun; mode="λ")
	# collecting the hold-up time functions of every edge as function of time t for system `sys` and the squared pressure solution functions `p2fun`. This function should be used, if parameters of the system are to be changes, e.g. column length or diameter, but the structure of the system is the same (same grape, same unknown pressures/flows)
	p_func = pressure_functions(sys, p2fun; mode=mode) #!!! mode "λ" "κ" !!!
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
	holdup_time_path(sys, numpaths; mode="λ")

Calculates the hold-up times of the `numpaths` paths as functions of time t of the system of capillaries `sys`.

The hold-up time over edge `i => j` is calculated as
	
# Arguments
* `sys`: System structure of the capillary system for which the flow balance is set up.
* `numpaths`: Number of the different paths between inlet and outlets of system `sys`.
* `mode`: Mode for flow equations to use flow permeabilities λ (`mode = λ`; default) or flow restrictions κ (`mode = κ`)
"""
function holdup_time_path(sys, num_paths; mode="λ")
	tM = holdup_time_functions(sys; mode=mode)
	paths = all_paths(sys.g, num_paths)[2]
	
	tMp = Array{Function}(undef, num_paths)
	for i=1:num_paths
		i_paths = GasChromatographySystems.index_parameter(sys.g, paths[i])
		f(t) = sum([tM[x](t) for x in i_paths])
		tMp[i] = f
	end
	return tMp
end

function holdup_time_path(sys, p2fun, num_paths; mode="λ")
	# collecting the hold-up time functions of every path as function of time t for system `sys` and the squared pressure solution functions `p2fun`. This function should be used, if parameters of the system are to be changes, e.g. column length or diameter, but the structure of the system is the same (same grape, same unknown pressures/flows)
	tM = holdup_time_functions(sys, p2fun; mode=mode)
	paths = all_paths(sys.g, num_paths)[2]
	
	tMp = Array{Function}(undef, num_paths)
	for i=1:num_paths
		i_paths = GasChromatographySystems.index_parameter(sys.g, paths[i])
		f(t) = sum([tM[x](t) for x in i_paths])
		tMp[i] = f
	end
	return tMp
end