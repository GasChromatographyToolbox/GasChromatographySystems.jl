# functions for thermal modulator

# this splicing function is with separating focussed and unfocussed peak segments
# area as an input is needed 'AR'
# A_focussed calculated in relation to it
# A_focussed is the complete area during a modulation period
# not focussed segment is already included in the focussed segment, assuming it will be focussed in the 2nd modulation in a multi stage modulation
function slicing(pl, PM, ratio, shift, par::GasChromatographySimulator.Parameters; nτ=6, τ₀=zeros(length(pl.τR)), abstol=1e-8, reltol=1e-8, alg=OwrenZen5())
	tR = pl.tR
	τR = pl.τR
	AR = pl.A
	tcold = PM*ratio
	thot = PM*(1-ratio)
	init_t_start = (fld.(tR .+ shift .- nτ.*τR, PM)).*PM .- shift # start time of the peaks, rounded down to multiple of PM
	init_t_end = (fld.(tR .+ shift .+ nτ.*τR, PM)).*PM .- shift # end time of the peaks, rounded down to multiple of PM (rounding up leads to an additional slice)
 	n_slice = round.(Int, (init_t_end .- init_t_start)./PM .+ 1) # number of slices for every substance
	println("slicing(): init_t_start=$(init_t_start)s, init_t_end=$(init_t_end)s, n_slice=$(n_slice).")
	sub_TM_focussed = Array{GasChromatographySimulator.Substance}(undef, sum(n_slice))
	A_focussed = Array{Float64}(undef, sum(n_slice))
	#A_unfocussed = Array{Float64}(undef, sum(n_slice))
	g(x,p) = 1/sqrt(2*π*p[2]^2)*exp(-(x-p[1])^2/(2*p[2]^2))
	ii = 1
	Name = Array{String}(undef, sum(n_slice))
	CAS = Array{String}(undef, sum(n_slice))
	Ann_focussed = Array{String}(undef, sum(n_slice))
	t0_foc = Array{Float64}(undef, sum(n_slice))
	for i=1:length(n_slice)
		for j=1:n_slice[i]
			if n_slice[i] == 1 # no slicing, peak fits completly inside the modulation periode
				t₀ = tR[i]
                A_focussed[ii] = AR[i]
			else
				t₀ = init_t_start[i]+(j-1)*PM # initial start time
                # Integrals:
                p = [tR[i], τR[i]]
                # approximated integrals
                #prob_focussed = IntegralProblem(g, init_t_start[i]+(j-1)*PM, init_t_start[i]+(j-1)*PM+tcold, p)
                #prob_unfocussed = IntegralProblem(g, init_t_start[i]+(j-1)*PM+tcold, init_t_start[i]+(j-1)*PM+tcold+thot, p)
                # all focussed:
                domain = (init_t_start[i]+(j-1)*PM, init_t_start[i]+j*PM)
                    println("domain: $(domain), i=$(i), j=$(j), type: $(typeof(domain))")
                prob_focussed = IntegralProblem(g, domain, p)
                    #println("prob_focussed: $(prob_focussed)")
                A_focussed[ii] = solve(prob_focussed, QuadGKJL(); reltol = reltol, abstol = abstol).u * AR[i]
                #A_unfocussed[ii] = solve(prob_unfocussed, QuadGKJL(); reltol = 1e-18, abstol = 1e-30).u * AR[i]
			end
            println("A_focussed[ii] = $(A_focussed[ii])")
			CAS_par = [par.sub[i].CAS for i in 1:length(par.sub)]
			i_sub = findfirst(pl.CAS[i] .== CAS_par)
			sub_TM_focussed[ii] = GasChromatographySimulator.Substance(par.sub[i_sub].name, par.sub[i_sub].CAS, par.sub[i_sub].Tchar, par.sub[i_sub].θchar, par.sub[i_sub].ΔCp, par.sub[i_sub].φ₀, "s$(j)_"*pl.Annotations[i], par.sub[i_sub].Cag, t₀, τ₀[i])
			# Areas in the same order as sub_TM_focussed
			Name[ii] = sub_TM_focussed[ii].name
			CAS[ii] = sub_TM_focussed[ii].CAS
			Ann_focussed[ii] = sub_TM_focussed[ii].ann
			t0_foc[ii] = t₀
			ii = ii + 1
		end
	end
	newopt = GasChromatographySimulator.Options(alg, abstol, reltol, par.opt.Tcontrol, par.opt.odesys, par.opt.ng, par.opt.vis, par.opt.control, par.opt.k_th)
	newpar_focussed = GasChromatographySimulator.Parameters(par.col, par.prog, sub_TM_focussed, newopt)

	#df_A_foc = DataFrame(Name=Name, CAS=CAS, Annotations=Ann_focussed, A=A_focussed+A_unfocussed, t0=t0_foc)
	df_A_foc = DataFrame(Name=Name, CAS=CAS, Annotations=Ann_focussed, A=A_focussed, t0=t0_foc)
	# t0 is the start time of the slice
	# A the area of the sliced peak (including the section during hot-jet [source of error])
	
	return newpar_focussed, df_A_foc
end

function simplifiedTM(T, par, df_A, PM, ratio, shift, Thot;)
	# rectangular function is assumed -> does this work with a changed definition (start modulation with the hot jet)
	tcold = PM*ratio
	thot = PM*(1-ratio)
	sort_df_A = sort(df_A, :t0)
	#println("simplifiedTM(): sort_df_A.t₀ = $(sort_df_A.t0)s.")
	t₀ = sort_df_A.t0 # no need to again calculate the time as it was already calculated with the slicing in df_A.t0
	A = sort_df_A.A
	No = [parse(Int, split(sort_df_A.Annotations[x], " ")[end]) for x in 1:length(sort_df_A.Annotations)]
	Name = sort_df_A.Name
	CAS = sort_df_A.CAS
    number_of_modulation = cld.(t₀ .+ shift, PM)
	tR = (number_of_modulation + 1) .* PM .- shift .- thot # the next hot jet is set as tR
	#println("simplifiedTM(): t₀ = $(t₀)s, tR = $(tR)s.")

	# using par.prog.T_itp can result in wrong temperatures at tR, because of rounding errors for Float64 in `mod()`-function inside the `therm_mod()`-function.
	# take the temperature/temperature program defined for the module and add Thot 
	T_itp = if typeof(T) <: Number
		gf(x) = zero(x).*ones(2)
		GasChromatographySimulator.temperature_interpolation([0.0, 3600.0], [T + Thot, T + Thot], gf, par.col.L)
	else
		GasChromatographySimulator.temperature_interpolation(T.time_steps, T.temp_steps .+ Thot, T.gf, par.col.L)
	end
	TR = T_itp(par.col.L, tR) .- 273.15
	kR = Array{Float64}(undef, length(tR))
	uR = Array{Float64}(undef, length(tR))
	τ₀ = Array{Float64}(undef, length(tR))
	for i=1:length(tR)
		i_sub = findfirst(Name[i] .== df_A.Name)
		kR[i] = GasChromatographySimulator.retention_factor(par.col.L, tR[i], T_itp, par.col.d, par.col.df, par.sub[i_sub].Tchar, par.sub[i_sub].θchar, par.sub[i_sub].ΔCp, par.sub[i_sub].φ₀)
		rM = GasChromatographySimulator.mobile_phase_residency(par.col.L, tR[i], T_itp, par.prog.Fpin_itp, par.prog.pout_itp, par.col.L, par.col.d, par.col.gas; ng=par.opt.ng, vis=par.opt.vis, control=par.opt.control)
		uR[i] = 1/(rM*(1+kR[i]))
		τ₀[i] = par.sub[i_sub].τ₀
	end
	τR = par.col.L./uR
	σR = uR./τR
	Annotations = sort_df_A.Annotations
	pl = sort!(DataFrame(No=No, Name=Name, CAS=CAS, tR=tR, τR=τR, TR=TR, σR=σR, uR=uR, kR=kR, Annotations=Annotations, A=A), [:tR])
	Res = Array{Float64}(undef, length(tR))
	Δs = Array{Float64}(undef, length(tR))
	for i=1:(length(tR)-1)
		Res[i] = (pl.tR[i+1] - pl.tR[i])/(2*(pl.τR[i+1] + pl.τR[i]))
        Δs[i] = (pl.tR[i+1] - pl.tR[i])/(pl.τR[i+1] - pl.τR[i]) * log(pl.τR[i+1]/pl.τR[i])
	end
	pl[!, :Res] = Res
	pl[!, :Δs] = Δs

	sol = Array{NamedTuple{(:t, :u), Tuple{Vector{Float64}, Vector{Tuple{Float64, Float64}}}}}(undef, length(tR))
	for i=1:length(tR)
		sol[i] = (t = [0.0, par.col.L], u = [(t₀[i], τ₀[i]), (tR[i], τR[i]^2)])
	end
	
	return select(pl, [:No, :Name, :CAS, :tR, :τR, :TR, :σR, :uR, :kR, :Res, :Δs, :Annotations, :A]), sol
end

# thermal modulator segment 
function simulate_ModuleTM(segment_par, segment_module, prev_peaklist; nτ=6, τ₀_focus=zeros(length(segment_par.sub)), refocus=false, kwargsTM...) # segment_par = par_sys[i_par[j]] , segment_module = sys.modules[i_par[j]], prev_peaklist = peaklists_[j-1]
	if refocus == true
		τ₀=τ₀_focus
	else
		τ₀=prev_peaklist.τR # PM?
	end
	new_segment_par, df_A = slicing(prev_peaklist, segment_module.PM, segment_module.ratio, segment_module.shift, segment_par; nτ=nτ, τ₀=τ₀, abstol=segment_module.opt.abstol, reltol=segment_module.opt.reltol, alg=segment_module.opt.alg)
	if segment_module.opt.alg == "simplifiedTM"
		peaklist, solutions = simplifiedTM(segment_module.T, new_segment_par, df_A, segment_module.PM, segment_module.ratio, segment_module.shift, segment_module.Thot)
	else # not tested
		sol = Array{Any}(undef, length(new_segment_par.sub))
		for i_sub=1:length(new_segment_par.sub)
			dt = segment_module.opt.dtinit
			sol[i_sub] = GasChromatographySimulator.solving_odesystem_r(new_segment_par.col, new_segment_par.prog, new_segment_par.sub[i_sub], new_segment_par.opt; kwargsTM..., dt=dt)
			tR = sol[i_sub].u[end][1]
			while (fld(tR + segment_module.shift, segment_module.PM) > fld(new_segment_par.sub[i_sub].t₀ + segment_module.shift, segment_module.PM) && dt > eps()) || (sol[i_sub].retcode != ReturnCode.Success && dt > eps()) # solute elutes not in the same modulation periode or the solving failed
				dt = dt/10 # reduce initial step-width
				dtmax = dt*1000
				@warn "Retention time $(tR) surpasses modulation period, t₀ = $(new_segment_par.sub[i_sub].t₀), PM = $(segment_module.PM). Initial step-width dt is decreased ($(dt))."
				sol[i_sub] = GasChromatographySimulator.solving_odesystem_r(new_segment_par.col, new_segment_par.prog, new_segment_par.sub[i_sub], new_segment_par.opt; kwargsTM..., dt=dt, dtmax=dtmax)
				tR = sol[i_sub].u[end][1]
			end
								#if isnan(tR) # does not work
								#	@warn "Simulation result is NaN. Approximate modulator."
								#	sol[i_sub] = approximate_modulator(new_par_sys[i_par[j]], df_A, sys.modules[i_par[j]].PM, sys.modules[i_par[j]].ratio, sys.modules[i_par[j]].shift)[2][i_sub]
								#end
		end
		peaklist = GasChromatographySimulator.peaklist(sol, new_segment_par)
		GasChromatographySystems.add_A_to_pl!(peaklist, df_A)
		solutions = sol
	end
	if maximum(peaklist.τR) > segment_module.PM
		return @warn "Peak width of focussed peaks $(maximum(peaklist.τR)) > modulation period $(segment_module.PM). Simulation is aborted. alg=$(segment_module.opt.alg), T=$(segment_module.T), peaklist=$(peaklist)."
	end
	return new_segment_par, peaklist, solutions
end