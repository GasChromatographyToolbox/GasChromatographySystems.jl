# definition of specific example systems

# begin - specific systems
function SeriesSystem(Ls, ds, dfs, sps, TPs, F, pin, pout; opt=GasChromatographySystems.Options(), kwargs...)
	# add test for correct lengths of input
	# ? make two versions
	# 1. defining flow over the columns (calculate pin)
	# 2. defining inlet pressure (calculate F)
	n = length(Ls)
	g = SimpleDiGraph(n+1)
	for i=1:n
		add_edge!(g, i, i+1) 
	end
	# common time steps
	#com_timesteps = []
	#for i=1:length(TPs)
	#	if typeof(TPs[i]) <: GasChromatographySystems.TemperatureProgram
	#		com_timesteps = GasChromatographySimulator.common_time_steps(com_timesteps, TPs[i].time_steps)
	#	end
	#end
	#if isempty(com_timesteps)
	#	com_timesteps = [0.0, 36000.0]
	#end
	# pressure points
	pp = Array{GasChromatographySystems.PressurePoint}(undef, n+1)
	#pins = pin*1000.0.*ones(length(com_timesteps))
	#nans = NaN.*ones(length(com_timesteps))
	if pout == 0.0
		pouts = eps(Float64)
	else 
		pouts = pout#*1000.0.*ones(length(com_timesteps))
	end
	pp[1] = GasChromatographySystems.PressurePoint("p1", pin) # inlet
	for i=2:n
		pp[i] = GasChromatographySystems.PressurePoint("p$(i)", NaN) #
	end
	pp[end] = GasChromatographySystems.PressurePoint("p$(n+1)", pouts) # outlet
	# modules
	modules = Array{GasChromatographySystems.AbstractModule}(undef, n)
	for i=1:n
		if i==1
			modules[i] = GasChromatographySystems.ModuleColumn("$(i) -> $(i+1)", Ls[i], ds[i]*1e-3, dfs[i]*1e-6, sps[i], TPs[i], F/60e6; kwargs...)
		else
			modules[i] = GasChromatographySystems.ModuleColumn("$(i) -> $(i+1)", Ls[i], ds[i]*1e-3, dfs[i]*1e-6, sps[i], TPs[i], NaN; kwargs...)
		end
	end
	# system
	sys = GasChromatographySystems.update_system(GasChromatographySystems.System(g, pp, modules, opt))

	# add test for the defined pressures and flows
	return sys
end

function SeriesSystem(; Ls = [10.0, 5.0, 2.0, 1.0], ds = [0.53, 0.32, 0.25, 0.1], dfs = [0.53, 0.32, 0.25, 0.1], sps = ["Rxi17SilMS", "Rxi17SilMS", "Rxi17SilMS", "Rxi17SilMS"], TPs = [default_TP(), default_TP(), default_TP(), default_TP()], F = NaN, pin = default_PP(), pout = 0.0, opt=GasChromatographySystems.Options(), kwargs...)
	sys = SeriesSystem(Ls, ds, dfs, sps, TPs, F, pin, pout; opt=opt, kwargs...)
	return sys
end

#example_SeriesSystem() = SeriesSystem([10.0, 5.0, 2.0, 1.0], [0.53, 0.32, 0.25, 0.1], [0.53, 0.32, 0.25, 0.1], ["Rxi17SilMS", "Rxi17SilMS", "Rxi17SilMS", "Rxi17SilMS"], [default_TP(), default_TP(), default_TP(), default_TP()], NaN, 300.0, 0.0)

function SplitSystem(Ls, ds, dfs, sps, TPs, Fs, pin, pout1, pout2; opt=GasChromatographySystems.Options(), kwargs...)
	g = SimpleDiGraph(4)
	add_edge!(g, 1, 2) # Inj -> GC column -> Split point
	add_edge!(g, 2, 3) # Split point -> TL column -> Det 1
	add_edge!(g, 2, 4) # Split point -> TL column -> Det 2
	# common time steps
	#com_timesteps = []
	#for i=1:length(TPs)
	#	if typeof(TPs[i]) <: GasChromatographySystems.TemperatureProgram
	#		com_timesteps = GasChromatographySimulator.common_time_steps(com_timesteps, TPs[i].time_steps)
	#	end
	#end
	#if isempty(com_timesteps)
	#	com_timesteps = [0.0, 36000.0]
	#end
	# pressure points
	pp = Array{GasChromatographySystems.PressurePoint}(undef, nv(g))
	#pins = pin*1000.0.*ones(length(com_timesteps))
	#nans = NaN.*ones(length(com_timesteps))
	if pout1 == 0.0
		pout1s = eps(Float64)#.*ones(length(com_timesteps))
	else 
		pout1s = pout1#*1000.0.*ones(length(com_timesteps))
	end
	if pout2 == 0.0
		pout2s = eps(Float64)#.*ones(length(com_timesteps))
	else 
		pout2s = pout2#*1000.0.*ones(length(com_timesteps))
	end
	pp[1] = GasChromatographySystems.PressurePoint("p₁", pin) # inlet 
	pp[2] = GasChromatographySystems.PressurePoint("p₂", NaN) # 
	pp[3] = GasChromatographySystems.PressurePoint("p₃", pout1s) # outlet 1 
	pp[4] = GasChromatographySystems.PressurePoint("p₄", pout2s) # outlet 2
	# modules
	modules = Array{GasChromatographySystems.AbstractModule}(undef, ne(g))
	modules[1] = GasChromatographySystems.ModuleColumn("1 -> 2", Ls[1], ds[1]*1e-3, dfs[1]*1e-6, sps[1], TPs[1], Fs[1]/60e6; kwargs...)
	modules[2] = GasChromatographySystems.ModuleColumn("2 -> 3", Ls[2], ds[2]*1e-3, dfs[2]*1e-6, sps[2], TPs[2], Fs[2]/60e6; kwargs...)
	modules[3] = GasChromatographySystems.ModuleColumn("2 -> 4", Ls[3], ds[3]*1e-3, dfs[3]*1e-6, sps[3], TPs[3], Fs[3]/60e6; kwargs...)
	# system
	sys = GasChromatographySystems.update_system(GasChromatographySystems.System(g, pp, modules, opt))

	# add test for the defined pressures and flows
	return sys
end

function SplitSystem(; Ls = [10.0, 1.0, 5.0], ds = [0.25, 0.1, 0.25], dfs = [0.25, 0.0, 0.0], sps = ["Rxi17SilMS", "", ""], TPs = [default_TP(), 300.0, 300.0], Fs = [1.0, NaN, NaN], pin = NaN, pout1 = 0.0, pout2 = 101300.0, opt=GasChromatographySystems.Options(), kwargs...)
	sys = SplitSystem(Ls, ds, dfs, sps, TPs, Fs, pin, pout1, pout2; opt=opt, kwargs...)
	return sys
end

# definition GCxGC system with thermal modulator
function GCxGC_TM(L1, d1, df1, sp1, TP1, L2, d2, df2, sp2, TP2, LTL, dTL, dfTL, spTL, TPTL, LM::Array{Float64,1}, dM, dfM, spM, shift, PM, ratioM, HotM, ColdM, TPM, F, pin, pout; opt=GasChromatographySystems.Options(), optTM=ModuleTMOptions(), optCol=ModuleColumnOptions())

	TPs = [TP1, TP2, TPM]
	
	g = SimpleDiGraph(9)
	add_edge!(g, 1, 2) # 1st-D GC
	add_edge!(g, 2, 3) # modulator
	add_edge!(g, 3, 4) # hot/cold 1
	add_edge!(g, 4, 5) # modulator
	add_edge!(g, 5, 6) # hot/cold 2 
	add_edge!(g, 6, 7) # modulator
	add_edge!(g, 7, 8) # 2nd-D GC
	add_edge!(g, 8, 9) # TL
	
	# common time steps
	#com_timesteps = []
	#for i=1:length(TPs)
	#	if typeof(TPs[i]) <: TemperatureProgram
	#		com_timesteps = GasChromatographySimulator.common_time_steps(com_timesteps, TPs[i].time_steps)
	#	end
	#end
	#if isempty(com_timesteps)
	#	com_timesteps = [0.0, 36000.0]
	#end
	
	# pressure points
	#if length(pin) == 1
	#	pins = pin*1000.0.*ones(length(com_timesteps))
	#else
	#end
	#nans = NaN.*ones(length(com_timesteps))
	if pout == 0.0
		pouts = eps(Float64)#.*ones(length(com_timesteps))
	else 
		pouts = pout#*1000.0.*ones(length(com_timesteps))
	end
	pp = Array{PressurePoint}(undef, nv(g))
	pp[1] = PressurePoint("p1", pin) # inlet 
	for i=2:(nv(g)-1)
		pp[i] = PressurePoint("p$(i)", NaN)
	end
	pp[end] = PressurePoint("p$(nv(g))", pouts) # outlet
	# modules
	modules = Array{AbstractModule}(undef, ne(g))
	modules[1] = ModuleColumn("GC column 1", L1, d1*1e-3, df1*1e-6, sp1, TP1, F/60e6, optCol)
	modules[2] = ModuleColumn("mod in", LM[1], dM*1e-3, dfM*1e-6, spM, TPM, optCol)
	modules[3] = ModuleTM("TM1", LM[2], dM*1e-3, dfM*1e-6, spM, TPM, shift, PM, ratioM, HotM, ColdM, NaN, optTM)
	modules[4] = ModuleColumn("mod loop", LM[3], dM*1e-3, dfM*1e-6, spM, TPM, optCol)
	modules[5] = ModuleTM("TM2", LM[4], dM*1e-3, dfM*1e-6, spM, TPM, shift, PM, ratioM, HotM, ColdM, NaN, optTM)
	modules[6] = ModuleColumn("mod out", LM[5], dM*1e-3, dfM*1e-6, spM, TPM, NaN, optCol)
	modules[7] = ModuleColumn("GC column 2", L2, d2*1e-3, df2*1e-6, sp2, TP2, NaN, optCol)
	modules[8] = ModuleColumn("TL", LTL, dTL*1e-3, dfTL*1e-6, spTL, TPTL, NaN, optCol)
	# system
	sys = update_system(System(g, pp, modules, opt))
	return sys
end

function GCxGC_TM(; L1 = 30.0, d1 = 0.25, df1 = 0.25, sp1 = "ZB1ms", TP1 = default_TP(), L2 = 2.0, d2 = 0.1, df2 = 0.1, sp2 = "Stabilwax", TP2 = default_TP(), LTL = 0.25, dTL = 0.1, dfTL = 0.1, spTL = "Stabilwax", TPTL = 280.0, LM = [0.30, 0.01, 0.90, 0.01, 0.30], dM = 0.1, dfM = 0.1, spM = "Stabilwax", shift = 0.0, PM = 4.0, ratioM = 0.9125, HotM = 30.0, ColdM = -120.0, TPM = default_TP(), F = 0.8, pin = NaN, pout = 0.0, opt=GasChromatographySystems.Options(), optTM=ModuleTMOptions(), optCol=ModuleColumnOptions())
	sys = GCxGC_TM(L1, d1, df1, sp1, TP1, L2, d2, df2, sp2, TP2, LTL, dTL, dfTL, spTL, TPTL, LM::Array{Float64,1}, dM, dfM, spM, shift, PM, ratioM, HotM, ColdM, TPM, F, pin, pout; opt=opt, optTM=optTM, optCol=optCol)
	return sys
end
# end - specific systems
