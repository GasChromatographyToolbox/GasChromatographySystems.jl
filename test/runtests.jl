using Test, CSV, DataFrames, GasChromatographySystems

@testset "example systems" begin
    # define some example systems
    ex_series = GasChromatographySystems.SeriesSystem(sps = ["SLB5ms", "SPB50", "Wax", "Wax"])
    ex_series_ = GasChromatographySystems.SeriesSystem(sps = ["SLB5ms", "SPB50", "Wax", "Wax"]; abstol=1e-9)
    @test ex_series_.modules[1].opt.abstol*10.0 == ex_series.modules[1].opt.abstol
    ex_split = GasChromatographySystems.SplitSystem(sps = ["SLB5ms", "SPB50", "Wax"])
    @test GasChromatographySystems.ne(ex_split.g) == 3
    ex_GCxGC_TM = GasChromatographySystems.GCxGC_TM(sp1 = "SLB5ms", sp2 = "Wax", spTL = "Wax", spM = "Wax")
    @test ex_GCxGC_TM.options.gas == "He"
    # run simulations on these systems
    # data for soultes
    db_file = string(@__DIR__, "/data/Database_test.csv")
    db_dataframe = DataFrame(CSV.File(db_file, header=1, silencewarnings=true))
	insertcols!(db_dataframe, 1, :No => collect(1:length(db_dataframe.Name)))
    selected_solutes = ["5-Nonanol", "Undecane", "2-Nonanol"]
    # graph to parameters
    sol_ex_series = GasChromatographySystems.solve_balance(ex_series)
    p2fun_series = GasChromatographySystems.build_pressure_squared_functions(ex_series, sol_ex_series)
    par_series = GasChromatographySystems.graph_to_parameters(ex_series, p2fun_series, db_dataframe, selected_solutes)
    @test par_series[1].col.sp == ex_series.modules[1].sp
    # regression: NaN junction pressures and default injection times (GCSim 0.6)
    @test all(isfinite, par_series[1].prog.Fpin_steps)
    @test all(isfinite, par_series[1].prog.pout_steps)
    @test all(s -> iszero(s.t₀) && iszero(s.τ₀), par_series[1].sub)

    sol_ex_split = GasChromatographySystems.solve_balance(ex_split)
    p2fun_split = GasChromatographySystems.build_pressure_squared_functions(ex_split, sol_ex_split)
    par_split = GasChromatographySystems.graph_to_parameters(ex_split, p2fun_split, db_dataframe, selected_solutes)
    @test par_split[2].col.sp == ex_split.modules[2].sp 

    sol_ex_GCxGC_TM = GasChromatographySystems.solve_balance(ex_GCxGC_TM)
    p2fun_GCxGC_TM = GasChromatographySystems.build_pressure_squared_functions(ex_GCxGC_TM, sol_ex_GCxGC_TM)
    par_GCxGC_TM = GasChromatographySystems.graph_to_parameters(ex_GCxGC_TM, p2fun_GCxGC_TM, db_dataframe, selected_solutes)
    @test par_GCxGC_TM[4].col.sp == ex_GCxGC_TM.modules[4].sp 
end

@testset "system with temperature gradient" begin
    # temperature program without gradient
    TP = GasChromatographySystems.TemperatureProgram(GasChromatographySystems.GasChromatographySimulator.conventional_program([40.0, 1.0, 5.0, 200.0, 2.0, 15.0, 300.0, 3.0])...)
    # temperature program with gradient
    ΔT = [0.0, 40.0, 80.0]
	x0 = [0.0, 0.0, 0.0]
	L0 = [2.0, 2.0, 2.0]
	alpha = [0.0, 3.0, 6.0]
	a_gf = [ΔT x0 L0 alpha]
	gf(x) = GasChromatographySystems.GasChromatographySimulator.gradient(x, a_gf)
	TP_grad = GasChromatographySystems.TemperatureProgram([0.0, 2000.0, 3000.0], [40.0, 160.0, 260.0], gf, a_gf)
    series_grad = GasChromatographySystems.SeriesSystem([10.0, 2.0], [0.25, 0.25], [0.25, 0.25], ["SLB5ms", "SLB5ms"], [TP, TP_grad], 1.0, NaN, 0.0; name="SeriesSystem", opt=GasChromatographySystems.Options())
    
    @test series_grad.modules[1].T.time_steps == series_grad.modules[2].T.time_steps
    @test series_grad.modules[1].opt.ng == true
    @test series_grad.modules[2].opt.ng == false
    @test series_grad.modules[2].T.gf(2.0)[end] == -ΔT[end]
end

@testset "temporal ModuleColumn (reverse-capable)" begin
    db_file = string(@__DIR__, "/data/Database_test.csv")
    db_dataframe = DataFrame(CSV.File(db_file, header=1, silencewarnings=true))
    insertcols!(db_dataframe, 1, :No => collect(1:length(db_dataframe.Name)))
    selected = ["Undecane"]

    # Regression: default domain stays :spatial and matches today's tR path
    sys_pos = GasChromatographySystems.SeriesSystem(
        [10.0], [0.25], [0.25], ["SLB5ms"],
        [GasChromatographySystems.TemperatureProgram([0.0, 600.0], [40.0, 200.0])],
        NaN, 150000.0, 101300.0;
        name="pos_reg", opt=GasChromatographySystems.Options(control="Pressure"),
    )
    @test sys_pos.modules[1].opt.domain === :spatial
    sol_bal = GasChromatographySystems.solve_balance(sys_pos)
    p2 = GasChromatographySystems.build_pressure_squared_functions(sys_pos, sol_bal)
    par = GasChromatographySystems.graph_to_parameters(sys_pos, p2, db_dataframe, selected)
    @test par[1].opt.domain === :spatial
    paths = GasChromatographySystems.all_paths(sys_pos.g, 1)[2]
    _, pls_pos, _, _ = GasChromatographySystems.simulate_along_paths(sys_pos, p2, paths, par)
    tR_pos = pls_pos[1][1].tR[1]
    @test isfinite(tR_pos)

    # Same graph with domain=:temporal should give a comparable forward elution time
    sys_td = GasChromatographySystems.SeriesSystem(
        [10.0], [0.25], [0.25], ["SLB5ms"],
        [GasChromatographySystems.TemperatureProgram([0.0, 600.0], [40.0, 200.0])],
        NaN, 150000.0, 101300.0;
        name="temporal_fwd", opt=GasChromatographySystems.Options(control="Pressure"),
        domain=:temporal,
    )
    @test sys_td.modules[1].opt.domain === :temporal
    sol_bal_t = GasChromatographySystems.solve_balance(sys_td)
    p2_t = GasChromatographySystems.build_pressure_squared_functions(sys_td, sol_bal_t)
    par_t = GasChromatographySystems.graph_to_parameters(sys_td, p2_t, db_dataframe, selected)
    @test par_t[1].opt.domain === :temporal
    _, pls_t, _, _ = GasChromatographySystems.simulate_along_paths(sys_td, p2_t, paths, par_t)
    @test isfinite(pls_t[1][1].tR[1])
    @test abs(pls_t[1][1].tR[1] - tR_pos) / tR_pos < 0.05
    @test hasproperty(pls_t[1][1], :exit_end)
    @test pls_t[1][1].exit_end[1] ≈ sys_td.modules[1].L atol=1e-6

    # Cold forward load then hot reverse flush on one short column (two phase programs)
    Ltrap = 0.5
    opt_td = GasChromatographySystems.ModuleColumnOptions(domain=:temporal)
    g = SimpleDiGraph(2)
    add_edge!(g, 1, 2)
    # Load: cold, pin > pout
    TP_load = GasChromatographySystems.TemperatureProgram([0.0, 30.0], [40.0, 40.0])
    pp_load = [
        GasChromatographySystems.PressurePoint("pin", 150000.0),
        GasChromatographySystems.PressurePoint("pout", 101300.0),
    ]
    mod_load = GasChromatographySystems.ModuleColumn("trap", Ltrap, 0.25e-3, 0.25e-6, "SLB5ms", TP_load, opt_td)
    sys_load = GasChromatographySystems.update_system(GasChromatographySystems.System("trap_load", g, pp_load, [mod_load], GasChromatographySystems.Options(control="Pressure")))
    sol_load = GasChromatographySystems.solve_balance(sys_load)
    p2_load = GasChromatographySystems.build_pressure_squared_functions(sys_load, sol_load)
    par_load = GasChromatographySystems.graph_to_parameters(sys_load, p2_load, db_dataframe, selected)[1]
    _, pl_load, sol_load_td = GasChromatographySystems.simulate_ModuleColumn(par_load, [0.0], [0.0])
    # Band should still be inside after a short cold load
    @test isnan(pl_load.tR[1]) || pl_load.exit_end[1] ≈ Ltrap
    x_park = sol_load_td[1].u[end][1]
    σ_park = sqrt(max(sol_load_td[1].u[end][2], 0.0))
    @test 0.0 < x_park < Ltrap || isapprox(x_park, Ltrap; atol=1e-8)

    # Reverse: hot, pin < pout; continue from parked (x, σ)
    TP_rel = GasChromatographySystems.TemperatureProgram([0.0, 120.0], [250.0, 250.0])
    pp_rel = [
        GasChromatographySystems.PressurePoint("pin", 101300.0),
        GasChromatographySystems.PressurePoint("pout", 150000.0),
    ]
    mod_rel = GasChromatographySystems.ModuleColumn("trap", Ltrap, 0.25e-3, 0.25e-6, "SLB5ms", TP_rel, opt_td)
    sys_rel = GasChromatographySystems.update_system(GasChromatographySystems.System("trap_rel", g, pp_rel, [mod_rel], GasChromatographySystems.Options(control="Pressure")))
    @test GasChromatographySystems.path_possible(sys_rel, GasChromatographySystems.build_pressure_squared_functions(sys_rel, GasChromatographySystems.solve_balance(sys_rel)), collect(edges(g)))
    sol_rel_bal = GasChromatographySystems.solve_balance(sys_rel)
    p2_rel = GasChromatographySystems.build_pressure_squared_functions(sys_rel, sol_rel_bal)
    par_rel0 = GasChromatographySystems.graph_to_parameters(sys_rel, p2_rel, db_dataframe, selected)[1]
    par_rel = GasChromatographySystems.change_initial_spatial(par_rel0, [0.0], [x_park], [σ_park])
    pl_rel, sol_rel = GasChromatographySystems.GasChromatographySimulator.simulate(par_rel)
    @test isfinite(pl_rel.tR[1])
    @test pl_rel.exit_end[1] ≈ 0.0 atol=1e-6
    @test sol_rel[1].u[end][1] ≈ 0.0 atol=1e-5
end

    println("Test run successful.")
