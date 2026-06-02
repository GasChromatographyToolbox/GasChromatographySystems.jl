using Test, CSV, DataFrames, Graphs, GasChromatographySystems

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

@testset "ValveProgram, ModuleValveOptions, ModuleValve" begin
    GCS = GasChromatographySystems

    @testset "ValveProgram" begin
        vp_manual = GCS.ValveProgram([5.0, 10.0], [true, false])
        @test vp_manual.time_steps == [5.0, 10.0]
        @test vp_manual.state_steps == [true, false]
        @test GCS.valve_state(vp_manual, 0.0)
        @test !GCS.valve_state(vp_manual, 5.0)
        @test !GCS.valve_state(vp_manual, 15.0)

        vp_def = GCS.default_ValveProgram()
        @test vp_def.time_steps == [0.0, 1800.0]
        @test vp_def.state_steps == [true, false]

        vp_per = GCS.ValveProgram(10.0, 2.0, 30.0)
        @test vp_per.time_steps == [2.0, 8.0, 2.0, 8.0, 2.0, 8.0]
        @test vp_per.state_steps == [false, true, false, true, false, true]
        @test sum(vp_per.time_steps) ≈ 30.0
        @test !GCS.valve_state(vp_per, 1.0)
        @test GCS.valve_state(vp_per, 2.0)
        @test GCS.valve_state(vp_per, 5.0)
        @test !GCS.valve_state(vp_per, 10.0)

        vp_inv = GCS.ValveProgram(10.0, 2.0, 10.0; inverted=true)
        @test GCS.valve_state(vp_inv, 1.0)
        @test !GCS.valve_state(vp_inv, 2.0)

        vp_long = GCS.default_periodic_ValveProgram()
        @test vp_long isa GCS.PeriodicValveProgram
        @test vp_long.mp == 10.0 && vp_long.t_closed == 2.0 && vp_long.t_end == 1800.0
        vp_long_exp = GCS.expand_valve_program(vp_long)
        @test sum(vp_long_exp.time_steps) ≈ 1800.0
        @test length(vp_long_exp.time_steps) == 360

        @test_throws ErrorException GCS.ValveProgram([1.0, 2.0], [true])
        @test_throws ErrorException GCS.ValveProgram(0.0, 2.0, 10.0)
        @test_throws ErrorException GCS.ValveProgram(10.0, 12.0, 10.0)
    end

    @testset "PeriodicValveProgram" begin
        pvp = GCS.PeriodicValveProgram(10.0, 2.0, 30.0)
        vp_exp = GCS.expand_valve_program(pvp)
        @test vp_exp.time_steps == [2.0, 8.0, 2.0, 8.0, 2.0, 8.0]
        for t in (0.0, 1.0, 2.0, 5.0, 10.0, 22.0, 29.0, 30.0, 35.0)
            @test GCS.valve_state(pvp, t) == GCS.valve_state(vp_exp, t)
        end
        pvp_inv = GCS.PeriodicValveProgram(10.0, 2.0, 10.0; inverted=true)
        @test GCS.valve_state(pvp_inv, 1.0)
        @test !GCS.valve_state(pvp_inv, 2.0)
        @test_throws ErrorException GCS.PeriodicValveProgram(0.0, 2.0, 10.0)
    end

    @testset "ModuleValveOptions" begin
        opt_def = GCS.ModuleValveOptions()
        @test opt_def.ng == true
        opt_ng = GCS.ModuleValveOptions(; ng=false)
        @test opt_ng.ng == false
    end

    @testset "ModuleValve" begin
        vp = GCS.ValveProgram([2.0, 8.0], [false, true])
        opt = GCS.ModuleValveOptions(; ng=true)
        T = 25.0

        v_full = GCS.ModuleValve("v1", 0.05, 1e-3, eps(), T, vp, 1.5, opt)
        @test v_full isa GCS.AbstractModule
        @test v_full.name == "v1"
        @test v_full.L == 0.05
        @test v_full.d_open == 1e-3
        @test v_full.d_closed == eps()
        @test v_full.T == T
        @test v_full.state === vp
        @test v_full.F == 1.5
        @test v_full.opt === opt

        v_nan = GCS.ModuleValve("v2", 0.05, 1e-3, eps(), T, vp, opt)
        @test isnan(v_nan.F)

        v_kw = GCS.ModuleValve("v3", 0.05, 1e-3, eps(), T, vp; ng=false)
        @test isnan(v_kw.F)
        @test v_kw.opt.ng == false

        v_short = GCS.ModuleValve("v4", T, vp; ng=true)
        @test v_short.L == 0.01
        @test v_short.d_open == 0.001
        @test v_short.d_closed == eps(Float64)
        @test isnan(v_short.F)
        @test v_short.opt.ng == true

        sys_empty = GCS.System("", Graphs.SimpleDiGraph(0), GCS.PressurePoint[], GCS.AbstractModule[], GCS.Options())
        _, temp_steps_const, _, _, _ = GCS.module_temperature(v_short, sys_empty)
        @test all(==(T), temp_steps_const)

        v_prog = GCS.ModuleValve("v5", GCS.default_TP(), vp; ng=true)
        _, temp_steps_prog, _, _, _ = GCS.module_temperature(v_prog, sys_empty)
        @test temp_steps_prog == v_prog.T.temp_steps
    end
end

@testset "Program synchronization (match_programs, update_system)" begin
    GCS = GasChromatographySystems

    """Tee with mismatched `TemperatureProgram`, `ValveProgram`, and `PressureProgram` grids."""
    function tee_mismatched_programs()
        g = SimpleDiGraph(4)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 4, 2)
        TP_col = GCS.TemperatureProgram([10.0, 20.0, 5.0], [40.0, 120.0, 200.0])
        VP = GCS.ValveProgram(10.0, 2.0, 30.0)
        TP_valve = GCS.TemperatureProgram([100.0, 50.0], [80.0, 160.0])
        pp = [
            GCS.PressurePoint("p1", 3.0e5),
            GCS.PressurePoint("p2", NaN),
            GCS.PressurePoint("p3", 1.013e5),
            GCS.PressurePoint("p4", GCS.PressureProgram([5.0, 25.0], [3.1e5, 3.2e5])),
        ]
        modules = GCS.AbstractModule[
            GCS.ModuleColumn("c12", 1.0, 0.25e-3, 0.25e-6, "Test", TP_col),
            GCS.ModuleColumn("c23", 0.5, 0.1e-3, 0.1e-6, "Test", GCS.default_TP()),
            GCS.ModuleValve("v42", 0.01, 1e-3, eps(), TP_valve, VP),
        ]
        GCS.System("tee_sync", g, pp, modules, GCS.Options())
    end

    sys = tee_mismatched_programs()
    vp_orig = sys.modules[3].state

    @test sys.modules[1].T.time_steps != vp_orig.time_steps
    @test sys.modules[2].T.time_steps != vp_orig.time_steps

    com = GCS.common_timesteps(sys)
    @test !isempty(com)
    @test length(com) > length(sys.modules[1].T.time_steps)

    com_mp, _, _, _, _, i_tempprog = GCS.match_programs(sys)
    @test com_mp == com
    @test 3 in i_tempprog
    @test GCS.index_modules_with_valve_program(sys) == [3]

    sys2 = GCS.update_system(sys)
    com2 = GCS.common_timesteps(sys2)
    @test com2 == com
    @test sys2.modules[3].state == vp_orig
    @test sys2.modules[1].T.time_steps == com2
    @test sys2.modules[2].T.time_steps == com2
    @test sys2.modules[3].T.time_steps == com2
    @test length(sys2.modules[1].T.temp_steps) == length(com2)
    @test sys2.pressurepoints[1].P == sys.pressurepoints[1].P
    @test sys2.pressurepoints[4].P.time_steps == com2

    # constant valve temperature: state program unchanged; column T synchronized
    g2 = SimpleDiGraph(3)
    add_edge!(g2, 1, 2)
    add_edge!(g2, 2, 3)
    VP2 = GCS.ValveProgram([3.0, 7.0], [false, true])
    TP2 = GCS.TemperatureProgram([50.0, 10.0], [30.0, 250.0])
    sys_c = GCS.System(
        "line",
        g2,
        [
            GCS.PressurePoint("in", 2.0e5),
            GCS.PressurePoint("mid", NaN),
            GCS.PressurePoint("out", 1.0e5),
        ],
        GCS.AbstractModule[
            GCS.ModuleColumn("col", 2.0, 0.25e-3, 0.25e-6, "Test", TP2),
            GCS.ModuleValve("valve", 0.02, 1e-3, eps(), 42.0, VP2),
        ],
        GCS.Options(),
    )
    sys_c2 = GCS.update_system(sys_c)
    com_c = GCS.common_timesteps(sys_c2)
    @test sys_c2.modules[2].T == 42.0
    @test sys_c2.modules[2].state == VP2
    @test sys_c2.modules[1].T.time_steps == com_c

    sol = GCS.solve_balance(sys2)
    @test length(sol) == 1
    p2fun = GCS.build_pressure_squared_functions(sys2, sol)
    FF = GCS.flow_functions(sys2, p2fun)
    @test isfinite(FF[3](5.0))

    VP_dense = GCS.PeriodicValveProgram(1.0 / 3, 0.1, 1800.0)
    @test VP_dense.mp ≈ 1.0 / 3
    g3 = SimpleDiGraph(4)
    add_edge!(g3, 1, 2)
    add_edge!(g3, 2, 3)
    add_edge!(g3, 4, 2)
    sys_dense = GCS.System(
        "dense_vp",
        g3,
        [
            GCS.PressurePoint("p1", 3.0e5),
            GCS.PressurePoint("p2", NaN),
            GCS.PressurePoint("p3", 1.013e5),
            GCS.PressurePoint("p4", GCS.default_PP()),
        ],
        GCS.AbstractModule[
            GCS.ModuleColumn("c12", 1.0, 0.25e-3, 0.25e-6, "Test", GCS.default_TP()),
            GCS.ModuleColumn("c23", 0.5, 0.1e-3, 0.1e-6, "Test", GCS.default_TP()),
            GCS.ModuleValve("v42", 0.01, 1e-3, eps(), 25.0, VP_dense),
        ],
        GCS.Options(),
    )
    n_vp = length(GCS.expand_valve_program(VP_dense).time_steps)
    @test length(GCS.common_timesteps(sys_dense)) < n_vp
    sys_dense2 = GCS.update_system(sys_dense)
    @test sys_dense2.modules[3].state === VP_dense
    @test sys_dense2.modules[3].state isa GCS.PeriodicValveProgram
    @test length(GCS.expand_valve_program(VP_dense).time_steps) > 100
end

@testset "Valve junction slicing" begin
    GCS = GasChromatographySystems
    VP = GCS.PeriodicValveProgram(10.0, 2.0, 100.0)
    @test GCS.valve_slicing_schedule(VP) == (mp=10.0, t_closed=2.0, phase_shift=0.0)
    @test GCS.t_start_next_open_window(1.0, 10.0, 2.0, 0.0) ≈ 0.0
    @test GCS.t_start_next_open_window(12.5, 10.0, 2.0, 0.0) ≈ 10.0
    @test !GCS.valve_state_varies(VP, 0.0, 1.0)
    @test GCS.valve_state_varies(VP, 0.0, 15.0)
    const_vp = GCS.ValveProgram([100.0], [true])
    @test GCS.valve_slicing_schedule(const_vp) === nothing
    pl = DataFrame(
        Name=["A"],
        CAS=["64-17-5"],
        tR=[12.5],
        τR=[0.5],
        Annotations=[""],
        A=[1.0],
    )
    g = SimpleDiGraph(4)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 4, 2)
    sys = GCS.System(
        "tee",
        g,
        fill(GCS.PressurePoint("p", 1.0e5), 4),
        GCS.AbstractModule[
            GCS.ModuleColumn("c12", 1.0, 0.25e-3, 0.25e-6, "Test", GCS.default_TP()),
            GCS.ModuleColumn("c23", 0.5, 0.1e-3, 0.1e-6, "Test", GCS.default_TP()),
            GCS.ModuleValve("v42", 0.01, 1e-3, eps(), 25.0, VP),
        ],
        GCS.Options(),
    )
    @test length(GCS.incident_valve_modules(sys, 2)) == 1
    @test GCS.edges_along_path_in_order(sys.g, collect(edges(g))[1:2]) == [1, 2]

    # Regression checks for junction slicing quality:
    # 1) area conservation, 2) monotonic downstream tR, 3) slice annotation/count.
    db_file = string(@__DIR__, "/data/Database_test.csv")
    db_dataframe = DataFrame(CSV.File(db_file, header=1, silencewarnings=true))
    insertcols!(db_dataframe, 1, :No => collect(1:length(db_dataframe.Name)))
    selected = [String(db_dataframe.Name[1])]
    sys_col = GCS.SeriesSystem([2.0], [0.25], [0.25], ["SLB5ms"], [GCS.default_TP()], 1.0, NaN, 0.0; name="SeriesSystem", opt=GCS.Options())
    sol_col = GCS.solve_balance(sys_col)
    p2fun_col = GCS.build_pressure_squared_functions(sys_col, sol_col)
    par_col = GCS.graph_to_parameters(sys_col, p2fun_col, db_dataframe, selected)[1]
    cas = par_col.sub[1].CAS
    pl_wide = DataFrame(
        Name=[selected[1]],
        CAS=[cas],
        tR=[12.5],
        τR=[2.0],         # wide enough to span multiple valve periods with nτ=6
        Annotations=["src_"],
        A=[1.0],
    )

    new_par, pl_out, _ = GCS.simulate_valve_junction(par_col, sys.modules[3], pl_wide; nτ=6)
    @test length(new_par.sub) == 3
    @test length(pl_out.tR) == 3
    @test all(ann -> startswith(ann, "v"), pl_out.Annotations)
    @test issorted(pl_out.tR)
    @test all(diff(pl_out.tR) .>= 0.0)
    @test isapprox(sum(pl_out.A), sum(pl_wide.A); rtol=1e-8, atol=1e-10)
end

@testset "change_initial finite-row guard" begin
    GCS = GasChromatographySystems
    db_file = string(@__DIR__, "/data/Database_test.csv")
    db_dataframe = DataFrame(CSV.File(db_file, header=1, silencewarnings=true))
    insertcols!(db_dataframe, 1, :No => collect(1:length(db_dataframe.Name)))
    selected_solutes = ["5-Nonanol"]

    sys = GCS.SeriesSystem([2.0], [0.25], [0.25], ["SLB5ms"], [GCS.default_TP()], 1.0, NaN, 0.0; name="SeriesSystem", opt=GCS.Options())
    sol = GCS.solve_balance(sys)
    p2fun = GCS.build_pressure_squared_functions(sys, sol)
    par = GCS.graph_to_parameters(sys, p2fun, db_dataframe, selected_solutes)[1]

    cas = par.sub[1].CAS
    pl_mixed = DataFrame(
        CAS=[cas, cas],
        tR=[1.23, NaN],
        τR=[0.11, NaN],
        Annotations=["ok", "bad"],
    )
    new_par = GCS.change_initial(par, pl_mixed)
    @test length(new_par.sub) == 1
    @test isfinite(new_par.sub[1].t₀)
    @test isfinite(new_par.sub[1].τ₀)

    pl_bad = DataFrame(
        CAS=[cas],
        tR=[NaN],
        τR=[NaN],
        Annotations=["bad"],
    )
    err = try
        GCS.change_initial(par, pl_bad)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("no finite peaks to pass downstream", sprint(showerror, err))
end

@testset "Chromatographic paths (exclude ModuleValve)" begin
    GCS = GasChromatographySystems
    g = SimpleDiGraph(4)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 4, 2)
    modules = GCS.AbstractModule[
        GCS.ModuleColumn("c12", 1.0, 0.25e-3, 0.25e-6, "Test", GCS.default_TP()),
        GCS.ModuleColumn("c23", 0.5, 0.1e-3, 0.1e-6, "Test", GCS.default_TP()),
        GCS.ModuleValve("v42", 0.01, 1e-3, eps(), 25.0, GCS.default_ValveProgram()),
    ]
    Eg = collect(edges(g))
    path_cols = [Eg[1], Eg[2]]
    path_with_valve = [Eg[3], Eg[2]]
    @test GCS.path_is_chromatographic(g, modules, path_cols)
    @test !GCS.path_is_chromatographic(g, modules, path_with_valve)
    _, Ep = GCS.all_paths(g, modules)
    @test !isempty(Ep)
    for ep in Ep
        @test GCS.path_is_chromatographic(g, modules, ep)
        @test all(e -> e != Eg[3], ep)
    end
end

println("Test run successful.")