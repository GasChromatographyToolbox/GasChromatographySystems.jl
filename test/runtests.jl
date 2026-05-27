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
        @test sum(vp_long.time_steps) ≈ 1800.0
        @test length(vp_long.time_steps) == 360
        @test vp_long.time_steps[1:2] == [2.0, 8.0]
        @test vp_long.state_steps[1:2] == [false, true]

        @test_throws ErrorException GCS.ValveProgram([1.0, 2.0], [true])
        @test_throws ErrorException GCS.ValveProgram(0.0, 2.0, 10.0)
        @test_throws ErrorException GCS.ValveProgram(10.0, 12.0, 10.0)
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
    end
end

println("Test run successful.")