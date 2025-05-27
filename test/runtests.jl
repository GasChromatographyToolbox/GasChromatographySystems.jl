using Test, CSV, DataFrames, GasChromatographySystems

@testset "example systems" begin
    # define some example systems
    ex_series = GasChromatographySystems.SeriesSystem(sps = ["SLB5ms", "SPB50", "Wax", "Wax"])
    ex_series_ = GasChromatographySystems.SeriesSystem(sps = ["SLB5ms", "SPB50", "Wax", "Wax"]; abstol=1e-9)
    @test ex_series_.modules[1].opt.abstol*10.0 == ex_series.modules[1].opt.abstol
    ex_split = GasChromatographySystems.SplitSystem(sps = ["SLB5ms", "SPB50", "Wax"])
    @test GasChromatographySystems.ne(ex_split.g) == 3
    #ex_GCxGC_TM_simp = GasChromatographySystems.GCxGC_TM_simp(sp1 = "SLB5ms", sp2 = "Wax")
    #@test isnan(ex_GCxGC_TM_simp.pressurepoints[2].pressure_steps[1])
    #ex_GCxGC_FM_simp = GasChromatographySystems.GCxGC_FM_simp(sp1 = "SLB5ms", sp2 = "Wax")
    #@test ex_GCxGC_FM_simp.options.mobile_phase == "He"
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
    sol_ex_split = GasChromatographySystems.solve_balance(ex_split)
    p2fun_split = GasChromatographySystems.build_pressure_squared_functions(ex_split, sol_ex_split)
    par_split = GasChromatographySystems.graph_to_parameters(ex_split, p2fun_split, db_dataframe, selected_solutes)
    @test par_split[2].col.sp == ex_split.modules[2].sp 
    #par_GCxGC_TM_simp = GasChromatographySystems.graph_to_parameters(ex_GCxGC_TM_simp, db_dataframe, selected_solutes)
    #@test par_GCxGC_TM_simp[1].col.sp == ex_GCxGC_TM_simp.modules[1].stationary_phase 
    #par_GCxGC_FM_simp = GasChromatographySystems.graph_to_parameters(ex_GCxGC_FM_simp, db_dataframe, selected_solutes) 
    #@test par_GCxGC_FM_simp[1].col.sp == ex_GCxGC_FM_simp.modules[1].stationary_phase
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

    println("Test run successful.")