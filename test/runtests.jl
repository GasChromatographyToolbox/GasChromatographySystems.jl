using Test, CSV, DataFrames, GasChromatographySystems

#@testset "example systems" begin
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
#end

println("Test run successful.")