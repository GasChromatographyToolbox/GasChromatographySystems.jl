using Test, GasChromatographySystems

@testset "example systems" begin
    # define some example systems
    ex_series = GasChromatographySystems.example_SeriesSystem()
    @test ex_series.modules[1].length == 10.0 
    ex_split = GasChromatographySystems.example_SplitSystem()
    @test GasChromatographySystems.ne(ex_split.g) == 3
    ex_GCxGC_TM = GasChromatographySystems.example_GCxGC_TM_simp()
    @test isnan(ex_GCxGC_TM.pressurepoints[2].pressure_steps[1])
    ex_GCxGC_FM = GasChromatographySystems.example_GCxGC_FM_simp()
    @test ex_GCxGC_FM.options.mobile_phase == "He"
end

println("Test run successful.")