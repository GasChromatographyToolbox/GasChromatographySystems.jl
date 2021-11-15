# test the functionality of GasChromatographySystems
using Pkg
Pkg.activate("/Users/janleppert/Documents/GitHub/GasChromatographySystems")
using GasChromatographySystems

# Definition of a simple GC-system

# Options:
Option = GasChromatographySystems.Options("He", OwrenZen5(), 10.0^-6, 10.0.^-3, "inlet", true)

# Modules:
TL1 = GasChromatographySystems.Transferline(0.5, 0.25e-3, 0.25e-6, "Wax", 300.0)

a_gf = [[20.0, 20.0, 40.0, 40.0] [0.0, 0.0, 0.0, 0.0] [10.0, 10.0, 10.0, 10.0] [-8, -8, -6, -6]]
gf(x) = GasChromatographySimulator.gradient(x, a_gf; Tcontrol=Option.Tcontrol)

TP = GasChromatographySystems.Temperature_Program([0.0, 60.0, 300.0, 60.0], 
                                                    [40.0, 40.0, 280.0, 280.0],
                                                    gf,
                                                    a_gf)
a_d = [0.25e-3]
d(x) = GasChromatographySimulator.gradient(x, a_d)
a_df = [0.25e-6]
df(x) = GasChromatographySimulator.gradient(x, a_df)
Column = GasChromatographySystems.Column(10.0, d, a_d, df, a_df, "Wax", TP)

# Pressures:
PP1 = GasChromatographySystems.Pressure_Point([0.0, 60.0, 150.0, 150.0, 60.0],
                                                [200.0, 200.0, 250.0, 280.0, 280.0].*1000.0.+101300)
PP2 = GasChromatographySystems.Pressure_Point([0.0, 500.0],
                                                [101300.0, 101300.0])                               
# Combine for GC-system:
GCsys1 = PP1, TL1, Column, PP2
GCsys2 = PP1, Column, PP2
GCsys3 = TL1, PP1, Column, PP2

# test the GC-system:
GasChromatographySystems.test_of_GCsys(GCsys1) # should be 'true'
GasChromatographySystems.test_of_GCsys(GCsys2) # should be 'true'
GasChromatographySystems.test_of_GCsys(GCsys3) # should error

GasChromatographySystems.modules_with_timesteps(GCsys1) == [1, 3, 4]
GasChromatographySystems.modules_with_diameter(GCsys1) == [2, 3]
GasChromatographySystems.modules_with_film_thickness(GCsys1) == [2, 3]
GasChromatographySystems.modules_index(GCsys1) == [2, 3]
GasChromatographySystems.pressure_points_index(GCsys1) == [1, 4]

# definition of solutes
solutes = ["C10", "C11", "C12"]

# load the Database
db = DataFrame(CSV.File("Database_test.csv"))

# translate the GC-system into GasChromatographySimulator.Parameters
par1 = GasChromatographySystems.initilize_parameters(GCsys1, Option, solutes, db)
# run the Simulation
par, sol = GasChromatographySystems.linear_GC_system_simulation(GCsys1, Option, solutes, db)

# compare results for GCsys2 and direct Simulation with
# GasChromatographySimulator:
par2, sol2 = GasChromatographySystems.linear_GC_system_simulation(GCsys2, Option, solutes, db)
tR2_C10 = sol2[1][1].u[end][1]

opt = GasChromatographySimulator.Options(Option.alg, Option.abstol, Option.reltol, Option.Tcontrol, Option.odesys)
sys = GasChromatographySimulator.System(Column.length, Column.diameter, Column.a_diameter, Column.film_thickness, Column.a_film_thickness, Column.stationary_phase, Option.mobile_phase)
prog = GasChromatographySimulator.constructor_Program(GasChromatographySystems.new_time_steps(GCsys2), 
                                                        GasChromatographySystems.new_temperature_steps(GCsys2, Option)[:,1], 
                                                        GasChromatographySystems.new_pressure_steps(GCsys2)[:,1], 
                                                        GasChromatographySystems.new_pressure_steps(GCsys2)[:,2], 
                                                        GasChromatographySystems.new_gradient_parameter_steps(GCsys2[2], GasChromatographySystems.new_time_steps(GCsys2))[:,1], 
                                                        GasChromatographySystems.new_gradient_parameter_steps(GCsys2[2], GasChromatographySystems.new_time_steps(GCsys2))[:,2], 
                                                        GasChromatographySystems.new_gradient_parameter_steps(GCsys2[2], GasChromatographySystems.new_time_steps(GCsys2))[:,3], 
                                                        GasChromatographySystems.new_gradient_parameter_steps(GCsys2[2], GasChromatographySystems.new_time_steps(GCsys2))[:,4], 
                                                        opt.Tcontrol, 
                                                        sys.L)
sub = GasChromatographySimulator.load_solute_database(db, sys.sp, sys.gas, solutes, zeros(length(solutes)), zeros(length(solutes)))
par2_Simulator = GasChromatographySimulator.Parameters(sys, prog, sub, opt)

sol2_Simulator = GasChromatographySimulator.solve_system_multithreads(par2_Simulator)
tR2_C10_Simulator = sol2_Simulator[1].u[end][1]

tR2_C10 == tR2_C10_Simulator