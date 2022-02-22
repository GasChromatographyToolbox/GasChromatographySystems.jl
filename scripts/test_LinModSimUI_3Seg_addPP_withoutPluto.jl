# test the more or less same code as in LinModSimUI_3Seg_addPP.jl
# which results in an error at teh execution of the simulation (linear_GC_system_simulation())

using Plots, GasChromatographySystems, GasChromatographyTools

using GasChromatographySimulator
# Options
gas = "He"
alg = OwrenZen5()#"OwrenZen5"
abs = -6
rel = -4
Tcontrol = "inlet"
odesys = true
Option = GasChromatographySystems.Options(gas, alg, 10.0^abs, 10.0.^rel, Tcontrol, odesys);

# Columns 
sp1 = "Rxi17SilMS"
sp2 = "Rxi17SilMS"
sp3 = ""
L1 = 0.16
L2 = 3.77
L3 = 0.395
d1 = 0.11
d2 = 0.11
d3 = 0.1
df1 = 0.1
df2 = 0.1
df3 = 0.0
T1 = 300.0
T3 = 300.0

# Temperature Program
tsteps_str = "0 2 15 10 60 40 13 30"
Tsteps_str = "70 70 135 149 235 235 295 295"
ΔT_str = "30 30 52 46 30 30 52 55"
x₀_str = "0 0 0 0 0 0 0 0"
α_str = "-4 -3 -4 -8 -11 -13 -13 -14"

# Pressure Program
p_tsteps_str  = "0 2 15 70 40 13 30"
pin_str = "300 300 350 420 420 470 470"
pcon_str = "200 200 200 200 200 200 200"
pout_str = "0 0 0 0 0 0 0"

tsteps = parse.(Float64,split(tsteps_str))
Tsteps = parse.(Float64,split(Tsteps_str))
ΔT = parse.(Float64,split(ΔT_str))
x₀ = parse.(Float64,split(x₀_str))
L₀ = L2.*ones(length(ΔT))
α = parse.(Float64,split(α_str))
a = [ΔT x₀ L₀ α]
grad_func(x) = GasChromatographySimulator.gradient(x, a; Tcontrol=Option.Tcontrol)
TP = GasChromatographySystems.Temperature_Program(tsteps, Tsteps, grad_func, a)

TL1 = GasChromatographySystems.Transferline(L1, d1*1e-3, df1*1e-6, sp1, T1)
a_d_2 = [d2*1e-3]
d_2(x) = GasChromatographySimulator.gradient(x, a_d_2)
a_df_2 = [df2*1e-6]
df_2(x) = GasChromatographySimulator.gradient(x, a_df_2)
GC  = GasChromatographySystems.Column(L2, d_2, a_d_2, df_2, a_df_2, sp2, TP)
TL2 = GasChromatographySystems.Transferline(L3, d3*1e-3, df3*1e-6, sp3, T3)

p_tsteps = parse.(Float64,split(p_tsteps_str))
pin_steps = parse.(Float64,split(pin_str))
pcon_steps = parse.(Float64,split(pcon_str))
pout_steps = parse.(Float64,split(pout_str))
PPin = GasChromatographySystems.Pressure_Point(p_tsteps, pin_steps.*1000.0.+101300.0)
PPcon = GasChromatographySystems.Pressure_Point(p_tsteps, pcon_steps.*1000.0.+101300.0)
PPout = GasChromatographySystems.Pressure_Point(p_tsteps, pout_steps.*1000.0)

GCsys = PPin, TL1, GC, PPcon, TL2, PPout

GasChromatographySystems.test_of_GCsys(GCsys)

db = DataFrame(CSV.File(string(pwd(),"/data/Database_append.csv")))
com_solutes = GasChromatographySystems.common_solutes(db, GCsys)

function PAH_selection(com_solutes)
    PAH = ["Naphthalin", "Acenaphthylene", "Acenaphthene", "Fluorene", "Phenanthrene", "Anthracene", "Fluoranthene", "Pyrene", "Benz[a]anthracene", "Chrysene", "Benzo[b]fluoranthene", "Benzo[k]fluoranthene", "Benzo[a]pyrene", "Dibenzo[a,h]anthracene", "Indeno[1,2,3-cd]pyrene", "Benzo[ghi]perylene"]
    PAHs = String[]
    for i=1:length(com_solutes)
        if !isa(findfirst(occursin.(com_solutes[i], PAH)), Nothing)
            push!(PAHs, com_solutes[i])
        end
    end
    return PAHs
end

solutes = PAH_selection(com_solutes)

# run the simulation
par, sol = GasChromatographySystems.linear_GC_system_simulation(GCsys, Option, solutes, db)
