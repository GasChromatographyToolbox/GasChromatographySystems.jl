### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ be16a106-f3c0-11ed-2fa4-7daccc96b7c2
begin
	import Pkg
    # activate the shared project environment
	Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()

	using CSV, DataFrames
	using Plots, CairoMakie, GraphMakie
	using Graphs, NetworkLayout, Symbolics
	using GasChromatographySimulator
	using PlutoUI
	include(joinpath(dirname(pwd()), "src", "GasChromatographySystems.jl"))
	TableOfContents()
end

# ╔═╡ 2460130a-1b09-4dde-99a6-2bc609f6df28
md"""
# Simplified GCxGC simulation
** for the poster of Tillman Brehmer presented at the 20th GCxGC Symposium in Canmore, Canada, 2023.**

Using isothermal determined retention parameters (RetentionData database) and parameters calculated from LSER data (UFZ database and Poole database) and comparing GCxGC simulations of both retention datasets.
"""

# ╔═╡ 361ad6f5-f412-45a2-8619-89a0dbcecff5
md"""
## Definition of the GCxGC system
"""

# ╔═╡ 50d291e9-a5d8-4897-a4be-6b57b0670bdf
begin
	tsteps_, Tsteps_ = GasChromatographySimulator.conventional_program([40.0, 2.0, 10.0, 340.0, 5.0]) 
	TP1 = GasChromatographySystems.TemperatureProgram(tsteps_, Tsteps_)
end

# ╔═╡ 1c3291a6-d0e7-4f68-808b-e39b16157a53
phase = ["Rxi5SilMS", "Rxi17SilMS"]

# ╔═╡ 7bc5ff7b-a6aa-4898-b8e2-70375ce16670
sys = GasChromatographySystems.GCxGC_TM_simp(20.0, 0.25, 0.5, phase[1], TP1, 2.0, 0.25, 0.25, phase[2], TP1, 1.4, NaN, 101.3; opt=GasChromatographySystems.Options(ng=true))

# ╔═╡ bd9279f2-4b24-49d2-b7ae-edebdfedcca6
GasChromatographySystems.plot_graph_with_flow(sys, 0; lay=SquareGrid(cols=3), node_size=90, nlabels_fontsize=20, elabels_fontsize=20, elabels_distance=28)

# ╔═╡ dd5d4628-ab5f-4096-8799-73c2ff0c8fbf
GasChromatographySystems.plot_graph_with_flow(sys, sum(tsteps_); lay=SquareGrid(cols=3), node_size=90, nlabels_fontsize=20, elabels_fontsize=20, elabels_distance=28)

# ╔═╡ c136c0f9-e243-4480-8314-12e58eb5c96f
GasChromatographySystems.plot_pressure_over_time(sys)

# ╔═╡ a55e113e-d93a-48ba-a50a-bbb6f137d390
md"""
## Retention data
"""

# ╔═╡ 8f4a9352-7e59-48ca-b35a-fd40237b8083
begin
	path_RetentionData = "/Users/janleppert/Documents/sciebo/Konferenz/GCxGC Symposium Kanada 2023/Poster Tillman/GCSim_database_nonflag.csv"
	paths_Poole = ["/Users/janleppert/Documents/sciebo/Konferenz/GCxGC Symposium Kanada 2023/Poster Tillman/Database_Poole_Rxi5SilMS.csv", "/Users/janleppert/Documents/sciebo/Konferenz/GCxGC Symposium Kanada 2023/Poster Tillman/Database_Poole_Rxi17SilMS.csv"]
	paths_UFZ = ["/Users/janleppert/Documents/sciebo/Konferenz/GCxGC Symposium Kanada 2023/Poster Tillman/Database_UFZ_Rxi5SilMS.csv", "/Users/janleppert/Documents/sciebo/Konferenz/GCxGC Symposium Kanada 2023/Poster Tillman/Database_UFZ_Rxi17SilMS.csv"]

	db_RetentionData = DataFrame(CSV.File(path_RetentionData, header=1, silencewarnings=true))
	filter!([:phi0] => x -> x == 0.001, db_RetentionData)
	db_Poole = DataFrame()
	db_UFZ = DataFrame()
	for i=1:2
		db_ = DataFrame(CSV.File(paths_Poole[i], header=1, silencewarnings=true))
		db_[!, :Phase] = fill(phase[i], length(db_.Name))
		db_[!, :phi0] = fill(0.001, length(db_.Name))
		db_[!, :Source] = fill("Poole", length(db_.Name))
		append!(db_Poole, db_)
		db__ = DataFrame(CSV.File(paths_UFZ[i], header=1, silencewarnings=true))
		db__[!, :Phase] = fill(phase[i], length(db__.Name))
		db__[!, :phi0] = fill(0.001, length(db__.Name))
		db__[!, :Source] = fill("UFZ", length(db__.Name))
		append!(db_UFZ, db__)
	end
	db = [db_RetentionData, db_Poole, db_UFZ]
end

# ╔═╡ 23b419b8-43c1-4ecd-8e06-9d5202e40d18
begin
	solutes = Array{DataFrame}(undef, 3)
	for i=1:3
		solutes[i] = GasChromatographySystems.common_solutes(db[i], sys)
	end
	solutes
end

# ╔═╡ fa180992-b2d6-4603-90bf-0eff6376f5dd
begin
	common_CAS = GasChromatographySimulator.common(GasChromatographySimulator.common(solutes[1].CAS, solutes[2].CAS), solutes[3].CAS)
	common_solutes = Array{String,1}(undef, length(common_CAS))
	for i=1:length(common_CAS)
		common_solutes[i] = solutes[1].Name[findfirst(common_CAS[i].==solutes[1].CAS)]
	end
	common_solutes
end

# ╔═╡ 3ef47737-51b0-45fc-9147-dc06dcb68ce1
md"""
## Parameters 
"""

# ╔═╡ db1ad723-7627-4831-9b47-ab380a74f239
begin
	par = Array{Array{GasChromatographySimulator.Parameters,1}}(undef, 3)
	for i=1:3
		par[i] = GasChromatographySystems.graph_to_parameters(sys, db[i], common_solutes)
	end
	par
end

# ╔═╡ 419bcc59-74b6-4f91-bdbd-e680104d0fed
par[3][1]

# ╔═╡ 4f9d42f3-e617-4ded-a884-b8c45f03e376
par[3][2]

# ╔═╡ 311690d7-556a-46c2-9e40-8676f1e73b2f
md"""
## Simulations
"""

# ╔═╡ 0022c392-5a5a-400e-9da5-5f6d0fd70549
paths = GasChromatographySystems.all_paths(sys.g, 1)[2]

# ╔═╡ 549a03ec-1d89-460c-acff-2812b160c377
begin
	pl = Array{Array{Array{DataFrame,1},1}}(undef, 3)
	for i=1:3
		pl[i] = GasChromatographySystems.simulate_along_paths(sys, paths, par[i]; refocus=trues(ne(sys.g)), τ₀_focus=0.1.*ones(length(par[i][1].sub)))[2]
	end
	pl
end

# ╔═╡ 447c50d6-527d-44cc-a49d-72ca6428d01e
pl_GCxGC_UFZ = GasChromatographySystems.peaklist_GCxGC(pl[3][1][1], pl[3][1][2])

# ╔═╡ 455fb9bc-04c6-4fbf-97e0-a47eb27cf8f1
md"""
## different solute groups
"""

# ╔═╡ d5463097-925d-424d-b791-cf4987561422
begin
	PCBs = ["PCB 28", "PCB 52", "PCB 101", "PCB 138", "PCB 153", "PCB 180"]
	Phenones = ["Propiophenone", "Valerophenone", "Hexanophenone", "Heptanophenone", "Octanophenone"]
	PAHs = ["Naphthalene", "Acenaphthylene", "Acenaphthene", "Fluorene", "Phenanthrene", "Anthracene", "Fluoranthene", "Pyrene", "Benz[a]anthracene", "Chrysene", "Benzo[b]fluoranthene", "Benzo[k]fluoranthene", "Benzo[a]pyrene", "Indeno[1,2,3-cd]pyrene", "Dibenz[a,h]anthracene", "Benzo[ghi]perylene"]
	Alkanes = ["Octane", "Nonane", "Decane", "Undecane", "Dodecane", "Tridecane", "Tetradecane", "Pentadecane", "Hexadecane", "Heptadecane", "Octadecane", "Nonadecane", "Eicosane"]
	FAMEs = ["Methyl butyrate", "Methyl hexanoate", "Methyl octanoate", "Methyl decanoate", "Methyl undecanoate", "Methyl laurate", "Methyl tridecanoate", "Methyl myristate", "Methyl pentadecanoate", "Methyl palmitate", "Methyl heptadecanoate", "Methyl stearate", "cis-9-Oleic acid methyl ester", "Methyl linoleate"]
	Alcohols = ["Pentanol", "Hexan-1-ol", "Heptanol", "Octan-1-ol", "Nonanol", "Decanol", "Undecanol", "Dodecanol"]
	solute_groups = [PCBs, Phenones, PAHs, Alkanes, FAMEs, Alcohols]
end

# ╔═╡ 14d95732-117e-4e47-a337-b7481d2b7783
md"""
### Parameters
"""

# ╔═╡ 4c37a861-475b-4940-9eda-6d345fa2c082
begin
	par_groups = Array{Array{GasChromatographySimulator.Parameters,1}}(undef, 3, length(solute_groups))
	for j=1:length(solute_groups)
		for i=1:3
			par_groups[i,j] = GasChromatographySystems.graph_to_parameters(sys, db[i], solute_groups[j])
		end
	end
end

# ╔═╡ a50674d3-9338-474f-97e4-b8498a550933
begin
	pl_groups = Array{Array{Array{DataFrame,1},1}}(undef, 3, length(solute_groups))
	for j=1:length(solute_groups)
		for i=1:3
			pl_groups[i,j] = GasChromatographySystems.simulate_along_paths(sys, paths, par_groups[i,j]; refocus=trues(ne(sys.g)), τ₀_focus=0.1.*ones(length(par_groups[i,j][1].sub)))[2]
		end
	end
	pl_groups
end

# ╔═╡ 5c97c65f-bfb1-44ab-9ddb-55020f14784a
par_groups[3,6]

# ╔═╡ 9126d196-d1a5-423e-bfe9-9d3df112e7be
par_groups[1,6]

# ╔═╡ 6c33c1bb-e6df-47e9-a690-f22f9ecd0a2f
begin
	par_Alcohol_RetentionData_D1 = GasChromatographySimulator.Parameters(par_groups[1,6][1].col, par_groups[1,6][1].prog, par_groups[1,6][1].sub, par_groups[1,6][1].opt)
	par_Alcohol_RetentionData_D2 = GasChromatographySimulator.Parameters(par_groups[1,6][2].col, par_groups[1,6][2].prog, par_groups[1,6][2].sub, par_groups[1,6][2].opt)
	par_Alcohol_UFZ_D1 = GasChromatographySimulator.Parameters(par_groups[3,6][1].col, par_groups[3,6][1].prog, par_groups[3,6][1].sub[1:8], par_groups[3,6][1].opt)
	par_Alcohol_UFZ_D2 = GasChromatographySimulator.Parameters(par_groups[3,6][2].col, par_groups[3,6][2].prog, par_groups[3,6][2].sub[1:8], par_groups[3,6][2].opt)
	par_Alcohol_Poole_D1 = GasChromatographySimulator.Parameters(par_groups[2,6][1].col, par_groups[2,6][1].prog, par_groups[2,6][1].sub, par_groups[2,6][1].opt)
	par_Alcohol_Poole_D2 = GasChromatographySimulator.Parameters(par_groups[2,6][2].col, par_groups[2,6][2].prog, par_groups[2,6][2].sub, par_groups[2,6][2].opt)
end

# ╔═╡ cf2794a8-9e39-4035-9651-db610b2ae654
begin
	par_Alkane_RetentionData_D1 = par_groups[1,4][1]
	par_Alkane_RetentionData_D2 = par_groups[1,4][2]
	par_Alkane_Poole_D1 = par_groups[2,4][1]
	par_Alkane_Poole_D2 = par_groups[2,4][2]
	par_Alkane_UFZ_D1 = par_groups[3,4][1]
	par_Alkane_UFZ_D2 = par_groups[3,4][2]
end

# ╔═╡ 6ba883bb-4681-4358-bbd6-73d2abc77f32
begin
	par_FAME_RetentionData_D1 = GasChromatographySimulator.Parameters(par_groups[1,5][1].col, par_groups[1,5][1].prog, par_groups[1,5][1].sub[1:14], par_groups[1,5][1].opt)
	par_FAME_RetentionData_D2 = GasChromatographySimulator.Parameters(par_groups[1,5][2].col, par_groups[1,5][2].prog, par_groups[1,5][2].sub[1:14], par_groups[1,5][2].opt)
	par_FAME_UFZ_D1 = GasChromatographySimulator.Parameters(par_groups[3,5][1].col, par_groups[3,5][1].prog, par_groups[3,5][1].sub[1:14], par_groups[3,5][1].opt)
	par_FAME_UFZ_D2 = GasChromatographySimulator.Parameters(par_groups[3,5][2].col, par_groups[3,5][2].prog, par_groups[3,5][2].sub[1:14], par_groups[3,5][2].opt)
	par_FAME_Poole_D1 = GasChromatographySimulator.Parameters(par_groups[2,5][1].col, par_groups[2,5][1].prog, par_groups[2,5][1].sub[1:3], par_groups[2,5][1].opt)
	par_FAME_Poole_D2 = GasChromatographySimulator.Parameters(par_groups[2,5][2].col, par_groups[2,5][2].prog, par_groups[2,5][2].sub[1:3], par_groups[2,5][2].opt)
end

# ╔═╡ fcba3fb7-1778-47c3-89a9-b68fe8078f81
par_Phenone_RetentionData_D1 = GasChromatographySimulator.Parameters(par_groups[1,2][1].col, par_groups[1,2][1].prog, par_groups[1,2][1].sub, par_groups[1,2][1].opt)

# ╔═╡ 0e71379e-b682-4fb6-80e9-33c79df879e7
par_Phenone_RetentionData_D2 = GasChromatographySimulator.Parameters(par_groups[1,2][2].col, par_groups[1,2][2].prog, par_groups[1,2][2].sub, par_groups[1,2][2].opt)

# ╔═╡ 32e274da-e0ec-4112-ac5f-5d952f079b25
par_Phenone_UFZ_D1 = GasChromatographySimulator.Parameters(par_groups[3,2][1].col, par_groups[3,2][1].prog, par_groups[3,2][1].sub, par_groups[3,2][1].opt)

# ╔═╡ 7938448f-9dc3-4019-8fdc-06440d525b06
par_Phenone_UFZ_D2 = GasChromatographySimulator.Parameters(par_groups[3,2][2].col, par_groups[3,2][2].prog, par_groups[3,2][2].sub, par_groups[3,2][2].opt)

# ╔═╡ 6bad5fde-5efd-4aa2-a78a-ec6c32002aa8
par_PAH_RetentionData_D1 = GasChromatographySimulator.Parameters(par_groups[1,3][1].col, par_groups[1,3][1].prog, [par_groups[1,3][1].sub[1:10]; par_groups[1,3][1].sub[12:14]], par_groups[1,3][1].opt)

# ╔═╡ 4d976b8f-55d9-4e37-9b42-58d130d5c0d3
par_PAH_RetentionData_D2 = GasChromatographySimulator.Parameters(par_groups[1,3][2].col, par_groups[1,3][2].prog, [par_groups[1,3][2].sub[1:10]; par_groups[1,3][2].sub[12:14]], par_groups[1,3][2].opt)

# ╔═╡ 94cf7ea2-eb1c-46bc-9625-6fcaccc7f771
par_PAH_UFZ_D1 = GasChromatographySimulator.Parameters(par_groups[3,3][1].col, par_groups[3,3][1].prog, par_groups[3,3][1].sub[1:13], par_groups[3,3][1].opt)

# ╔═╡ f31a0c40-908e-47a9-a3dd-fac6be8904bc
par_PAH_UFZ_D2 = GasChromatographySimulator.Parameters(par_groups[3,3][2].col, par_groups[3,3][2].prog, par_groups[3,3][2].sub[1:13], par_groups[3,3][2].opt)

# ╔═╡ 8d9ab99d-091b-44fc-905f-3d21f8cb3f00
par_PAH_Poole_D1 = GasChromatographySimulator.Parameters(par_groups[2,3][1].col, par_groups[2,3][1].prog, par_groups[2,3][1].sub[1:8], par_groups[2,3][1].opt)

# ╔═╡ 8f7672ac-c8e0-4583-853e-1b4bc9471e47
par_PAH_Poole_D2 = GasChromatographySimulator.Parameters(par_groups[2,3][2].col, par_groups[2,3][2].prog, par_groups[2,3][2].sub[1:8], par_groups[2,3][2].opt)

# ╔═╡ 8f4d9af7-0c94-4246-8509-de2b8f3c03e3
par_PCB_RetentionData_D1 = GasChromatographySimulator.Parameters(par_groups[1,1][1].col, par_groups[1,1][1].prog, par_groups[1,1][1].sub[1:6], par_groups[1,1][1].opt)

# ╔═╡ 9e3ab0b6-e0b1-4732-a34e-1900635d45e5
par_PCB_RetentionData_D2 = GasChromatographySimulator.Parameters(par_groups[1,1][2].col, par_groups[1,1][2].prog, par_groups[1,1][2].sub[1:6], par_groups[1,1][2].opt)

# ╔═╡ c865cde9-9c7b-49bd-bb85-89677e002a2b
par_PCB_UFZ_D1 = GasChromatographySimulator.Parameters(par_groups[3,1][1].col, par_groups[3,1][1].prog, par_groups[3,1][1].sub[1:6], par_groups[3,1][1].opt)

# ╔═╡ 7b3f5215-5a5c-43a4-8519-6cf7ab2941ae
par_PCB_UFZ_D2 = GasChromatographySimulator.Parameters(par_groups[3,1][2].col, par_groups[3,1][2].prog, par_groups[3,1][2].sub[1:6], par_groups[3,1][2].opt)

# ╔═╡ 03b3633f-a1d6-4951-83b5-6c2332183e2c
md"""
### Simulations
"""

# ╔═╡ 30c01f19-1892-4e7a-9fe2-7b66cfbc37a8
md"""
### PCB
"""

# ╔═╡ 58237852-30da-477c-a10f-155407f47cd2
pl_PCB_RetentionData = GasChromatographySystems.simulate_along_paths(sys, paths, [par_PCB_RetentionData_D1, par_PCB_RetentionData_D2]; refocus=trues(ne(sys.g)), τ₀_focus=0.1.*ones(length(par_PCB_RetentionData_D1.sub)))[2]

# ╔═╡ d7d78054-d1c7-4d88-971d-cd809aae18dc
pl_PCB_UFZ = GasChromatographySystems.simulate_along_paths(sys, paths, [par_PCB_UFZ_D1, par_PCB_UFZ_D2]; refocus=trues(ne(sys.g)), τ₀_focus=0.1.*ones(length(par_PCB_UFZ_D1.sub)))[2]

# ╔═╡ 351cb0fe-fcaa-46fd-a6e3-caef23e75e47
pl_GCxGC_PCB_RetentionData = GasChromatographySystems.peaklist_GCxGC(pl_PCB_RetentionData[1][1], pl_PCB_RetentionData[1][2])

# ╔═╡ e8a09fee-1f13-4eae-b070-4225945db950
pl_GCxGC_PCB_UFZ = GasChromatographySystems.peaklist_GCxGC(pl_PCB_UFZ[1][1], pl_PCB_UFZ[1][2])

# ╔═╡ 2d817ece-a5c7-4f72-a949-2248fd70e9f3
md"""
### PAH
"""

# ╔═╡ 53bf9e63-fbdd-480a-9a3f-50c1dea25871
pl_PAH_RetentionData = GasChromatographySystems.simulate_along_paths(sys, paths, [par_PAH_RetentionData_D1, par_PAH_RetentionData_D2]; refocus=trues(ne(sys.g)), τ₀_focus=0.1.*ones(length(par_PAH_RetentionData_D1.sub)))[2]

# ╔═╡ 38dccb75-80eb-4f1d-a859-608da5e1fbac
pl_PAH_UFZ = GasChromatographySystems.simulate_along_paths(sys, paths, [par_PAH_UFZ_D1, par_PAH_UFZ_D2]; refocus=trues(ne(sys.g)), τ₀_focus=0.1.*ones(length(par_PAH_UFZ_D1.sub)))[2]

# ╔═╡ 02b49171-0068-4a5d-8002-c2b2fe266bf8
pl_PAH_Poole = GasChromatographySystems.simulate_along_paths(sys, paths, [par_PAH_Poole_D1, par_PAH_Poole_D2]; refocus=trues(ne(sys.g)), τ₀_focus=0.1.*ones(length(par_PAH_Poole_D1.sub)))[2]

# ╔═╡ 72decd94-304d-474b-b793-71010be7d2c6
pl_GCxGC_PAH_RetentionData = GasChromatographySystems.peaklist_GCxGC(pl_PAH_RetentionData[1][1], pl_PAH_RetentionData[1][2])

# ╔═╡ 8d0d0823-44e2-4f65-8c81-f039ecfc4b0d
pl_GCxGC_PAH_UFZ = GasChromatographySystems.peaklist_GCxGC(pl_PAH_UFZ[1][1], pl_PAH_UFZ[1][2])

# ╔═╡ 4fef2f8d-e91f-4e30-b183-0196aad9ec17
pl_GCxGC_PAH_Poole = GasChromatographySystems.peaklist_GCxGC(pl_PAH_Poole[1][1], pl_PAH_Poole[1][2])

# ╔═╡ 34c4a88e-78d2-4232-9f35-51139df8b2d3
md"""
### Phenone
"""

# ╔═╡ 1d1d772a-60d6-4951-a0c0-7e0045d2e107
pl_Phenone_RetentionData = GasChromatographySystems.simulate_along_paths(sys, paths, [par_Phenone_RetentionData_D1, par_Phenone_RetentionData_D2]; refocus=trues(ne(sys.g)), τ₀_focus=0.1.*ones(length(par_Phenone_RetentionData_D1.sub)))[2]

# ╔═╡ 686fc196-3a39-4e98-8607-edbb1445e930
pl_Phenone_UFZ = GasChromatographySystems.simulate_along_paths(sys, paths, [par_Phenone_UFZ_D1, par_Phenone_UFZ_D2]; refocus=trues(ne(sys.g)), τ₀_focus=0.1.*ones(length(par_Phenone_UFZ_D1.sub)))[2]

# ╔═╡ 544a02a1-25db-47f1-a6e1-940e11457d2b
pl_GCxGC_Phenone_RetentionData = GasChromatographySystems.peaklist_GCxGC(pl_Phenone_RetentionData[1][1], pl_Phenone_RetentionData[1][2])

# ╔═╡ 761dc3e3-36f8-4c36-b7d8-ffb17d9a16bb
pl_GCxGC_Phenone_UFZ = GasChromatographySystems.peaklist_GCxGC(pl_Phenone_UFZ[1][1], pl_Phenone_UFZ[1][2])

# ╔═╡ ad448ab0-3789-479e-8df8-4cda0d40d8a3
md"""
### Alkane
"""

# ╔═╡ 75548590-01e8-4d13-a03a-bb2352ce8a88
pl_Alkane_RetentionData = GasChromatographySystems.simulate_along_paths(sys, paths, [par_Alkane_RetentionData_D1, par_Alkane_RetentionData_D2]; refocus=trues(ne(sys.g)), τ₀_focus=0.1.*ones(length(par_Alkane_RetentionData_D1.sub)))[2]

# ╔═╡ e3d31326-4c69-412e-b9e2-6d7a0c685b46
pl_Alkane_Poole = GasChromatographySystems.simulate_along_paths(sys, paths, [par_Alkane_Poole_D1, par_Alkane_Poole_D2]; refocus=trues(ne(sys.g)), τ₀_focus=0.1.*ones(length(par_Alkane_Poole_D1.sub)))[2]

# ╔═╡ 7620a6a0-71c4-484a-8d77-925a2bd1b0ae
pl_Alkane_UFZ = GasChromatographySystems.simulate_along_paths(sys, paths, [par_Alkane_UFZ_D1, par_Alkane_UFZ_D2]; refocus=trues(ne(sys.g)), τ₀_focus=0.1.*ones(length(par_Alkane_UFZ_D1.sub)))[2]

# ╔═╡ 6be685cf-22fc-4723-ac78-c87c5c283e1f
pl_GCxGC_Alkane_RetentionData = GasChromatographySystems.peaklist_GCxGC(pl_Alkane_RetentionData[1][1], pl_Alkane_RetentionData[1][2])

# ╔═╡ bb174895-d065-4d6f-a7f0-d88c1f06ee44
pl_GCxGC_Alkane_Poole = GasChromatographySystems.peaklist_GCxGC(pl_Alkane_Poole[1][1], pl_Alkane_Poole[1][2])

# ╔═╡ 1d771f43-9e96-473f-ab03-ebe646f60aa1
pl_GCxGC_Alkane_UFZ = GasChromatographySystems.peaklist_GCxGC(pl_Alkane_UFZ[1][1], pl_Alkane_UFZ[1][2])

# ╔═╡ 6fb1f245-3285-4662-ac0d-57f9384dac6d
md"""
### FAME
"""

# ╔═╡ 6a678c63-4aa9-43b4-bdfd-5119dc14dea1
begin
	pl_FAME_RetentionData = GasChromatographySystems.simulate_along_paths(sys, paths, [par_FAME_RetentionData_D1, par_FAME_RetentionData_D2]; refocus=trues(ne(sys.g)), τ₀_focus=0.1.*ones(length(par_FAME_RetentionData_D1.sub)))[2]
	pl_FAME_UFZ = GasChromatographySystems.simulate_along_paths(sys, paths, [par_FAME_UFZ_D1, par_FAME_UFZ_D2]; refocus=trues(ne(sys.g)), τ₀_focus=0.1.*ones(length(par_FAME_UFZ_D1.sub)))[2]
	pl_FAME_Poole = GasChromatographySystems.simulate_along_paths(sys, paths, [par_FAME_Poole_D1, par_FAME_Poole_D2]; refocus=trues(ne(sys.g)), τ₀_focus=0.1.*ones(length(par_FAME_Poole_D1.sub)))[2]
end

# ╔═╡ ccceaebe-1438-4484-94bb-bf6c0c182a70
pl_GCxGC_FAME_RetentionData = GasChromatographySystems.peaklist_GCxGC(pl_FAME_RetentionData[1][1], pl_FAME_RetentionData[1][2])

# ╔═╡ d6966984-023f-4a73-a778-c0c66349fd69
pl_GCxGC_FAME_UFZ = GasChromatographySystems.peaklist_GCxGC(pl_FAME_UFZ[1][1], pl_FAME_UFZ[1][2])

# ╔═╡ a7e09a23-d501-41bc-969a-ae48fc64d28d
pl_GCxGC_FAME_Poole = GasChromatographySystems.peaklist_GCxGC(pl_FAME_Poole[1][1], pl_FAME_Poole[1][2])

# ╔═╡ 88381961-abf2-4fdc-ae14-ea6e60f9e292
md"""
### Alcohol
"""

# ╔═╡ 9a4a6e9c-9dcf-4f66-b947-d42cd1f861e9
begin
	pl_Alcohol_RetentionData = GasChromatographySystems.simulate_along_paths(sys, paths, [par_Alcohol_RetentionData_D1, par_Alcohol_RetentionData_D2]; refocus=trues(ne(sys.g)), τ₀_focus=0.1.*ones(length(par_Alcohol_RetentionData_D1.sub)))[2]
	pl_Alcohol_UFZ = GasChromatographySystems.simulate_along_paths(sys, paths, [par_Alcohol_UFZ_D1, par_Alcohol_UFZ_D2]; refocus=trues(ne(sys.g)), τ₀_focus=0.1.*ones(length(par_Alcohol_UFZ_D1.sub)))[2]
	pl_Alcohol_Poole = GasChromatographySystems.simulate_along_paths(sys, paths, [par_Alcohol_Poole_D1, par_Alcohol_Poole_D2]; refocus=trues(ne(sys.g)), τ₀_focus=0.1.*ones(length(par_Alcohol_Poole_D1.sub)))[2]
end

# ╔═╡ d8114bb8-ae36-42d3-82ee-8178ac8f09ba
pl_GCxGC_Alcohol_RetentionData = GasChromatographySystems.peaklist_GCxGC(pl_Alcohol_RetentionData[1][1], pl_Alcohol_RetentionData[1][2])

# ╔═╡ 2669af8b-0449-4753-9fe9-7de05e2bc13e
pl_GCxGC_Alcohol_UFZ = GasChromatographySystems.peaklist_GCxGC(pl_Alcohol_UFZ[1][1], pl_Alcohol_UFZ[1][2])

# ╔═╡ bf44b9c6-fdd2-45ee-9ea6-5cf6af03a88e
pl_GCxGC_Alcohol_Poole = GasChromatographySystems.peaklist_GCxGC(pl_Alcohol_Poole[1][1], pl_Alcohol_Poole[1][2])

# ╔═╡ 4afeb72c-bf72-487c-a1af-f51752316a50
md"""
### Combine Chromatograms
"""

# ╔═╡ 6d44e4a0-745e-486e-a8d7-750150d57687
begin
	pl_RetentionData = DataFrame()
	append!(pl_RetentionData, pl_GCxGC_PCB_RetentionData)
	append!(pl_RetentionData, pl_GCxGC_PAH_RetentionData)
	append!(pl_RetentionData, pl_GCxGC_Phenone_RetentionData)
	append!(pl_RetentionData, pl_GCxGC_Alkane_RetentionData)
	append!(pl_RetentionData, pl_GCxGC_FAME_RetentionData)
	append!(pl_RetentionData, pl_GCxGC_Alcohol_RetentionData)
end

# ╔═╡ b2c43093-7c09-4bf3-9211-20f73c65006c
#savefig(p_GCxGC_RetentionData, "chrom_GCxGC_RetentionData_vs_LSER-UFZ_cross.svg")

# ╔═╡ f72c5e04-4bff-40ab-a430-8158f78a3f16
Plots.heatmap([0.0, 2099.0, 2100.0], [0.0, 0.0, 15.0], [0 0 0;0 1 0; 0 0 0], c=:jet1)

# ╔═╡ 913e485f-762f-495d-a35e-633e63c22707
md"""
# End
"""

# ╔═╡ bb9db28c-4bcc-4b0a-aac2-422531d7b203
function peaklist_GCxGC_with_categories(pl_D1, pl_D2, db)
	pl = GasChromatographySystems.peaklist_GCxGC(pl_D1, pl_D2)
	cat = Array{String}(undef, length(pl.Name))
	for i=1:length(pl.Name)
		ii = findfirst(pl.Name[i].==db.Name)
		cat[i] = ""
		for j=1:length(findall(occursin.("Cat", names(db))))
			if !ismissing(db[ii, Symbol("Cat_$(j)")])
				if j==1
					cat[i] = cat[i]*db[ii, Symbol("Cat_$(j)")]
				else
					cat[i] = cat[i]*", "*db[ii, Symbol("Cat_$(j)")]
				end
			end
		end
	end
	pl[!, :Cat] = cat
	return pl
end

# ╔═╡ adba7751-d82b-4136-acad-e34dafcc41e6
pl_GCxGC = peaklist_GCxGC_with_categories(pl[1][1][1], pl[1][1][2], db[1])

# ╔═╡ c9a523b7-ded8-45b5-bee9-1c08d03fd521
pl_GCxGC_Poole = peaklist_GCxGC_with_categories(pl[2][1][1], pl[2][1][2], db[1])

# ╔═╡ c07421e4-7d14-4e83-a0d3-08890d8337fa
function plot_GCxGC(pl_GCxGC, sys; categories = String[])
	t1end = sum(GasChromatographySystems.common_timesteps(sys))
	t2end = ceil(maximum(pl_GCxGC.tR2.-pl_GCxGC.tR1))
	p_gcxgc = Plots.scatter(pl_GCxGC.tR1, pl_GCxGC.tR2.-pl_GCxGC.tR1, xlims=(0.0,t1end), ylims=(0.0, t2end), label="", xlabel=string(sys.modules[1].stationary_phase, ", t¹ in s"), ylabel=string(sys.modules[2].stationary_phase, ", t² in s"))
	if !isempty(categories)
		for i=1:length(categories) 
			pl_f = filter([:Cat] => x -> occursin(categories[i], x), pl_GCxGC)
			Plots.scatter!(p_gcxgc, pl_f.tR1, pl_f.tR2.-pl_f.tR1, label=categories[i], legend=:bottomright)
		end
	end
	return p_gcxgc
end

# ╔═╡ 26e02657-9cf3-4244-8dd1-b157245bd900
begin
	plotly()
	p_GCxGC = plot_GCxGC(pl_GCxGC, sys; categories = ["FAME", "aromatic", "PAH", "PCB", "alkanes", "alcohols", "alkanones", "phenones"])
	#Plots.plot!(p_GCxGC, legend=:topleft)
	p_GCxGC
end

# ╔═╡ 8dca230b-6184-4766-bcb8-5fc13659e91a
function plot_GCxGC_contour_(pl_GCxGC; split_ratio=ones(length(pl_GCxGC.Name)), mod_τ1=1.0)
	t¹ = 0.9*minimum(pl_GCxGC.tR1):(1.1*maximum(pl_GCxGC.tR1)-0.9*minimum(pl_GCxGC.tR1))/1000:1.1*maximum(pl_GCxGC.tR1)
	t² = 0.9*minimum(pl_GCxGC.tR2.-pl_GCxGC.tR1):(1.1*maximum(pl_GCxGC.tR2.-pl_GCxGC.tR1)-0.9*minimum(pl_GCxGC.tR2.-pl_GCxGC.tR1))/1000:1.1*maximum(pl_GCxGC.tR2.-pl_GCxGC.tR1)
	chrom2D = Array{Float64}(undef, length(t¹), length(t²), length(pl_GCxGC.Name))
	for j=1:length(t²)
		for i=1:length(t¹)
			for k=1:length(pl_GCxGC.Name)
				chrom2D[i,j,k] = split_ratio[k]*exp(-((t¹[i]-pl_GCxGC.tR1[k])^2/(2*(mod_τ1*pl_GCxGC.τR1[k])^2)+(t²[j]-(pl_GCxGC.tR2[k]-pl_GCxGC.tR1[k]))^2/(2*pl_GCxGC.τR2[k]^2)))
			end
		end
	end
	#p_2D = Plots.contourf(t¹, t², sum(chrom2D, dims=3)[:,:,1]', fill=true, legend=false, xlabel="tR1 in s", ylabel="tR2 in s", c=:turbo, levels=40)
	return t¹, t², sum(chrom2D, dims=3)[:,:,1]'
end

# ╔═╡ c5aebbf9-5fa6-4871-a811-9e041bc9ebb7
t¹, t², chrom = plot_GCxGC_contour_(pl_GCxGC; mod_τ1=1.5)

# ╔═╡ a018f1ef-b49b-4bfb-9391-c38be7bbfb5e
begin
	plotly()
	p_GCxGC_ = Plots.heatmap(t¹, t², chrom.^(1//3), legend=false, xlabel="tR1 in s", ylabel="tR2 in s", c=:jet1, dpi=500)#, xlims=(0.0, 1700.0), ylims=(0.0, 5.0), grid=false)
	#gui()
end

# ╔═╡ 19bffce9-397a-4a07-ac53-9ddb6500b0c5
begin
	Plots.scatter!(p_GCxGC_, pl_GCxGC_Poole.tR1, pl_GCxGC_Poole.tR2.-pl_GCxGC_Poole.tR1)
	Plots.scatter!(p_GCxGC_, pl_GCxGC_UFZ.tR1, pl_GCxGC_UFZ.tR2.-pl_GCxGC_UFZ.tR1)
	p_GCxGC_
end

# ╔═╡ acc3dbc8-f66a-4773-a83d-9a66beaf02fe
t¹_PCB, t²_PCB, chrom_PCB = plot_GCxGC_contour_(pl_GCxGC_PCB_RetentionData; mod_τ1=1.5)

# ╔═╡ a12957c6-c0ec-46b6-aae4-1f992eb25037
begin
	plotly()
	p_GCxGC_PCB = Plots.heatmap(t¹_PCB, t²_PCB, chrom_PCB.^(1//1), legend=false, xlabel="tR1 in s", ylabel="tR2 in s", c=:jet1, dpi=500)#, xlims=(0.0, 1700.0), ylims=(0.0, 5.0), grid=false)
	#gui()
	Plots.scatter!(p_GCxGC_PCB, pl_GCxGC_PCB_UFZ.tR1, pl_GCxGC_PCB_UFZ.tR2.-pl_GCxGC_PCB_UFZ.tR1)
	p_GCxGC_PCB
end

# ╔═╡ 618cea77-6724-48fe-a1d9-fa2be2428ba2
t¹_PAH, t²_PAH, chrom_PAH = plot_GCxGC_contour_(pl_GCxGC_PAH_RetentionData; mod_τ1=1.5)

# ╔═╡ bd3824cd-5dc9-42e8-aaef-6b0750fa0898
begin
	plotly()
	p_GCxGC_PAH = Plots.heatmap(t¹_PAH, t²_PAH, chrom_PAH.^(1//1), legend=false, xlabel="tR1 in s", ylabel="tR2 in s", c=:jet1, dpi=500)#, xlims=(0.0, 1700.0), ylims=(0.0, 5.0), grid=false)
	#gui()
	Plots.scatter!(p_GCxGC_PAH, pl_GCxGC_PAH_UFZ.tR1, pl_GCxGC_PAH_UFZ.tR2.-pl_GCxGC_PAH_UFZ.tR1)
	Plots.scatter!(p_GCxGC_PAH, pl_GCxGC_PAH_Poole.tR1, pl_GCxGC_PAH_Poole.tR2.-pl_GCxGC_PAH_Poole.tR1)
	p_GCxGC_PAH
end

# ╔═╡ b77d4a3b-9933-453c-b318-01f24658ee58
t¹_Phenone, t²_Phenone, chrom_Phenone = plot_GCxGC_contour_(pl_GCxGC_Phenone_RetentionData; mod_τ1=1.5)

# ╔═╡ 25262675-b410-4f1e-9937-6cd18886bdca
begin
	plotly()
	p_GCxGC_Phenone = Plots.heatmap(t¹_Phenone, t²_Phenone, chrom_Phenone.^(1//1), legend=false, xlabel="tR1 in s", ylabel="tR2 in s", c=:jet1, dpi=500)#, xlims=(0.0, 1700.0), ylims=(0.0, 5.0), grid=false)
	#gui()
	Plots.scatter!(p_GCxGC_Phenone, pl_GCxGC_Phenone_UFZ.tR1, pl_GCxGC_Phenone_UFZ.tR2.-pl_GCxGC_Phenone_UFZ.tR1)
	p_GCxGC_Phenone
end

# ╔═╡ e0846d06-6d97-4ece-8b11-e4b4476892bb
t¹_Alkane, t²_Alkane, chrom_Alkane = plot_GCxGC_contour_(pl_GCxGC_Alkane_RetentionData; mod_τ1=1.5)

# ╔═╡ 2b214b26-06b0-41b8-b72b-0ef7779309f2
begin
	plotly()
	p_GCxGC_Alkane = Plots.heatmap(t¹_Alkane, t²_Alkane, chrom_Alkane.^(1//1), legend=false, xlabel="tR1 in s", ylabel="tR2 in s", c=:jet1, dpi=500)#, xlims=(0.0, 1700.0), ylims=(0.0, 5.0), grid=false)
	#gui()
	Plots.scatter!(p_GCxGC_Alkane, pl_GCxGC_Alkane_UFZ.tR1, pl_GCxGC_Alkane_UFZ.tR2.-pl_GCxGC_Alkane_UFZ.tR1)
	Plots.scatter!(p_GCxGC_Alkane, pl_GCxGC_Alkane_Poole.tR1, pl_GCxGC_Alkane_Poole.tR2.-pl_GCxGC_Alkane_Poole.tR1)
	p_GCxGC_Alkane
end

# ╔═╡ 4c404315-70c1-42cf-a944-043e5f03c12a
t¹_FAME, t²_FAME, chrom_FAME = plot_GCxGC_contour_(pl_GCxGC_FAME_RetentionData; mod_τ1=1.5)

# ╔═╡ a2e7c37b-39d1-40f8-afb9-9064d790e655
begin
	plotly()
	p_GCxGC_FAME = Plots.heatmap(t¹_FAME, t²_FAME, chrom_FAME.^(1//1), legend=false, xlabel="tR1 in s", ylabel="tR2 in s", c=:jet1, dpi=500)#, xlims=(0.0, 1700.0), ylims=(0.0, 5.0), grid=false)
	#gui()
	Plots.scatter!(p_GCxGC_FAME, pl_GCxGC_FAME_UFZ.tR1, pl_GCxGC_FAME_UFZ.tR2.-pl_GCxGC_FAME_UFZ.tR1)
	Plots.scatter!(p_GCxGC_FAME, pl_GCxGC_FAME_Poole.tR1, pl_GCxGC_FAME_Poole.tR2.-pl_GCxGC_FAME_Poole.tR1)
	p_GCxGC_FAME
end

# ╔═╡ dfbcd7c8-ac62-4c6c-96ee-3220345f5577
t¹_Alcohol, t²_Alcohol, chrom_Alcohol = plot_GCxGC_contour_(pl_GCxGC_Alcohol_RetentionData; mod_τ1=1.5)

# ╔═╡ 668080ca-6abc-4920-bfaa-babfb3f15ed9
begin
	plotly()
	p_GCxGC_Alcohol = Plots.heatmap(t¹_Alcohol, t²_Alcohol, chrom_Alcohol.^(1//1), legend=false, xlabel="tR1 in s", ylabel="tR2 in s", c=:jet1, dpi=500)#, xlims=(0.0, 1700.0), ylims=(0.0, 5.0), grid=false)
	#gui()
	Plots.scatter!(p_GCxGC_Alcohol, pl_GCxGC_Alcohol_UFZ.tR1, pl_GCxGC_Alcohol_UFZ.tR2.-pl_GCxGC_Alcohol_UFZ.tR1)
	Plots.scatter!(p_GCxGC_Alcohol, pl_GCxGC_Alcohol_Poole.tR1, pl_GCxGC_Alcohol_Poole.tR2.-pl_GCxGC_Alcohol_Poole.tR1)
	p_GCxGC_Alcohol
end

# ╔═╡ 4ef5ae52-688c-4632-b947-67197b70d75f
t¹_RetentionData, t²_RetentionData, chrom_RetentionData = plot_GCxGC_contour_(pl_RetentionData; mod_τ1=3.0)

# ╔═╡ 2f477231-8cc0-4df8-a8b6-cd672c0fcfe5
begin
	#plotly()
	gr()
	shape = :cross
	size = 4
	p_GCxGC_RetentionData = Plots.heatmap([0.0, 2099.0, 2100.0], [0.0, 0.0, 15.0], [0 0 0;0 0 0; 0 0 0], legend=false, c=:jet1)
	Plots.heatmap!(p_GCxGC_RetentionData, t¹_RetentionData, t²_RetentionData, chrom_RetentionData.^(1//1), xlabel="tR1 in s", ylabel="tR2 in s", c=:jet1, dpi=500, xlims=(0.0, 2000.0), ylims=(0.0, 15.0), grid=false)
	#gui()
	Plots.scatter!(p_GCxGC_RetentionData, pl_GCxGC_PCB_UFZ.tR1, pl_GCxGC_PCB_UFZ.tR2.-pl_GCxGC_PCB_UFZ.tR1, label="PCB", palette=:Set1_6, markershape=shape, msize=size)#:Dark2_6)
	Plots.scatter!(p_GCxGC_RetentionData, pl_GCxGC_PAH_UFZ.tR1, pl_GCxGC_PAH_UFZ.tR2.-pl_GCxGC_PAH_UFZ.tR1, label="PAH", markershape=shape, msize=size, c=:maroon1)
	Plots.scatter!(p_GCxGC_RetentionData, pl_GCxGC_Phenone_UFZ.tR1, pl_GCxGC_Phenone_UFZ.tR2.-pl_GCxGC_Phenone_UFZ.tR1, label="Phenone", markershape=shape, msize=size)
	Plots.scatter!(p_GCxGC_RetentionData, pl_GCxGC_Alkane_UFZ.tR1, pl_GCxGC_Alkane_UFZ.tR2.-pl_GCxGC_Alkane_UFZ.tR1, label="Alkane", markershape=shape, msize=size)
	Plots.scatter!(p_GCxGC_RetentionData, pl_GCxGC_FAME_UFZ.tR1, pl_GCxGC_FAME_UFZ.tR2.-pl_GCxGC_FAME_UFZ.tR1, label="FAME", markershape=shape, msize=size)
	Plots.scatter!(p_GCxGC_RetentionData, pl_GCxGC_Alcohol_UFZ.tR1, pl_GCxGC_Alcohol_UFZ.tR2.-pl_GCxGC_Alcohol_UFZ.tR1, label="Alcohol", markershape=shape, msize=size)

	p_GCxGC_RetentionData
end

# ╔═╡ ba8976fd-efb9-4eac-938e-f29e956231a0
size(solute_groups)

# ╔═╡ Cell order:
# ╠═be16a106-f3c0-11ed-2fa4-7daccc96b7c2
# ╠═2460130a-1b09-4dde-99a6-2bc609f6df28
# ╠═361ad6f5-f412-45a2-8619-89a0dbcecff5
# ╠═50d291e9-a5d8-4897-a4be-6b57b0670bdf
# ╠═1c3291a6-d0e7-4f68-808b-e39b16157a53
# ╠═7bc5ff7b-a6aa-4898-b8e2-70375ce16670
# ╠═bd9279f2-4b24-49d2-b7ae-edebdfedcca6
# ╠═dd5d4628-ab5f-4096-8799-73c2ff0c8fbf
# ╠═c136c0f9-e243-4480-8314-12e58eb5c96f
# ╠═a55e113e-d93a-48ba-a50a-bbb6f137d390
# ╠═8f4a9352-7e59-48ca-b35a-fd40237b8083
# ╠═23b419b8-43c1-4ecd-8e06-9d5202e40d18
# ╠═fa180992-b2d6-4603-90bf-0eff6376f5dd
# ╠═3ef47737-51b0-45fc-9147-dc06dcb68ce1
# ╠═db1ad723-7627-4831-9b47-ab380a74f239
# ╠═419bcc59-74b6-4f91-bdbd-e680104d0fed
# ╠═4f9d42f3-e617-4ded-a884-b8c45f03e376
# ╠═311690d7-556a-46c2-9e40-8676f1e73b2f
# ╠═0022c392-5a5a-400e-9da5-5f6d0fd70549
# ╠═549a03ec-1d89-460c-acff-2812b160c377
# ╠═adba7751-d82b-4136-acad-e34dafcc41e6
# ╠═c9a523b7-ded8-45b5-bee9-1c08d03fd521
# ╠═447c50d6-527d-44cc-a49d-72ca6428d01e
# ╠═26e02657-9cf3-4244-8dd1-b157245bd900
# ╠═c5aebbf9-5fa6-4871-a811-9e041bc9ebb7
# ╠═a018f1ef-b49b-4bfb-9391-c38be7bbfb5e
# ╠═19bffce9-397a-4a07-ac53-9ddb6500b0c5
# ╠═455fb9bc-04c6-4fbf-97e0-a47eb27cf8f1
# ╠═d5463097-925d-424d-b791-cf4987561422
# ╠═ba8976fd-efb9-4eac-938e-f29e956231a0
# ╠═14d95732-117e-4e47-a337-b7481d2b7783
# ╠═4c37a861-475b-4940-9eda-6d345fa2c082
# ╠═a50674d3-9338-474f-97e4-b8498a550933
# ╠═5c97c65f-bfb1-44ab-9ddb-55020f14784a
# ╠═9126d196-d1a5-423e-bfe9-9d3df112e7be
# ╠═6c33c1bb-e6df-47e9-a690-f22f9ecd0a2f
# ╠═cf2794a8-9e39-4035-9651-db610b2ae654
# ╠═6ba883bb-4681-4358-bbd6-73d2abc77f32
# ╠═fcba3fb7-1778-47c3-89a9-b68fe8078f81
# ╠═0e71379e-b682-4fb6-80e9-33c79df879e7
# ╠═32e274da-e0ec-4112-ac5f-5d952f079b25
# ╠═7938448f-9dc3-4019-8fdc-06440d525b06
# ╠═6bad5fde-5efd-4aa2-a78a-ec6c32002aa8
# ╠═4d976b8f-55d9-4e37-9b42-58d130d5c0d3
# ╠═94cf7ea2-eb1c-46bc-9625-6fcaccc7f771
# ╠═f31a0c40-908e-47a9-a3dd-fac6be8904bc
# ╠═8d9ab99d-091b-44fc-905f-3d21f8cb3f00
# ╠═8f7672ac-c8e0-4583-853e-1b4bc9471e47
# ╠═8f4d9af7-0c94-4246-8509-de2b8f3c03e3
# ╠═9e3ab0b6-e0b1-4732-a34e-1900635d45e5
# ╠═c865cde9-9c7b-49bd-bb85-89677e002a2b
# ╠═7b3f5215-5a5c-43a4-8519-6cf7ab2941ae
# ╠═03b3633f-a1d6-4951-83b5-6c2332183e2c
# ╠═30c01f19-1892-4e7a-9fe2-7b66cfbc37a8
# ╠═58237852-30da-477c-a10f-155407f47cd2
# ╠═d7d78054-d1c7-4d88-971d-cd809aae18dc
# ╠═351cb0fe-fcaa-46fd-a6e3-caef23e75e47
# ╠═e8a09fee-1f13-4eae-b070-4225945db950
# ╠═acc3dbc8-f66a-4773-a83d-9a66beaf02fe
# ╠═a12957c6-c0ec-46b6-aae4-1f992eb25037
# ╠═2d817ece-a5c7-4f72-a949-2248fd70e9f3
# ╠═53bf9e63-fbdd-480a-9a3f-50c1dea25871
# ╠═38dccb75-80eb-4f1d-a859-608da5e1fbac
# ╠═02b49171-0068-4a5d-8002-c2b2fe266bf8
# ╠═72decd94-304d-474b-b793-71010be7d2c6
# ╠═8d0d0823-44e2-4f65-8c81-f039ecfc4b0d
# ╠═4fef2f8d-e91f-4e30-b183-0196aad9ec17
# ╠═618cea77-6724-48fe-a1d9-fa2be2428ba2
# ╠═bd3824cd-5dc9-42e8-aaef-6b0750fa0898
# ╠═34c4a88e-78d2-4232-9f35-51139df8b2d3
# ╠═1d1d772a-60d6-4951-a0c0-7e0045d2e107
# ╠═686fc196-3a39-4e98-8607-edbb1445e930
# ╠═544a02a1-25db-47f1-a6e1-940e11457d2b
# ╠═761dc3e3-36f8-4c36-b7d8-ffb17d9a16bb
# ╠═b77d4a3b-9933-453c-b318-01f24658ee58
# ╠═25262675-b410-4f1e-9937-6cd18886bdca
# ╠═ad448ab0-3789-479e-8df8-4cda0d40d8a3
# ╠═75548590-01e8-4d13-a03a-bb2352ce8a88
# ╠═e3d31326-4c69-412e-b9e2-6d7a0c685b46
# ╠═7620a6a0-71c4-484a-8d77-925a2bd1b0ae
# ╠═6be685cf-22fc-4723-ac78-c87c5c283e1f
# ╠═bb174895-d065-4d6f-a7f0-d88c1f06ee44
# ╠═1d771f43-9e96-473f-ab03-ebe646f60aa1
# ╠═e0846d06-6d97-4ece-8b11-e4b4476892bb
# ╠═2b214b26-06b0-41b8-b72b-0ef7779309f2
# ╠═6fb1f245-3285-4662-ac0d-57f9384dac6d
# ╠═6a678c63-4aa9-43b4-bdfd-5119dc14dea1
# ╠═ccceaebe-1438-4484-94bb-bf6c0c182a70
# ╠═d6966984-023f-4a73-a778-c0c66349fd69
# ╠═a7e09a23-d501-41bc-969a-ae48fc64d28d
# ╠═4c404315-70c1-42cf-a944-043e5f03c12a
# ╠═a2e7c37b-39d1-40f8-afb9-9064d790e655
# ╠═88381961-abf2-4fdc-ae14-ea6e60f9e292
# ╠═9a4a6e9c-9dcf-4f66-b947-d42cd1f861e9
# ╠═d8114bb8-ae36-42d3-82ee-8178ac8f09ba
# ╠═2669af8b-0449-4753-9fe9-7de05e2bc13e
# ╠═bf44b9c6-fdd2-45ee-9ea6-5cf6af03a88e
# ╠═dfbcd7c8-ac62-4c6c-96ee-3220345f5577
# ╠═668080ca-6abc-4920-bfaa-babfb3f15ed9
# ╠═4afeb72c-bf72-487c-a1af-f51752316a50
# ╠═6d44e4a0-745e-486e-a8d7-750150d57687
# ╠═4ef5ae52-688c-4632-b947-67197b70d75f
# ╠═2f477231-8cc0-4df8-a8b6-cd672c0fcfe5
# ╠═b2c43093-7c09-4bf3-9211-20f73c65006c
# ╠═f72c5e04-4bff-40ab-a430-8158f78a3f16
# ╠═913e485f-762f-495d-a35e-633e63c22707
# ╠═bb9db28c-4bcc-4b0a-aac2-422531d7b203
# ╠═c07421e4-7d14-4e83-a0d3-08890d8337fa
# ╠═8dca230b-6184-4766-bcb8-5fc13659e91a
