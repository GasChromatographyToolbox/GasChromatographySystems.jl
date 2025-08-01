### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 9d7383a1-d34f-4b0c-977a-fd53919ce93d
begin
    import Pkg
    # activate a temporary environment
    Pkg.activate(mktempdir())
    Pkg.add([
        Pkg.PackageSpec(name="GasChromatographySystems", version="v0.2.6"),
		Pkg.PackageSpec(name="GasChromatographySimulator", version="v0.5.4"),
		Pkg.PackageSpec(name="HypertextLiteral", version="v0.9.5"),
		Pkg.PackageSpec(name="OrdinaryDiffEq", version="v6.95.1"),
		Pkg.PackageSpec(name="Plots", version="v1.40.11"),
		Pkg.PackageSpec(name="PlutoUI", version="v0.7.61"),
		Pkg.PackageSpec(name="UrlDownload", version="v1.0.1"),
    ])
    using CSV, DataFrames
	using GasChromatographySystems
	using GasChromatographySimulator
	using HypertextLiteral
	using OrdinaryDiffEq
	using Plots
	using PlutoUI
	using UrlDownload

	TableOfContents(depth=4)
end

# ╔═╡ 98217474-a16f-406a-83a7-17fee89c951a
html"""<style>
main {
    max-width: 80%;
    margin-left: 1%;
    margin-right: 10% !important;
}
"""

# ╔═╡ 22091a27-80e1-4e98-abe7-4b9652cb832c
md"""
# GCxGC with thermal modulator

Zoex 2 Setup: One Cold-/Hot-Jet used twice by looping the 2nd D column.
"""

# ╔═╡ 370950e7-7c8e-4503-a771-01227f56874d
function UI_column(;title="Columns", def=(29.5, 0.3, 0.005, 0.9, 0.005, 0.545, 0.235, 0.25, 0.25, 0.25, 0.25), sp=["Rxi5SilMS", "Rxi17SilMS"], sp_def=("Rxi5SilMS", "Rxi17SilMS"))
	PlutoUI.combine() do Child
		@htl("""
		<table>
			<caption> $(title) </capiton>
			<tr>
				<th></th>
				<th>1st D</th>
				<th>mod in</th>
				<th>1st spot</th>
				<th>loop</th>
				<th>2nd spot</th>
				<th>2nd D</th>
				<th>TL</th>
			</tr>
			<tr>
				<td><i>L</i> [m]</td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[1])))</center></td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[2])))</center></td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[3])))</center></td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[4])))</center></td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[5])))</center></td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[6])))</center></td>
				<td><center>$(Child(NumberField(0.0:0.001:100.0; default=def[7])))</center></td>
			</tr>
			<tr>
				<td><i>d</i> [mm]</td>
				<td><center>$(Child(NumberField(0.0:0.001:1.0; default=def[8])))</center></td>
				<td colspan="7"><center> &larrb; $(Child(NumberField(0.0:0.001:1.0; default=def[9]))) &rarrb; </center></td>
			</tr>
			<tr>
				<td><i>d</i><sub>f</sub> [µm]</td>
				<td><center>$(Child(NumberField(0.0:0.001:1.0; default=def[10])))</center></td>
				<td colspan="7"><center> &larrb; $(Child(NumberField(0.0:0.001:1.0; default=def[11]))) &rarrb; </center></td>
			</tr>
			<tr>
				<td>stat. phase</td>
				<td><center>$(Child(Select(sp; default=sp_def[1])))</center></td>
				<td colspan="7"><center> &larrb; $(Child(Select(sp; default=sp_def[2]))) &rarrb; </center></td>
			</tr>
		</table>
		""")
	end
end

# ╔═╡ 43f92dc5-7911-40d9-978c-54f0873736a9
md"""
### Column settings
$(@bind col_values confirm(UI_column(title="", 
										def=(29.74, 0.3, 0.005, 0.9, 0.005, 0.53, 0.24, 0.25, 0.1, 0.25, 0.1), 
										sp=["ZB1ms", "Stabilwax"], 
										sp_def=("ZB1ms", "Stabilwax")
									)))
"""

# ╔═╡ ef6820c9-baba-4185-ac1a-5187cb9aa736
md"""
### Flow and outlet
$(@bind flow_mode confirm(Select(["flow", "pressure"]; default="pressure")))
"""

# ╔═╡ 5ccefb78-a43d-45c0-bf85-b98997e3f0fc
function UI_flows(; title="Flows", def=(0.8, 0.0))
	PlutoUI.combine() do Child
		@htl("""
		<table>
			<caption> Flows </capiton>
			<tr>
				<td> <i>F</i> [mL/min] </td>
				<td><center> $(Child(NumberField(0.0:0.01:10.0; default=def[1]))) </center></td>
			</tr>
			<tr>
				<td> <i>p</i><sub>out</sub> [kPa(g)] </td>
				<td><center> $(Child(NumberField(0.0:0.01:1000.0; default=def[2]))) </center></td>
			</tr>
		</table>
		""")
	end
end

# ╔═╡ b1524be6-14f7-4f33-9fed-390d755b1f99
begin
	if flow_mode == "flow"
		md"""
		$(@bind flow_values confirm(UI_flows(title="")))
		"""
	else
		@htl("""<i>p</i><sub>out</sub> [kPa(g)]: $(@bind flow_values confirm(NumberField(0.0:0.01:1000.0; default=0.0)))""")
	end
end

# ╔═╡ 2514f7a9-2e58-4716-8b61-3a3c4e35360d
md"""
### Temperature Program
Number of ramps: $(@bind n_ramp confirm(NumberField(0:1:100; default=1)))
"""

# ╔═╡ 7aebd4d7-c190-4ea8-a09d-b33c19fbb990
function UI_TP(n_ramp; def=[[40.0, 2.0]; reduce(vcat, [[i*5.0, i*300.0/n_ramp, i] for i=1:n_ramp])], def_TTL=300.0)
	PlutoUI.combine() do Child
		@htl("""
		<table>
			<tr>
				<th>ramp [°C/min]</th>
				<th><i>T</i> [°C]</th>
				<th><i>t</i><sub>hold</sub> [min]</th>
			</tr>
			<tr>
				<td></td>
				<td><center>$(Child(NumberField(0.0:0.1:400.0; default=def[1])))</center></td>
				<td><center>$(Child(NumberField(0.0:0.1:100.0; default=def[2])))</center></td>
			</tr>
			$([
				@htl("
					<tr>
						<td><center>$(Child(NumberField(0.0:0.1:100.0; default=def[3+3*(i-1)])))</center></td>
						<td><center>$(Child(NumberField(0.0:0.1:400.0; default=def[4+3*(i-1)])))</center></td>
						<td><center>$(Child(NumberField(0.0:0.1:100.0; default=def[5+3*(i-1)])))</center></td>
					</tr>
				")
				for i=1:n_ramp
			])
		</table>
		<table>
			<tr>
				<td><i>T</i><sub>TL</sub> [°C]: $(Child(NumberField(0.0:0.1:400.0; default=def_TTL)))</td>
			</tr>
		</table>
		""")
	end
end

# ╔═╡ cf17e634-98ab-4fe0-bf9f-4f87d426bbc4
function UI_TP_pressure(n_ramp; def=[[40.0, 2.0, 10.0]; reduce(vcat, [[i*5.0, i*300.0/n_ramp, i, i*20.0] for i=1:n_ramp])], def_TTL=300.0)
	PlutoUI.combine() do Child
		@htl("""
		<table>
			<tr>
				<th>ramp [°C/min]</th>
				<th><i>T</i> [°C]</th>
				<th><i>t</i><sub>hold</sub> [min]</th>
				<th><i>p</i><sub>in</sub> [kPa(g)]</th>
			</tr>
			<tr>
				<td></td>
				<td><center>$(Child(NumberField(0.0:0.1:400.0; default=def[1])))</center></td>
				<td><center>$(Child(NumberField(0.0:0.1:100.0; default=def[2])))</center></td>
				<td><center>$(Child(NumberField(0.0:0.01:1000.0; default=def[3])))</center></td>
			</tr>
			$([
				@htl("
					<tr>
						<td><center>$(Child(NumberField(0.0:0.1:100.0; default=def[4+4*(i-1)])))</center></td>
						<td><center>$(Child(NumberField(0.0:0.1:400.0; default=def[5+4*(i-1)])))</center></td>
						<td><center>$(Child(NumberField(0.0:0.1:100.0; default=def[6+4*(i-1)])))</center></td>
						<td><center>$(Child(NumberField(0.0:0.01:1000.0; default=def[7+4*(i-1)])))</center></td>
					</tr>
				")
				for i=1:n_ramp
			])
		</table>
		<table>
			<tr>
				<td><i>T</i><sub>TL</sub> [°C]: $(Child(NumberField(0.0:0.1:400.0; default=def_TTL)))</td>
			</tr>
		</table>
		""")
	end
end

# ╔═╡ 53c91d41-18d8-4ad2-8dfa-65df326746c4
begin
	if flow_mode == "flow"
		def_TP = [[50.0, 2.0]; reduce(vcat, [[i*3.0, i*225.0/n_ramp, 5.0] for i=1:n_ramp])]
		md"""
		$(@bind TP_values confirm(UI_TP(n_ramp; def=def_TP, def_TTL=225.0)))
		"""
	else
		def_TP = [[50.0, 2.0, 160.3]; reduce(vcat, [[i*3.0, i*225.0/n_ramp, 5.0, 274.96] for i=1:n_ramp])]
		md"""
		$(@bind TP_values confirm(UI_TP_pressure(n_ramp; def=def_TP, def_TTL=250.0)))
		"""
	end
end

# ╔═╡ 62173d6c-7e66-4952-a0e5-919230be481e
function UI_modulator(;title="Modulator", def=(4.0, 3.5, 0.35, 25.0, -130.0))
	PlutoUI.combine() do Child
		@htl("""
		<table>
			<caption> $(title) </capiton>
			<tr>
				<td> <i>t</i><sub>PM</sub> [s] </td>
				<td><center> $(Child(NumberField(0.0:0.000001:10.0; default=def[1]))) </center></td>
			</tr>
			<tr>
				<td> shift [s] </td>
				<td><center> $(Child(NumberField(0.0:0.01:100.0; default=def[2]))) </center></td>
			</tr>
			<tr>
				<td> <i>t</i><sub>hot</sub> [s] </td>
				<td><center> $(Child(NumberField(0.0:0.01:100.0; default=def[3]))) </center></td>
			</tr>
			<tr>
				<td> <i>T</i><sub>hot</sub> [°C] </td>
				<td><center> $(Child(NumberField(-100.0:0.1:200.0; default=def[4]))) </center></td>
			</tr>
			<tr>
				<td> <i>T</i><sub>cold</sub> [°C] </td>
				<td><center> $(Child(NumberField(-200.0:0.1:200.0; default=def[5]))) </center></td>
			</tr>
		</table>
		""")
	end
end

# ╔═╡ 54af6817-92bd-4094-9c12-1515c44b0f5d
md"""
### Modulator settings

$(@bind mod_values confirm(UI_modulator(title="", def=(3.000236, 0.13, 0.35, 25.0, -80.0))))
"""

# ╔═╡ 91c8369f-b5cb-4f06-bcb5-740859c63718
md"""
### System
"""

# ╔═╡ a2b5fb96-9738-4a2e-8c18-fc830fef295e
function setup_GCxGC_TM(col_values, flow_mode, flow_values, n_ramp, TP_values, mod_values; opt_values=("simplifiedTM", "absolut"))
 
	Tcold_abs = if opt_values[2] == "relativ"
		false # fixed relativ temperature difference at modulation point related to temperature program
	else
		true
	end
	optTM = if opt_values[1] == "uniformTM"
		GasChromatographySystems.ModuleTMOptions(Tcold_abs=Tcold_abs, ng=true, alg=OrdinaryDiffEq.Vern9(), tflank=20, sflank=Inf, dtinit=0.5e-7, abstol=1e-12, reltol=1e-10)
	else
		GasChromatographySystems.ModuleTMOptions(Tcold_abs=Tcold_abs, ng=true, alg="simplifiedTM", tflank=Inf, sflank=Inf, dtinit=0.5e-7, abstol=1e-14, reltol=1e-12)
	end
	

	# column settings
	L1 = col_values[1]
	d1 = col_values[8] # 9
	df1 = col_values[10] # 11
	sp1 = col_values[12] # 13
	L2 = col_values[6] # 7
	d2 = col_values[9] # 10
	df2 = col_values[11] # 12
	sp2 = col_values[13] # 14
	LTL = col_values[7] # 8
	dTL = d2
	dfTL = df2
	spTL = sp2
	LM = [col_values[2], col_values[3], col_values[4], col_values[5]]
	dM = d2
	dfM = df2
	spM = sp2
	# flow, pressure and temperature program
	if flow_mode == "pressure"
		time_steps, temp_steps = GasChromatographySimulator.conventional_program([TP_values...][Not([[3+i*4 for i=0:n_ramp]; length(TP_values)])])
		TP = GasChromatographySystems.TemperatureProgram(time_steps, temp_steps) #length of TP depends on n_ramp
		PP = GasChromatographySystems.PressureProgram(time_steps, [TP_values[reduce(vcat,[[3+i*4, 3+i*4] for i=0:n_ramp])].*1000.0.+101300.0...]) # also depends on n_ramp
		pout = flow_values
		F = NaN
	elseif flow_mode == "flow"
		time_steps, temp_steps = GasChromatographySimulator.conventional_program([TP_values[i] for i=1:(length(TP_values)-1)])
		TP = GasChromatographySystems.TemperatureProgram(time_steps, temp_steps.*TP_factor) #length of TP depends on n_ramp
		PP = NaN
		pout = flow_values[2]
		F = flow_values[1]
	end
	# temperature transfer line
	TPTL = TP_values[end]
	# modulator
	PM = mod_values[1]
	shift = mod_values[2]
	ratioM = (mod_values[1]-mod_values[3])/mod_values[1]
	HotM = mod_values[4]
	ColdM = mod_values[5]
	
	sys = GasChromatographySystems.GCxGC_TM(L1, d1, df1, sp1, TP, L2, d2, df2, sp2, TP, LTL, dTL, dfTL, spTL, TPTL, LM, dM, dfM, spM, shift, PM, ratioM, HotM, ColdM, TP, F, PP, pout; name="GCxGC_TM", opt=GasChromatographySystems.Options(), optTM=optTM, optCol=GasChromatographySystems.ModuleColumnOptions(ng=true))
	return sys
end

# ╔═╡ 863fad4d-d3cf-4d37-aa67-11ddcf3c5827
sys = setup_GCxGC_TM(col_values, flow_mode, flow_values, n_ramp, TP_values, mod_values)

# ╔═╡ 30b954f1-3ca3-433e-841f-4c7a1e8f3191
p2fun = GasChromatographySystems.build_pressure_squared_functions(sys, GasChromatographySystems.solve_balance(sys))

# ╔═╡ b1ac781a-6f39-4da4-86de-cf2fb22d076b
md"""
### Substances
"""

# ╔═╡ bddfa563-b0c1-40fb-a0bf-e2fc991f3112
begin #switch later to url
	db = DataFrame(CSV.File("/Users/janleppert/Documents/GitHub/GCxGC_TM/data/exp_pro/GCxGC_simulation/Messungen für RetentionParameterEstimator für GCxGC/Database_m1.csv")) 
	insertcols!(db, 1, :No => collect(1:length(db.Name)))
	db
end

# ╔═╡ a720dc6c-0e35-4e35-b7d3-f9c231457eaa
begin
	selected_solutes = GasChromatographySystems.common_solutes(db, sys).Name
	md"""
	#### Select Substances
	$(@bind solute_values confirm(MultiSelect(selected_solutes; default=selected_solutes[1:end])))
	"""
end

# ╔═╡ 0a1d4088-a5c4-4984-b6bb-23e9a018caa8
md"""
## Simulation
"""

# ╔═╡ 9397d018-8e3c-489b-9250-c3e74bc4681a
par = GasChromatographySystems.graph_to_parameters(sys, p2fun, db, solute_values) # offline version needed -> GCSim.jl - line 693 

# ╔═╡ f7fbf4e5-31c8-4abb-a799-9f460e6187e2
begin
	plotly()
	# temperature program
	p_TP = Plots.plot(cumsum(par[1].prog.time_steps), par[1].prog.T_itp(0.0, cumsum(par[1].prog.time_steps)).-273.15, label="", xlabel="time in s", ylabel="temperature in °C", title="Oven temperature")
	# temperature program at modulator
	p_TPmod = Plots.plot(0.0:0.01:(sys.modules[3].PM*1000.0), par[3].prog.T_itp.(sys.modules[3].L/2, 0.0:0.01:(sys.modules[3].PM*1000.0)).-273.115, title="temperature at modulator point", xlabel="time in s", ylabel="temperature in °C", label="")
	Plots.plot!(p_TPmod, 0.0:0.01:(sys.modules[3].PM*1000.0), par[1].prog.T_itp(0.0, 0.0:0.01:(sys.modules[3].PM*1000.0)).-273.15, label="")
	# temperature modulation in time 
	p_Mod = Plots.plot(0.0:0.01:sys.modules[3].PM, par[3].prog.T_itp.(sys.modules[3].L/2, 0.0:0.01:sys.modules[3].PM).-273.115, title="temperature at modulator point", xlabel="time in s", ylabel="temperature in °C", label="")
	Plots.plot!(p_Mod, [-sys.modules[3].shift, -sys.modules[3].shift], [par[3].prog.T_itp(sys.modules[3].L/2,-sys.modules[3].shift), par[3].prog.T_itp(sys.modules[3].L/2,(1+sys.modules[3].ratio)/2*sys.modules[3].PM-sys.modules[3].shift)].-273.15, c=:black, label="")
	Plots.plot!(p_Mod, [sys.modules[3].ratio*sys.modules[3].PM-sys.modules[3].shift, sys.modules[3].ratio*sys.modules[3].PM-sys.modules[3].shift], [par[3].prog.T_itp(sys.modules[3].L/2,sys.modules[3].ratio*sys.modules[3].PM-sys.modules[3].shift-eps()), par[3].prog.T_itp(sys.modules[3].L/2,(1+sys.modules[3].ratio)/2*sys.modules[3].PM-sys.modules[3].shift)].-273.15, c=:black, linestyle=:dash, label="")
	Plots.plot!(p_Mod, [sys.modules[3].PM-sys.modules[3].shift, sys.modules[3].PM-sys.modules[3].shift], [par[3].prog.T_itp(sys.modules[3].L/2,sys.modules[3].PM-sys.modules[3].shift), par[3].prog.T_itp(sys.modules[3].L/2,(1+sys.modules[3].ratio)/2*sys.modules[3].PM-sys.modules[3].shift)].-273.15, c=:black, label="")
	# temperature modulation in space
	p_Modx = Plots.plot(xlabel="position x in m", ylabel="temperature in °C", title="temperature at modulator point", legend=false)
	Plots.plot!(p_Modx, 0.0:sys.modules[3].L/100.0:sys.modules[3].L, par[3].prog.T_itp.(0.0:sys.modules[3].L/100.0:sys.modules[3].L, sys.modules[3].ratio*sys.modules[3].PM-sys.modules[3].shift).-273.15, label="t=$(sys.modules[3].ratio*sys.modules[3].PM-sys.modules[3].shift)s")
	Plots.plot!(p_Modx, 0.0:sys.modules[3].L/100.0:sys.modules[3].L, par[3].prog.T_itp.(0.0:sys.modules[3].L/100.0:sys.modules[3].L, (1+3*sys.modules[3].ratio)/4*sys.modules[3].PM-sys.modules[3].shift).-273.15, label="t=$((1+3*sys.modules[3].ratio)/4*sys.modules[3].PM-sys.modules[3].shift)s")
	Plots.plot!(p_Modx, 0.0:sys.modules[3].L/100.0:sys.modules[3].L, par[3].prog.T_itp.(0.0:sys.modules[3].L/100.0:sys.modules[3].L, (1+sys.modules[3].ratio)/2*sys.modules[3].PM-sys.modules[3].shift).-273.15, label="t=$((1+sys.modules[3].ratio)/2*sys.modules[3].PM-sys.modules[3].shift)s")
	Plots.plot!(p_Modx, 0.0:sys.modules[3].L/100.0:sys.modules[3].L, par[3].prog.T_itp.(0.0:sys.modules[3].L/100.0:sys.modules[3].L, 3.95-sys.modules[3].shift).-273.15, label="t=$(3.95-sys.modules[3].shift)s")
	Plots.plot!(p_Modx, 0.0:sys.modules[3].L/100.0:sys.modules[3].L, par[3].prog.T_itp.(0.0:sys.modules[3].L/100.0:sys.modules[3].L, sys.modules[3].PM-sys.modules[3].shift).-273.15, label="t=$(sys.modules[3].PM-sys.modules[3].shift)s")
	# flow
	p_F = GasChromatographySystems.plot_flow_over_time(sys, p2fun; dt=10.0)
	Plots.plot(p_F, legend=false)
	# pressure
	p_p = GasChromatographySystems.plot_pressure_over_time(sys, p2fun; dt=30.0)
	Plots.plot(p_p, legend=false)
	
	Plots.plot(p_TP, p_TPmod, p_Mod, p_Modx, p_F, p_p, size=(1000,600), layout = (3, 2), legend=false)
end

# ╔═╡ 7675a0d6-9f32-4882-b77c-6f46fc52d5b2
sim = GasChromatographySystems.simulate_along_paths(sys, p2fun, GasChromatographySystems.all_paths(sys.g, 1)[2], par; nτ=5)

# ╔═╡ b520173c-710a-466a-ac3f-0f87d0268c72
md"""
## Peaklist
"""

# ╔═╡ 8bd4e0a8-8979-4859-b34c-38c0db4cc53b
pl_GCxGC = GasChromatographySystems.peaklist_GCxGC(sim[2][1][end], mod_values[1])

# ╔═╡ ff60df98-8a93-454a-a855-6412e1ddad9e
begin # rework to export sim-results
	column = ["column:", string(col_values)]
	flow = ["flow_mode:", string(flow_mode), "flow:", string(flow_values)]
	temp = ["n_ramp:", string(n_ramp), "TempProg:", string(TP_values)]
	modulator = ["modulator", string(mod_values)]
	io = IOBuffer()
	CSV.write(io, DataFrame[], header=column)
	CSV.write(io, DataFrame[], header=flow, append=true)
	CSV.write(io, DataFrame[], header=temp, append=true)
	CSV.write(io, DataFrame[], header=modulator, append=true)
	CSV.write(io, pl_GCxGC[!,[:Name, :CAS, :tR1, :tR2, :τR1, :τR2]], append=true, writeheader=true)
	#export_str_ = export_str(opt_values, col_values, prog_values, peaklist)
	md"""
	#### Export peak list
	Filename: $(@bind result_filename TextField((20,1); default="Result.csv"))
	"""
end

# ╔═╡ 5b328cc9-a7c0-4fc1-ad74-db56b00ee70e
md"""
$(DownloadButton(io.data, result_filename))
"""

# ╔═╡ d8ec5371-c512-4bef-9778-ebf317159ab1
md"""
## Chromatograms

Elution outside the temperature program (assuming keeping the end temperature constant indefinitly) are shown outside the heatmap.
"""

# ╔═╡ fc7baa4e-e6b8-40f0-94b6-d7d8df5c6052
slice_mat, t_D1, t_D2, c_slices, t_, chrom_sliced_sum  = GasChromatographySystems.chrom2d(sim[2][1][end], sys, mod_values[1]);

# ╔═╡ d9b6ae38-2904-453c-bfbf-f537c2c2f390
begin
	plotly()
	x1 = floor(minimum(pl_GCxGC.tR1*0.99))
	x2 = ceil(maximum(pl_GCxGC.tR1*1.01))
	y1 = floor(minimum(pl_GCxGC.tR2*0.99))
	y2 = ceil(maximum(pl_GCxGC.tR2*1.01))
	p_heatmap = Plots.heatmap(t_D1[1:end-1], t_D2[1], slice_mat', levels=10, xlabel="tR1 in s", ylabel="tR2 in s", c=:jet1)
	Plots.scatter!(p_heatmap, pl_GCxGC.tR1, pl_GCxGC.tR2, msize=3, shape=:diamond, c=:black, markeralpha=1)
	Plots.plot!(p_heatmap, xlims=(x1, x2), ylims=(y1, y2), size=(1000,400))
end

# ╔═╡ 40462520-dad4-4464-be60-07fddcb0b063
md"""
# End
"""

# ╔═╡ Cell order:
# ╟─9d7383a1-d34f-4b0c-977a-fd53919ce93d
# ╠═98217474-a16f-406a-83a7-17fee89c951a
# ╟─22091a27-80e1-4e98-abe7-4b9652cb832c
# ╟─370950e7-7c8e-4503-a771-01227f56874d
# ╟─43f92dc5-7911-40d9-978c-54f0873736a9
# ╟─ef6820c9-baba-4185-ac1a-5187cb9aa736
# ╟─5ccefb78-a43d-45c0-bf85-b98997e3f0fc
# ╟─b1524be6-14f7-4f33-9fed-390d755b1f99
# ╟─2514f7a9-2e58-4716-8b61-3a3c4e35360d
# ╟─53c91d41-18d8-4ad2-8dfa-65df326746c4
# ╟─7aebd4d7-c190-4ea8-a09d-b33c19fbb990
# ╟─cf17e634-98ab-4fe0-bf9f-4f87d426bbc4
# ╟─54af6817-92bd-4094-9c12-1515c44b0f5d
# ╟─62173d6c-7e66-4952-a0e5-919230be481e
# ╟─91c8369f-b5cb-4f06-bcb5-740859c63718
# ╟─a2b5fb96-9738-4a2e-8c18-fc830fef295e
# ╠═863fad4d-d3cf-4d37-aa67-11ddcf3c5827
# ╟─30b954f1-3ca3-433e-841f-4c7a1e8f3191
# ╟─f7fbf4e5-31c8-4abb-a799-9f460e6187e2
# ╟─b1ac781a-6f39-4da4-86de-cf2fb22d076b
# ╠═bddfa563-b0c1-40fb-a0bf-e2fc991f3112
# ╟─a720dc6c-0e35-4e35-b7d3-f9c231457eaa
# ╟─0a1d4088-a5c4-4984-b6bb-23e9a018caa8
# ╟─9397d018-8e3c-489b-9250-c3e74bc4681a
# ╟─7675a0d6-9f32-4882-b77c-6f46fc52d5b2
# ╟─b520173c-710a-466a-ac3f-0f87d0268c72
# ╟─8bd4e0a8-8979-4859-b34c-38c0db4cc53b
# ╟─ff60df98-8a93-454a-a855-6412e1ddad9e
# ╟─5b328cc9-a7c0-4fc1-ad74-db56b00ee70e
# ╟─d8ec5371-c512-4bef-9778-ebf317159ab1
# ╟─fc7baa4e-e6b8-40f0-94b6-d7d8df5c6052
# ╟─d9b6ae38-2904-453c-bfbf-f537c2c2f390
# ╠═40462520-dad4-4464-be60-07fddcb0b063
