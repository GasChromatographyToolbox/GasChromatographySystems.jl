### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 1cce655a-1ba6-11ec-1bf4-edfded0254c8
begin
	import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
	#using DataFrames, CSV, Interpolations, Plots, QuadGK, DifferentialEquations, ForwardDiff, Intervals
	using Plots
	using GasChromatographySystems, GasChromatographyTools
	using HypertextLiteral
	using PlutoUI
	TableOfContents()
end

# ╔═╡ ed577478-23b5-4c76-95ef-8de32307225e
begin
	# display options
	
	#gr()
	plotly()
	
	html"""
	<style>
	  main {
		max-width: 900px;
	  }
	</style>
	"""
end

# ╔═╡ f5075b30-4b3f-46b8-bb9a-6c8c01d6af9f
md"""
# Simulation of a GC-system with 3 segments
"""

# ╔═╡ f82f4dfc-4bb5-47a0-9a7a-f4faf86360f7
md"""
## Definition of Parameters
"""

# ╔═╡ c5ec5bcb-4bcb-49c7-8ec3-8f4d1317fcbe
function UI_Options()
	PlutoUI.combine() do Child
		@htl("""
		<h3>Option settings</h3>
		
		mobile phase: $(
			Child(Select(["H₂", "He", "N₂"]; default="He"))
		) ODE algorithm: $(
			Child(Select([OwrenZen3(), OwrenZen4(), OwrenZen5()]; default=OwrenZen5()))
		) abstol: 1e $(
			Child(NumberField(-10:1:-3; default=-8))
		) reltol: 1e $(
			Child(NumberField(-8:1:-2; default=-5))
		) Tcontrol: $(
			Child(Select(["inlet", "outlet"]; default="inlet"))
		) ODE-system: $(
			Child(CheckBox(default=true))
		)
		""")
	end
end

# ╔═╡ 6d694987-d98f-4087-b0f3-515d4332f9af
@bind opt_values confirm(UI_Options())

# ╔═╡ 75d2dcdc-25fc-404e-aa62-2b75170a7372
Option = GasChromatographySystems.Options(opt_values[1], opt_values[2], 10.0^opt_values[3], 10.0.^opt_values[4], opt_values[5], opt_values[6])

# ╔═╡ 531977ba-c3d7-49ed-ab64-05b04c8895a2
md"""
### GC-system

$(LocalResource(string(pwd(),"/pic/PressurePoint.png")))
$(LocalResource(string(pwd(),"/pic/Transferline.png")))
$(LocalResource(string(pwd(),"/pic/Column.png")))
$(LocalResource(string(pwd(),"/pic/PressurePoint.png")))
$(LocalResource(string(pwd(),"/pic/Transferline.png")))
$(LocalResource(string(pwd(),"/pic/PressurePoint.png")))
"""

# ╔═╡ c45726d2-2920-4ac8-839f-ff876e8183ac
function UI_System()
	PlutoUI.combine() do Child
		@htl("""
		<h3>System settings</h3>
		<table>
			<tr>
				<th>Parameter</th>
				<th>Column 1</th>
				<th>Column 2</th>
				<th>Column 3</th>
			</tr>
			<tr>
				<td>stat. phase</td>
				<td><center>$(Child(Select(["", "Rxi17SilMS", "SLB5ms", "SPB50", "Wax", "FS5ms", "DB5ms", "Rxi5MS"]; default="DB5ms")))</center></td>
				<td><center>$(Child(Select(["", "Rxi17SilMS", "SLB5ms", "SPB50", "Wax", "FS5ms", "DB5ms", "Rxi5MS"]; default="DB5ms")))</center></td>
				<td><center>$(Child(Select(["", "Rxi17SilMS", "SLB5ms", "SPB50", "Wax", "FS5ms", "DB5ms", "Rxi5MS"]; default="DB5ms")))</center></td>
			</tr>
			<tr>
				<td>L [m]</td>
				<td><center>$(Child(NumberField(0.10:0.01:10.00; default=0.15)))</center></td>
				<td><center>$(Child(NumberField(0.10:0.01:10.00; default=2.0)))</center></td>
				<td><center>$(Child(NumberField(0.10:0.01:10.00; default=2.0)))</center></td>
			</tr>
			<tr>
				<td>d [mm]</td>
				<td><center>$(Child(NumberField(0.01:0.01:0.50; default=0.10)))</center></td>
				<td><center>$(Child(NumberField(0.01:0.01:0.50; default=0.10)))</center></td>
				<td><center>$(Child(NumberField(0.01:0.01:0.50; default=0.10)))</center></td>
			</tr>
			<tr>
				<td>df [µm]</td>
				<td><center>$(Child(NumberField(0.00:0.01:0.50; default=0.10)))</center></td>
				<td><center>$(Child(NumberField(0.00:0.01:0.50; default=0.10)))</center></td>
				<td><center>$(Child(NumberField(0.00:0.01:0.50; default=0.10)))</center></td>
			</tr>
			<tr>
				<td>T [°C]</td>
				<td><center>$(Child(NumberField(30.0:1.0:400.00; default=300.0)))</center></td>
				<td><center>Temp. Prog.</center></td>
				<td><center>$(Child(NumberField(30.0:1.0:400.00; default=300.0)))</center></td>
			</tr>
		</table>
		""")
	end
end

# ╔═╡ 10ad0aea-a5c1-4422-a632-e95a9cb1f7ef
@bind sys_values confirm(UI_System())

# ╔═╡ 6ae3059f-6cb3-417a-b6ab-78a3e44f74bd
function UI_Program()
	PlutoUI.combine() do Child
		@htl("""
		<h3>Temperature program</h3>
		<table>
			<tr>
				<th>Parameter</th>
				<th>Value</th>
			</tr>
			<tr>
				<td>time steps [s]</td>
				<td><center>$(Child(TextField((40,1); default="0 10 60 120")))</center></td>		
			</tr>
			<tr>
				<td>temperature steps [°C]</td>
				<td><center>$(Child(TextField((40,1); default="40 40 350 350")))</center></td>
			</tr>
			<tr>
				<td>ΔT steps [°C]</td>
				<td><center>$(Child(TextField((40,1); default="40 40 50 50")))</center></td>
			</tr>
			<tr>
				<td>α steps [-]</td>
				<td><center>$(Child(TextField((40,1); default="-5 -5 -5 -5")))</center></td>
			</tr>
		</table>
		<h3>Pressure Program</h3>
		<table>
			<tr>
				<th>Parameter</th>
				<th>Value</th>
			</tr>
			<tr>
				<td>time steps [s]</td>
				<td><center>$(Child(TextField((40,1); default="0 10 60 120")))</center></td>
			</tr>
			<tr>
				<td>inlet pressure [kPa(g)]</td>
				<td><center>$(Child(TextField((40,1); default="400 400 500 500")))</center></td>
			</tr>
			<tr>
				<td>connector pressure [kPa(g)]</td>
				<td><center>$(Child(TextField((40,1); default="300 300 300 300")))</center></td>
			</tr>
			<tr>
				<td>outlet pressure [kPa(g)]</td>
				<td>$(Child(Select(["atmospheric", "vacuum"]; default="atmospheric")))</td>
			</tr>
		</table>
		""")
	end
end

# ╔═╡ c1b31764-00e3-45fc-91e9-eeeaa2d09f75
@bind prog_values confirm(UI_Program())

# ╔═╡ e42bf4b4-11c5-4db8-b018-52d355322e0d
begin
	L2 = sys_values[5]
	tsteps = parse.(Float64,split(prog_values[1]))
	Tsteps = parse.(Float64,split(prog_values[2]))
	ΔT = parse.(Float64,split(prog_values[3]))
	x₀ = zeros(length(ΔT))
	L₀ = L2.*ones(length(ΔT))
	α = parse.(Float64,split(prog_values[4]))
	a = [ΔT x₀ L₀ α]
	grad_func(x) = GasChromatographySimulator.gradient(x, a; Tcontrol=Option.Tcontrol)

	TP = GasChromatographySystems.Temperature_Program(tsteps, Tsteps, grad_func, a)
end

# ╔═╡ f041e68e-d629-4f36-902c-7b77be067b1f
begin
	TL1 = GasChromatographySystems.Transferline(sys_values[4], sys_values[7]*1e-3, sys_values[10]*1e-6, sys_values[1], sys_values[13])
	a_d_2 = [sys_values[8]*1e-3]
	d_2(x) = GasChromatographySimulator.gradient(x, a_d_2)
	a_df_2 = [sys_values[11]*1e-6]
	df_2(x) = GasChromatographySimulator.gradient(x, a_df_2)
	GC  = GasChromatographySystems.Column(sys_values[5], d_2, a_d_2, df_2, a_df_2, sys_values[2], TP)
	TL2 = GasChromatographySystems.Transferline(sys_values[6], sys_values[9]*1e-3, sys_values[12]*1e-6, sys_values[3], sys_values[14])
end

# ╔═╡ 87ce0821-e017-4600-8e3e-b24759482a15
begin
	p_tsteps = parse.(Float64,split(prog_values[5]))
	pin_steps = parse.(Float64,split(prog_values[6]))
	pcon_steps = parse.(Float64,split(prog_values[7]))
	if prog_values[8] == "atmospheric"
		pout_steps = 101300.0.*ones(length(p_tsteps))
	elseif prog_values[8] == "vacuum"
		pout_steps = zeros(length(p_tsteps))
	end
	PPin = GasChromatographySystems.Pressure_Point(p_tsteps, pin_steps.*1000.0.+101300.0)
	PPcon = GasChromatographySystems.Pressure_Point(p_tsteps, pcon_steps.*1000.0.+101300.0)
	PPout = GasChromatographySystems.Pressure_Point(p_tsteps, pout_steps)
end

# ╔═╡ 08142bda-450f-42fe-aece-6f0bd5e28717
GCsys = PPin, TL1, GC, PPcon, TL2, PPout

# ╔═╡ ed14d2ca-9a09-41fa-a67a-8190161646b6
begin
	# Definition of functions
	function p_str(pvalue, n)
		# constructes a string with the 'pvalue' repeated 'n' times separated by " "
		str = ""
		for i=1:n-1
			str = string(str,pvalue," ")
		end
		str = string(str,pvalue)
		return str
	end
	
	function alkane_selection(com_solutes)
		alkanes = String[]
		regex = r"C[1-9]+"
		for i=1:length(com_solutes)
			if occursin(regex, com_solutes[i])
				push!(alkanes, com_solutes[i])
			end
		end
		return alkanes
	end
	
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
	
	function plot_chromatogram(peaklist, Δtsteps; labelpeaks=true, neg=false)
		if peaklist.tR[end]<sum(Δtsteps)
			tend = sum(Δtsteps)
		else
			tend = 1.1*peaklist.tR[end]
		end
		if neg==true
			sgn = -1.0
		else
			sgn = 1.0
		end
		gauss(t, tR, τR) = 1/sqrt(2*π*τR^2)*exp(-(t-tR)^2/(2*τR^2))
		gauss_sum(t) = sum(gauss.(t, peaklist.tR, peaklist.τR))
		time = 0:tend/10000:tend
		chrom = gauss_sum.(time)
		pc = plot(time, sgn.*chrom,
					title="Chromatogram",
					legend=false,
					grid=false,
					xlims=(0,tend),
					size=(800,500))
		# peak annotations
		if labelpeaks==true
			peaklabel = Array{Plots.PlotText}(undef, size(peaklist)[1])
			for i=1:size(peaklist)[1]
				peaklabel[i] = Plots.text(string(i), 8)
			end
			scatter!(peaklist.tR, sgn.*gauss_sum.(peaklist.tR).*1.01, 
						series_annotations=peaklabel,
						markersize=0,
						m=:vline)
		end
		xlabel!("time in s")

		return pc, time, chrom
	end
	
	function option_str(Option)
		fn_o = fieldnames(typeof(Option))
		opt_str_array = String[]
		for i=1:length(fn_o)
			push!(opt_str_array, string(fn_o[i], ": ", getfield(Option, fn_o[i])))
		end
		opt_str = string(join(opt_str_array, ", "), "\n")
		return opt_str
	end
	
	function pl_to_string(df)
		# translates peaklist Dataframe to a string for download as csv file 
		a, b = size(df)
		dfstrarray = Array{String}(undef, a+1, b)
		for i=1:a+1
			for j=1:b
				if i==1 && j<b
					dfstrarray[i,j] = string(names(df)[j],",")
				elseif i==1 && j==b
					dfstrarray[i,j] = string(names(df)[j])
				elseif j==b
					dfstrarray[i,j] = string(df[!, j][i-1])
				else 
					dfstrarray[i,j] = string(df[!, j][i-1], ",")
				end
			end
		end
		str=string(convert(Matrix,dfstrarray))
		str1 = replace(str, ";" => "\n")
		str2 = replace(str1, "[" => "")
		str3 = replace(str2, "]" => "")
		str4 = replace(str3, "\"" => "")
		return str4
	end
	
	function df_to_string(df)
		# translates Dataframe to a string for download as csv file 
		a, b = size(df)
		dfstrarray = Array{String}(undef, a+1, b)
		for i=1:a+1
			for j=1:b
				if i==1 && j<b
					dfstrarray[i,j] = string(names(df)[j],",")
				elseif i==1 && j==b
					dfstrarray[i,j] = string(names(df)[j])
				else
					dfstrarray[i,j] = string(df[!, j][i-1])
				end
			end
		end
		str=string(convert(Matrix,dfstrarray))
		str1 = replace(str, ";" => "\n")
		str2 = replace(str1, "[" => "")
		str3 = replace(str2, "]" => "")
		str4 = replace(str3, "\"" => "")
		return str4
	end
	
	function settings(GCsys)
		# export the settings of a GC-system 'GCsys' into a dataframe
		# includes only the moduls 'Pressure_Point', 'Transferline' and 'Column' (together with 
		# 'Temperature_Program')
		modul_str = String[]
		parameter_str = String[]
		tp_parameter_str = String[]
		values_str = String[]
		for i=1:length(GCsys)
			n = length(fieldnames(typeof(GCsys[i])))
			if typeof(GCsys[i])==Pressure_Point || typeof(GCsys[i])==Transferline
				push!(modul_str, string(typeof(GCsys[i]),","))
				push!(parameter_str, string(fieldnames(typeof(GCsys[i]))[1],","))
				push!(tp_parameter_str, ",")
				push!(values_str, string(getfield(GCsys[i], fieldnames(typeof(GCsys[i]))[1])))
				for j=2:n
					push!(modul_str, ",")
					push!(parameter_str, string(fieldnames(typeof(GCsys[i]))[j],","))
					push!(tp_parameter_str, ",")
					push!(values_str, string(getfield(GCsys[i], fieldnames(typeof(GCsys[i]))[j])))
				end
			elseif typeof(GCsys[i])==Column
				col_modul, col_param, col_tp_param, col_values = column_settings(GCsys[i])
				for j=1:length(col_modul)
					push!(modul_str, col_modul[j])
					push!(parameter_str, col_param[j])
					push!(tp_parameter_str, col_tp_param[j])
					push!(values_str, col_values[j])
				end
			end
		end
		df_settings = DataFrame(Modul=modul_str, Parameter=parameter_str, TP_Parameter=tp_parameter_str, Values=values_str)
		return df_settings
	end
	
	function column_settings(Column_Modul)
		modul_str = String[]
		parameter_str = String[]
		tp_parameter_str = String[]
		values_str = String[]
		fn = fieldnames(typeof(Column_Modul))
		m = length(fn)
		for i=1:m # loop over the fieldnames of Column
			if !isa(getfield(Column_Modul, fn[i]), Function)
				if i==1 # first field of Column (length -> one value)
					push!(modul_str, string(typeof(Column_Modul), ","))
					push!(parameter_str, string(fn[i], ","))
					push!(tp_parameter_str, ",")
					push!(values_str, string(getfield(Column_Modul, fn[i])))
				else
					#push!(modul_str, ",")
					#push!(parameter_str, string(fn[i], ","))
					if fn[i]==:temperature_program
						# value is a structure, additional rows needed
						tp_fn = fieldnames(typeof(Column_Modul.temperature_program))
						n = length(tp_fn)
						for j=1:n # loop over the fieldnames of Temperature_Program
							if !isa(getfield(Column_Modul.temperature_program, tp_fn[j]), Function)
								if j==1
									push!(modul_str, ",")
									push!(parameter_str, ",")
									push!(tp_parameter_str, string(tp_fn[j], ","))
									push!(values_str, string(getfield(Column_Modul.temperature_program, tp_fn[j])))
								else
									if isa(getfield(GCsys[3].temperature_program, tp_fn[j]), Matrix)
										# test for the size of the values, if 2d array, add more rows
										o = size(getfield(Column_Modul.temperature_program, tp_fn[j]))[1]
										for k=1:o
											push!(modul_str, ",")
											push!(parameter_str, ",")
											if k==1
												push!(tp_parameter_str, string(tp_fn[j], ","))
											else
												push!(tp_parameter_str, ",")
											end
											push!(values_str, string(getfield(Column_Modul.temperature_program, tp_fn[j])[k,:]))
										end
									else
										push!(modul_str, ",")
										push!(parameter_str, ",")
										push!(tp_parameter_str, string(tp_fn[j], ","))
										push!(values_str, string(getfield(Column_Modul.temperature_program, tp_fn[j])))
									end									
								end
							else # fieldname is a function, this is not exported into the string
							end
						end

					else
						push!(modul_str, ",")
						push!(parameter_str, string(fn[i], ","))
						push!(tp_parameter_str, ",")
						push!(values_str, string(getfield(Column_Modul, fn[i])))
					end
				end
			else # fieldname is a function, this is not exported into the string
			end
		end
		return modul_str, parameter_str, tp_parameter_str, values_str
	end
end

# ╔═╡ bb5c36f1-5b47-4055-b10c-0b1873c4a5c8
GasChromatographySystems.test_of_GCsys(GCsys)

# ╔═╡ 0cac5ac9-31f7-4bf0-b71b-b8e691e8d5ea
md"""
### Solutes

Selection of solutes: $(@bind selection Select(["all", "Alkanes", "PAH", "manually"]; default="Alkanes"))
"""

# ╔═╡ fd69382e-5e7e-4409-9754-f9177d8f35c4
begin
	db = DataFrame(CSV.File(string(pwd(),"/../data/Database_append.csv")))
	com_solutes = GasChromatographySystems.common_solutes(db, GCsys)
	if selection=="manually"
		md"""
		manually select from the following $(length(com_solutes)) solutes: 

		$(@bind solutes MultiSelect(com_solutes; default=com_solutes[1:2]))
		"""
	elseif selection=="all"
		md"""
		all $(length(com_solutes)) solutes are selected: 

		$(@bind solutes MultiSelect(com_solutes; default=com_solutes[1:end]))
		"""
	elseif selection=="Alkanes"
		solutes = alkane_selection(com_solutes)
		if length(solutes)==0
			md"""
			No data for alkanes on selected stationary phases available.
			"""
		end
	elseif selection=="PAH"
		solutes = PAH_selection(com_solutes)
		if length(solutes)==0
			md"""
			No data for PAHs on selected stationary phases available.
			"""
		end
	end
end

# ╔═╡ 6fef3cba-8917-480a-bdf2-a3ae472c02e3
md"""
selected solutes:
$(embed_display(solutes))
"""

# ╔═╡ 0542086b-0f06-40e5-8ea6-9e788a4ded7f
md"""
## Show Temperature, Flow and Pressure
"""

# ╔═╡ f781ca29-9d94-4a82-9ab8-41046b3453dc
md"""
### Temperature
select temperature plot: $(@bind Tplot Select(["T(x,t)", "T(x)", "T(t)"]; default="T(x,t)"))
"""

# ╔═╡ 6114231c-3f8a-4e87-8a19-629e4f03e900
begin
	L = GasChromatographySystems.length_vector(GCsys)
	if Tplot=="T(x)"
		md"""
		select time: $(@bind tselect Slider(0:sum(tsteps), show_value=true))s
		"""
	elseif Tplot=="T(t)"
		md"""
		select time: $(@bind xselect Slider(0:sum(L)/100:sum(L), show_value=true))m
		"""
	end
end

# ╔═╡ e2b9d8bc-038b-4d09-8558-27bf660c9882
begin
	if Tplot=="T(x,t)"
		md"""
		**_Temperature T(x,t)_**
		
		$(embed_display(GasChromatographySystems.temperature_plot(GCsys, Option, Tplot)[1]))
		"""
	elseif Tplot=="T(x)"
		md"""
		**_Temperature T(x)_**
		
		$(embed_display(GasChromatographySystems.temperature_plot(GCsys, Option, Tplot, t₀=tselect)[1]))
		"""
	elseif Tplot=="T(t)"
		md"""
		**_Temperature T(t)_**
		$(embed_display(GasChromatographySystems.temperature_plot(GCsys, Option, Tplot, x₀=xselect)[1]))
		"""
	end
end

# ╔═╡ b4b3ce39-9d6b-49d2-9b3e-2852d0401c9a
md"""
### Flow

**_Flow F(t)_**

$(embed_display(GasChromatographySystems.flow_plot(GCsys, Option)[1]))
"""

# ╔═╡ adaf8dfb-7c62-476a-9bb9-fb34ff949abf
md"""
### Pressure
select pressure plot: $(@bind pplot Select(["p(x,t)", "p(x)", "p(t)"]; default="p(x,t)"))
"""

# ╔═╡ da3f7878-99d7-41fd-9294-57664cfd1c10
begin
	if pplot=="p(x)"
		md"""
		select time: $(@bind tselectp Slider(0:sum(tsteps), show_value=true))s
		"""
	elseif pplot=="p(t)"
		md"""
		select time: $(@bind xselectp Slider(0:sum(L)/100:sum(L), show_value=true))m
		"""
	end
end

# ╔═╡ a4aadd32-84e8-4d5d-a54d-877da644e2bc
begin
	if pplot=="p(x,t)"
		md"""
		**_Pressure p(x,t)_**
		
		$(embed_display(GasChromatographySystems.pressure_plot(GCsys, Option, pplot)[1]))
		"""
	elseif pplot=="p(x)"
		md"""
		**_Pressure p(x)_**
		
		$(embed_display(GasChromatographySystems.pressure_plot(GCsys, Option, pplot, t₀=tselectp)[1]))
		"""
	elseif pplot=="p(t)"
		md"""
		**_Pressure p(t)_**
		$(embed_display(GasChromatographySystems.pressure_plot(GCsys, Option, pplot, x₀=xselectp)[1]))
		"""
	end
end

# ╔═╡ bbe1f0ae-3675-4b3b-b012-c25304c38689
md"""
## Run the Simulation
Simulate: $(@bind go CheckBox(default=false))
"""

# ╔═╡ f382e2e5-38ae-4558-99ff-74213b488c1b
begin
	#if go==true
		par, pl, sol = GasChromatographySystems.linear_GC_system_simulation(GCsys, Option, solutes, db)
	#end
end;

# ╔═╡ af016139-05bc-4623-b72e-1c4726024c45
md"""
## Results
"""

# ╔═╡ 16ef205f-5129-4bdc-8bb7-84b1282c18ef
begin
	pl1 = pl[1]
	pl2 = pl[2]
	pl3 = pl[3]
	peaklist_combi = DataFrame(Name=pl3.Name, tR_D=pl3.tR, τR_D=pl3.τR, TR_C=pl2.TR, σR_D=pl3.σR, uR_D=pl3.uR, kR_C=pl2.kR, Res_D=pl3.Res)
	md"""
	### Peaklist
	Index 'D' ... values at the end of the GC-system
	
	Index 'C' ... values at the end of the column with the temperature program
	
	$(embed_display(peaklist_combi))
	"""
end

# ╔═╡ 7b0ede51-4a76-4360-b09d-23b002ae6652
begin
	p_chrom1, time1, chrom1 = plot_chromatogram(pl1, TP.timesteps)
	p_chrom2, time2, chrom2 = plot_chromatogram(pl2, TP.timesteps)
	p_chrom3, time3, chrom3 = plot_chromatogram(pl3, TP.timesteps)
	md"""
	### Chromatogram
	
	Chromatogram at the end of the first transferline:
	
	$(embed_display(p_chrom1))
	
	Chromatogram at the end of the column:
	
	$(embed_display(p_chrom2))
	
	Chromatogram at the end of the second transferline:
	
	$(embed_display(p_chrom3))
	"""
end

# ╔═╡ 3f52ec2f-cac1-4c60-b589-f52f4b388273
begin
	p_Res = plot(pl1.Name, pl1.Res, ylabel="resolution", label="TL1")
	plot!(p_Res, pl2.Name, pl2.Res, label="GC1")
	plot!(p_Res, pl3.Name, pl3.Res, label="TL2")
	p_Res
end

# ╔═╡ 59e8e493-4d40-42cd-bc2e-c3dc415aac45
begin
	export_str = string("Options: \n", option_str(Option), "GC-system settings: \n", df_to_string(settings(GCsys)), "\n Peaklist: \n", pl_to_string(peaklist_combi))
	md"""
	## Export Results
	Filename: $(@bind result_filename TextField((20,1); default="Result.txt"))
	"""
end

# ╔═╡ 15ee4fd0-7e13-4e3c-ba54-9e6601279a68
md"""
$(DownloadButton(export_str, result_filename))
"""

# ╔═╡ 60738259-20ea-4c69-be9f-aadddbc05b8f
# add additional graphs, e.g. T(x), T(t), τ(x), σ(t) or u(x) of the solutes

# ╔═╡ 233612d8-ab4d-462e-b825-5106222b1059
md"""
# End
"""

# ╔═╡ Cell order:
# ╠═1cce655a-1ba6-11ec-1bf4-edfded0254c8
# ╟─ed577478-23b5-4c76-95ef-8de32307225e
# ╟─ed14d2ca-9a09-41fa-a67a-8190161646b6
# ╟─f5075b30-4b3f-46b8-bb9a-6c8c01d6af9f
# ╟─f82f4dfc-4bb5-47a0-9a7a-f4faf86360f7
# ╟─6d694987-d98f-4087-b0f3-515d4332f9af
# ╟─75d2dcdc-25fc-404e-aa62-2b75170a7372
# ╟─c5ec5bcb-4bcb-49c7-8ec3-8f4d1317fcbe
# ╟─531977ba-c3d7-49ed-ab64-05b04c8895a2
# ╟─10ad0aea-a5c1-4422-a632-e95a9cb1f7ef
# ╟─f041e68e-d629-4f36-902c-7b77be067b1f
# ╟─c45726d2-2920-4ac8-839f-ff876e8183ac
# ╟─c1b31764-00e3-45fc-91e9-eeeaa2d09f75
# ╟─6ae3059f-6cb3-417a-b6ab-78a3e44f74bd
# ╟─e42bf4b4-11c5-4db8-b018-52d355322e0d
# ╟─87ce0821-e017-4600-8e3e-b24759482a15
# ╟─08142bda-450f-42fe-aece-6f0bd5e28717
# ╟─bb5c36f1-5b47-4055-b10c-0b1873c4a5c8
# ╟─0cac5ac9-31f7-4bf0-b71b-b8e691e8d5ea
# ╟─fd69382e-5e7e-4409-9754-f9177d8f35c4
# ╟─6fef3cba-8917-480a-bdf2-a3ae472c02e3
# ╟─0542086b-0f06-40e5-8ea6-9e788a4ded7f
# ╟─f781ca29-9d94-4a82-9ab8-41046b3453dc
# ╟─6114231c-3f8a-4e87-8a19-629e4f03e900
# ╟─e2b9d8bc-038b-4d09-8558-27bf660c9882
# ╟─b4b3ce39-9d6b-49d2-9b3e-2852d0401c9a
# ╟─adaf8dfb-7c62-476a-9bb9-fb34ff949abf
# ╟─da3f7878-99d7-41fd-9294-57664cfd1c10
# ╟─a4aadd32-84e8-4d5d-a54d-877da644e2bc
# ╟─bbe1f0ae-3675-4b3b-b012-c25304c38689
# ╟─f382e2e5-38ae-4558-99ff-74213b488c1b
# ╟─af016139-05bc-4623-b72e-1c4726024c45
# ╠═16ef205f-5129-4bdc-8bb7-84b1282c18ef
# ╟─7b0ede51-4a76-4360-b09d-23b002ae6652
# ╠═3f52ec2f-cac1-4c60-b589-f52f4b388273
# ╟─59e8e493-4d40-42cd-bc2e-c3dc415aac45
# ╟─15ee4fd0-7e13-4e3c-ba54-9e6601279a68
# ╠═60738259-20ea-4c69-be9f-aadddbc05b8f
# ╟─233612d8-ab4d-462e-b825-5106222b1059
