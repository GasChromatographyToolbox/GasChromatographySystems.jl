### A Pluto.jl notebook ###
# v0.20.20

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

# ╔═╡ e50d12ba-5e53-11f1-aa44-a152c26b54bb
begin
	import Pkg
	# careful: this is _not_ a reproducible environment
	# activate the shared project environment
	Pkg.activate(Base.current_project())
	# instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()
	using DataFrames
	using GasChromatographySystems
	using GasChromatographySimulator
	using Plots
	using UrlDownload
	using PlutoUI
end

# ╔═╡ 8b8d2ec9-0cc3-4537-ab3a-a1f364af72a8
TableOfContents(depth=4)

# ╔═╡ 8d49ebdd-0b20-41c9-8794-dc605aa61e91
md"""
# Two Columns in Serie
Combining two Columns in serie with different stationary phases and different length ratios to adjust the selectivity. No additional control of the pressure between the two columns.
"""

# ╔═╡ 07bcea72-f32d-4c05-9c98-996cc7efd038
md"""
## Setup
"""

# ╔═╡ c3a5bce8-b532-429f-9dae-7ea3526829a8
begin
	Ltotal = 30.0 # m
    d = 0.25 # mm
    df = 0.25 # µm
    sp1 = "ZB1ms"
	sp2 = "Stabilwax"
    TP = GasChromatographySystems.TemperatureProgram([30.0, 1.0, 20.0, 300.0, 1.0])
	F = 1.0 # mL/min
    pout = 0.0 # Pa
	Lratio_default = 0.5
end

# ╔═╡ 800e9468-ee5e-4ced-afe7-711bdd74946d
sys = GasChromatographySystems.SeriesSystem([Ltotal*Lratio_default, Ltotal*(1-Lratio_default)], [d, d], [df, df], [sp1, sp2], [TP, TP], F, NaN, pout; name="SeriesSystem", opt=GasChromatographySystems.Options(control="Pressure"))

# ╔═╡ 7c1bad0c-1ca8-4d67-9310-55b361276893
sys_reverse = GasChromatographySystems.SeriesSystem([Ltotal*Lratio_default, Ltotal*(1-Lratio_default)], [d, d], [df, df], [sp2, sp1], [TP, TP], F, NaN, pout; name="SeriesSystem", opt=GasChromatographySystems.Options(control="Pressure"))

# ╔═╡ 9aad392d-788a-4e3e-aa22-b61f8132fef6
solution = GasChromatographySystems.solve_balance(sys; mode="λ")

# ╔═╡ b28decc2-8067-4fb6-891a-e084c9ed07bc
p2fun = GasChromatographySystems.build_pressure_squared_functions(sys, solution; mode="λ")

# ╔═╡ 8de3f128-6bda-441e-9f7b-c34700ada555
md"""
## Graph to Parameters
"""

# ╔═╡ a6029cb9-a58d-4ba9-9bdd-d729fece1fc8
begin
	db = DataFrame(urldownload("https://raw.githubusercontent.com/GasChromatographyToolbox/GasChromatographySystems.jl/refs/heads/main/data/Database_GCxGC-TM.csv"))
	insertcols!(db, 1, :No => collect(1:length(db.Name)))
	db
end

# ╔═╡ f0ec2f20-4699-415e-bc84-023ed11cb4dd
selected_solutes = GasChromatographySystems.common_solutes(db, sys).Name

# ╔═╡ f3baf9fd-5719-47f5-8f08-20a883bd8a7e
par_default = GasChromatographySystems.graph_to_parameters(sys, p2fun, db, selected_solutes; interp=true, dt=1, mode="λ")

# ╔═╡ c5e2f032-e6fc-44c5-8f4d-3332ad26bdcb
par_reverse_default = GasChromatographySystems.graph_to_parameters(sys_reverse, p2fun, db, selected_solutes; interp=true, dt=1, mode="λ")

# ╔═╡ d2225282-6399-41d3-b922-de59461495c1
function update_L_in_parameters(par_, Ltotal, Lratio)
	col1 = GasChromatographySimulator.Column(Ltotal*Lratio, par_[1].col.d, par_[1].col.df, par_[1].col.sp, par_[1].col.gas)
	col2 = GasChromatographySimulator.Column(Ltotal*(1-Lratio), par_[2].col.d, par_[2].col.df, par_[2].col.sp, par_[2].col.gas)
	par1 = GasChromatographySimulator.Parameters(col1, par_[1].prog, par_[1].sub, par_[1].opt)
	par2 = GasChromatographySimulator.Parameters(col2, par_[2].prog, par_[2].sub, par_[2].opt)
	par = [par1, par2]
	return par
end

# ╔═╡ 829c7f3a-60c0-4886-b127-cc07933abcf4
md"""
## Simulate
"""

# ╔═╡ 9065ecec-a02b-419b-9e2a-ae2c74b961f6
vertex_paths, edge_paths = GasChromatographySystems.all_paths(sys.g, sys.modules)

# ╔═╡ 35254488-6bf0-4df9-9fe4-97bc4abca2e6
sim_default = GasChromatographySystems.simulate_along_paths(sys, p2fun, edge_paths, par_default)

# ╔═╡ abfcc64c-0c45-4935-a90e-23c104b47376
sim_reverse_default = GasChromatographySystems.simulate_along_paths(sys_reverse, p2fun, edge_paths, par_reverse_default)

# ╔═╡ 5a36a232-7b13-4688-8ef2-7e6fcd14b182
@bind Lratio Slider(0.01:0.01:0.99; default=Lratio_default, show_value=true)

# ╔═╡ 7b775cad-3f3f-4a9f-b1b0-2435cf1edd94
par = update_L_in_parameters(par_default, Ltotal, Lratio)

# ╔═╡ 7a047cef-1b87-44df-bf4e-1a9932df236e
sim = GasChromatographySystems.simulate_along_paths(sys, p2fun, edge_paths, par)

# ╔═╡ 67c87a76-de8a-4556-bbdd-0dc1f7b655aa
par_reverse = update_L_in_parameters(par_reverse_default, Ltotal, Lratio)

# ╔═╡ 549b07b9-6199-459a-9ea8-b9a9e37a6dc8
sim_reverse = GasChromatographySystems.simulate_along_paths(sys_reverse, p2fun, edge_paths, par_reverse)

# ╔═╡ 749c9f0f-fb05-49b1-9b82-b38c5741a190
let
	chrom_default = GasChromatographySimulator.plot_chromatogram(sim_default[2][1][end], (0.0, 1000.0), annotation=false)
	chrom = GasChromatographySimulator.plot_chromatogram(sim[2][1][end], (0.0, 1000.0), annotation=false)
	plot(chrom_default[1], chrom[2], chrom[3], c=:red)
end

# ╔═╡ 9765efe5-6c51-4417-8e1b-2f9231208d1d
let
	chrom_reverse_default = GasChromatographySimulator.plot_chromatogram(sim_reverse_default[2][1][end], (0.0, 1000.0), annotation=false)
	chrom_reverse = GasChromatographySimulator.plot_chromatogram(sim_reverse[2][1][end], (0.0, 1000.0), annotation=false)
	plot(chrom_reverse_default[1], chrom_reverse[2], chrom_reverse[3], c=:red)
end

# ╔═╡ Cell order:
# ╠═e50d12ba-5e53-11f1-aa44-a152c26b54bb
# ╠═8b8d2ec9-0cc3-4537-ab3a-a1f364af72a8
# ╠═8d49ebdd-0b20-41c9-8794-dc605aa61e91
# ╠═07bcea72-f32d-4c05-9c98-996cc7efd038
# ╠═c3a5bce8-b532-429f-9dae-7ea3526829a8
# ╠═800e9468-ee5e-4ced-afe7-711bdd74946d
# ╠═7c1bad0c-1ca8-4d67-9310-55b361276893
# ╠═9aad392d-788a-4e3e-aa22-b61f8132fef6
# ╠═b28decc2-8067-4fb6-891a-e084c9ed07bc
# ╠═8de3f128-6bda-441e-9f7b-c34700ada555
# ╠═a6029cb9-a58d-4ba9-9bdd-d729fece1fc8
# ╠═f0ec2f20-4699-415e-bc84-023ed11cb4dd
# ╠═f3baf9fd-5719-47f5-8f08-20a883bd8a7e
# ╠═c5e2f032-e6fc-44c5-8f4d-3332ad26bdcb
# ╠═7b775cad-3f3f-4a9f-b1b0-2435cf1edd94
# ╠═67c87a76-de8a-4556-bbdd-0dc1f7b655aa
# ╟─d2225282-6399-41d3-b922-de59461495c1
# ╠═829c7f3a-60c0-4886-b127-cc07933abcf4
# ╠═9065ecec-a02b-419b-9e2a-ae2c74b961f6
# ╠═7a047cef-1b87-44df-bf4e-1a9932df236e
# ╠═35254488-6bf0-4df9-9fe4-97bc4abca2e6
# ╠═549b07b9-6199-459a-9ea8-b9a9e37a6dc8
# ╠═abfcc64c-0c45-4935-a90e-23c104b47376
# ╠═5a36a232-7b13-4688-8ef2-7e6fcd14b182
# ╟─749c9f0f-fb05-49b1-9b82-b38c5741a190
# ╟─9765efe5-6c51-4417-8e1b-2f9231208d1d
