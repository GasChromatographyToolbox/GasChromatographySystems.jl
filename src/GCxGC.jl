# GCxGC peaklist / chromatogram / fit helpers

function add_A_to_pl!(pl, df_A)
	# add the areas in the dataframe df_A to the peaklist, same CAS and Annotations are used to identify the correct solute
#	sort_A_foc = Array{Float64}(undef, length(df_A.A_foc))
#	sort_A_unfoc = Array{Float64}(undef, length(df_A.A_foc))
	sort_A = Array{Float64}(undef, length(df_A.A))
	for i=1:length(df_A.A)
		ii = common_index(pl, df_A.CAS[i], join(split(df_A.Annotations[i], "_")[1:end-1], "_"))
		#ii = common_index(pl, df_A.CAS[i], split(df_A.Annotations[i], ", ")[1])
		sort_A[ii] = df_A.A[i]
	end
	pl[!, :A] = sort_A
	return pl
end

# identfy the common index of a substances CAS number and Annotations in a peaklist 
function common_index(pl, CAS, Annotation)
	ii_CAS = findall(CAS.==pl.CAS)
	ii_ann = findall(occursin.(Annotation, pl.Annotations))
	ii = intersect(ii_CAS, ii_ann)[1]
	return ii
end

## further functions for evaluation of GCxGC_TM results
function peaklist_GCxGC(pl_end, pl_1D, PM)#, shift)
	# estimated from the gaussian peak fit to apex of projected retention times
	# here the shift is not known
	#fit_D1 = fit_gauss_D1(pl_end, pl_1D, PM, shift)
	#fit_D2 = fit_gauss_D2(pl_end, PM, shift)
	fit_D1 = fit_gauss_D1(pl_end, pl_1D, PM)
	fit_D2 = fit_gauss_D2(pl_end, PM)
	tR1 = Array{Float64}(undef, length(fit_D1.Name))
	tR2 = Array{Float64}(undef, length(fit_D1.Name))
	for i=1:length(fit_D1.Name)
		ii = findfirst(fit_D1.Name[i].==fit_D2.Name)
		tR1[i] = fit_D1.fits[i].param[1]
		tR2[i] = fit_D2.fits[ii].param[1]
	end
	return DataFrame(Name=fit_D1.Name, tR1=tR1, tR2=tR2)
end

"""
    peaklist_GCxGC(pl_end, PM; digits=6)

Calculate retention times and peak widths for compounds in comprehensive two-dimensional gas chromatography (GC×GC).

This function processes a peak list from a GC×GC separation to calculate:
1. First dimension retention times (tR1) and peak widths (τR1)
2. Second dimension retention times (tR2) and peak widths (τR2)
3. Weighted averages of these parameters based on peak heights

# Arguments
- `pl_end`: DataFrame containing the peak list with columns:
  - `Name`: Compound names
  - `CAS`: CAS numbers
  - `tR`: Retention times
  - `τR`: Peak widths
  - `A`: Peak areas
- `PM`: Modulation period in seconds
- `digits`: Number of digits for rounding (default: 6)

# Returns
A DataFrame containing:
- `Name`: Unique compound names
- `CAS`: CAS numbers
- `tR1`: Weighted mean retention time in first dimension
- `tR2`: Weighted mean retention time in second dimension
- `τR1`: Peak width in first dimension (weighted standard deviation)
- `τR2`: Peak width in second dimension (weighted mean)
- `tR1s`: Array of individual retention times in first dimension
- `tR2s`: Array of individual retention times in second dimension
- `τR2s`: Array of individual peak widths in second dimension
- `heights`: Array of peak heights used for weighting

# Notes
- First dimension retention times (tR1) are calculated using `mod_base_time`
- Second dimension retention times (tR2) are calculated using `mod_time`
- Peak widths are calculated using weighted statistics:
  - τR1: Weighted standard deviation of tR1 values
  - τR2: Weighted mean of individual τR values
- Weights are based on peak heights to give more importance to larger peaks
- The function handles multiple peaks of the same compound by combining their parameters
"""
function peaklist_GCxGC(pl_end, PM; digits=6)
	# estimated from the weighted (height) means of projected retention times
	# here the shift is not known and therefore ignored
	heights = GasChromatographySystems.heights_of_peaks(pl_end)
	tR1 = mod_base_time.(pl_end.tR, PM; digits=digits)
	tR2 = mod_time.(pl_end.tR, PM; digits=digits)
	Name = unique(pl_end.Name)
	nsub = length(Name)
	CAS = pl_end.CAS[[findfirst(Name[i].==pl_end.Name) for i=1:length(Name)]]
	
	sort_heights = Array{Array{Float64}}(undef, nsub)
	sort_tR1 = Array{Array{Float64}}(undef, nsub)
	sort_tR2 = Array{Array{Float64}}(undef, nsub)
	mean_tR1 = Array{Float64}(undef, nsub)
	mean_tR2 = Array{Float64}(undef, nsub)
	sort_τR2 = Array{Array{Float64}}(undef, nsub)
	mean_τR2 = Array{Float64}(undef, nsub)
	τR1 = Array{Float64}(undef, nsub)  # standard deviation (peak width) in 1st dimension
	for i=1:nsub
		ii_name = findall(Name[i].==pl_end.Name)
		sort_heights[i] = heights[ii_name]
		sort_tR1[i] = tR1[ii_name]
		sort_tR2[i] = tR2[ii_name]
		mean_tR1[i] = sum(sort_tR1[i].*sort_heights[i])/sum(sort_heights[i])
		mean_tR2[i] = sum(sort_tR2[i].*sort_heights[i])/sum(sort_heights[i])
		sort_τR2[i] = pl_end.τR[ii_name]
		mean_τR2[i] = sum(sort_τR2[i].*sort_heights[i])/sum(sort_heights[i])
		
		# Calculate weighted variance and standard deviation for 1st dimension
		τR1[i] = sqrt(sum(sort_heights[i] .* (sort_tR1[i] .- mean_tR1[i]).^2) / sum(sort_heights[i]))
		# for second dimension this value is not used as it is smaller than the widths of each peak 
		#τR2[i] = sqrt(sum(sort_heights[i] .* (sort_tR2[i] .- mean_tR2[i]).^2) / sum(sort_heights[i]))
	end
	return DataFrame(Name=Name, CAS=CAS, tR1=mean_tR1, tR2=mean_tR2, τR1=τR1, τR2=mean_τR2, tR1s=sort_tR1, tR2s=sort_tR2, τR2s=sort_τR2, heights=sort_heights)
end

function heights_of_peaks(pl)
	heights = Array{Float64}(undef, length(pl.tR))
 	for i=1:length(pl.tR)
		heights[i] = (GasChromatographySimulator.chromatogram([pl.tR[i]], [pl.tR[i]], [pl.τR[i]])[1].*pl.A[i])
	end
	return heights
end

function fit_envelope(pl_D2, pl_D1)
	heights = heights_of_peaks(pl_D2)
	tR = pl_D2.tR
	Name = unique(pl_D2.Name)
	nsub = length(Name)

	@. model_g(x,p) = p[3]/sqrt(2*π*p[2]^2)*exp(-(x-p[1])^2/(2*p[2]^2))
	
	sort_heights = Array{Array{Float64}}(undef, nsub)
	sort_tR = Array{Array{Float64}}(undef, nsub)
	fits = Array{LsqFit.LsqFitResult}(undef, nsub)
	for i=1:nsub
		ii_name = findall(Name[i].==pl_D2.Name)
		ii = findfirst(Name[i].==pl_D1.Name)
		sort_heights[i] = heights[ii_name]
		sort_tR[i] = tR[ii_name]
		fits[i] = curve_fit(model_g, sort_tR[i], sort_heights[i], [pl_D1.tR[ii], pl_D1.τR[ii], 1.0])
	end
	return DataFrame(Name=Name, tRs=sort_tR, heights=sort_heights, fits=fits)
end

function fit_gauss_D1(pl_D2, pl_D1, PM)
	# here the shift is not known and therefore ignored
	heights = heights_of_peaks(pl_D2)
	# test this change replace fld() with e.g. mod_time():
	tR = fld.(pl_D2.tR, PM).*PM
	#tR = mod_base_time.(pl_D2.tR, PM)
	Name = unique(pl_D2.Name)
	nsub = length(Name)

	@. model_g(x,p) = p[3]/sqrt(2*π*p[2]^2)*exp(-(x-p[1])^2/(2*p[2]^2))
	
	sort_heights = Array{Array{Float64}}(undef, nsub)
	sort_tR = Array{Array{Float64}}(undef, nsub)
	fits = Array{LsqFit.LsqFitResult}(undef, nsub)
	for i=1:nsub
		ii_name = findall(Name[i].==pl_D2.Name)
		ii = findfirst(Name[i].==pl_D1.Name)
		sort_heights[i] = heights[ii_name]
		sort_tR[i] = tR[ii_name]
		fits[i] = curve_fit(model_g, sort_tR[i], sort_heights[i], [pl_D1.tR[ii], pl_D1.τR[ii], 1.0])
	end
	return DataFrame(Name=Name, tRs=sort_tR, heights=sort_heights, fits=fits)
end

function fit_gauss_D2(pl_D2, PM)
	# here the shift is not known and therefore ignored
	heights = heights_of_peaks(pl_D2)
	# test this change replace fld() with e.g. mod_time():
	tR = pl_D2.tR .- (fld.(pl_D2.tR, PM).*PM)
	#tR = mod_time.(pl_D2.tR, PM)
	Name = unique(pl_D2.Name)
	nsub = length(Name)

	@. model_g(x,p) = p[3]/sqrt(2*π*p[2]^2)*exp(-(x-p[1])^2/(2*p[2]^2))
	
	sort_heights = Array{Array{Float64}}(undef, nsub)
	sort_tR = Array{Array{Float64}}(undef, nsub)
	fits = Array{LsqFit.LsqFitResult}(undef, nsub)
	for i=1:nsub
		ii_name = findall(Name[i].==pl_D2.Name)
		#ii = findfirst(Name[i].==pl_D1.Name)
		sort_heights[i] = heights[ii_name]
		sort_tR[i] = tR[ii_name]
		mean_tR = sum(sort_tR[i])/length(sort_tR[i])
		fits[i] = curve_fit(model_g, sort_tR[i], sort_heights[i], [mean_tR, mean_tR/10, 1.0])
	end
	return DataFrame(Name=Name, tRs=sort_tR, heights=sort_heights, fits=fits)
end

function chrom(pl; nτ=6)
	tstart = Array{Float64}(undef, length(pl.Name))
	tend = Array{Float64}(undef, length(pl.Name))
	t = Array{Array{Float64}}(undef, length(pl.Name))
	c = Array{Array{Float64}}(undef, length(pl.Name))
	for i=1:length(pl.Name)
		tstart[i] = pl.tR[i] - nτ * pl.τR[i]
		tend[i] = pl.tR[i] + nτ * pl.τR[i]
		t[i] = collect(tstart[i]:(2*nτ*pl.τR[i]/100):tend[i])
		c[i] = GasChromatographySimulator.chromatogram(t[i], [pl.tR[i]], [pl.τR[i]])*pl.A[i]
	end
	t1 = minimum(tstart)
	t2 = maximum(tend)
	dt = 2*nτ*minimum(pl.τR)/100
	t_sum = collect(t1:dt:t2)
	c_sum = fill(0.0, length(t_sum))
	for i=1:length(pl.Name)
		c_ = GasChromatographySimulator.chromatogram(t_sum, [pl.tR[i]], [pl.τR[i]])*pl.A[i]
		c_sum = c_sum .+ c_
	end

	names = unique(pl.Name)
	
	p_chrom = Plots.plot(t_sum, c_sum, xlabel="time in s", label="Chromatogram")
	for i=1:length(pl.Name)
		i_names = findfirst(pl.Name[i].==names)
		#if i > length(names)
		#	lbl = ""
		#else
			lbl = pl.Name[i]
		#end
		Plots.plot!(p_chrom, t[i], c[i], label=lbl, color=i_names+1)
	end
	Plots.plot!(p_chrom, ylims=(-0.02*maximum(c_sum), 1.02*maximum(c_sum)), xlims=(minimum(t_sum), maximum(t_sum)))
	return p_chrom, t_sum, c_sum, t, c 
end

function chrom_marked(pl, PM, ratio, shift; nτ=6)
	p_chrom, t_sum, c_sum, t, c = chrom(pl; nτ=nτ)
	max_y = maximum(c_sum)
	# add modulation period
	# test this change replace fld() with e.g. mod_time():
	n = unique(fld.(t_sum, PM + shift))
	#n = unique(mod_number.(t_sum, shift, PM, ratio))
	for i=1:length(n)
		Plots.plot!(p_chrom, [n[i]*PM-shift, n[i]*PM-shift], [0.0, 1.1*maximum(c_sum)], c=:black, label="")
		Plots.plot!(p_chrom, [n[i]*PM-shift+PM*ratio, n[i]*PM-shift+PM*ratio], [0.0, 1.1*maximum(c_sum)], c=:black, linestyle=:dash, label="")
	end
	Plots.plot!(p_chrom, ylims=(-0.02*maximum(c_sum), 1.02*maximum(c_sum)), xlims=(minimum(t_sum), maximum(t_sum)))
	return p_chrom, t_sum, c_sum, t, c
end

function collect_chrom(pl_array, sys; markings=true)
	# collect chromatograms for all segments
	# chromatograms of segments, which are followed by a ModuleTM should be marked with the modulation period (coldjet, hotjet)
	# add titles (module name)
	p = Array{Any}(undef, length(pl_array))
	for i=1:length(pl_array)
		if markings == true && i < length(pl_array)
			if typeof(sys.modules[i+1]) == GasChromatographySystems.ModuleTM
				p[i] = chrom_marked(pl_array[i], sys.modules[i+1].PM, sys.modules[i+1].ratio, sys.modules[i+1].shift)[1]
				
			else
				p[i] = chrom(pl_array[i])[1]
			end
		else
			p[i] = chrom(pl_array[i])[1]
		end
		Plots.plot!(p[i], title=sys.modules[i].name)
	end
	return p
end

#=function chrom_slicing(t, c, PM, shift) # correctly account for the shift!!!
	n = Int.(fld.(collect(t).-shift, PM)) # number of the slices
	slices = Array{Array{Float64, 1}, 1}(undef, length(unique(n)))
	t_D2 = Array{Array{Float64, 1}, 1}(undef, length(unique(n)))
	for i=1:length(unique(n))
		i1 = findfirst(unique(n)[i].==n)
		i2 = findlast(unique(n)[i].==n)
		slices[i] = c[i1:i2]
		t_D2[i] = t[i1:i2] .- unique(n)[i] * PM
	end
	t_D1 = (0.0:PM:t[end]).- shift# shift?
	return slices, t_D1, t_D2
end=#

"""
    chrom2d(pl_final, sys, PM)

Generate a 2D chromatogram from a peak list by slicing the chromatogram into modulation periods.

# Arguments
- `pl_final`: DataFrame containing the peak list with columns for retention times (tR), 
              peak widths (τR), and peak areas (A)
- `sys`: System structure containing the GC system configuration
- `PM`: Modulation period in seconds

# Returns
- `slice_mat`: 2D matrix containing the sliced chromatogram data, where:
  - rows represent modulation periods
  - columns represent time points within each period
- `t_D1`: Array of first dimension retention times (start of each modulation period)
- `t_D2`: Array of arrays containing second dimension retention times for each slice
- `slices`: Array of arrays containing the chromatogram data for each slice
- `t`: Time points used for chromatogram generation
- `chrom_sum`: Summed chromatogram before slicing

# Notes
- The function generates a chromatogram by summing individual peak contributions
- Time points are generated with 0.01s intervals
- The chromatogram is sliced into modulation periods of length PM
- For incomplete slices at the end, the matrix is padded with zeros
- Shift and ratio parameters are not used as they are not needed for 2D chromatogram generation
"""
function chrom2d(pl_final, sys, PM)
	t = 0.0:0.01:sum(GasChromatographySystems.common_timesteps(sys))
	chrom = Array{Array{Float64}}(undef, length(pl_final.tR))
	for i=1:length(pl_final.tR)
		chrom[i] = GasChromatographySimulator.chromatogram(collect(t), [pl_final.tR[i]], [pl_final.τR[i]]).*pl_final.A[i]
	end
	chrom_sum = chrom[1]
	for i=2:length(chrom)
		chrom_sum = chrom_sum .+ chrom[i]
	end

	t_D1 = (0.0:PM:t[end])
	
	# replace fld() with e.g. mod_number():
	# a shift and ratio are unknown here (e.g. in ChromSpace only PM is defined at 2D chromatogram generation) 
	n = Int.(floor.(t./PM))
	slices = Array{Array{Float64, 1}, 1}(undef, length(unique(n)))
	t_D2 = Array{Array{Float64, 1}, 1}(undef, length(unique(n)))
	for i=1:length(unique(n))
		i1 = findfirst(unique(n)[i].==n)
		i2 = findlast(unique(n)[i].==n)
		slices[i] = chrom_sum[i1:i2]
		t_D2[i] = t[i1:i2] .- (unique(n)[i] * PM )
	end

	# expected number of complete slices
	expected_slices = Int(floor(t[end]/PM))
	slice_mat = Array{Float64}(undef, expected_slices, length(t_D2[1]))
	for j=1:length(t_D2[1])
		for i=1:expected_slices
			 if i <= length(slices) && j <= length(slices[i])
                slice_mat[i,j] = slices[i][j]
            else
                slice_mat[i,j] = 0.0  # default value
            end
		end
	end
	return slice_mat, t_D1, t_D2, slices, t, chrom_sum
end

function comparison_meas_sim(meas, pl_sim)
	Name = meas.Name
	tR1_meas = meas.tR1
	tR2_meas = meas.tR2
	index = [findfirst(meas.Name[x].==pl_sim.Name) for x in 1:length(meas.Name)]
	tR1_sim = Array{Union{Missing,Float64}}(undef, length(meas.Name))
	tR2_sim = Array{Union{Missing,Float64}}(undef, length(meas.Name))
	for i=1:length(meas.Name)
		if isnothing(index[i])
			tR1_sim[i] = missing
			tR2_sim[i] = missing
		else
			tR1_sim[i] = pl_sim.tR1[index[i]]
			tR2_sim[i] = pl_sim.tR2[index[i]]
		end
	end
	comp = DataFrame(Name=Name, tR1_meas=tR1_meas, tR1_sim=tR1_sim, ΔtR1=tR1_meas.-tR1_sim, relΔtR1_percent=(tR1_meas.-tR1_sim)./tR1_meas.*100.0, tR2_meas=tR2_meas, tR2_sim=tR2_sim, ΔtR2=tR2_meas.-tR2_sim, relΔtR2_percent=(tR2_meas.-tR2_sim)./tR2_meas.*100.0)
	#for i=1:length(comp.Name)
	#	ii = findfirst(comp.Name[i].==meas.Name)
	#	if ismissing(meas.tR1[ii])
	#		comp[i, :tR1_meas] = NaN
	#		comp[i, :tR2_meas] = NaN
	#		comp[i, :ΔtR1] = NaN
	#		comp[i, :ΔtR2] = NaN
	#	else
	#		comp[i, :tR1_meas] = meas.tR1[ii]
	#		comp[i, :tR2_meas] = meas.tR2[ii]
	#		comp[i, :ΔtR1] = meas.tR1[ii] - comp[i, :tR1_sim]
	#		comp[i, :ΔtR2] = meas.tR2[ii] - comp[i, :tR2_sim]
	#	end
	#end
	return comp
end

function chrom_slicing(t1, c, PM)
	# test this change replace fld() with e.g. mod_time():
	n = Int.(fld.(collect(t1), PM)) # number of the slices
	#n = mod_number.(collect(t1), 0.0, PM, 0.5) # number of the slices
	slices = Array{Array{Float64, 1}, 1}(undef, length(unique(n)))
	t_D2 = Array{Array{Float64, 1}, 1}(undef, length(unique(n)))
	for i=1:length(unique(n))
		i1 = findfirst(unique(n)[i].==n)
		i2 = findlast(unique(n)[i].==n)
		slices[i] = c[i1:i2]
		# test this change replace fld() with e.g. mod_time():
		#t_D2[i] = t1[i1:i2] .- unique(n)[i] * PM
		t_D2[i] = mod_time.(t1[i1:i2], PM)
	end
	t_D1 = 0.0:PM:t1[end]
	return slices, t_D1, t_D2
end

function check_area(pl)
	names = unique(pl.Name)
	area = Array{Float64}(undef, length(names))
	for i=1:length(names)
		area[i] = sum(filter([:Name] => x -> x == names[i], pl).A)
	end
	return DataFrame(Name=names, sum_A = area)
end

function check_peakwidths(pl; τ_threshold = 0.5)
	pl_f = filter([:τR] => x -> x > τ_threshold, pl)
	names = unique(pl_f.Name)
	return names
end

function check_duration_modulation(pl_array, par, PM, ratio)
	Δts = duration_in_module(pl_array, par)
	ok_TM1 = Array{Bool}(undef, length(Δts[3].Δt))
	for i=1:length(Δts[3].Δt)
		if (Δts[3].Δt[i] > ratio*PM) && (Δts[3].Δt[i] < PM)
			ok_TM1[i] = true
		else
			ok_TM1[i] = false
		end
	end
	ok_TM2 = Array{Bool}(undef, length(Δts[5].Δt))
	
	for i=1:length(Δts[5].Δt)
		#if (Δts[5].Δt[i] > ratio*PM) && (Δts[5].Δt[i] < PM) # change this condition
		if (Δts[5].Δt[i] < PM)
			ok_TM2[i] = true
		else
			ok_TM2[i] = false
		end
	end
	index_TM1 = findall(ok_TM1.==false)
	index_TM2 = findall(ok_TM2.==false)
	TM1 = DataFrame(Index=index_TM1, Name=Δts[3].Name[index_TM1], CAS=Δts[3].CAS[index_TM1], tR=pl_array[3].tR[index_TM1], τR=pl_array[3].τR[index_TM1], Δt=Δts[3].Δt[index_TM1], Annotations=Δts[3].Annotations[index_TM1])
	TM2 = DataFrame(Index=index_TM2, Name=Δts[5].Name[index_TM2], CAS=Δts[5].CAS[index_TM2], tR=pl_array[5].tR[index_TM2], τR=pl_array[5].τR[index_TM2], Δt=Δts[5].Δt[index_TM2], Annotations=Δts[5].Annotations[index_TM2])
	return TM1, TM2
end

function duration_in_module(pl_array, par)
	Δts = Array{DataFrame}(undef, length(pl_array))
	for i=1:length(pl_array)
		CAS_par = [par[i].sub[x].CAS for x in 1:length(par[i].sub)]
		ann_par = [par[i].sub[x].ann for x in 1:length(par[i].sub)]

		Δt = Array{Float64}(undef, length(CAS_par))
		name = Array{String}(undef, length(CAS_par))
		for j=1:length(CAS_par)
			jj = GasChromatographySystems.common_index(pl_array[i], CAS_par[j], ann_par[j])
			Δt[j] = pl_array[i].tR[jj] - par[i].sub[j].t₀
			name[j] = pl_array[i].Name[jj]
		end
		df_Δt = DataFrame(Name=name, CAS=CAS_par, Δt=Δt, Annotations=ann_par)
		Δts[i] = df_Δt
	end
	return Δts
end

