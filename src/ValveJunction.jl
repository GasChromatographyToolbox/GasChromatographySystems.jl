# Peak splitting at junction vertices where a ModuleValve meets the chromatographic path.

"""
    ValveSlicingSchedule

Named tuple returned by [`valve_slicing_schedule`](@ref):

- `mp`: valve modulation period (s)
- `t_closed`: closed-phase duration per period (s); open duration is `mp - t_closed`
- `phase_shift`: time offset of the closed/open pattern (s), aligned with [`PeriodicValveProgram`](@ref) `t_start` when applicable
"""
const ValveSlicingSchedule = NamedTuple{(:mp, :t_closed, :phase_shift)}

"""
    edges_along_path_in_order(g, path)

Return global edge indices for `path` in path order (one index per edge in `path`).
"""
function edges_along_path_in_order(
	g::Graphs.AbstractGraph,
	path::AbstractVector{<:Graphs.AbstractEdge},
)::Vector{Int}
	E = collect(edges(g))
	idx = Vector{Int}(undef, length(path))
	for (k, e) in pairs(path)
		i = findfirst(==(e), E)
		i === nothing && throw(ArgumentError("edge $e not in graph"))
		idx[k] = i
	end
	return idx
end

"""
    incident_valve_modules(sys, vertex)

All `ModuleValve` instances on edges incident to `vertex`.
"""
function incident_valve_modules(sys::GasChromatographySystems.System, vertex::Integer)
	E = collect(edges(sys.g))
	valves = ModuleValve[]
	for i in eachindex(E)
		e = E[i]
		(src(e) == vertex || dst(e) == vertex) || continue
		mod = sys.modules[i]
		mod isa ModuleValve && push!(valves, mod)
	end
	return valves
end

"""
    valve_modulation_period(vp::PeriodicValveProgram)

Modulation period `mp` in seconds.
"""
valve_modulation_period(vp::PeriodicValveProgram) = vp.mp

"""Fraction of period `mp` spent closed (for [`slicing`](@ref) / [`mod_number`](@ref) only)."""
_slicing_closed_fraction(t_closed::Real, mp::Real) = Float64(t_closed) / Float64(mp)

"""
    valve_slicing_schedule(vp)

Build a [`ValveSlicingSchedule`](@ref) `(mp, t_closed, phase_shift)` for junction peak slicing.

**DPM / package convention** ([`valve_state`](@ref): `true` = open, `false` = closed):

- **Closed** (`σ = 0`): mod line 4→2 nearly blocked; tee `p₂` set mainly by columns (series-like limit).
- **Open** (`σ = 1`): mod line conducts (`d_open`); `p₂` couples to programmed `p₄` — forward flow on
  2→3 and the split at the tee change with time (already in [`flow_functions`](@ref) / `path_possible`).

Junction transport slices peaks by closed/open periods and retimes slices to the next **open** window
(see [`simplified_valve_junction`](@ref)). That is a lumped timing model, not a re-simulation of column 1.

Returns `nothing` when the schedule is effectively constant over the run.
"""
function valve_slicing_schedule(vp::PeriodicValveProgram)
	mp = vp.mp
	tol = eps(Float64) * 100
	if mp <= tol
		return nothing
	end
	if vp.t_closed <= tol || open_dur(mp, vp.t_closed) <= tol
		return nothing
	end
	return (mp=mp, t_closed=vp.t_closed, phase_shift=vp.t_start)
end

function valve_slicing_schedule(vp::ValveProgram)
	tol = eps(Float64) * 100
	length(vp.state_steps) == 1 && return nothing
	unique(vp.state_steps) == 1 && return nothing
	mp = sum(vp.time_steps)
	mp <= tol && return nothing
	t_closed = sum(vp.time_steps[i] for i in eachindex(vp.time_steps) if !vp.state_steps[i])
	if t_closed <= tol || t_closed >= mp - tol
		return nothing
	end
	return (mp=mp, t_closed=t_closed, phase_shift=0.0)
end

"""
    valve_state_varies(vp, t_lo, t_hi)

Return `true` if [`valve_state`](@ref) is not constant on `[t_lo, t_hi]`.
"""
function valve_state_varies(vp::AbstractValveProgram, t_lo::Real, t_hi::Real)
	t_lo = Float64(t_lo)
	t_hi = Float64(t_hi)
	if t_hi <= t_lo
		return false
	end
	period = vp isa PeriodicValveProgram ? vp.mp : sum(vp.time_steps)
	n = max(3, ceil(Int, (t_hi - t_lo) / max(period / 4, eps(Float64))))
	ts = range(t_lo, t_hi; length=n)
	states = [valve_state(vp, t) for t in ts]
	return !allequal(states)
end

function valve_state_varies_for_peaklist(vp::AbstractValveProgram, pl; nτ::Integer=6)
	τmax = isempty(pl.τR) ? 0.0 : maximum(pl.τR)
	t_lo = minimum(pl.tR) - nτ * τmax
	t_hi = maximum(pl.tR) + nτ * τmax
	return valve_state_varies(vp, t_lo, t_hi)
end

"""
    slice_peaks_by_valve(pl, mp, t_closed, phase_shift, par; nτ=6, τ₀=zeros(length(pl.τR)), kwargs...)

Slice peaks into one row per valve period (Gaussian area split; annotations `v1_`, `v2_`, …).

Uses the same period grid as [`slicing`](@ref), with `t_closed` / `mp` defining the closed fraction of
each period (`true` = open in [`valve_state`](@ref); closed intervals are the complement).
"""
function slice_peaks_by_valve(
	pl,
	mp::Real,
	t_closed::Real,
	phase_shift::Real,
	par::GasChromatographySimulator.Parameters;
	nτ=6,
	τ₀=zeros(length(pl.τR)),
	kwargs...,
)
	closed_fraction = _slicing_closed_fraction(t_closed, mp)
	return slicing(
		pl, mp, closed_fraction, phase_shift, par;
		nτ=nτ, τ₀=τ₀, ann_prefix="v", kwargs...,
	)
end

"""
    t_start_next_open_window(t_slice, mp, t_closed, phase_shift)

Absolute time (s) when the next **open** valve phase starts, for a slice that begins at `t_slice`.

Each period `mp` consists of `t_closed` s closed then `(mp - t_closed)` s open (see [`valve_state`](@ref)).
"""
function t_start_next_open_window(t_slice::Real, mp::Real, t_closed::Real, phase_shift::Real)
	t_open = open_dur(Float64(mp), Float64(t_closed))
	closed_fraction = _slicing_closed_fraction(t_closed, mp)
	period_index = mod_number(t_slice, phase_shift, mp, closed_fraction)
	return period_index * mp - (t_closed - phase_shift) - t_open
end

"""
    simplified_valve_junction(prev_peaklist, df_A, mp, t_closed, phase_shift)

Assign each valve slice a downstream injection time `tR` at the start of the next **open** window.

# What this does / does not do

- **Does:** After [`slice_peaks_by_valve`](@ref), retime slices leaving the upstream column so the next
  [`simulate_ModuleColumn`](@ref) sees bands that enter column 2 when the mod line is **open** (DPM:
  mod port active, tee coupled to `p₄`). Areas come from `df_A`; `τR` is copied from the parent peak (CAS).
- **Does not:** Change the already simulated upstream column. Time-varying flows on 1→2 and 2→3 come from
  the hydraulic solve (`p₂(t)`, `σ(t)`). Tee split fractions are deferred (workplan §5.6).

Use `inverted=true` on [`PeriodicValveProgram`](@ref) if hardware open/closed labels are swapped.
"""
function simplified_valve_junction(prev_peaklist, df_A, mp::Real, t_closed::Real, phase_shift::Real)
	sort_df_A = sort(df_A, :t0)
	t_slice = sort_df_A.t0
	A = sort_df_A.A
	Name = sort_df_A.Name
	CAS = sort_df_A.CAS
	Ann = sort_df_A.Annotations
	t_release_open = t_start_next_open_window.(t_slice, mp, t_closed, phase_shift)
	τR = Vector{Float64}(undef, length(t_release_open))
	for i in eachindex(t_release_open)
		ii = findfirst(==(CAS[i]), prev_peaklist.CAS)
		τR[i] = ii === nothing ? 0.0 : prev_peaklist.τR[ii]
	end
	pl = DataFrame(Name=Name, CAS=CAS, tR=t_release_open, τR=τR, Annotations=Ann, A=A)
	return pl, nothing
end

"""
    simulate_valve_junction(segment_par, valve::ModuleValve, prev_peaklist; nτ=6)

Split peaks at a path junction where an incident [`ModuleValve`](@ref) modulates the tee.

Runs [`slice_peaks_by_valve`](@ref) then [`simplified_valve_junction`](@ref). Constant `state` over the
peak width passes through unchanged.
"""
function simulate_valve_junction(
	segment_par::GasChromatographySimulator.Parameters,
	valve::ModuleValve,
	prev_peaklist;
	nτ=6,
)
	vp = valve.state
	if !valve_state_varies_for_peaklist(vp, prev_peaklist; nτ=nτ)
		return segment_par, prev_peaklist, nothing
	end
	sched = valve_slicing_schedule(vp)
	sched === nothing && return segment_par, prev_peaklist, nothing
	@assert sched isa ValveSlicingSchedule
	(; mp, t_closed, phase_shift) = sched
	τ₀ = prev_peaklist.τR
	new_segment_par, df_A = slice_peaks_by_valve(
		prev_peaklist, mp, t_closed, phase_shift, segment_par; nτ=nτ, τ₀=τ₀,
	)
	peaklist, solutions = simplified_valve_junction(prev_peaklist, df_A, mp, t_closed, phase_shift)
	return new_segment_par, peaklist, solutions
end

"""
    apply_valve_junctions_at_vertex(sys, vertex, segment_par, peaklist; nτ=6)

Apply every incident valve at `vertex` in sequence (typical tee: one valve).
"""
function apply_valve_junctions_at_vertex(
	sys::GasChromatographySystems.System,
	vertex::Integer,
	segment_par::GasChromatographySimulator.Parameters,
	peaklist;
	nτ=6,
)
	par = segment_par
	pl = peaklist
	sols = Any[]
	for valve in incident_valve_modules(sys, vertex)
		par, pl, sol = simulate_valve_junction(par, valve, pl; nτ=nτ)
		sol !== nothing && push!(sols, sol)
	end
	return par, pl, isempty(sols) ? nothing : sols
end
