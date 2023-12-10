
# return variable inds_amp is empty if no state pair was pruned.
function getΔcm(
    ms::Vector{Vector{T}},
    coherence_state_pairs::Vector{Tuple{Int,Int}},
    intensities::Vector{T};
    intensity_tol::T = zero(T), # by default, disable pruncing based on intensity.
    ) where T
    
    c_states_prune = coherence_state_pairs
    inds_amp = Vector{Int}(undef, 0)
    if intensity_tol > zero(T) && isfinite(intensity_tol)

        inds_amp = findall(xx->(xx>intensity_tol), intensities)
        c_states_prune = coherence_state_pairs[inds_amp]
    end

    # assemble coherence drop dual vectors.
    c_m_r = collect( ms[r] for (r,s) in c_states_prune )
    c_m_s = collect( ms[s] for (r,s) in c_states_prune )
    Δc_m = collect( c_m_r[j] - c_m_s[j] for j in eachindex(c_m_r))

    return Δc_m, inds_amp, c_states_prune
end

function getΔcm(A::SpinSystem{T}; intensity_tol::T = zero(T)) where T
    return getΔcm(
        A.partial_quantum_numbers,
        A.coherence_state_pairs,
        A.intensities;
        intensity_tol = intensity_tol,
    )
end

# TODO check and return gracefully.
# αs and Ωs must not contain singlet groups.
function partitionresonances(
    spin_systems::Vector{SpinSystem{T}},
    N_spins_sys::Vector{Int},
    config::SHConfig;
    ME::Vector{Vector{Vector{Int}}} = Vector{Vector{Vector{Int}}}(undef, 0),
    ) where {T <: AbstractFloat}

    relative_α_threshold, total_α_threshold = config.relative_α_threshold, config.total_α_threshold
    max_dev, acceptance_factor = config.max_deviation_from_mean, config.acceptance_factor
    
    N_systems = length(spin_systems)
    @assert N_systems == length(N_spins_sys)

    part_inds_set = Vector{Vector{Vector{Int}}}(undef, N_systems)
    Δc_set = Vector{Vector{Vector{T}}}(undef, N_systems)

    as = Vector{Vector{T}}(undef, N_systems)
    Fs = Vector{Vector{T}}(undef, N_systems)

    Δc_bar = Vector{Vector{Vector{T}}}(undef, N_systems)
    c_states = Vector{Vector{Tuple{Int64, Int64}}}(undef, N_systems)

    for i = 1:N_systems

        # ## prune resonance components that have a low intensity.
        sp = spin_systems[i]

        # debug.
        αs = copy(sp.intensities)
        Ωs = copy(sp.frequencies)
        ms = sp.partial_quantum_numbers
        c_states_prune = sp.coherence_state_pairs
        c_m_r = collect( ms[r] for (r,s) in c_states_prune )
        c_m_s = collect( ms[s] for (r,s) in c_states_prune )
        Δc_m = collect( c_m_r[j] - c_m_s[j] for j in eachindex(c_m_r))
    
        # ## prune resonance components that have a low intensity.
        # Δc_m, inds_amp, c_states_prune = getΔcm(
        #     sp;
        #     intensity_tol = relative_α_threshold*maximum(sp.intensities),
        # )

        # αs_i_prune = copy(sp.intensities)
        # Ωs_i_prune = copy(sp.frequencies)
        # if !isempty(inds_amp)
        #     αs_i_prune = sp.intensities[inds_amp]
        #     Ωs_i_prune = sp.frequencies[inds_amp]
        # end

        ## reduce the cardinality of Δc_m if there is magnetic equivalence in this spin system.
        if !isempty(ME)
            if !isempty(ME[i])
                Δc_m = reduceΔc(Δc_m, ME[i], N_spins_sys[i])
            end
        end

        # remove components.
        #@show length(αs), αs
        αs, Ωs, Δc = processcomponents(αs, Ωs, Δc_m, relative_α_threshold)
        #@show length(αs), αs

        # @show Δc_m
        # @show Δc
        # println()

        # rescale  before clustering.
        #Δc_m_rescaled = scalebydim(Δc_m)

        P = getpartition(Δc, max_dev, acceptance_factor)

        #Δc_centroids = assemblecentroids(MaxAmplitudeCentroids(), part_inds, Δc_m, αs_i_prune)
        P2, as[i], Fs[i], Δc_set[i], Δc_bar[i] = postprocesscomponents(
            P, αs, Ωs, Δc;
            total_α_threshold = total_α_threshold,
        )

        # debug.
        P2_flat = collect( Iterators.flatten(P2))
        #@show length(P2_flat), length(as[i]), length(Δc_bar[i]), length(Fs[i])
        @assert norm(unique(sort(P2_flat)) - collect(1:length(as[i]))) < eps(T)*10        
        #end debug/

        keep_inds = collect( Iterators.flatten(P2))
        c_states[i] = c_states_prune[keep_inds]
        
        part_inds_set[i] = P2
    end

    return as, Fs, part_inds_set, Δc_set, Δc_bar, c_states
end

# check for neligible intensities, duplicates in coherence order differences, and prune enough components to get approximately the target intensity_err.
function processcomponents(
    αs0::Vector{T},
    Ωs0::Vector{T},
    Δc_m0::Vector{Vector{T}},
    relative_α_threshold::T;
    zero_tol = eps(T)*100,
    ) where T <: AbstractFloat

    total_intensity = sum(αs0)

    # # remove components with negligible intensity.
    inds = findall(xx->xx>zero_tol, αs0)
    αs = αs0[inds]
    Ωs = Ωs0[inds]
    Δc_m = Δc_m0[inds]
    
    # # detect duplicates, including frequency.
    level_trait = SL.UseSLDistance(zero_tol)
    center_trait = SL.UseScore(SL.UseMaximum(), αs)

    # augment with frequency.
    X = collect([Δc_m[n]; Ωs[n]] for n in eachindex(Ωs))

    Xc, _, partition = SL.reducepts(
        level_trait,
        center_trait,
        SL.geteuclideanmetric(),
        X,
    )

    # intensity: sum up the duplicates.
    xs0 = collect(
        sum( αs[i] for i in partition[k] )
        for k in eachindex(partition)
    )

    # frequency: take the mean among the duplicates.
    Fs0 = collect(
        Statistics.mean(Ωs[i] for i in partition[k] )
        for k in eachindex(partition)
    )

    # coherence order differenes: remove the last dimension, which was Ωs[n].
    Δc0 = collect(Xc[n][begin:end-1] for n in eachindex(Xc)) # the last dim of Xc is Ωs.

    if relative_α_threshold == zero(T)
        # no further discarding of components required.
        return xs0, Fs0, Δc0
    end

    # # prune enough low intensity components such that xs_intensity = sum(xs) is approximately relative_α_threshold*total_intensity.
    # use a binary search for this.

    f = pp->prunecomponents(pp, xs0, total_intensity)[1]

    intensity_tol, iters_ran = SL.binarysearch(
        f,
        relative_α_threshold,
        zero(T),
        relative_α_threshold*total_intensity;
        atol = relative_α_threshold/10,
        max_iters = 100,
    )
    
    # get the processed components.
    rel_err, inds = prunecomponents(intensity_tol, xs0, total_intensity)
    xs = xs0[inds]
    Fs = Fs0[inds]
    Δc = Δc0[inds]
    #@show rel_err, iters_ran

    return xs, Fs, Δc
end

# for use with Sl.binarysearch().
function prunecomponents(intensity_tol, xs0, nominal_val)
    inds = findall(xx->xx>intensity_tol, xs0)
    
    candidate_val = sum(xs0[i] for i in inds)
    return abs(candidate_val-nominal_val)/nominal_val, inds
end

# given a length r, find the radius of the D-dimensional ball with radius h, such that:
# such that sqrt(r^2 + r^2 + ... + r^2) = sqrt(h). The LHS sum has D times.
function getradius(r::T, D::Int)::T where T
    
    # based on:
    # h^2 = D*r^2 # strive to get the D-dim distance with the given 1D length.
    
    h = sqrt(D*r^2)
    return h
end

# use SL.
function getpartition(
    X::Vector{Vector{T}},
    max_dev::T,
    acceptance_factor::T,
    ) where T <: AbstractFloat

    metric = SL.geteuclideanmetric()
    level_config = SL.UseMaxDeviation(max_dev)

    P, _ = SL.iteratedsl(
        level_config,
        metric,
        X;
        acceptance_factor = acceptance_factor,
        max_iter = 100,
    )

    return P
end

# further reduce components.
function postprocesscomponents(
    P::Vector{Vector{Int}},
    xs::Vector{T}, # intensity
    Fs::Vector{T}, # frequency
    X::Vector{Vector{T}}; # Δc
    total_α_threshold::T = zero(T)
    ) where T <: AbstractFloat

    @assert zero(T) <= total_α_threshold < one(T)

    # for the frequency and Δc, we chose the component with the largest intensity.
    inds_bar = SL.getcenters(
        SL.UseScore(SL.UseMaximum(), xs),
        collect(1:length(X)),
        P,
    )
    #Fs_bar = Fs[inds_bar]
    X_bar = X[inds_bar]

    # for the intensity, we use the total intensity of the part.
    xs_bar = collect(
        sum( xs[i] for i in P[k] )
        for k in eachindex(P)
    )

    if total_α_threshold == zero(T)
        return P, xs, Fs, X, X_bar
    end

    inds = sortperm(xs_bar)
    s = xs_bar[inds]
    Z = sum(xs)
    st_ind = findfirst(xx->xx>(Z*total_α_threshold), cumsum(s))
    @assert !isnothing(st_ind)
    st_ind = max(1, st_ind-1) # make sure we stay within the 0.01.

    P2 = P[inds][st_ind:end]
    keep_inds = collect( Iterators.flatten(P2))

    xs2 = xs[keep_inds]
    Fs2 = Fs[keep_inds]
    X2 = X[keep_inds]

    # remap P2 with respect to the indices of x2, Fs2, X2.
    j = 0
    for k in eachindex(P2)
        for i in eachindex(P2[k])
            j += 1
            P2[k][i] = j
        end
    end

    X2_bar = SL.getcenters(
        SL.UseScore(SL.UseMaximum(), xs2),
        X2,
        P2,
    )

    return P2, xs2, Fs2, X2, X2_bar
end

