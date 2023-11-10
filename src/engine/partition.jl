
# return variable inds_amp is empty if no state pair was pruned.
function getΔcm(
    ms::Vector{Vector{T}},
    coherence_state_pairs::Vector{Tuple{Int,Int}},
    intensities::Vector{T};
    intensity_tol::T = zero(T),
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
    N_spins_sys::Vector{Int};
    ME::Vector{Vector{Vector{Int}}} = Vector{Vector{Vector{Int}}}(undef, 0),
    relative_α_threshold::T = convert(T, 0.01),
    ) where {T <: AbstractFloat}

    N_systems = length(spin_systems)
    @assert N_systems == length(N_spins_sys)

    part_inds_set = Vector{Vector{Vector{Int}}}(undef, N_systems)
    Δc_m_set = Vector{Vector{Vector{T}}}(undef, N_systems)

    as = Vector{Vector{T}}(undef, N_systems)
    Fs = Vector{Vector{T}}(undef, N_systems)

    Δc_bar = Vector{Vector{Vector{T}}}(undef, N_systems)
    c_states = Vector{Vector{Tuple{Int64, Int64}}}(undef, N_systems)

    for i = 1:N_systems

        # ## prune resonance components that have a low intensity.
        sp = spin_systems[i]

        ## prune resonance components that have a low intensity.
        Δc_m, inds_amp, c_states_prune = getΔcm(
            sp;
            intensity_tol = relative_α_threshold*maximum(sp.intensities),
        )

        αs_i_prune = copy(sp.intensities)
        Ωs_i_prune = copy(sp.frequencies)
        if !isempty(inds_amp)
            αs_i_prune = sp.intensities[inds_amp]
            Ωs_i_prune = sp.frequencies[inds_amp]
        end

        ## reduce the cardinality of Δc_m if there is magnetic equivalence in this spin system.
        if !isempty(ME)
            if !isempty(ME[i])
                Δc_m = reduceΔc(Δc_m, ME[i], N_spins_sys[i])
            end
        end

        # rescale  before clustering.
        Δc_m_rescaled = scalebydim(Δc_m)

        part_inds = getpartition(
            Δc_m_rescaled,
        )

        # won't yield Δc_bar such that they approximately sum to -1.
        #Δc_centroids = assemblecentroids(WeightedCentroids(), part_inds, Δc_m, αs_i_prune)
        
        #Δc_centroids = assemblecentroids(AveragedCentroids(), part_inds, Δc_m, αs_i_prune)
        
        # safest option.
        Δc_centroids = assemblecentroids(MaxAmplitudeCentroids(), part_inds, Δc_m, αs_i_prune)

        as[i] = αs_i_prune
        Fs[i] = Ωs_i_prune

        part_inds_set[i] = part_inds
        Δc_m_set[i] = Δc_m

        Δc_bar[i] = Δc_centroids
        c_states[i] = c_states_prune
    end

    return as, Fs, part_inds_set, Δc_m_set, Δc_bar, c_states
end

# rescale before partitioning, since sum(Δc_m[l]) is approximately -1, as length(Δc_m[l]) (number of chemical shifts) increases, the values in Δc_m[l] tend to take on smaller values.  This affects the partitioning algorithm's required stopping duality gap tolderance. Instead of changing the tolerance, do this rescaling procedure.
function scalebydim(Δc_m::Vector{Vector{T}})::Vector{Vector{T}} where T
    out = similar(Δc_m)
    for i in eachindex(Δc_m)
        out[i] = Δc_m[i] .* length(Δc_m[i]) 
    end

    return out
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
    Δc_m_rescaled::Vector{Vector{T}};
    nuc_factor::T = convert(T, 1.5),
    tol_radius_1D::T = convert(T, 0.1), # strictly between 0 and 1. The lower, the better the approximation, but would a larger partition (i.e. more resonance groups).
    ) where T <: AbstractFloat

    X = Δc_m_rescaled

    N = length(X)
    D = length(X[begin])
    N_nuc = D

    metricfunc = (xx,yy)->norm(xx-yy)
    distance_set, partition_set = SL.runsinglelinkage(
        X,
        metricfunc;
        early_stop_distance = convert(T, Inf),
    )
    default_size = floor(Int, N_nuc *nuc_factor)

    h = getradius(tol_radius_1D, D)
    ind = findfirst(xx->xx>h, distance_set)
    if isnothing(ind)
        ind = max(1, length(X) - default_size)
    else
        ind -= 1
    end

    h_size = length(X) - ind
    if h_size < default_size
        #
        h_size = default_size
    end

    ind = max(1, length(X) - default_size)

    ind2 = findfirst(xx->xx>h, cumsum(distance_set))
    if !isnothing(ind2)
        if ind2 > 1
            ind2 -= 1
        else
            ind2 = nothing
        end
    end
    
    if !isnothing(ind2)
        if ind2 > ind
            # ind2 points to a smaller partition than ind, and is safe to use, since cumsum < h.
            ind = ind2
        end
    end

    part_inds = partition_set[ind]

    return part_inds
end

# use the Δc that has the largest amplitude.
function assemblecentroids(
    ::MaxAmplitudeCentroids,
    G::Vector{Vector{Int}},
    Δc_m::Vector{Vector{T}},
    αs::Vector{T},
    )::Vector{Vector{T},
    } where T <: AbstractFloat
    
    centroids = Vector{Vector{T}}(undef, length(G))
    for k in eachindex(G)
        inds = G[k]
        val, j = findmax(αs[inds])

        centroids[k] = Δc_m[inds][j]
        
        #@assert isapprox(sum(centroids[k]), -1; atol = 0.1)
    end
    
    return centroids
end


# simple average to get centroid.
function assemblecentroids(
    ::AveragedCentroids,
    G::Vector{Vector{Int}},
    Δc_m::Vector{Vector{T}},
    args...
    )::Vector{Vector{T}} where T <: AbstractFloat
    
    centroids = Vector{Vector{T}}(undef, length(G))
    for k in eachindex(G)
        inds = G[k]
        centroids[k] = sum(Δc_m[inds]) ./ length(inds)
    end
    
    return centroids
end
