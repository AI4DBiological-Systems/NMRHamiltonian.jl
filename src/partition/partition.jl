
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
"""
αs and Ωs must not contain singlet groups.
"""
function partitionresonances(
    name::String,
    partition_configs::Vector{PartitionAlgorithmType},
    #coherence_state_pairs_sys,
    #ms_sys,
   # αs::Vector{Vector{T}},
    #Ωs::Vector{Vector{T}},
    spin_systems::Vector{SpinSystem{T}},
    N_spins_sys::Vector{Int};
    ME::Vector{Vector{Vector{Int}}} = Vector{Vector{Vector{Int}}}(undef, 0),
    relative_α_thresholds::Vector{T} = ones(T, length(N_spins_sys)) .* 0.01
    ) where {T <: AbstractFloat}

    N_systems = length(spin_systems)
    @assert N_systems == length(N_spins_sys) == length(relative_α_thresholds)

    part_inds_set = Vector{Vector{Vector{Int}}}(undef, N_systems)
    Δc_m_set = Vector{Vector{Vector{T}}}(undef, N_systems)

    as = Vector{Vector{T}}(undef, N_systems)
    Fs = Vector{Vector{T}}(undef, N_systems)

    Δc_bar = Vector{Vector{Vector{T}}}(undef, N_systems)
    c_states = Vector{Vector{Tuple{Int64, Int64}}}(undef, N_systems)

    search_results = Vector{PartitionSearchRecord{T}}(undef, N_systems)

    for i = 1:N_systems

        # ## prune resonance components that have a low intensity.
        # αs_i_prune = copy(spin_systems[i].intensities)
        # Ωs_i_prune = copy(spin_systems[i].frequencies)

        # c_states_prune = spin_systems[i].coherence_state_pairs
        # if relative_α_thresholds[i] > 0

        #     α_tol = relative_α_thresholds[i]*maximum(αs[i])
        #     inds_amp = findall(xx->(xx>α_tol), αs[i])

        #     c_states_prune = (coherence_state_pairs_sys[i])[inds_amp]
        #     αs_i_prune = αs[i][inds_amp]
        #     Ωs_i_prune = Ωs[i][inds_amp]
        # end        

        sp = spin_systems[i]

        ## prune resonance components that have a low intensity.
        Δc_m, inds_amp, c_states_prune = getΔcm(
            sp;
            intensity_tol = relative_α_thresholds[i]*maximum(sp.intensities),
        )

        αs_i_prune = copy(sp.intensities)
        Ωs_i_prune = copy(sp.frequencies)
        if !isempty(inds_amp)
            αs_i_prune = sp.intensities[inds_amp]
            Ωs_i_prune = sp.frequencies[inds_amp]
        end

        # # ## to see Δc_0, before magnetic equivalence is accounted for.
        # @show Δc_m
        # # #inds_c_bar = [[1891], [1622], [56], [277], [1970], [2981]] # isoleucine. 700
        # # #inds_c_bar = [[442], [481], [436], [329], [99], [228], [40], [74], [330]] # alpha-glucose. 700
        # # inds_c_bar = [[42], [15]] # ethanol 700.

        # # for i1 in eachindex(inds_c_bar)
        # #     @show Δc_m[i1]
        # # end
        # @show αs_i_prune
        # @show Ωs_i_prune
        # @assert 1==2

        ## reduce the cardinality of Δc_m if there is magnetic equivalence in this spin system.
        if !isempty(ME)
            if !isempty(ME[i])
                Δc_m = reduceΔc(Δc_m, ME[i], N_spins_sys[i])
            end
        end
        
        # for i1 in eachindex(inds_c_bar)
        #     @show Δc_m[i1]
        # end
        # @assert 1==2
        
        # rescale  before clustering.
        Δc_m_rescaled = scalebydim(Δc_m)
        part_inds, Gs, θs, γs, chosen_ind_γ, connectivity_label,
        connectivity_parameter = getpartitioningfunc(
            partition_configs[i], 
            name,
            i,
            Δc_m_rescaled,
            αs_i_prune,
        )

        # won't yield Δc_bar such that they approximately sum to -1.
        #Δc_centroids = assemblecentroids(WeightedCentroids(), part_inds, Δc_m, αs_i_prune)
        
        #Δc_centroids = assemblecentroids(AveragedCentroids(), part_inds, Δc_m, αs_i_prune)
        
        # safest option.
        Δc_centroids = assemblecentroids(MaxAmplitudeCentroids(), part_inds, Δc_m, αs_i_prune)

        part_inds = Gs[chosen_ind_γ]

        ## update.
        search_results[i] = PartitionSearchRecord(Gs, γs, chosen_ind_γ, θs, connectivity_label, connectivity_parameter)

        as[i] = αs_i_prune
        Fs[i] = Ωs_i_prune

        part_inds_set[i] = part_inds
        Δc_m_set[i] = Δc_m

        Δc_bar[i] = Δc_centroids
        c_states[i] = c_states_prune
    end

    return as, Fs, part_inds_set, Δc_m_set, Δc_bar, c_states, search_results
end

# rescale before partitioning, since sum(Δc_m[l]) is approximately -1, as length(Δc_m[l]) (number of chemical shifts) increases, the values in Δc_m[l] tend to take on smaller values.  This affects the partitioning algorithm's required stopping duality gap tolderance. Instead of changing the tolerance, do this rescaling procedure.
function scalebydim(Δc_m::Vector{Vector{T}})::Vector{Vector{T}} where T
    out = similar(Δc_m)
    for i in eachindex(Δc_m)
        out[i] = Δc_m[i] .* length(Δc_m[i]) 
    end

    return out
end

# mutates z with  PartitionAlgorithmType. Make part_rec an element of Vector ANY. to save trouble.
function getpartitioningfunc(
    part_config::PartitionAlgorithmType,
    name::String,
    i::Int,
    Δc_m_rescaled::Vector{Vector{T}},
    αs_i_prune::Vector{T},
    ) where T <: AbstractFloat

    # create the configs for this input state: (name, i, Δc_m_rescaled, αs_i_prune).
    γ_default = part_config.getdefaultγfunc(name, i, Δc_m_rescaled, αs_i_prune)
    θ_default = part_config.getdefaultθfunc(name, i, Δc_m_rescaled, αs_i_prune)
    
    graph_config = part_config.getgraphconfigfunc(name, i, Δc_m_rescaled, αs_i_prune)
    optim_config = part_config.getoptimconfigfunc(name, i, Δc_m_rescaled, αs_i_prune)
    assignment_config = part_config.getassignmentconfigfunc(name, i, Δc_m_rescaled, αs_i_prune)

    config_θ = part_config.getsearchθconfigfunc(name, i, Δc_m_rescaled, αs_i_prune)
    config_γ = part_config.getsearchγconfigfunc(name, i, Δc_m_rescaled, αs_i_prune)

    if part_config.verbose_kernel || part_config.report_cost
        println("Molecule: $name, spin system: $i, number of resonance components: $(length(Δc_m_rescaled)).")
    end

    config_γ_send = config_γ
    if isfinite(γ_default)
        config_γ_send = nothing
    end

    config_θ_send = config_θ
    if isfinite(θ_default)
        config_θ_send = nothing
    end

    # run convex clustering to get partition of Δc_m_rescaled.
    Gs, θs, γs, connectivity_label, connectivity_parameter = getpartitionθγ!(
        graph_config,
        Δc_m_rescaled,
        αs_i_prune,
        optim_config,
        assignment_config,
        config_γ_send,
        config_θ_send;
        γ_default = γ_default,
        θ_default = θ_default,
        store_trace = part_config.store_trace,
        report_cost = part_config.report_cost,
        verbose_kernel = part_config.verbose_kernel,
    )

    # if the final iteration is less than the number of nuclei, then we try to look for another solution in Gs.

    # default case: the final solution is at the end of the trace buffer. If no trace is stored, the only entry of Gs is the final solution.
    part_inds = Gs[end]

    # need to refine the partition solution: sometimes the size of the partition is smaller than the number of nuclei in the spin system that already had magnetic equivalence taken care of. This might give rise to significant restriction in expression for the resultant model.
    N_nuclei = length(Δc_m_rescaled[begin])

    chosen_ind_γ = length(Gs)
    if length(Gs[end]) < N_nuclei
        preferred_inds = findall(xx->(length(xx)>N_nuclei), Gs)
        if !isempty(preferred_inds)
            # take the smallest partition, which is the most closest to the end.
            chosen_ind_γ = maximum(preferred_inds)
            part_inds = Gs[chosen_ind_γ]
        end
    end

    # # store configurations and search results into record.
    # part_rec.iter_θ = iter_θ # 0 if no search was done.
    # part_rec.iter_γ = iter_γ # 0 if no search was done.
    # part_rec.optim_config = optim_config
    # part_rec.assignment_config = assignment_config
    # part_rec.graph_config = graph_config
    # part_rec.config_θ = config_θ
    # part_rec.config_γ = config_γ

    #rec_γ = SearchγRecord(iter_γ, config_γ)
    #rec_θ = SearchθRecord(iter_θ, config_θ)
    #return part_inds, Δc_centroids, rec_θ, rec_γ # type unstable since config_γ or config_θ are uions of types.

    if part_config.verbose_kernel || part_config.report_cost
        println("")
    end

    return part_inds, Gs, θs, γs, chosen_ind_γ,
    connectivity_label, connectivity_parameter
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


#= # weighted average with weights from αs.
function assemblecentroids(
    ::WeightedCentroids,
    G::Vector{Vector{Int}},
    Δc_m::Vector{Vector{T}},
    αs::Vector{T},
    )::Vector{Vector{T}} where T <: AbstractFloat
    
    centroids = Vector{Vector{T}}(undef, length(G))
    for k in eachindex(G)
        inds = G[k]
        centroids[k] = computeweightedcentroid(αs[inds], Δc_m[inds])
    end
    
    return centroids


# make c into convex weights a, then return the weighted centroid.
function computeweightedcentroid(c::Vector{T}, X::Vector{Vector{T}})::Vector{T} where T <: AbstractFloat
    @assert length(c) == length(X)

    sum_c = sum(c)
    a = c ./ sum_c

    return sum( a .* X ) ./ length(X)
end
end =#


############ for testing.

"""
function checkNgroups(As::Vector{SHType{T}}) where T

Checks if the number of Δc_bar resonance groups for each spin system is at least as many of the number of nuclei (after accounting for magnetic equivalence) in that system.
Returns `valid_flag::Bool`, `n::Int`, `i::Int`
If `valid_flag` is `false`, then `(n,i)` is the (compound index, spin system index) of the first offending spin system.

Set `γ_base` in `createsearchγconfigs` to a smaller number to start with more resonance groups.

For some compounds, this could make sense. Such as GISSMO's DSS; it has 6 nuclei in the non-0 ppm system, but often gets only 3 resonance groups for the range of model construction parameters we tested. It means the compound is more constrained than magnetic equivalence, and warrants future study. For this reason, failure of this test does not mean a bug.
"""
function checkNgroups(As::Vector{SHType{T}}; verbose = false) where T
    for n in eachindex(As)
        c_bar = As[n].Δc_bar

        for i in eachindex(c_bar)
            N_Nuclei = length(c_bar[i][begin])
            N_non_singlet_groups = length(c_bar[i])

            if verbose
                @show (n,i), N_Nuclei, N_non_singlet_groups
            end

            if N_non_singlet_groups < N_Nuclei
                return false, n, i
            end
        end
    end

    return true, 0, 0
end