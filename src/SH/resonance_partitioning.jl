

function allabssmaller(x, threshold)
    if maximum(abs.(x)) > threshold
        return false
    end

    return true
end

function sumabssmaller(x, threshold)
    if sum(abs.(x)) > threshold
        return false
    end

    return true
end

function prunecombocoherences!(A, α_relative_threshold, tol_coherence, Δc_partition_radius)

    for i = 1:length(A.Δc_m_compound)

        keep_flags = collect( allabssmaller(A.Δc_m_compound[i][l],1+tol_coherence) for l = 1:length(A.Δc_m_compound[i]) )
        #keep_flags = collect( sumabssmaller(A.Δc_m_compound[i][l],1+tol_coherence) for l = 1:length(A.Δc_m_compound[i]) )

        A.αs[i] = A.αs[i][keep_flags]
        A.Ωs[i] = A.Ωs[i][keep_flags]
        A.Δc_m_compound[i] = A.Δc_m_compound[i][keep_flags]

        # remake Δc_bar.
        α_tol = α_relative_threshold*maximum(A.αs[i])
        part_inds, Δc_centroids = partitionresonancesbyneighbors(A.Δc_m_compound[i],
            A.αs[i], α_tol; radius = Δc_partition_radius)

        A.Δc_bar[i] = Δc_centroids
        A.part_inds_compound[i] = part_inds
    end

    return nothing
end

function prunecombocoherencesbar!(A, α_relative_threshold, tol_coherence, Δc_partition_radius)

    for i = 1:length(A.Δc_m_compound)

        keep_flags = collect( allabssmaller(A.Δc_m_compound[i][l],1+tol_coherence) for l = 1:length(A.Δc_m_compound[i]) )
        #keep_flags = collect( sumabssmaller(A.Δc_m_compound[i][l],1+tol_coherence) for l = 1:length(A.Δc_m_compound[i]) )

        A.αs[i] = A.αs[i][keep_flags]
        A.Ωs[i] = A.Ωs[i][keep_flags]
        A.Δc_m_compound[i] = A.Δc_m_compound[i][keep_flags]

        # remake Δc_bar.
        α_tol = α_relative_threshold*maximum(A.αs[i])
        part_inds, Δc_centroids = partitionresonancesbyneighbors(A.Δc_m_compound[i],
            A.αs[i], α_tol; radius = Δc_partition_radius)

        #
        keep_flags = collect( sumabssmaller(Δc_centroids[l],1+tol_coherence) for l = 1:length(Δc_centroids) )

        A.Δc_bar[i] = Δc_centroids[keep_flags]
        A.part_inds_compound[i] = part_inds[keep_flags]
    end

    return nothing
end

function prunecombocoherencesbar2!(Δc_m_compound,
    αs, Ωs, Δc_bar, part_inds_compound;
    α_relative_threshold = 0.05, tol_coherence = 1e-2, Δc_partition_radius = 0.3)

    for i = 1:length(Δc_m_compound)

        keep_flags = collect( allabssmaller(Δc_m_compound[i][l],1+tol_coherence) for l = 1:length(Δc_m_compound[i]) )
        #keep_flags = collect( sumabssmaller(A.Δc_m_compound[i][l],1+tol_coherence) for l = 1:length(A.Δc_m_compound[i]) )

        αs[i] = αs[i][keep_flags]
        Ωs[i] = Ωs[i][keep_flags]
        Δc_m_compound[i] = Δc_m_compound[i][keep_flags]

        # remake Δc_bar.
        α_tol = α_relative_threshold*maximum(αs[i])
        part_inds, Δc_centroids = partitionresonancesbyneighbors(Δc_m_compound[i],
            αs[i], α_tol; radius = Δc_partition_radius)

        #
        keep_flags = collect( sumabssmaller(Δc_centroids[l],1+tol_coherence) for l = 1:length(Δc_centroids) )

        Δc_bar[i] = Δc_centroids[keep_flags]
        part_inds_compound[i] = part_inds[keep_flags]
    end

    return nothing
end

#####
"""
Uses the NearestNeighbors.jl search library to determine which resonances should be grouped together.
"""
function partitionresonancesbyneighbors(x_in::Vector{Vector{T}},
    amplitudes_in::Vector{T},
    threshold_amplitude::T;
    radius = 1e-1,
    Minkowski_parameter = 3.5) where T <: Real

    @assert length(x_in) == length(amplitudes_in)
    @assert !isempty(x_in)

    # take largest amplitude, set as new partition centre.
    amplitudes = copy(amplitudes_in)
    x = copy(x_in)
    inds_buffer = collect(1:length(x))

    out_inds = Vector{Vector{Int}}(undef, 0)

    x_maxs = Vector{Vector{T}}(undef, 0)

    max_amplitude, max_ind = findmax(amplitudes)
    while !isempty(amplitudes) || max_amplitude < threshold_amplitude

        # set up tree.
        X = Matrix{T}(undef, length(x[1]), length(x))
        for r = 1:length(x)
            X[:,r] = x[r]
        end

        balltree = NearestNeighbors.BallTree(X,
            NearestNeighbors.Minkowski(Minkowski_parameter);
            reorder = false)

        inds = NearestNeighbors.inrange(balltree, x[max_ind], radius, true)
        push!(x_maxs, x[max_ind])

        # book keep.
        push!(out_inds, inds_buffer[inds])

        deleteat!(amplitudes, inds)
        deleteat!(x, inds)
        deleteat!(inds_buffer, inds)

        if isempty(amplitudes)
            return out_inds, x_maxs
        end
        max_amplitude, max_ind = findmax(amplitudes)
    end

    return out_inds, x_maxs
end

"""
αs and Ωs must not contain singlet groups.
"""
function partitionresonances(coherence_state_pairs_sys, ms_sys,
    αs::Vector{Vector{T}}, Ωs::Vector{Vector{T}},
    N_spins_sys::Vector{Int};
    α_relative_threshold = 0.05,
    Δc_partition_radius = 1e-1,
    ME::Vector{Vector{Vector{Int}}} = Vector{Vector{Vector{Int}}}(undef, 0),
    simple_coherence_atol::T = -1.2) where T

    N_groups = length(αs)

    part_inds_set = Vector{Vector{Vector{Int}}}(undef, N_groups)
    Δc_m_set = Vector{Vector{Vector{T}}}(undef, N_groups)

    # as = Vector{Vector{Vector{T}}}(undef, N_groups)
    # Fs = Vector{Vector{Vector{T}}}(undef, N_groups)
    as = Vector{Vector{T}}(undef, N_groups)
    Fs = Vector{Vector{T}}(undef, N_groups)

    N_systems = length(coherence_state_pairs_sys)
    @assert N_groups == N_systems == length(N_spins_sys)
    Δc_bar = Vector{Vector{Vector{T}}}(undef, N_systems)


    for i = 1:N_systems

        α_tol = α_relative_threshold*maximum(αs[i])
        inds_amp = findall(xx->(xx>α_tol), αs[i])

        c_states_prune = (coherence_state_pairs_sys[i])[inds_amp]
        αs_i_prune = αs[i][inds_amp]
        Ωs_i_prune = Ωs[i][inds_amp]

        c_m_r = collect( ms_sys[i][r] for (r,s) in c_states_prune )
        c_m_s = collect( ms_sys[i][s] for (r,s) in c_states_prune )
        Δc_m = collect( c_m_r[j] - c_m_s[j] for j = 1:length(c_m_r))

        # println("Δc_m = ", Δc_m)
        # println("ME[i] = ", ME[i])
        # println("N_spins_sys[i] = ", N_spins_sys[i])
        if !isempty(ME)
            if !isempty(ME[i])
                Δc_m = reduceΔc(Δc_m, ME[i], N_spins_sys[i])
            end
        end


        # # check if we should discard the non-simple coherences.
        # if simple_coherence_atol > 0
        #
        #     inds = findsimplecoherences(Δc_m, atol = simple_coherence_atol)
        #
        #     αs_i_prune = αs_i_prune[inds]
        #     Ωs_i_prune = Ωs_i_prune[inds]
        #     Δc_m = Δc_m[inds]
        # end

        part_inds, Δc_centroids = partitionresonancesbyneighbors(Δc_m,
            αs_i_prune, α_tol; radius = Δc_partition_radius)



        # partition_size = length(part_inds)
        # as[i] = Vector{Vector{T}}(undef, partition_size)
        # Fs[i] = Vector{Vector{T}}(undef, partition_size)
        # for m = 1:partition_size
        #
        #     as[i][m] = αs_i_prune[part_inds[m]]
        #     Fs[i][m] = Ωs_i_prune[part_inds[m]]
        # end
        as[i] = αs_i_prune
        Fs[i] = Ωs_i_prune

        part_inds_set[i] = part_inds
        Δc_m_set[i] = Δc_m

        Δc_bar[i] = Δc_centroids
    end



    # for i = 1:N_systems # over elements in a spin group.

    #     N_partition_elements = length(part_inds_compound[i])
    #     Δc_bar[i] = Vector{Vector{T}}(undef, N_partition_elements)

    #     for k = 1:N_partition_elements

    #         inds = part_inds_compound[i][k]


    #         Δc_bar[i][k] = Statistics.mean( Δc_m_compound[i][inds] )

    #         # # weighted mean.
    #         α = as[i][inds]
    #         Ω = Ωs[i][inds]
    #         # tmp = Δc_m_compound[i][inds]
    #         # Δc_bar[i][k] = sum(tmp[l] .* α[l]) / sum(α)
    #     end
    # end

    return as, Fs, part_inds_set, Δc_m_set, Δc_bar
end




function setupmixtureSH(target_names::Vector{String},
    fs::T, SW::T, ν_0ppm::T,
    Phys::Vector{PhysicalParamsType{T}};
    config_path = "",
    tol_coherence = 1e-2,
    α_relative_threshold = 0.05,
    Δc_partition_radius = 1e-1,
    simple_coherence_atol::T = -1.2,
    prune_combo_Δc_bar_flag::Bool = true) where {T <: Real, SST}

    ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

    N_compounds = length(target_names)
    As = Vector{SHType{T}}(undef, N_compounds)

    for n = 1:N_compounds

        # ME = Vector{Vector{Vector{Int}}}(undef, 0)
        # if !isempty(MEs)
        #     ME = MEs[n]
        # end

        αs, Ωs, part_inds_compound, Δc_m_compound, Δc_bar, N_spins_sys,
            αs_singlets, Ωs_singlets = setupcompoundSH(target_names[n],
            Phys[n].J_inds_sys_local,
            Phys[n].J_vals_sys,
            Phys[n].cs_sys,
            Phys[n].H_inds_singlets,
            Phys[n].cs_singlets,
            ppm2hzfunc,
            fs, SW, ν_0ppm;
            config_path = config_path,
            ME = Phys[n].ME,
            tol_coherence = tol_coherence,
            α_relative_threshold = α_relative_threshold,
            Δc_partition_radius = Δc_partition_radius,
            simple_coherence_atol = simple_coherence_atol,
            prune_combo_Δc_bar_flag = prune_combo_Δc_bar_flag)

        As[n] = SHType(αs, Ωs, Δc_m_compound, part_inds_compound,
            Δc_bar, N_spins_sys, αs_singlets, Ωs_singlets, fs, SW, ν_0ppm)
    end

    return As
end




function setupcompoundSH(name,
    J_inds_sys_local,
    J_vals_sys,
    cs_sys,
    H_inds_singlets,
    cs_singlets,
    ppm2hzfunc,
    fs::T, SW::T, ν_0ppm::T;
    ME::Vector{Vector{Vector{Int}}} = Vector{Vector{Vector{Int}}}(undef, 0),
    config_path::String = "",
    tol_coherence = 1e-2,
    α_relative_threshold = 0.05,
    Δc_partition_radius = 1e-1,
    simple_coherence_atol::T = -1.2,
    prune_combo_Δc_bar_flag = true) where T <: Real

    if ispath(config_path)

        # TODO add error-handling if name is not found in the dictionary or filename does not exist.
        db_dict = JSON.parsefile(config_path)
        dict = db_dict[name]

        tol_coherence = dict["coherence tolerance"]
        α_relative_threshold = dict["relative amplitude threshold"]
        Δc_partition_radius = dict["maximum Δc deviation"]
        simple_coherence_atol = dict["simple coherence absolute tolerance"]
    end


    ### SH.
    # #TODO add error-handling if name is not found in the dictionary, or filename does not exist.
    # load_path = joinpath(base_path, dict_compound_to_filename[name]["file name"])
    # H_IDs, H_css, J_IDs, J_vals = loadcouplinginfojson(load_path)

    # J_inds_sys, J_inds_sys_local, J_IDs_sys, J_vals_sys, H_inds_sys,
    #     cs_sys, H_inds_singlets, cs_singlets, H_inds, J_inds,
    #     g = setupcsJ(H_IDs, H_css, J_IDs, J_vals)

    N_spins_singlet = length.(H_inds_singlets)

    N_spins_sys = collect( length(cs_sys[m]) for m = 1:length(cs_sys) )
    intermediates_sys = prepcouplingalgorithm(N_spins_sys)

    # css_sys becomes cs_sys
    # cs_singlets_compact becomes cs_singlets.

    αs_inp, Ωs_inp, coherence_mat_sys, eig_vals_sys, Q_sys,
    coherence_state_pairs_sys, H_sys, states_sys, ms_sys,
    M_sys = setupcompoundspectrum!(cs_sys,
        J_vals_sys, J_inds_sys_local, intermediates_sys,
        ppm2hzfunc, cs_singlets, N_spins_singlet, fs, SW;
        tol_coherence = tol_coherence)

    k = findfirst(xx->length(xx)==1, αs_inp)
    αs_spin_sys = copy(αs_inp)
    Ωs_spin_sys = copy(Ωs_inp)
    αs_singlets = Vector{T}(undef, 0)
    Ωs_singlets = Vector{T}(undef, 0)
    if typeof(k) == Int
        # the case there are singlet groups in this compound.
        αs_spin_sys = αs_spin_sys[1:k-1]
        Ωs_spin_sys = Ωs_spin_sys[1:k-1]

        αs_singlets = collect( αs_inp[l][1] for l = k:length(αs_inp) )
        Ωs_singlets = collect( Ωs_inp[l][1] for l = k:length(Ωs_inp) )
    end

    # spectra eval func for spin systems.
    N_spins_sys = collect( length(cs_sys[i]) for i = 1:length(cs_sys))

    αs, Ωs, part_inds_compound,
    Δc_m_compound, Δc_bar = partitionresonances(coherence_state_pairs_sys,
    ms_sys, αs_spin_sys, Ωs_spin_sys, N_spins_sys;
    ME = ME,
    α_relative_threshold = α_relative_threshold,
    Δc_partition_radius = Δc_partition_radius,
    simple_coherence_atol = simple_coherence_atol)

    # f = uu->evalcLcompoundviapartitions(uu, d,
    # αs, Ωs, κs_λ, κs_β, λ0, Δc_m_compound, part_inds_compound)

    if prune_combo_Δc_bar_flag

        # The following pruning strategy gave dubious spectra shape for L-Leucine.
        #prunecombocoherences!(As[1], α_relative_threshold, tol_coherence, Δc_partition_radius)

        #prunecombocoherencesbar!(As[n], α_relative_threshold, tol_coherence, Δc_partition_radius)
        prunecombocoherencesbar2!(Δc_m_compound,
            αs, Ωs, Δc_bar, part_inds_compound;
            α_relative_threshold = α_relative_threshold, tol_coherence = tol_coherence, Δc_partition_radius = Δc_partition_radius)
    end

    return αs, Ωs, part_inds_compound, Δc_m_compound, Δc_bar, N_spins_sys,
    αs_singlets, Ωs_singlets
end
