


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


function fitproxies(As::Vector{SHType{T}},
    dummy_SSFID::SST,
    λ0::T;
    names::Vector{String} = Vector{String}(undef, 0),
    config_path::String = "",
    Δcs_max_scalar_default = 0.2,
    κ_λ_lb_default = 0.5,
    κ_λ_ub_default = 2.5,
    u_min = Inf,
    u_max = Inf,
    Δr_default = 1.0,
    Δκ_λ_default = 0.05) where {T,SST}

    config_dict = Dict()
    if ispath(config_path)

        # TODO add error-handling if name is not found in the dictionary or filename does not exist.
        config_dict = JSON.parsefile(config_path)
    end

    cores = Vector{FIDModelType{T,SST}}(undef, length(As))

    for n = 1:length(As)
        A = As[n]

        # fit surrogate, save into `core`.
        cores[n] = fitproxy(dummy_SSFID, A, λ0, config_dict;
        compound_name = names[n],
        Δcs_max_scalar_default = Δcs_max_scalar_default,
        κ_λ_lb_default = κ_λ_lb_default,
        κ_λ_ub_default = κ_λ_ub_default,
        u_min = u_min,
        u_max = u_max,
        Δr_default = Δr_default,
        Δκ_λ_default = Δκ_λ_default)
    end

    return cores
end


function fitproxy(dummy_SSFID::SST,
    A::SHType{T},
    λ0::T,
    config_dict;
    compound_name::String = "",
    Δcs_max_scalar_default = 0.2,
    κ_λ_lb_default = 0.5,
    κ_λ_ub_default = 2.5,
    u_min = Inf,
    u_max = Inf,
    Δr_default = 1.0,
    Δκ_λ_default = 0.05)::FIDModelType{T,SST} where {T,SST}

    hz2ppmfunc = uu->(uu - A.ν_0ppm)*A.SW/A.fs
    ppm2hzfunc = pp->(A.ν_0ppm + pp*A.fs/A.SW)

    # allocate `core` data structure.
    N_singlets = length(A.αs_singlets)
    κs_λ_singlets = ones(T, N_singlets)
    κs_β_singlets = zeros(T, N_singlets)
    d_singlets = zeros(T, N_singlets)

    #N_β_vars_sys = A.N_spins_sys # no equivalence used.
    N_β_vars_sys::Vector{Int} = collect( length(A.Δc_bar[i][1]) for i = 1:length(A.Δc_bar) )

    # proxy placeholder.
    #qs = Vector{Vector{Function}}(undef, length(N_spins_sys))

    SSFID_obj = setupSSFIDparams(dummy_SSFID, A.part_inds_compound, N_β_vars_sys)

    # prepare configuration parameters.
    κ_λ_lb = κ_λ_lb_default
    κ_λ_ub = κ_λ_ub_default
    Δκ_λ = Δκ_λ_default
    Δcs_max::Vector{T} = collect( Δcs_max_scalar_default for i = 1:length(A.N_spins_sys))
    Δr = Δr_default

    d_max::Vector{T} = ppm2hzfunc.(Δcs_max) .- ppm2hzfunc(zero(T))

    if !isempty(config_dict) && !isempty(compound_name)
        dict = config_dict[compound_name] # TODO graceful error-handle.

        κ_λ_lb = dict["κ_λ lower bound"]
        κ_λ_ub = dict["κ_λ upper bound"]
        Δκ_λ = dict["κ_λ surrogate sampling step size"]
        if length(dict["Δcs_max"]) == length(A.N_spins_sys)
            Δcs_max = dict["Δcs_max"]
        # else
        #     println("Warning: problem an entry's Δcs_max value in config file. Using default scalar value for Δcs_max")
        end
        Δr = dict["radians surrogate sampling step size"]

        d_max = collect( ppm2hzfunc(Δcs_max[i])-ppm2hzfunc(0.0) for i = 1:length(Δcs_max) )
    end

    # threshold and partition the resonance components.
    if !isfinite(u_min) || !isfinite(u_max)
        Ωs_ppm = hz2ppmfunc.( combinevectors(A.Ωs) ./ (2*π) )
        min_ppm = minimum(Ωs_ppm) - 0.5
        max_ppm = maximum(Ωs_ppm) + 0.5
        u_min = ppm2hzfunc(min_ppm)
        u_max = ppm2hzfunc(max_ppm)
    end

    qs = setupcompoundpartitionitp(d_max,
        SSFID_obj.κs_β,
        A.Δc_bar,
        A.part_inds_compound,
        A.αs, A.Ωs,
        λ0, u_min, u_max;
        κ_λ_lb = κ_λ_lb,
        κ_λ_ub = κ_λ_ub,
        Δr = Δr,
        Δκ_λ = Δκ_λ)

    core = FIDModelType(qs, SSFID_obj, κs_λ_singlets, κs_β_singlets, d_singlets,
        Δcs_max, λ0)

    return core
end

function setupmixtureSH(target_names::Vector{String},
    H_params_path::String,
    dict_compound_to_filename,
    fs::T, SW::T, ν_0ppm::T;
    MEs::Vector{Vector{Vector{Vector{Int}}}} = Vector{Vector{Vector{Vector{Int}}}}(undef, 0),
    config_path = "",
    tol_coherence = 1e-2,
    α_relative_threshold = 0.05,
    Δc_partition_radius = 1e-1,
    simple_coherence_atol::T = -1.2) where {T <: Real, SST}

    ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

    N_compounds = length(target_names)
    As = Vector{SHType{T}}(undef, N_compounds)

    for n = 1:N_compounds

        ME = Vector{Vector{Vector{Int}}}(undef, 0)
        if !isempty(MEs)
            ME = MEs[n]
        end

        αs, Ωs, part_inds_compound, Δc_m_compound, Δc_bar, N_spins_sys,
            αs_singlets, Ωs_singlets = setupcompoundSH(target_names[n],
            H_params_path, dict_compound_to_filename, ppm2hzfunc,
            fs, SW, ν_0ppm;
            config_path = config_path,
            ME = ME,
            tol_coherence = tol_coherence,
            α_relative_threshold = α_relative_threshold,
            Δc_partition_radius = Δc_partition_radius,
            simple_coherence_atol = simple_coherence_atol)

        As[n] = SHType(αs, Ωs, Δc_m_compound, part_inds_compound,
            Δc_bar, N_spins_sys, αs_singlets, Ωs_singlets, fs, SW, ν_0ppm)
    end

    return As
end




function setupcompoundSH(name, base_path, dict_compound_to_filename,
    ppm2hzfunc,
    fs::T, SW::T, ν_0ppm::T;
    ME::Vector{Vector{Vector{Int}}} = Vector{Vector{Vector{Int}}}(undef, 0),
    config_path::String = "",
    tol_coherence = 1e-2,
    α_relative_threshold = 0.05,
    Δc_partition_radius = 1e-1,
    simple_coherence_atol::T = -1.2) where T <: Real

    # TODO add error-handling if name is not found in the dictionary, or filename does not exist.
    load_path = joinpath(base_path, dict_compound_to_filename[name]["file name"])
    H_IDs, H_css, J_IDs, J_vals = loadcouplinginfojson(load_path)

    if ispath(config_path)

        # TODO add error-handling if name is not found in the dictionary or filename does not exist.
        db_dict = JSON.parsefile(config_path)
        dict = db_dict[name]

        tol_coherence = dict["coherence tolerance"]
        α_relative_threshold = dict["relative amplitude threshold"]
        Δc_partition_radius = dict["maximum Δc deviation"]
        simple_coherence_atol = dict["simple coherence absolute tolerance"]
    end


    # SH.
    # p_cs_sys, css_sys, cs_singlets, J_vals_sys, J_IDs_sys, intermediates_sys,
    # cs_LUT, cs_singlets_compact, N_spins_singlet, css, p_compound,
    # cs_compound = fetchSHparameters(H_IDs, H_css, J_IDs, J_vals)
    J_inds_sys, J_inds_sys_local, J_IDs_sys, J_vals_sys, H_inds_sys,
        cs_sys, H_inds_singlets, cs_singlets, H_inds, J_inds,
        g = setupcsJ(H_IDs, H_css, J_IDs, J_vals)

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

    return αs, Ωs, part_inds_compound, Δc_m_compound, Δc_bar, N_spins_sys,
    αs_singlets, Ωs_singlets
end

function setupSSFIDparams(dummy_SSFID::SpinSysParamsType1{T}, part_inds_compound, N_β_vars_sys)::SpinSysParamsType1{T} where T
    L = length(part_inds_compound)

    κs_λ = ones(T, L)
    κs_β = collect( zeros(T, N_β_vars_sys[i]) for i = 1:length(N_β_vars_sys))
    #d = rand(length(αs))
    d = zeros(T, L)

    return constructorSSFID(dummy_SSFID, κs_λ, κs_β, d)
end

function setupSSFIDparams(dummy_SSFID::SpinSysParamsType2{T}, part_inds_compound::Vector{Vector{Vector{Int}}}, N_β_vars_sys)::SpinSysParamsType2{T} where T

    N_sys = length(part_inds_compound)
    κs_λ = Vector{Vector{T}}(undef, N_sys)
    d = Vector{Vector{T}}(undef, N_sys)

    for i = 1:length(d)
        N_partition_elements = length(part_inds_compound[i])

        κs_λ[i] = ones(T, N_partition_elements)
        d[i] = zeros(T, N_partition_elements)
    end

    κs_β = collect( zeros(T, N_β_vars_sys[i]) for i = 1:length(N_β_vars_sys))

    return constructorSSFID(dummy_SSFID, κs_λ, κs_β, d)
end
