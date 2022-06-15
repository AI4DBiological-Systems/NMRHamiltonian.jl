function prepcouplingalgorithm(N_spins::Int)
    # for Hamiltonian.
    Id = getsingleId()
    Ix = getsingleIx()
    Iy_no_im = getsingleIynoim()
    Iz = getsingleIz()

    # for testing or coherence.
    Ip = getsingleI⁺()
    Im = getsingleI⁻()
    Iy = getsingleIy()

    Ims = collect( Im for j = 1:N_spins )
    Im_full = Kronecker.kroneckersum(Ims...)

    Izs = collect( Iz for j = 1:N_spins )
    Iz_full = Kronecker.kroneckersum(Izs...)

    Ixs = collect( Ix for j = 1:N_spins )
    Ix_full = Kronecker.kroneckersum(Ixs...)

    Iys_no_im = collect( Iy_no_im for j = 1:N_spins )
    Iys_no_im_full = Kronecker.kroneckersum(Iys_no_im...)

    return Id, Ix, Iy_no_im, Iz, Ip, Im, Iy,
        Im_full, Iz_full, Ix_full, Iys_no_im_full
end

function prepcouplingalgorithm(N_protons_group::Vector{Int})
    N_groups = length(N_protons_group)

    out_array = Vector{Any}(undef, N_groups)
    for i = 1:N_groups
        out_array[i] = prepcouplingalgorithm(N_protons_group[i])
    end

    out = tuple(out_array...)
    return out
end

function setupcompoundspectrum!( css_sys::Vector{Vector{T}},
    J_vals_sys,
    J_inds_sys,
    intermediates_sys,
    ppm2hzfunc,
    cs_singlets,
    N_spins_singlet::Vector{Int},
    fs::T,
    SW::T;
    tol_coherence = 1e-2) where T <: Real

    @assert length(J_vals_sys) == length(J_inds_sys) == length(css_sys)
    #
    N_sys = length(css_sys)
    N_singlets = length(cs_singlets)

    N_groups = N_sys + N_singlets

    # non-singlet objects.
    states_sys = Vector{Vector{Int}}(undef, N_sys)
    coherence_state_pairs_sys = Vector{Vector{Tuple{Int,Int}}}(undef, N_sys)
    eig_vals_sys = Vector{Vector{T}}(undef, N_sys)
    Q_sys = Vector{Vector{Vector{T}}}(undef, N_sys)
    H_sys = Vector{Matrix{T}}(undef, N_sys)
    coherence_mat_sys = Vector{Matrix{T}}(undef, N_sys)
    ms_sys = Vector{Vector{Vector{T}}}(undef, N_sys)
    M_sys = Vector{Vector{T}}(undef, N_sys)

    # output buffers.
    αs = Vector{Vector{T}}(undef, N_groups)
    Ωs = Vector{Vector{T}}(undef, N_groups)

    # Spin systems.
    for i = 1:N_sys

        αs[i], Ωs[i], coherence_mat_sys[i], eig_vals_sys[i], Q_sys[i], coherence_state_pairs_sys[i],
        H_sys[i], H1, H2, M_sys[i] = evalSCalgorithm!(css_sys[i],
        J_vals_sys[i], J_inds_sys[i],
        intermediates_sys[i],
        ppm2hzfunc;
        tol_coherence = tol_coherence)

        # normalize intensity according to number of spins.
        N_spins = length(css_sys[i])
        normalizetoNspins!(αs[i], N_spins)

        Id = getsingleId()
        ms_sys[i] = computequantumnumbers(Q_sys[i], Id)

        tmp1 = collect( coherence_state_pairs_sys[i][j][1] for j = 1:length(coherence_state_pairs_sys[i]))
        tmp2 = collect( coherence_state_pairs_sys[i][j][2] for j = 1:length(coherence_state_pairs_sys[i]))
        states_sys[i] = unique([tmp1; tmp2])
    end

    # singlets. evalRayleighproxy!() updates singlets for ΩS.
    for i = 1:N_singlets
        αs[i+N_sys] = [ N_spins_singlet[i] ]
        Ωs[i+N_sys] = [ ppm2hzfunc.(cs_singlets[i]) .* (2*π) ]
    end

    return αs, Ωs, coherence_mat_sys, eig_vals_sys, Q_sys,
    coherence_state_pairs_sys, H_sys, states_sys, ms_sys, M_sys
end

function setupapproximateΩmolecule!( css_sys::Vector{Vector{T}},
                                    p_cs_sys::Vector{Vector{T}},
                                    J_vals_sys,
                                    J_IDs_sys,
                                    intermediates_sys,
                                    cs_LUT,
                                    ppm2hzfunc,
                                    cs_singlets_compact,
                                    N_spins_singlet::Vector{Int},
                                    N_spins_sys::Vector{Int},
                                    fs::T,
                                    SW::T;
                                    tol_coherence = 1e-2,
                                    prune_a_tol_factor = 0.05,
                                    fuse_radius = ppm2hzfunc(0.001) - ppm2hzfunc(0.0),
                                    Minkowski_parameter = 3.5,
                                    update_α_flag = false) where T <: Real

    #
    N_sys = length(p_cs_sys)
    N_singlets = length(cs_singlets_compact)

    N_groups = N_sys + N_singlets

    # objects from full solve.
    updateΩFSfuncs = Vector{Function}(undef, N_sys)

    # common objects.
    statess = Vector{Vector{Int}}(undef,N_sys)
    state_labelss = Vector{Vector{Tuple{Int,Int}}}(undef,N_sys)
    eig_vals = Vector{Vector{T}}(undef, N_groups)

    # output buffers.
    αs = Vector{Vector{T}}(undef, N_groups)
    Ωs = Vector{Vector{T}}(undef, N_groups)

    # Spin systems.
    for i = 1:N_sys

        # sanity-check.
        @assert N_spins_sys[i] == length(css_sys[i])

        a0, F0, statess[i], state_labelss[i], H, H1, H2, eig_vals0,
            gs = setupanchorspectrum(p_cs_sys[i], css_sys[i],
                                J_vals_sys[i], J_IDs_sys[i],
                                intermediates_sys[i],
                                cs_LUT[i], ppm2hzfunc;
                                tol_coherence = tol_coherence,
                                prune_a_tol_factor = prune_a_tol_factor,
                                fuse_radius = fuse_radius,
                                Minkowski_parameter = Minkowski_parameter)
        # set up buffers.
        eig_vals[i] = copy(eig_vals0)
        Ωs[i] = copy(F0)
        αs[i] = copy(a0)

        # fullsolve.
        updateΩFSfuncs[i] = setupapproximateΩfullsolve!(eig_vals[i],
            αs[i], Ωs[i], statess[i], state_labelss[i],
            H, gs, fs, SW; update_α_flag = update_α_flag)

    end

    # singlets. evalRayleighproxy!() updates singlets for ΩS.
    for i = 1:N_singlets
        αs[i+N_sys] = [ N_spins_singlet[i] ]
        Ωs[i+N_sys] = [ ppm2hzfunc.(cs_singlets_compact[i]) .* (2*π) ]
    end

    return αs, eig_vals, updateΩFSfuncs, Ωs,
            statess, state_labelss
end


# A version of setupapproximateΩrayleigh() that:
#   - J remains constant, taken from pivot J-coupling values.
#   - fixed coherence using eigenvectors from pivot p_cs.
#   - use full eigensolve to get Ω.
#   - option to use updated eigenvectors to update α.
function setupanchorspectrum( p_cs::Vector{T},
                            css::Vector{T},
                            J_vals::Vector{T},
                            J_IDs::Vector{Tuple{Int,Int}},
                            intermediates,
                            cs_LUT::Vector{Vector{Int}},
                            ppm2hzfunc;
                            tol_coherence = 1e-2,
                            prune_a_tol_factor = 0.05,
                            fuse_radius = ppm2hzfunc(0.001) - ppm2hzfunc(0.0),
                            Minkowski_parameter = 3.5) where T <: Real

    #
    N_spins = length(css)

    F0, a0, state_labels, H, H1, H2, eig_vals, gs,
       _ = getreducedcomponentsforspinsystem(p_cs, css,
                     J_vals, J_IDs, intermediates,
                     cs_LUT, ppm2hzfunc;
                     tol_coherence = tol_coherence,
                     prune_a_tol_factor = prune_a_tol_factor,
                     fuse_radius = fuse_radius,
                     Minkowski_parameter = Minkowski_parameter)

    # prepare intermediate objects and buffers.
    tmp1 = collect( state_labels[i][1] for i = 1:length(state_labels))
    tmp2 = collect( state_labels[i][2] for i = 1:length(state_labels))
    states::Vector{Int} = unique([tmp1; tmp2])

    return a0, F0, states, state_labels, H, H1, H2, eig_vals, gs
end



function getreducedcomponentsforspinsystem(p_cs::Vector{T},
            css::Vector{T},
            J_vals::Vector{T},
            J_IDs::Vector{Tuple{Int,Int}},
            intermediates,
            cs_LUT::Vector{Vector{Int}},
            ppm2hzfunc;
            tol_coherence = 1e-2,
            prune_a_tol_factor = 0.05,
            fuse_radius = ppm2hzfunc(0.001) - ppm2hzfunc(0.0),
            Minkowski_parameter = 3.5) where T <: Real

    #
    a0, F0, coherence_mat, eig_vals, Q, state_labels0,
    H, H1, H2, M_array = evalSCalgorithm!(p_cs, css,
    J_vals, J_IDs, intermediates,
    cs_LUT, ppm2hzfunc;
    tol_coherence = tol_coherence)

    # default: no pruning nor merging.
    F3 = F0
    a2 = a0
    state_labels = state_labels0

    if isfinite(prune_a_tol_factor)

        ### prune.
        prune_zero_tol = maximum(a0) * prune_a_tol_factor
        keep_flags = abs.(a0) .> prune_zero_tol
        a_pruned = a0[keep_flags]
        F_pruned = F0[keep_flags]
        state_labels_pruned = state_labels0[keep_flags]

        # # debug.
        # println("fuse_radius = ", fuse_radius)

        # # drop.
        # a2, F2 = dropcomponents(F, F, a; radius = fuse_radius*2*π,
        #                            Minkowski_parameter = Minkowski_parameter)

        ### merge.
        F2, a2 = approxΩsαs(F_pruned, F_pruned, a_pruned; radius = fuse_radius*2*π,
                                   Minkowski_parameter = Minkowski_parameter)

        # state labels.
        state_labels = findapproxstates(F_pruned, state_labels_pruned, F2)

        F3 = getΩfroms(eig_vals, state_labels)
    end

    # normalize intensity according to number of spins.
    N_spins = length(css)
    normalizetoNspins!(a2, N_spins)

    #TODO old possibly irrelevant note: normalize by # of hydrogens css_sys. hunt down the code for cs_singlet (probably at runtime), and do the same.

    return F3, a2, state_labels, H, H1, H2, eig_vals, Q, state_labels0, a0, F0, coherence_mat, F2
end

# #### proxy to the strong coupling.
# parsefunc converts low-dim p to dim(cs).
# - It is molecule and spin group specific.
# uses p_cs to modify cs.
# TODO: add ! to function name since this mutates cs.
function evalSCalgorithm!(  cs::Vector{T},
                            J_vals::Vector{T},
                            J_inds::Vector{Tuple{Int,Int}},
                            intermediates,
                            ppm2hzfunc;
                            tol_coherence::T = 1e-3) where T <: Real

    # set up.
    #pcstocs!(cs, p_cs, cs_LUT)

    Id, Ix, Iy_no_im, Iz, Ip, Im, Iy, Im_full, Iz_full, Ix_full,
    Iys_no_im_full = intermediates # see prepcouplingalgorithm(Int) for details.

    # strong coupling algorithm.
    ω0 = ppm2hzfunc.(cs) .* 2*π
    #H = getgenericHamiltonian(Id, Ix, Iy_no_im, Iz, ω0, J_vals, J_inds)
    H, H1, H2 = getgenericHamiltoniantest(Id, Ix, Iy_no_im, Iz, ω0, J_vals, J_inds)

    a, F, p, s, Q, coherence_labels, M_array = getaΩ(Iz_full, H,
        Iys_no_im_full, Ix_full; tol = tol_coherence)
    #println("hi")

    return a, F, p, s, Q, coherence_labels, H, H1, H2, M_array
end


# for now, assign p is Δcs, and hard code the conversion of
#   Δcs to ΔΩ0 as unit conversion for fullsolve.
# Here, updateΩfuncs takes Δcs instead of cs.
## TODO later: generic version where p reads into p_cs, then converted to css_sys.
#   the units of p_cs needs to be defined before calling this func.
function evalfullsolveproxy!(css_sys_buffer::Vector{Vector{T}},
                            cs_singlets_compact0::Vector{T},
                            #N_singlets::Int,
                            Ωs::Vector{Vector{T}},
                            p::Vector{T},
                            st_ind::Int,
                            cs_LUT::Vector{Vector{Vector{Int}}},
                            updateΩfuncs::Vector{Function},
                            fs::T,
                            SW::T,
                            ppm2hzfunc) where T <: Real
    ### checks.
    N_sys = length(css_sys_buffer)
    @assert length(updateΩfuncs) == N_sys == length(cs_LUT)

    N_singlets = length(cs_singlets_compact0)

    N_groups = N_sys + N_singlets
    @assert length(Ωs) == N_groups

    ### set up.
    fin_ind = st_ind -1

    ### multi-spin systems.
    for i = 1:N_sys

        ## parse to css_sys.
        N_cs_vars = length(cs_LUT[i])

        st_ind = fin_ind +1
        fin_ind = st_ind + N_cs_vars -1
        p_cs = p[st_ind:fin_ind]
        pcstocs!(css_sys_buffer[i], p_cs, cs_LUT[i])

        ## convert units.
        convertΔcstoΔω0!(css_sys_buffer[i], fs, SW)

        #updaterfuncs(ω0)
        updateΩfuncs[i](css_sys_buffer[i])
    end

    ### singlets.
    for i = 1:N_singlets

        ## parse.
        fin_ind += 1

        # ## convert units and update.
        # Ωs[i+N_sys][1] = p[fin_ind]*2*π*fs/SW

        # p is Δcs, and new cs = old cs + Δcs.
        tmp = cs_singlets_compact0[i] + p[fin_ind]

        ## update.
        Ωs[i+N_sys][1] = ppm2hzfunc(tmp) * 2*π
    end

    return fin_ind
end

function convertΔcstoΔω0!(x::Vector{T}, fs::T, SW::T) where T
    for i = 1:length(x)
        x[i] = x[i]*2*π*fs/SW
    end
end

function convertΔcstoΔω0(x::T, fs::T, SW::T)::T where T
    return x*2*π*fs/SW
end

### from eigen vectors to Zeeman states. See the lower case m_j^{(r)}
#   in sec. 18.2 (pag 468), Spin Dynamics.
function computequantumnumbers(basis::Vector{Vector{T}}, Id) where T

    N_spins = round(Int, log(2,length(basis)))
    I_jz_set = collect( getmultiIz(Id, j, N_spins) for j = 1:N_spins )

    m_basis = Vector{Vector{T}}(undef, length(basis))
    for r = 1:length(basis)
        m_basis[r] = collect( getzangularmomentum(basis[r], I_jz_set[j]) for j = 1:N_spins )
    end

    return m_basis
end
