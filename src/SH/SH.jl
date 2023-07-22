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

function computeSH(
    css_sys::Vector{Vector{T}},
    J_vals_sys,
    J_inds_sys,
    intermediates_sys,
    ppm2hzfunc,
    cs_singlets,
    N_spins_singlet::Vector{Int};
    coherence_tols::Vector{T} = collect( 1e-2 for _ in eachindex(css_sys)),
    normalize_αs::Vector{Bool} = trues(length(css_sys)), # not using BitVector since the current serialization/deserialization routine has issues with it.
    ) where T <: AbstractFloat

    N_sys = length(css_sys)
    @assert length(J_vals_sys) == length(J_inds_sys) == N_sys == length(coherence_tols) == length(normalize_αs)
    #
    
    N_singlets = length(cs_singlets)

    N_groups = N_sys + N_singlets

    # non-singlet objects.
    # states_sys = Vector{Vector{Int}}(undef, N_sys)
    # coherence_state_pairs_sys = Vector{Vector{Tuple{Int,Int}}}(undef, N_sys)
    # eig_vals_sys = Vector{Vector{T}}(undef, N_sys)
    # Q_sys = Vector{Vector{Vector{T}}}(undef, N_sys)
    # H_sys = Vector{Matrix{T}}(undef, N_sys)
    # coherence_mat_sys = Vector{Matrix{T}}(undef, N_sys)
    # ms_sys = Vector{Vector{Vector{T}}}(undef, N_sys)
    # M_sys = Vector{Vector{T}}(undef, N_sys)

    # output buffers.
    

    # # Spin systems.

    sp = Vector{SpinSystem{T}}(undef, N_sys)

    for i = 1:N_sys

        intensities, frequencies, coherence_mat, eigenvalues, eigenvectors,
        coherence_state_pairs,
        H, H1, H2, Ms = evalSCalgorithm(
            css_sys[i],
            J_vals_sys[i],
            J_inds_sys[i],
            intermediates_sys[i],
            ppm2hzfunc;
            coherence_tol = coherence_tols[i],
        )

        # normalize intensities according to number of spins.
        if normalize_αs[i]
            N_spins = length(css_sys[i])
            normalizetoNspins!(intensities, N_spins)
        end

        Id = getsingleId()
        ms = computequantumnumbers(eigenvectors, Id)

        # get the unique list of all states that appear in the coherence state paris.
        tmp1 = collect( coherence_state_pairs[j][1] for j in eachindex(coherence_state_pairs) )
        tmp2 = collect( coherence_state_pairs[j][2] for j in eachindex(coherence_state_pairs) )
        states = unique([tmp1; tmp2])

        # package.
        H_contributions = Vector{Matrix{T}}(undef, 2)
        H_contributions[begin] = H1
        H_contributions[begin+1] = H2
        sp[i] = SpinSystem(
            intensities,
            frequencies,
            Hamiltonian(
                H,
                H_contributions,
                eigenvalues,
                eigenvectors,
            ),
            coherence_mat,
            coherence_state_pairs,
            states,
            ms,
            Ms,
            coherence_tols[i],
        )
    end

    # # singlets. evalRayleighproxy!() updates singlets for ΩS.
    
    αs_singlets = Vector{T}(undef, N_singlets)
    Ωs_singlets = Vector{T}(undef, N_singlets)

    for i = 1:N_singlets
        αs_singlets[i] = N_spins_singlet[i]
        Ωs_singlets[i] = ppm2hzfunc.(cs_singlets[i]) .* (2*π)
    end

    # return αs, Ωs, coherence_mat_sys, eig_vals_sys, Q_sys,
    # coherence_state_pairs_sys, H_sys, states_sys, ms_sys, M_sys
    return MoleculeSpinSystem(sp, αs_singlets, Ωs_singlets)
end


# #### proxy to the strong coupling.
# parsefunc converts low-dim p to dim(cs).
# - It is molecule and spin group specific.
# uses p_cs to modify cs.
function evalSCalgorithm(  cs::Vector{T},
                            J_vals::Vector{T},
                            J_inds::Vector{Tuple{Int,Int}},
                            intermediates,
                            ppm2hzfunc;
                            coherence_tol::T = 1e-3) where T <: AbstractFloat

    # set up.
    #pcstocs!(cs, p_cs, cs_LUT)

    Id, Ix, Iy_no_im, Iz, Ip, Im, Iy, Im_full, Iz_full, Ix_full,
    Iys_no_im_full = intermediates # see prepcouplingalgorithm(Int) for details.

    # strong coupling algorithm.
    ω0 = ppm2hzfunc.(cs) .* 2*π
    #H = getgenericHamiltonian(Id, Ix, Iy_no_im, Iz, ω0, J_vals, J_inds)
    H, H1, H2 = getgenericHamiltoniantest(Id, Ix, Iy_no_im, Iz, ω0, J_vals, J_inds)

    a, F, p, s, Q, coherence_labels, M_array = getaΩ(Iz_full, H,
        Iys_no_im_full, Ix_full; tol = coherence_tol)
    #println("hi")

    return a, F, p, s, Q, coherence_labels, H, H1, H2, M_array
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

function computequantumnumbers(A::SpinSystem{T})::Vector{Vector{T}} where T
    Id = getsingleId()
    return computequantumnumbers(A.Q, Id)
end

