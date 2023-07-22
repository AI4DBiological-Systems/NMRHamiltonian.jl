"""
```
simulate(
    part_algs::Vector{Vector{PartitionAlgorithmType}},
    names::Vector{String},
    fs::T,
    SW::T,
    ν_0ppm::T,
    Phys::Vector{PhysicalParamsType{T}},
    configs::Vector{SHsConfigType{T}}
    ) where T <: AbstractFloat
```

Performs spin hamiltonian simulation and resonance group partitioning of the molecule entries in `names`.

### Inputs

- `part_algs::Vector{Vector{PartitionAlgorithmType}}`
Element `[n][i]` is the partition algorithm config molecule entry `n`, non-singlet spin system `i`.

- `names::Vector{String}`
List of molecule entries.

- `fs`$DOCSTRING_fs

- `SW`$DOCSTRING_SW

- `ν_0ppm`$DOCSTRING_ν_0ppm

- `Phys::Vector{PhysicalParamsType{T}}`
List of chemical shift and J-coupling information for each molecule entry.

- `configs::Vector{SHsConfigType{T}}`
List of spin Hamiltonian simulation configurations for each molecule entry.

### Outputs

- `As::SHType{T}`
A list of spin Hamiltonian result data structure variables. 

- `Rs::Vector{Vector{PartitionSearchRecord}}`
A list of iterative search results for the partitioning algorithm. This variable is only useful if an iterative search was done for regulariation parameter γ or/and kernel hyperparameter θ.

"""
function simulate(
    part_algs::Vector{Vector{PartitionAlgorithmType}},
    names::Vector{String},
    fs::T,
    SW::T,
    ν_0ppm::T,
    Phys::Vector{PhysicalParamsType{T}},
    configs::Vector{SHsConfigType{T}}
    ) where T <: AbstractFloat



    # set up.
    N_molecules = length(names)
    @assert length(part_algs) == length(configs) == length(Phys) == N_molecules

    ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

    # assemble.
    As = Vector{SHType{T}}(undef, N_molecules)
    Rs = Vector{Vector{PartitionSearchRecord{T}}}(undef, N_molecules)
    MSPs = Vector{MoleculeSpinSystem{T}}(undef, N_molecules)

    # check
    if length(unique(names)) != length(names)
        println("Error: Please supply a list of unique compound names. Exit without simulation.")
        return As, Rs, MSPs
    end
    
    # simulate.
    for n in eachindex(Phys)

        αs, Ωs, parts, Δc, Δc_bar, N_spins_sys,
            αs_singlets, Ωs_singlets, N_spins_singlet, search_results,
            MSPs[n] = setupmoleculeSH( part_algs[n],
            names[n],
            Phys[n].J_inds_sys_local,
            Phys[n].J_vals_sys,
            Phys[n].cs_sys,
            Phys[n].H_inds_singlets,
            Phys[n].cs_singlets,
            ppm2hzfunc,
            configs[n];
            ME = Phys[n].ME)

        # update with singlets.
        setupsingletspinsystem!(
            αs,
            Ωs,
            Δc,
            parts,
            Δc_bar,
            N_spins_sys,
            αs_singlets,
            Ωs_singlets,
            N_spins_singlet,
        )

        # package up.
        As[n] = SHType(
            αs,
            Ωs,
            Δc,
            parts,
            Δc_bar,
            N_spins_sys,
            #αs_singlets, Ωs_singlets,
            fs, SW, ν_0ppm
        )
        
        Rs[n] = search_results
    end

    return As, Rs, MSPs
end

function setupmoleculeSH(
    part_algs::Vector{PartitionAlgorithmType},
    name,
    J_inds_sys_local,
    J_vals_sys,
    cs_sys,
    H_inds_singlets,
    cs_singlets,
    ppm2hzfunc,
    config::SHsConfigType{T};
    ME::Vector{Vector{Vector{Int}}} = Vector{Vector{Vector{Int}}}(undef, 0)
    ) where T <: AbstractFloat

    ### spin Hamiltonian simulation.
    N_spins_singlet = length.(H_inds_singlets)

    N_spins_sys = collect( length(cs_sys[m]) for m in eachindex(cs_sys) )
    intermediates_sys = prepcouplingalgorithm(N_spins_sys)

    # αs_inp, Ωs_inp, coherence_mat_sys, eig_vals_sys, Q_sys,
    # coherence_state_pairs_sys, H_sys, states_sys, ms_sys, M_sys = computeSH(
    MSP = computeSH( # MSP stands for molecule spin system.
        cs_sys,
        J_vals_sys, J_inds_sys_local, intermediates_sys,
        ppm2hzfunc, cs_singlets, N_spins_singlet;
        coherence_tols = config.coherence_tols,
        normalize_αs = config.normalize_αs,
    )
    
    # # assemble intensities α and frequencies Ω for singlets and non-singlet spin systems.
    # k = findfirst(xx->length(xx)==1, MSP.αs_inp)
    # αs_spin_sys = copy(MSP.αs_inp)
    # Ωs_spin_sys = copy(MSP.Ωs_inp)
    # αs_singlets = Vector{T}(undef, 0)
    # Ωs_singlets = Vector{T}(undef, 0)
    # if typeof(k) == Int # check if there are singlet groups in this molecule.

    #     αs_spin_sys = αs_spin_sys[1:k-1]
    #     Ωs_spin_sys = Ωs_spin_sys[1:k-1]

    #     αs_singlets = collect( MSP.αs_inp[l][1] for l = k:length(MSP.αs_inp) )
    #     Ωs_singlets = collect( MSP.Ωs_inp[l][1] for l = k:length(MSP.Ωs_inp) )
    # end

    # partition resonance components into resonance groups.
    #N_spins_sys = collect( length(cs_sys[i]) for i in eachindex(cs_sys))

    αs, Ωs, parts,
    Δc, Δc_bar, c_states, search_results = partitionresonances(
        name,
        part_algs,
        # coherence_state_pairs_sys,
        # ms_sys,
        # αs_spin_sys,
        # Ωs_spin_sys,
        MSP.spin_systems,
        N_spins_sys;
        ME = ME,
        relative_α_thresholds = config.relative_α_thresholds)

    return αs, Ωs, parts, Δc, Δc_bar, N_spins_sys,
    MSP.singlet_intensities, MSP.singlet_frequencies, N_spins_singlet, search_results, MSP
end

