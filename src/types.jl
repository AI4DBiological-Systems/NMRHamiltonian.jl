##### confiuration

# molecule (multiple spin systems).
"""
```
@kwdef struct SHConfig{T <: Real}
    coherence_tol::T = convert(T, 0.01)
    relative_α_threshold::T = convert(T, 0.005) # initial intensity pruning.
    tol_radius_1D::T = convert(T, 0.1) # strictly between 0 and 1. The lower, the better the approximation, but would a larger partition (i.e. more resonance groups).
    max_deviation_from_mean::T = convert(T, 0.2)
    acceptance_factor::T = convert(T, 0.99)
    total_α_threshold::T = convert(T, 0.01) # final intensity pruning.
end
```

Configuration settings for the spin Hamiltonian simulation portion of `simulate()`.

The parameters here control the trade-off between the number of coherence difference features, the set of Δc's, and the number of resonance groups, the set of Δc_bar's, that are generated.
The more resonance groups, the better the approximation to the ideal model where each Δc is a resonance group.
The resultant model becomes more computationally intensive as the number of resonance groups increase.

≡≡≡≡≡≡≡≡ Details ≡≡≡≡≡≡≡≡

### `coherence_tol::T`
This controls the number of Δc features we get, by how far their sum deviate from `-1`.

### `relative_α_threshold::T`
This controls the number of Δc feature we get, by how intense their corresponding resonance intensity is relative to the maximum intensity within their local spin system.

### `max_deviation_from_mean::T`
This non-negative controls the number of resonance groups we get. The smaller, the more resonance groups, but better approximation to the preliminary model.

### `acceptance_factor::T`
This non-negative controls the number of resonance groups we get. Takes on value strictly between 0 to 1. Tends to create less resonance groups when it is close to 1, but takes longer to compute.

### `total_α_threshold::T`
This non-negative controls the number of resonance groups we get. This dictates the cumulative intensity threshold for resonance groups, as a proportionality factor of the total intensity of the spin system. The least intense groups that cumulatively add to less than this value is discarded, and αs and Ωs and Δc are updated according to the kept components. Takes on value from [0, 1). When set to zero, does not prune any resonance groups of low intensity during the post-processing step.


"""
@kwdef struct SHConfig{T <: Real}

    coherence_tol::T = convert(T, 0.01)
    relative_α_threshold::T = convert(T, 0.005) # initial intensity pruning.
    tol_radius_1D::T = convert(T, 0.1) # strictly between 0 and 1. The lower, the better the approximation, but would a larger partition (i.e. more resonance groups).
    max_deviation_from_mean::T = convert(T, 0.2)
    acceptance_factor::T = convert(T, 0.99)
    total_α_threshold::T = convert(T, 0.01) # final intensity pruning.
end

# multiple spin systems. Partition


abstract type CentroidTrait end
#struct WeightedCentroids <: CentroidTrait end # won't have sum(Δc_bar) ≈ -1.
struct AveragedCentroids <: CentroidTrait end
struct MaxAmplitudeCentroids <: CentroidTrait end # in used. the safest option, but yields non-differentiable clustering hyperparamets. Note that convex clustering hyperparameters aren't differetiable already.

######
"""
```
struct SHType{T}

    # resonance components for each spin systems (indexed by first/outer index).
    αs::Vector{Vector{T}}
    Ωs::Vector{Vector{T}}

    Δc::Vector{Vector{Vector{T}}}
    parts::Vector{Vector{Vector{Int}}}
    Δc_bar::Vector{Vector{Vector{T}}}
    cs_Δc::Vector{Vector{T}}

    # the number of spins for each spin system.
    N_spins_sys::Vector{Int}

    # misc input for the SH simulation.
    fs::T
    SW::T
    ν_0ppm::T
end
```

Container for the spin Hamiltonian simulation result of a molecule entry. All arrays follow the 1-indexing scheme.

Frequency information about the NMR data and spectrometer can be obtained from these information.
$DOCSTRING_spectrometer_freq
$DOCSTRING_hz2ppm_ppm2hz

Fields
≡≡≡≡≡≡≡≡
- `αs`$(DOCSTRING_SHType_αs("T"))
- `Ωs`$(DOCSTRING_SHType_Ωs("T"))
- `Δc`$(DOCSTRING_SHType_Δc("T"))
- `parts`$(DOCSTRING_SHType_parts())
- `Δc_bar`$(DOCSTRING_SHType_Δc_bar("T"))
- `cs_Δc`: the reference chemical shift values, one for each dimension of a Δc vector.
- `N_spins_sys`$(DOCSTRING_SHType_N_spins_sys())
- `αs_singlets`$(DOCSTRING_SHType_αs_singlets("T"))
- `Ωs_singlets`$(DOCSTRING_SHType_Ωs_singlets("T"))
- `fs`$(DOCSTRING_fs)
- `SW`$(DOCSTRING_SW)
- `ν_0ppm`$(DOCSTRING_ν_0ppm)
"""
struct SHType{T} # output of the SH simulation.

    # resonance components in spin systems.
    αs::Vector{Vector{T}}
    Ωs::Vector{Vector{T}}

    Δc::Vector{Vector{Vector{T}}} # every resonance frequency forms its own singleton part, i.e. the set to be partitioned.
    parts::Vector{Vector{Vector{Int}}}
    Δc_bar::Vector{Vector{Vector{T}}} # multiple resonance frequencies form a part.
    cs_Δc::Vector{Vector{T}}
    
    N_spins_sys::Vector{Int}

    # resonance components in singlets spin systems.
    #αs_singlets::Vector{T}
    #Ωs_singlets::Vector{T}

    # misc input for the SH simulation.
    fs::T
    SW::T
    ν_0ppm::T
end



###
"""
```
struct PhysicalParamsType{T}
- H_IDs::Vector{Int}
- H_inds_sys::Vector{Vector{Int}}
- cs_sys::Vector{Vector{T}}
- H_inds_singlets::Vector{Vector{Int}}
- cs_singlets::Vector{T}
- J_inds_sys::Vector{Vector{Tuple{Int,Int}}}
- J_inds_sys_local::Vector{Vector{Tuple{Int,Int}}}
- J_vals_sys::Vector{Vector{T}}
- ME::Vector{Vector{Vector{Int}}}
```

The chemical shift and J-coupling values with respect to nuclei labels for a molecule. Based on this information, an algorithm is used to create magnetic equivalence information in `ME`.

Details
≡≡≡≡≡≡≡≡

- `H_IDs`$(DOCSTRING_PhysicalParamsType_H_IDs)
- `H_inds_sys`$(DOCSTRING_PhysicalParamsType_H_inds_sys)
- `cs_sys`$(DOCSTRING_PhysicalParamsType_cs_sys("T"))
- `H_inds_singlets`$(DOCSTRING_PhysicalParamsType_H_inds_singlets)
- `cs_singlets`$(DOCSTRING_PhysicalParamsType_cs_singlet("T"))
- `J_inds_sys`$(DOCSTRING_PhysicalParamsType_J_inds_sys)
- `J_inds_sys_local`$(DOCSTRING_PhysicalParamsType_J_inds_sys_local)
- `J_vals_sys`$(DOCSTRING_PhysicalParamsType_J_vals_sys("T"))
- `ME`$(DOCSTRING_PhysicalParamsType_ME)
"""
struct PhysicalParamsType{T}

    H_IDs::Vector{Int}

    H_inds_sys::Vector{Vector{Int}}
    cs_sys::Vector{Vector{T}}

    H_inds_singlets::Vector{Vector{Int}}
    cs_singlets::Vector{T}

    J_inds_sys::Vector{Vector{Tuple{Int,Int}}}
    J_inds_sys_local::Vector{Vector{Tuple{Int,Int}}} # redundant?
    J_vals_sys::Vector{Vector{T}}

    ME::Vector{Vector{Vector{Int}}}
end

function getNnonsingletspinsystems(A::PhysicalParamsType)::Int
    N_spins_sys = collect( length(A.cs_sys[i]) for i in eachindex(A.cs_sys))
    return N_spins_sys
end

############# spin Hamiltonian. see SH.jl and front_end.jl.

struct Hamiltonian{T}
    matrix::Matrix{T} # sum(order). redundent representation, but keep it as a cache.
    contributions::Vector{Matrix{T}} # the i-th entry is the additive contribution to H from the i-th order terms. i.e., [1] are due to chem shifts only, [2] are due to J-coupling only.

    eigenvalues::Vector{T}
    eigenvectors::Vector{Vector{T}}
end

# output of computeSH(), no singlets.
struct SpinSystem{T}

    intensities::Vector{T}
    frequencies::Vector{T} # in radians.

    # The rest are not in SHtype. Does not include singlets.
    H::Hamiltonian{T}
    coherence_mat::Matrix{T}
    coherence_state_pairs::Vector{Tuple{Int,Int}}
    
    states::Vector{Int} # the unique states that contribute to the signal
    partial_quantum_numbers::Vector{Vector{T}} # [state pair][nuclei] partial quantum numbers. Sum to full quantum numbers.
    quantum_numbers::Vector{T} # full quantum numbers, e.g. coherence drop.  getorderofcoherence(). Redundant since should be sum of ms. keep for error-checking.

    coherence_tol::T
end

struct MoleculeSpinSystem{T}
    spin_systems::Vector{SpinSystem{T}}

    # non-singlets
    singlet_intensities::Vector{T}
    singlet_frequencies::Vector{T} # in radians.
end


###### internals, IO physical parameters.

struct CSJContainer{T}
    J_inds_sys::Vector{Vector{Tuple{Int,Int}}}
    J_inds_sys_local
    J_IDs_sys::Vector{Vector{Tuple{Int,Int}}}

    J_vals_sys::Vector{Vector{T}}
    H_inds_sys::Vector{Vector{Int}}

    cs_sys::Vector{Vector{T}}

    H_inds_singlets::Vector{Vector{Int64}}
    cs_singlets::Vector{T}
    H_inds::Vector{Int64}
    J_inds::Vector{Tuple{Int64, Int64}}
end