"""
```
simulate(
    Phys::Vector{PhysicalParamsType{T}},
    molecule_entries::Vector{String},
    fs::T,
    SW::T,
    ν_0ppm::T,
    config::SHConfig;
) where T <: AbstractFloat
```

Simulate the resonance components' intensities and frequencies, compute order of coherences for each component, and group the components into resonance groups.

### Inputs

- `Phys` -- A `Vector` of physical chemistry parameters, e.g. the output of getphysicalparameters().
- `molecule_entries` -- A `Vector` of compound entries.
- `fs` -- the sampling frequency in Hz for use in the simulation.
- `SW` -- the spectral window in ppm for use in the simulation.
- `ν_0ppm` -- the 0 ppm peak frequency in Hz in the spectrum, for use in the simulation.
- `config` -- a configuration file of type `SHConfig`.

### Outputs

- `As::Vector{SHType{Float64}}` -- the simulated resonance intensities and frequencies, sub-system order of coherences, and resonance groups.

- `MSPs` -- the spin Hamiltonian matrices and coherence-related quantities. For diagnostic purposes.

"""
function simulate(
    Phys::Vector{PhysicalParamsType{T}},
    molecule_entries::Vector{String},
    fs::T,
    SW::T,
    ν_0ppm::T,
    config::SHConfig;
    ) where T <: AbstractFloat

    # set up.
    N_molecules = length(molecule_entries)
    @assert length(Phys) == N_molecules

    ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

    # assemble.
    As = Vector{SHType{T}}(undef, N_molecules)
    #Rs = Vector{Vector{PartitionSearchRecord{T}}}(undef, N_molecules)
    MSPs = Vector{MoleculeSpinSystem{T}}(undef, N_molecules)

    # check
    if length(unique(molecule_entries)) != length(molecule_entries)
        println("Error: Please supply a list of unique compound molecule_entries. Exit without simulation.")
        return As, MSPs
    end
    
    # simulate.
    for n in eachindex(Phys)

        αs, Ωs, parts, Δc, Δc_bar, N_spins_sys,
        αs_singlets, Ωs_singlets, N_spins_singlet, 
        MSPs[n] = setupmoleculeSH( 
            Phys[n].J_inds_sys_local,
            Phys[n].J_vals_sys,
            Phys[n].cs_sys,
            Phys[n].H_inds_singlets,
            Phys[n].cs_singlets,
            ppm2hzfunc,
            config;
            ME = Phys[n].ME
        )

        # update with singlets.
        setupsingletspinsystem!(
            αs,
            Ωs,
            Δc,
            parts,
            Δc_bar,
            N_spins_sys, # this and the above mutates.

            # here, N_spins_singlet gets appended to N_spins_sys, as with αs_singlets to αs, etc.
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
        
        #Rs[n] = search_results
    end

    return As, MSPs
end

function setupmoleculeSH(
    J_inds_sys_local::Vector{Vector{Tuple{Int64, Int64}}},
    J_vals_sys::Vector{Vector{T}},
    cs_sys::Vector{Vector{T}},
    H_inds_singlets::Vector{Vector{Int64}},
    cs_singlets::Vector{T},
    ppm2hzfunc,
    config::SHConfig;
    ME::Vector{Vector{Vector{Int}}} = Vector{Vector{Vector{Int}}}(undef, 0)
    ) where T <: AbstractFloat

    ### spin Hamiltonian simulation.
    N_spins_singlet = length.(H_inds_singlets)

    N_spins_sys = collect( length(cs_sys[m]) for m in eachindex(cs_sys) )
    intermediates_sys = prepcouplingalgorithm(T, N_spins_sys)
    
    # αs_inp, Ωs_inp, coherence_mat_sys, eig_vals_sys, Q_sys,
    # coherence_state_pairs_sys, H_sys, states_sys, ms_sys, M_sys = computeSH(
    MSP = computeSH( # MSP stands for molecule spin system.
        cs_sys,
        J_vals_sys,
        J_inds_sys_local,
        intermediates_sys,
        ppm2hzfunc,
        cs_singlets,
        N_spins_singlet,
        config,
    )
    
    αs, Ωs, parts,
    Δc, Δc_bar, c_states = partitionresonances(
        MSP.spin_systems,
        N_spins_sys;
        ME = ME,
        relative_α_threshold = config.relative_α_threshold,
    )
    
    return αs, Ωs, parts, Δc, Δc_bar, N_spins_sys,
    MSP.singlet_intensities, MSP.singlet_frequencies, N_spins_singlet, MSP
end


####################### for convinence. Based on presets from real data.

function getpresetspectrometer(::Type{T}, tag) where T <: AbstractFloat

    fs = 14005.602240896402
    SW = 20.0041938620844
    ν_0ppm = 10656.011933076665

    if tag == "700"
        # machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
        fs = 14005.602240896402
        SW = 20.0041938620844
        ν_0ppm = 10656.011933076665

    elseif tag == "600"
        ## machine values from a 600 MHz experiment: bmse000915, methionine.
        fs = 9615.38461538462
        SW = 16.022093454391
        ν_0ppm = 6685.791496181313

    elseif tag == "900"
        ## machine values from a 900 MHz experiment.
        fs = 14423.0769230769
        SW = 16.0300195009073
        ν_0ppm = 10160.027322585376

    elseif tag == "500"
        fs = 6493.50649350649
        SW = 12.9911090156122
        ν_0ppm = 4035.6644246816795

    elseif tag == "400"

        ### 400 MHz, bmse000297, ethanol.
        fs = 4807.69230769231
        SW = 12.0152693165838
        ν_0ppm = 2884.905244600881
    end

    return convert(T, fs), convert(T, SW), convert(T, ν_0ppm)
end

function loadandsimulate(
    ::Type{T},
    machine_settings_tag::String,
    molecule_entries::Vector{String},
    H_params_path,
    molecule_mapping_file_path;
    config::SHConfig{T} = SHConfig{T}(
        coherence_tol = convert(T, 0.01),
        relative_α_threshold = convert(T, 0.005),
        tol_radius_1D = convert(T, 0.1), # strictly between 0 and 1. The lower, the better the approximation, but would a larger partition (i.e. more resonance groups).
        nuc_factor = convert(T, 1.5),
    ),
    unique_cs_atol = convert(T, 1e-6),
    ) where T <: AbstractFloat

    fs, SW, ν_0ppm = getpresetspectrometer(T, machine_settings_tag)

    return loadandsimulate(
        fs,
        SW,
        ν_0ppm,
        molecule_entries::Vector{String},
        H_params_path,
        molecule_mapping_file_path,
        config;
        unique_cs_atol = unique_cs_atol,
    )
end

function loadandsimulate(
    fs::T,
    SW::T,
    ν_0ppm::T,
    molecule_entries::Vector{String},
    H_params_path::String,
    molecule_mapping_file_path::String,
    config::SHConfig{T};
    unique_cs_atol = convert(T, 1e-6),
    )::Tuple{
        Vector{PhysicalParamsType{T}},
        Vector{SHType{T}},
        Vector{MoleculeSpinSystem{T}},
    } where T <: AbstractFloat

    Phys = getphysicalparameters(
        T,
        molecule_entries,
        H_params_path,
        molecule_mapping_file_path;
        unique_cs_atol = unique_cs_atol,
    )

    As, MSPs = simulate(
        Phys,
        molecule_entries,
        fs,
        SW,
        ν_0ppm,
        config,
    )

    return Phys, As, MSPs
end