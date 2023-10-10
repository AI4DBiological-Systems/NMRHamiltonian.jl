# SH simulation via auto-generated configuration from hyperparameters.

mutable struct SHConfig{T <: Real}

    coherence_tol::T # for use with defaultmixtureshsconfig(). Must be above 0.
    relative_α_threshold::T
    max_partition_size_offset::Int
    partition_rate::T
    #search_θ::Bool
    #θ_default::T
    #search_γ::Bool
    γ_base::T
    γ_rate::T
    max_iter_γ::Int
    fully_connected_convex_clustering::Bool
    max_connected_components_offset::Int
    starting_manual_knn::Int
    length_scale_base::T
    length_scale_rate::T
    min_dynamic_range::T
    cc_gap_tol::T
    cc_max_iters::Int
    assignment_zero_tol::T
    report_γ_cost::Bool
    verbose_kernel::Bool
end

"""
function SHConfig(
    ::Type{T};
    coherence_tol = 0.01,
    relative_α_threshold = 0.01,
    max_partition_size_offset::Integer = 0,
    partition_rate::Real = 2.0,
    #search_θ::Bool = true,
    #search_γ::Bool = true,
    #θ_default = 0.0,
    γ_base = 0.1,
    γ_rate = 1.05,
    max_iter_γ = 100,
    fully_connected_convex_clustering = false, #  overides all knn-related optional inputs
    max_connected_components_offset = -1,
    starting_manual_knn = 60,
    length_scale_base = 10.0,
    length_scale_rate = 0.7,
    min_dynamic_range = 0.95,
    cc_gap_tol = 1e-8,
    cc_max_iters = 300,
    assignment_zero_tol = 1e-3,
    report_γ_cost::Bool = true,
    verbose_kernel::Bool = true,
    )::SHConfig{T} where T
"""
function SHConfig(
    ::Type{T};
    coherence_tol = 0.01,
    relative_α_threshold = 0.01,
    max_partition_size_offset::Integer = 0,
    partition_rate::Real = 2.0,
    #search_θ::Bool = true,
    #search_γ::Bool = true,
    #θ_default = 0.0,
    γ_base = 0.1,
    γ_rate = 1.05,
    max_iter_γ = 100,
    fully_connected_convex_clustering = false, #  overides all knn-related optional inputs
    max_connected_components_offset = -1,
    starting_manual_knn = 60,
    length_scale_base = 10.0,
    length_scale_rate = 0.7,
    min_dynamic_range = 0.95,
    cc_gap_tol = 1e-8,
    cc_max_iters = 300,
    assignment_zero_tol = 1e-3,
    report_γ_cost::Bool = true,
    verbose_kernel::Bool = true,
    )::SHConfig{T} where T

    @assert coherence_tol > zero(T)

    return SHConfig(
        coherence_tol,
        relative_α_threshold,
        max_partition_size_offset,
        partition_rate,
        #search_θ,
        #θ_default,
        #search_γ,
        γ_base,
        γ_rate,
        max_iter_γ,
        fully_connected_convex_clustering,
        max_connected_components_offset,
        starting_manual_knn,
        length_scale_base,
        length_scale_rate,
        min_dynamic_range,
        cc_gap_tol,
        cc_max_iters,
        assignment_zero_tol,
        report_γ_cost,
        verbose_kernel,
    )
end

"""
```
simulate(
    Phys::Vector{PhysicalParamsType{T}},
    molecule_entries::Vector{String},
    fs::T,
    SW::T,
    ν_0ppm::T,
    config::SHConfig;

    θs::Vector{Vector{T}} = defaultspinsysbuffer(convert(T,Inf), Phys), # non-finite value means use search strategy.
    γs::Vector{Vector{T}} = defaultspinsysbuffer(convert(T,Inf), Phys),

    γ_bases = defaultspinsysbuffer(config.γ_base, Phys),
    γ_rates = defaultspinsysbuffer(config.γ_rate, Phys),
    γ_max_iters = defaultspinsysbuffer(config.max_iter_γ, Phys),

    length_scale_bases = defaultspinsysbuffer(config.length_scale_base, Phys),
    length_scale_rates = defaultspinsysbuffer(config.length_scale_rate, Phys),
    min_dynamic_ranges = defaultspinsysbuffer(config.min_dynamic_range, Phys),
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

### Optional inputs

- `θs` -- This is a nested array for manually specifying the weight kernel bandwidth parameter (positive) for the convex clustering weights. The kernel is a square exponential kernel, i.e., `exp(-θ*norm(X[i]-X[j])^2)` for the edge weight between data points `X[i]` and `X[j]`. θs[n][i] is the positive weight kernel bandwidth of your choice for the n-th compound entry, i-th spin system. Use the default (i.e., leave this optional input alone) if you want to use the automated search strategy for the kernel bandwidth.

- `γs` -- The sizing of this nested array is the same as θs, but for the sparsity-inducing parameter (positive, smaller means more resonance groups) for the convex clustering optimization problem. 

- `γ_bases` -- This needs to be assigned to `defaultspinsysbuffer(config.γ_base, Phys),`, but replace `config.γ_base` with your choice of the starting γ parameter. For use with the automated search strategy.

- `γ_rates` -- This needs to be assigned to `defaultspinsysbuffer(config.γ_rate, Phys),`, but replace `config.γ_base` with your choice of the rate of increase for the next candidate γ parameter. For use with the automated search strategy. We've `γ_next = γ_bases*γ_rate^(iter-1)`, where iter is the current search iteration. For use with the automated search strategy.

- `γ_max_iters` -- This is similar to θs, but for the maximum allowed iterations in the automated search strategy for γ.

- `length_scale_bases` and `length_scale_rates` are similar to `γ_bases` and `γ_rates`, but for length_scale = 1/(2*θ).

- `min_dynamic_ranges` -- The sizing of this nested array is the same as θs, but contain the stopping condition dynamic range value for which the kernel bandwidth automated search strategy terminates upon.

### Outputs

- `As::Vector{SHType{Float64}}` -- the simulated resonance intensities and frequencies, sub-system order of coherences, and resonance groups.

- `Rs` -- the partition tree used to generate the resonance groups. For diagnostic purposes.

- `MSPs` -- the spin Hamiltonian matrices and coherence-related quantities. For diagnostic purposes.

"""
function simulate(
    Phys::Vector{PhysicalParamsType{T}},
    molecule_entries::Vector{String},
    fs::T,
    SW::T,
    ν_0ppm::T,
    config::SHConfig;

    θs::Vector{Vector{T}} = defaultspinsysbuffer(convert(T,Inf), Phys), # non-finite value means use search strategy.
    γs::Vector{Vector{T}} = defaultspinsysbuffer(convert(T,Inf), Phys),

    γ_bases = defaultspinsysbuffer(config.γ_base, Phys),
    γ_rates = defaultspinsysbuffer(config.γ_rate, Phys),
    γ_max_iters = defaultspinsysbuffer(config.max_iter_γ, Phys),

    length_scale_bases = defaultspinsysbuffer(config.length_scale_base, Phys),
    length_scale_rates = defaultspinsysbuffer(config.length_scale_rate, Phys),
    min_dynamic_ranges = defaultspinsysbuffer(config.min_dynamic_range, Phys),
    ) where T <: AbstractFloat

    part_algs, mixture_sh_config = generateconfigs(
        Phys,
        molecule_entries,
        config;
        θs = θs,
        γs = γs,
        γ_bases = γ_bases,
        γ_rates = γ_rates,
        γ_max_iters = γ_max_iters,

        length_scale_bases = length_scale_bases,
        length_scale_rates = length_scale_rates,
        min_dynamic_ranges = min_dynamic_ranges,
    )
    
    As, Rs, MSPs = simulate(
        part_algs,
        molecule_entries,
        fs,
        SW,
        ν_0ppm,
        Phys,
        mixture_sh_config,
    )

    return As, Rs, MSPs
end

# creates a 2-level nested array, where each entry is assigned `content`.
function defaultspinsysbuffer(
    content::OT,
    Phys::Vector{PhysicalParamsType{T}},
    )::Vector{Vector{OT}} where {T,OT}

    out = collect(
        collect(
            content for _ in eachindex(Phys[n].cs_sys)
        ) for n in eachindex(Phys)
    )

    return out
end

# if θs values are not finite, then simulate() uses the search strategy.
function generateconfigs(
    Phys::Vector{PhysicalParamsType{T}},
    molecule_entries::Vector{String},
    config::SHConfig;

    θs::Vector{Vector{T}} = defaultspinsysbuffer(convert(T,Inf), Phys), # non-finite value means use search strategy.
    γs::Vector{Vector{T}} = defaultspinsysbuffer(convert(T,Inf), Phys),

    γ_bases = defaultspinsysbuffer(config.γ_base, Phys),
    γ_rates = defaultspinsysbuffer(config.γ_rate, Phys),
    γ_max_iters = defaultspinsysbuffer(config.max_iter_γ, Phys),

    length_scale_bases = defaultspinsysbuffer(config.length_scale_base, Phys),
    length_scale_rates = defaultspinsysbuffer(config.length_scale_rate, Phys),
    min_dynamic_ranges = defaultspinsysbuffer(config.min_dynamic_range, Phys),
    ) where T <: AbstractFloat

    # extract chemical shifts for spin systems and singlets.
    cs_sys_mixture, cs_singlets_mixture = extractcs(Phys)

    # if using default tolerances for coherence and intensity thresholding.
    mixture_sh_config = defaultmixtureshsconfig(
        cs_sys_mixture;
        coherence_tol = config.coherence_tol,
        relative_α_threshold = config.relative_α_threshold,
    )
   
    searchknnconfigfunc = defaultknnsearchconfig
    if config.fully_connected_convex_clustering
        searchknnconfigfunc = (nn, ii, cc, aa)->defaultknnsearchconfig(
            nn, ii, cc, aa;
            verbose = true,
            start_knn = length(cc),
        )
    end

    name_to_ind = nn->findfirst(xx->xx==nn, molecule_entries)

    getsearchθconfigfunc = (nn, ii, cc, aa)->createsearchθconfigs(
        nn, ii, cc, aa;
        length_scale_base = length_scale_bases[name_to_ind(nn)][ii],
        length_scale_rate = length_scale_rates[name_to_ind(nn)][ii],
        min_dynamic_range = min_dynamic_ranges[name_to_ind(nn)][ii],
    )
    searchγconfigfunc = (nn, ii, cc, aa)->createsearchγconfigs(
        nn, ii, cc, aa;
        max_partition_size_offset = config.max_partition_size_offset,
        partition_rate = config.partition_rate,
        γ_base = γ_bases[name_to_ind(nn)][ii],
        γ_rate = γ_rates[name_to_ind(nn)][ii],
        max_iter_γ = γ_max_iters[name_to_ind(nn)][ii],
    )

    getassignmentfunc = (nn, ii, cc, aa)->defeaultassignmentconfig(
        nn, ii, cc, aa;
        assignment_zero_tol = config.assignment_zero_tol,
    )

    part_algs = generatemixturepartitionalgorithm(
        molecule_entries,
        θs,
        γs,
        Phys;
        getgraphconfigfunc = searchknnconfigfunc,
        getsearchθconfigfunc = getsearchθconfigfunc,
        getsearchγconfigfunc = searchγconfigfunc,
        getassignmentfunc = getassignmentfunc,
        report_cost = config.report_γ_cost,
        verbose_kernel = config.verbose_kernel,
        gap_tol = config.cc_gap_tol,
        max_iters = config.cc_max_iters,
    )

    return part_algs, mixture_sh_config
    # As, Rs, MSPs = simulate(
    #     part_algs,
    #     molecule_entries,
    #     fs, SW, ν_0ppm,
    #     Phys,
    #     mixture_sh_config,
    # )
    #return As, Rs, MSPs
end
