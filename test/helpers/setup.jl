
# taken from /examples/simulate_SH.jl
function setupSH(
    molecule_entries::Vector{String},
    H_params_path,
    molecule_mapping_file_path;
    machine_settings_tag = "700",
    γ_base = 0.1, # larger means less parts in the partition.
    γ_rate = 1.1,
    max_iter_γ = 10,
    length_scale_base = 10.0,
    length_scale_rate = 0.7,
    θ_default = 1.0,
    relative_α_threshold = 0.01,
    )

    fs, SW, ν_0ppm = fetchsamplemachinesettings(machine_settings_tag)

    ## pull the sample coupling values into dictionary structures.
    
    #root_data_path = getdatapath() # coupling values data repository root path
    #H_params_path = joinpath(root_data_path, "coupling_info") # folder of coupling values. # replace with your own values in actual usage.

    ## doesn't produce partition tree due to numerical accuracy. 
    # cc_gap_tol = 1e-6
    # cc_max_iters = 300
    # assignment_zero_tol = 1e-6

    # These produce a aprtition tree over the searched γ values for many molecule entries.
    cc_gap_tol = 1e-8
    cc_max_iters = 300
    assignment_zero_tol = 1e-3

    max_partition_size_offset = 3
    min_dynamic_range = 0.95

    ### end inputs.

    Phys = NMRHamiltonian.getphysicalparameters(
        Float64,
        molecule_entries,
        H_params_path,
        molecule_mapping_file_path;
        unique_cs_atol = 1e-6)

    #
    cs_sys_mixture, cs_singlets_mixture = NMRHamiltonian.extractcs(Phys)

    # if using default tolerances for coherence and intensity thresholding.
    mixture_sh_config = NMRHamiltonian.defaultmixtureshsconfig(
        cs_sys_mixture;
        relative_α_threshold = relative_α_threshold,
    )
    mixture_parts_params = NMRHamiltonian.defaultmixturepartitionsparameters(cs_sys_mixture, θ_default)
    # # if loading from file.
    # file_filder = "./configs"
    # mixture_sh_config = NMRHamiltonian.loadmixtureshsconfig(
    #     cs_sys_mixture,
    #     joinpath(file_folder, "mixture_SH.json"),
    #     molecule_entries,
    #     Float64,
    # )
    # mixture_parts_params = NMRHamiltonian.loadmixturepartitionsparameters(
    #     cs_sys_mixture,
    #     joinpath(file_folder, "mixture_partition.json"),
    #     molecule_entries,
    #     2.0
    # )

    # # View imported chemical shifts.
    # import Dictionaries
    # cs_systems = Dictionaries.ArrayDictionary(Symbol.(molecule_entries), cs_sys_mixture)
    # cs_singlets = Dictionaries.ArrayDictionary(Symbol.(molecule_entries), cs_sys_mixture)

    # dict = Dictionaries.ArrayDictionary(Symbol.(molecule_entries), mixture_parts_params)
    #Cs = collect( As[n].Δc for n in eachindex(As) )

    ### simulate spin Hamiltonian, and partition resonance components into resonance groups.


    constantknnfunc, constantradiusfunc, θs,
        γs = NMRHamiltonian.setupconstantparameteroptions(
            molecule_entries, mixture_parts_params)

    # TODO: try every combination of the following 3 funcs. manual, no search, all default funcs.
    # verify manual actually passes the right values.
    # - add examples to post-process and explore data using TypedTables or Dictionaries.
    #   see edit_tools.jl

    # getgraphconfigfunc can be one of the following.:
    #   - defaultknnsearchconfig
    #   - defaultknnconfig
    #   - defaultradiusconfig
    #   - defaultradiussearchconfig
    #   - constantknnfunc
    #   - constantradiusfunc

    # set to negative large number to force max_connected_components to be 1.
    #max_connected_components_offset = -1
    #start_knn = 60

    searchknnconfigfunc = (nn, ii, cc, aa)->NMRHamiltonian.defaultknnsearchconfig(
        nn, ii, cc, aa;
        #start_knn = start_knn,
        #start_knn = round(Int, length(cc)*0.2),
        #start_knn = length(cc)-1, # fully-connected graph.
        start_knn = max(60, round(Int, length(cc)*0.05)),
        max_knn = max(60, round(Int, length(cc)*0.2)),
        max_connected_components_offset = -1,
        verbose = true,
    )

    getsearchθconfigfunc = (nn, ii, cc, aa)->NMRHamiltonian.createsearchθconfigs(
        nn, ii, cc, aa;
        length_scale_base = length_scale_base,
        length_scale_rate = length_scale_rate,
        min_dynamic_range =min_dynamic_range,
    )


    
    searchγconfigfunc = (nn, ii, cc, aa)->NMRHamiltonian.createsearchγconfigs(
        nn, ii, cc, aa;
        max_partition_size_offset = max_partition_size_offset,
        γ_base = γ_base,
        γ_rate = γ_rate,
        max_iter_γ = max_iter_γ,
    )


    getassignmentfunc = (nn, ii, cc, aa)->NMRHamiltonian.defeaultassignmentconfig(
        nn, ii, cc, aa;
        assignment_zero_tol = assignment_zero_tol,
    )

    part_algs = NMRHamiltonian.generatemixturepartitionalgorithm(
        molecule_entries,
        θs,
        γs,
        Phys;
        #getgraphconfigfunc = constantradiusfunc,
        #getgraphconfigfunc = NMRHamiltonian.defaultknnsearchconfig,
        getgraphconfigfunc = searchknnconfigfunc,
        getsearchθconfigfunc = getsearchθconfigfunc,
        #getsearchγconfigfunc = NMRHamiltonian.createsearchγconfigs,
        getsearchγconfigfunc = searchγconfigfunc,
        getassignmentfunc = getassignmentfunc,
        #store_trace = false,
        report_cost = true,
        verbose_kernel = true,

        gap_tol = cc_gap_tol,
        max_iters = cc_max_iters,
    )

    #@assert 1==2

    # import JLD
    # freq_tag = round(Int, fs/SW)
    # global save_prefix = "$(molecule_entries[1])_$freq_tag"
    # save_prefix = replace(save_prefix, " "=> "", "," => "") # remove spaces and commas.
    # getpartitioningfunc = (xx,yy,ii)->JLD.save("$(save_prefix)_$(ii).jld", "Δc_m", xx, "αs_i_prune", yy)

    println("Timing: simulate()")
    @time As, Rs, MSPs = NMRHamiltonian.simulate(
        part_algs,
        molecule_entries,
        fs, SW, ν_0ppm,
        Phys,
        mixture_sh_config)

    return Phys, As, Rs, MSPs
end