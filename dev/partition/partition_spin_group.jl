# test partitioning.

# run a.jl first.

# hi(1.23) # comment/uncomment to test Revise.
# ConvexClustering.bye(4.56) # comment/uncomment to test Revise.


Random.seed!(25)

include("../../src/NMRHamiltonian.jl")
import .NMRHamiltonian
#import NMRHamiltonian

include("./helpers/convexclustering_utils.jl")

# include("../../src/partition/partition_resonances.jl")
# include("../../src/partition/convex_clustering.jl")


# verbose, trace, and stopping condition configs.
verbose_subproblem = false
report_cost = true # want to see the objective score per θ run or γ run.
store_trace = true
verbose_kernel = true


function rungetpartitionθγ(Δc_m::Vector{Vector{T}},
    αs_i_prune::Vector{T},
    γ_test::T,
    θ_test::T;
    verbose_subproblem = false,
    report_cost = true, # want to see the objective score per θ run or γ run.
    store_trace = true,
    verbose_kernel = true,
    use_config_γ = true,
    use_config_θ = true) where T

    dimensions = length(Δc_m[1])
    N_components = length(Δc_m)

    ## makre sure enough connected components are made. knn should be a function of the number of points, N_components.
    # knn = min( round(Int, N_components/2), max(30, round(Int, N_components/10)))
    # graph_config, optim_config, assignment_config = getccdefaultconfigs(Δc_m;
    #     knn = knn, verbose_subproblem = verbose_subproblem)
    # @show graph_config.knn

    graph_config, optim_config, assignment_config = getccdefaultconfigs(Δc_m, αs_i_prune)

    ## hyperparameter config.
    config_θ, config_γ = getccθγdefaultconfigs(Δc_m, αs_i_prune)
    
    if !use_config_γ
        # use γ_default, don't search γ.
        config_γ = nothing
    end
    
    if !use_config_θ
        # use θ_default, don't search θ.
        config_θ = nothing
    end

    ## 
    part_inds, θs, γs = NMRHamiltonian.getpartitionθγ(Δc_m,
        αs_i_prune, optim_config, assignment_config, config_γ, config_θ, graph_config;
        γ_default = γ_test,
        θ_default = θ_test,
        store_trace = store_trace,
        report_cost = report_cost,
        verbose_kernel = verbose_kernel)
    @show graph_config.connectivity
    println()

    iter_θ = length(θs)
    iter_γ = length(γs)

    length_scale = θ2lengthscale( config_θ.getθfunc(iter_θ) )
    γ = config_γ.getγfunc(iter_γ)
    @show (length_scale, γ)

    @show (iter_θ, iter_γ)
    @show (θs[end], γs[end])

    num_parts = length(part_inds)
    @show (N_components, dimensions, num_parts)
    println()
    println("**************")

    return part_inds, θs, γs, graph_config, optim_config, assignment_config, config_θ, config_γ
end

function runbatchgetpartitionθγ(file_paths::Vector{String},
    γ_test::T,
    θ_test::T;
    verbose_subproblem = false,
    report_cost = true, # want to see the objective score per θ run or γ run.
    store_trace = true,
    verbose_kernel = true,
    use_config_γ = true,
    use_config_θ = true) where T

    for i in eachindex(file_paths)

        file_path = file_paths[i]
        dict = JLD.load(file_path)

        Δc_m = dict["Δc_m"]
        αs_i_prune = dict["αs_i_prune"]

        dimensions = length(Δc_m[1])
        filename = split(file_path, "/")[end]
        @show filename
        @show dimensions

        rungetpartitionθγ(Δc_m, αs_i_prune, γ_test, θ_test;
            verbose_subproblem = verbose_subproblem,
            report_cost = report_cost,
            store_trace = store_trace,
            verbose_kernel = verbose_kernel,
            use_config_γ = use_config_γ,
            use_config_θ = use_config_θ)
        
        # skip visualization for batch.
    end
end


# test inputs. no visualization.
γ_test = 0.5
θ_test = 3.0
use_config_γ = true
use_config_θ = true


base_path = "./JLD/large_N"
#base_path = "./JLD"

file_paths = joinpath.(base_path, readdir(base_path))
filter!(xx->xx[end-3:end]==".jld", file_paths)

### batch.
runbatchgetpartitionθγ(file_paths, γ_test, θ_test;
    verbose_subproblem = verbose_subproblem,
    report_cost = report_cost,
    store_trace = store_trace,
    verbose_kernel = verbose_kernel,
    use_config_γ = use_config_γ,
    use_config_θ = use_config_θ)

@assert 1==2

# TODO make mixture1 example work. then package runbatchgetpartitionθγ() 
# into a stand alone update for As, given a part_config::PartitionAlgorithmType, 
#   with verbose like runbatchgetpartitionθγ(), less parameters, cleaned up.
#   use as examples. then move to signalsimulator help for visualizing the svgs.

# then, SPSA + modelfit.

### single.
file_select = 3
file_path = file_paths[file_select]

dict = JLD.load(file_path)

Δc_m = dict["Δc_m"]
αs_i_prune = dict["αs_i_prune"]

dimensions = length(Δc_m[1])
@show filename
@show dimensions

part_inds, θs, γs, graph_config, optim_config, assignment_config,
config_θ, config_γ = rungetpartitionθγ(Δc_m, αs_i_prune, γ_test, θ_test;
    verbose_subproblem = verbose_subproblem,
    report_cost = report_cost,
    store_trace = store_trace,
    verbose_kernel = verbose_kernel,
    use_config_γ = use_config_γ,
    use_config_θ = use_config_θ)

Δc_centroids = NMRHamiltonian.assemblecentroids(part_inds, αs_i_prune, Δc_m)

### visualize.
P = ConvexClustering.applyassignment(part_inds, Δc_m)

d1 = 2
d2 = min(7, length(Δc_m[1]))
plot_handle = plot2Dclusters(P, d1, d2)
display(plot_handle)



## next, # part_inds, Δc_centroids = getpartitioningfunc(Δc_m, αs_i_prune, partition_config)
# Then generate publication simulation plots. save as svg.

# setup localregistry for Utilities.jl. make guide or video on how to do so.

