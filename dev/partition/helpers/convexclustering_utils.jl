
# no error-checking on inputs.
# P must be a standard 1-indexing scheme array.
function plot2Dclusters(P::Vector{Vector{Vector{T}}},
    d1::Int, d2::Int;
    legend_position = :outertopright,
    canvas_size = (1000,1000),
    line_width = 2,
    prefix_string = "part",
    x_label = "dimension $d1",
    y_label = "dimension $d2",
    title = "Points labelled with its parts index for dimensions $d2 vs. $d1") where T <: AbstractFloat

    x1 = collect( first(P)[n][d1] for n in eachindex(first(P)) )
    x2 = collect( first(P)[n][d2] for n in eachindex(first(P)) )

    plot_handle = Plots.plot(x1, x2,
        seriestype = :scatter,
        #labels = display_labels,
        label = "$prefix_string 1",
        title =  title,
        xlabel = x_label,
        ylabel = y_label,
        linewidth = line_width, legend = legend_position, size = canvas_size)

    for k = 2:length(P)
        x1 = collect( P[k][n][d1] for n in eachindex(P[k]) )
        x2 = collect( P[k][n][d2] for n in eachindex(P[k]) )

        Plots.plot!(x1, x2, seriestype = :scatter, label = "$prefix_string $k")
    end

    return plot_handle
end

#################### default config generation.
# heuristic procedure for simulation plots without data.


"""
function evalgeometricsequence(n::Int, a0::T, r::T)::T where T <: AbstractFloat

return a0*r^n
"""
function evalgeometricsequence(n::Int, a0::T, r::T)::T where T <: AbstractFloat
    return a0*r^n
end

### non-negative-valued kernels, used as weight function.

function evalSqExpkernel(x::Vector{T}, z::Vector{T}, θ::T)::T where T <: AbstractFloat
    return exp(-θ*norm(x-z)^2)
end


function lengthscale2θ(l::T)::T where T <: AbstractFloat
    return 1/(2*l^2)
end

function θ2lengthscale(θ::T)::T where T<: AbstractFloat
    return 1/sqrt(2*θ)
end



## variable hyperparameters
function getccθγdefaultconfigs(
    Δc_m::Vector{Vector{T}}, αs::Vector{T};
    length_scale_base::T = convert(T, 10.0),
    length_scale_rate::T = convert(T, 0.7),
    min_dynamic_range::T = convert(T, 0.95),
    length_scale_max_iters::Int = 1000,
    getθfunc::Function = nn::Int->convert(T, lengthscale2θ(evalgeometricsequence(nn-1, length_scale_base, length_scale_rate))),
    γ_base::T = convert(T, 0.1),
    γ_rate::T = convert(T, 1.05),
    max_partition_size::Int = length(Δc_m[1]) + 2, # stop searching once the size of the returned partition is less than `max_partition_size`.
    max_iter_γ::Int = 100,
    getγfunc::Function = nn::Int->convert(T, evalgeometricsequence(nn-1, γ_base, γ_rate)),
    )where T <: AbstractFloat

    # weight function hyperparameter search.
    config_θ = ConvexClustering.SearchθConfigType(length_scale_max_iters, min_dynamic_range, getθfunc)

    # regularization parameter search.    
    config_γ = ConvexClustering.SearchγConfigType(max_iter_γ, max_partition_size, getγfunc)

    return config_θ, config_γ
end

# kernelfuncmust be a positive-definite RKHS kernel that does not output negative numbers.  
function getccdefaultconfigs(
    Δc_m::Vector{Vector{T}}, αs::Vector{T};
    verbose_subproblem = false,
    metric::Distances.Metric = Distances.Euclidean(),
    kernelfunc::Function = evalSqExpkernel,
    gap_tol::T = convert(T, 1e-6),
    σ_base::T = convert(T, 0.4),
    σ_rate::T = convert(T, 1.05),
    updateσfunc::Function = nn::Int->convert(T, σ_base*σ_rate^(nn-1)),
    updateϵfunc::Function = nn::Int->convert(T, 1/(nn-1)^2),
    runoptimfunc::Function = (xx,ff,dff,gg_tol)->ConvexClustering.runOptimjl(xx, ff, dff, gg_tol; verbose = verbose_subproblem),
    cc_max_iters::Int = 200,
    assignment_zero_tol::T = convert(T, 1e-6))where T <: AbstractFloat

    # constant hyperparameters to the convex clustering problem.

    # # explicit radius-based contrsuction of graph.
    # D = length(first(Δc_m))
    # radius = convert(T, norm(ones(T,D),2)/(2*D))
    # connectivity = ConvexClustering.RadiusType(radius) # make an edge for all points within this radius of a given point i, cycle through all i in the point set to be partitioned. Takes a finite floating point number.
    
    ### k-nearest neighbour-based contrsuction of graph.
    N = length(Δc_m)
    knn = max(32, round(Int, N/10))
    knn = convert(Int, min(knn, N-1)) # cannot be more than the cardinality -1 of the set to be partitioned.
    connectivity = ConvexClustering.KNNType(knn) # make an edge for this number of nearest neighbours of a given point i, cycle through all i in the point set to be partitioned. Takes a positive integer.

    graph_config = CC.WeightedGraphConfigType(connectivity, metric, kernelfunc)

    ## convex clustering optimization algorithm configuration.
    optim_config = ConvexClustering.ALMConfigType(gap_tol;
        max_iters = cc_max_iters,
        runoptimfunc = runoptimfunc,
        updateσfunc = updateσfunc,
        updateϵfunc = updateϵfunc)

    ## assignment.
    assignment_config = ConvexClustering.AssignmentConfigType(metric, assignment_zero_tol)

    return graph_config, optim_config, assignment_config
end
