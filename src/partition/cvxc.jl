


# αs not used.
# no error-check on knn vs. length(Δc_m)
# mutates graph_config.connectivity if searching for knn or radius.
function getpartitionθγ!(
    graph_config::CC.WeightedGraphConfigType,
    Δc_m::Vector{Vector{T}},
    αs::Vector{T},
    optim_config::CC.ALMConfigType,
    assignment_config::CC.AssignmentConfigType,
    config_γ::Union{CC.SearchγConfigType, Nothing},
    config_θ::Union{CC.SearchθConfigType, Nothing};
    γ_default::T = one(T),
    θ_default::T = one(T),
    store_trace::Bool = false,
    report_cost::Bool = false,
    verbose_kernel::Bool = false,
    ) where T <: AbstractFloat

    @assert length(αs) == length(Δc_m)
    @assert length(Δc_m) > 1

    # println("pre")
    # @show graph_config.connectivity
    # @show typeof(graph_config.connectivity)
    # println()

    ### setup convex clustering problem.
    A, edge_pairs, w, θs = getselectedgraph(
        config_θ,
        Δc_m,
        graph_config;
        θ_default = θ_default,
        verbose_kernel = verbose_kernel,
    )

    connectivity_label = :nothing
    connectivity_parameter = :nothing

    # .connectivity can be of the following types: KNNType, RadiusType, KNNSearchType, RadiusSearchType. See the subtypes of ConvexClustering.AbstractConnectivityType.
    if typeof(graph_config.connectivity) <: ConvexClustering.KNNType || typeof(graph_config.connectivity) <: ConvexClustering.KNNSearchType
         
        connectivity_label = :knn
        connectivity_parameter = Symbol(string(graph_config.connectivity.knn))

        if verbose_kernel
            println("Connectivity: knn, ", graph_config.connectivity.knn)
        end
    elseif typeof(graph_config.connectivity) <: ConvexClustering.RadiusType || typeof(graph_config.connectivity) <: ConvexClustering.RadiusSearchType
        
        connectivity_label = :radius
        connectivity_parameter = Symbol(string(graph_config.connectivity.radius))

        if verbose_kernel
            println("Connectivity: radius, ", graph_config.connectivity.radius)
        end
    end

    ### run search.
    Gs, rets, γs = solveselectedproblem(
        config_γ,
        A,
        w,
        edge_pairs,
        optim_config, assignment_config;
        X_initial = zeros(T, size(A)),
        Z_initial = zeros(T, size(A,1), length(edge_pairs)),
        γ_default = γ_default,
        store_trace = store_trace,
        report_cost = report_cost,
    )
    
    return Gs, θs, γs, connectivity_label, connectivity_parameter
end

function solveselectedproblem(
    config_γ::CC.SearchγConfigType,
    A::Matrix{T},
    w::Vector{T},
    edge_pairs::Vector{Tuple{Int,Int}},
    optim_config::CC.ALMConfigType{T},
    assignment_config::CC.AssignmentConfigType;
    X_initial::Matrix{T} = zeros(T, size(A)),
    Z_initial::Matrix{T} = zeros(T, size(A,1), length(edge_pairs)),
    γ_default::T = one(T),
    store_trace::Bool = false,
    store_trace_assignments = true,
    report_cost::Bool = false) where T <: AbstractFloat

    # The γ field can be set to anything because it will get modified before use by config_γ.getγfunc() in searchγ().
    problem = ConvexClustering.ProblemType(A, NaN, w, edge_pairs)

    Gs, rets, γs = ConvexClustering.searchγ(
    X_initial, Z_initial, problem, optim_config, assignment_config, config_γ;
    store_trace = store_trace,
    store_trace_assignments = store_trace_assignments,
    report_cost = report_cost)

    return Gs, rets, γs
end

function solveselectedproblem(config_γ::Nothing,
    A::Matrix{T}, w::Vector{T}, edge_pairs::Vector{Tuple{Int,Int}},
    optim_config::CC.ALMConfigType{T},
    assignment_config::CC.AssignmentConfigType;
    X_initial::Matrix{T} = zeros(T, size(A)),
    Z_initial::Matrix{T} = zeros(T, size(A,1), length(edge_pairs)),
    γ_default::T = one(T),
    store_trace::Bool = false,
    store_trace_assignments = false,
    report_cost::Bool = false) where T <: AbstractFloat

    @assert isfinite(γ_default)
    @assert γ_default > zero(T)

    problem = ConvexClustering.ProblemType(A, γ_default, w, edge_pairs)

    G, ret = ConvexClustering.runconvexclustering(X_initial, Z_initial,
        problem, optim_config, assignment_config;
        store_trace = store_trace,
        report_cost = report_cost)

    Gs = Vector{Vector{Vector{Int}}}(undef, 1)
    Gs[begin] = G

    rets = Vector{ConvexClustering.ALMSolutionType{T}}(undef, 1)
    rets[begin] = ret

    return Gs, rets, [problem.γ;]
end

# search for θ according to the stopping condition and θ sequence in config_θ, then set up the problem.
function getselectedgraph(
    config_θ::CC.SearchθConfigType{T},
    Δc_m::Vector{Vector{T}},
    graph_config::CC.WeightedGraphConfigType;
    θ_default::T = one(T), # not used.
    verbose_kernel::Bool = false,
    ) where T <: AbstractFloat

    A, edge_pairs, w, A_neighbourhoods, θs = ConvexClustering.searchkernelparameters(
        T,
        Δc_m,
        config_θ,
        graph_config;
        verbose = verbose_kernel,
    )

    return A, edge_pairs, w, θs 
end

# use a supplied θ_default to set up the problem. No error-checking on the validity of θ with graph_config.kernelfunc(). To implement a check on θ, could put a θ validity check function in graph_config.
function getselectedgraph(
    ::Nothing,
    Δc_m::Vector{Vector{T}},
    graph_config::CC.WeightedGraphConfigType;
    θ_default::T = one(T),
    verbose_kernel::Bool = false,
    ) where T <: AbstractFloat

    ### setup convex clustering problem.
    #θ::T = convert(T, config_θ.getθfunc(1))

    A, edge_pairs, w, A_neighbourhoods = ConvexClustering.setupproblem(
        Δc_m,
        θ_default,
        graph_config,
    )
    
    return A, edge_pairs, w, [θ_default;]
end
