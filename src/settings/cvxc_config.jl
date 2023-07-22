
# k-nearest neighbour-based contrsuction of graph, with search.
# keep incrementing knn until the number of connected_components is equal or below max_connected_components.
function packageknnsearchconfig(
    start_knn::Int,
    max_knn::Int,
    max_connected_components::Int;
    metric::Distances.Metric = Distances.Euclidean(),
    kernelfunc::Function = evalSqExpkernel,
    verbose = false)

    connectivity = ConvexClustering.KNNSearchType(0; 
        start_knn = start_knn, verbose = verbose,
        searchfunc = (xx,mm)->ConvexClustering.searchknn(start_knn,
            mm, xx, max_connected_components;
            max_knn = max_knn,
            verbose = verbose)
        )

    return CC.WeightedGraphConfigType(connectivity, metric, kernelfunc)
end

# k-nearest neighbour-based contrsuction of graph, with search.
function defaultknnsearchconfig(
    name,
    i,
    Δc_m::Vector{Vector{T}},
    αs::Vector{T};
    start_knn::Int = max(60, round(Int, length(Δc_m)*0.05)),
    max_knn = max(60, round(Int, length(Δc_m)*0.2)),
    max_connected_components_offset::Int = -1, # set to negative large number to force max_connected_components to be 1.
    kernelfunc::Function = evalSqExpkernel,
    verbose = false,
    ) where T <: AbstractFloat
    
    @assert !isempty(Δc_m)
    N = length(Δc_m)
    D = length(first(Δc_m)) # number of effective/unique chemical shifts.

    max_connected_components = clamp(D + max_connected_components_offset, 1, N)

    # make sure the these are valid.
    start_knn = clamp(start_knn, 2, N-1)
    max_knn = convert(Int, min(max_knn, N-1))
    #@show max_knn

    return packageknnsearchconfig(
        start_knn,
        max_knn,
        max_connected_components;
        kernelfunc = kernelfunc,
        verbose = verbose)
end


################

function packageradiussearchconfig(
    max_connected_components::Int,
    ::Type{T};
    max_iters::Int = 100,
    metric::Distances.Metric = Distances.Euclidean(),
    kernelfunc::Function = evalSqExpkernel,
    verbose = false) where T

    # need ConvexClustering.RadiusSearchType() instead of RadiusSearchType since we're invoking a custom constructor func from ConvexClustering, not the default contructor of RadiusSearchType.
    connectivity = ConvexClustering.RadiusSearchType(-one(T); 
        max_iters = max_iters, verbose = verbose,
        searchfunc = (xx,mm)->ConvexClustering.searchradius(
            mm, xx, max_connected_components)
        )

    return CC.WeightedGraphConfigType(connectivity, metric, kernelfunc)
end

function defaultradiussearchconfig(
    name,
    i,
    Δc_m::Vector{Vector{T}},
    αs::Vector{T};
    max_iters::Int = 100,
    max_connected_components_offset::Int = -1, # set to negative large number to force max_connected_components to be 1.
    kernelfunc::Function = evalSqExpkernel,
    verbose = false
    ) where T <: AbstractFloat
    
    @assert !isempty(Δc_m)
    D = length(first(Δc_m)) # number of effective/unique chemical shifts.

    max_connected_components::Int = max(1, D + max_connected_components_offset)

    return packageradiussearchconfig(
        max_connected_components,
        T;
        max_iters = max_iters,
        kernelfunc = kernelfunc,
        verbose = verbose)
end


##################

function packageradiusconfig(
    radius::T;
    metric::Distances.Metric = Distances.Euclidean(),
    kernelfunc::Function = evalSqExpkernel,
    ) where T <: AbstractFloat

    connectivity = ConvexClustering.RadiusType(radius) # make an edge for all points within this radius of a given point i, cycle through all i in the point set to be partitioned. Takes a finite floating point number.

    return CC.WeightedGraphConfigType(connectivity, metric, kernelfunc)
end

function defaultradiusconfig(
    name,
    i,
    Δc_m::Vector{Vector{T}},
    αs::Vector{T};
    kernelfunc::Function = evalSqExpkernel,
    ) where T <: AbstractFloat

    @assert !isempty(Δc_m)
    D = length(first(Δc_m))

    radius = convert(T, norm(ones(T,D))/(2*D) )

    return packageradiusconfig(radius; kernelfunc = kernelfunc)
end


########

function packageknnconfig(
    knn::Int;
    metric::Distances.Metric = Distances.Euclidean(),
    kernelfunc::Function = evalSqExpkernel,
    )

    ### k-nearest neighbour-based contrsuction of graph.
    connectivity = ConvexClustering.KNNType(knn) # make an edge for this number of nearest neighbours of a given point i, cycle through all i in the point set to be partitioned. Takes a positive integer.

    return CC.WeightedGraphConfigType(connectivity, metric, kernelfunc)
end

function defaultknnconfig(
    name,
    i,
    Δc_m::Vector{Vector{T}},
    αs::Vector{T};
    kernelfunc::Function = evalSqExpkernel,
    ) where T <: AbstractFloat

    ### k-nearest neighbour-based contrsuction of graph.
    N = length(Δc_m)
    knn = max(64, round(Int, N/10))
    knn = convert(Int, min(knn, N-1)) # cannot be more than the cardinality -1 of the set to be partitioned.

    return packageknnconfig(knn; kernelfunc = kernelfunc)
end