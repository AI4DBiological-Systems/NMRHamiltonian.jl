###### l-infinite norm as a distance metric for NearestNeighbors.jl.
function evalmaxnorm(x, y)
    return maximum(abs.(x-y))
end

# Declare the metric
struct MaxNormMetric <: NearestNeighbors.Metric end
(::MaxNormMetric)(x,y) = evalmaxnorm(x,y)

"""
Uses the NearestNeighbors.jl search library to determine which resonances should be grouped together.
"""
function partitionresonancesbyneighbors(x_in::Vector{Vector{T}},
    amplitudes_in::Vector{T},
    threshold_amplitude::T;
    radius = 1e-1,
    ) where T <: AbstractFloat

    @assert length(x_in) == length(amplitudes_in)
    @assert !isempty(x_in)

    # take largest amplitude, set as new partition centre.
    amplitudes = copy(amplitudes_in)
    x = copy(x_in)
    inds_buffer = collect(1:length(x))

    out_inds = Vector{Vector{Int}}(undef, 0)

    x_maxs = Vector{Vector{T}}(undef, 0)

    max_amplitude, max_ind = findmax(amplitudes)
    while !isempty(amplitudes) || max_amplitude < threshold_amplitude

        # set up tree.
        X = Matrix{T}(undef, length(x[1]), length(x))
        for r in eachindex(x)
            X[:,r] = x[r]
        end

        balltree = NearestNeighbors.BallTree(X, MaxNormMetric(); reorder = false)


        inds = NearestNeighbors.inrange(balltree, x[max_ind], radius, true)
        push!(x_maxs, x[max_ind])

        # book keep.
        push!(out_inds, inds_buffer[inds])

        deleteat!(amplitudes, inds)
        deleteat!(x, inds)
        deleteat!(inds_buffer, inds)

        if isempty(amplitudes)
            return out_inds, x_maxs
        end
        max_amplitude, max_ind = findmax(amplitudes)
    end

    return out_inds, x_maxs
end

