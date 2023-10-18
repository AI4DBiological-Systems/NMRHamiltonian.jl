"""
uniqueinds(a_in::Vector{T}; atol::T = convert(T, 1e-6)) where T

Returns the unique values of `a_in`, with absolute tolerance `atol`, and the indices for each unique value.
"""
function uniqueinds(a_in::Vector{T}; atol::T = convert(T, 1e-6)) where T <: AbstractFloat

    if length(a_in) < 2
        inds = Vector{Vector{Int}}(undef, length(a_in))

        if length(a_in) == 1
            inds[1] = Vector{Int}(undef, 1)
            inds[1][1] = 1
        end

        return copy(a_in), inds
    end

    a = copy(a_in)
    inds_a = collect(1:length(a))

    b = Vector{T}(undef, 0)
    inds_b = Vector{Vector{Int}}(undef, 0)

    while !isempty(a)

        target = pop!(a)
        inds = findall(xx->isapprox(target, xx; atol = atol), a)

        target_ind = pop!(inds_a)

        # save to output.
        push!(b, target)
        push!(inds_b, [target_ind; inds_a[inds]])

        # delete saved.
        deleteat!(a, inds)
        deleteat!(inds_a, inds)
    end

    return b, inds_b
end

"""
keeptargetintegers(C::Vector{Vector{Int}}, search_list::Vector{Int})

keep an entry from C if the any of its integer array values appear in any entry of `search_list`.
"""
function keeptargetintegers(C::Vector{Vector{Int}}, search_list::Vector{Int})

    keep_flags = falses(length(C))
    for i = 1:length(C)

        if any(C[i][1] .== search_list)
            keep_flags[i] = true
        end
    end


    return C[keep_flags]
end

"""
getpairs(inds::Vector{T})

get all exhaustive pairwise combos without symmetry of the 1D array `inds`.
"""
function getpairs(inds::Vector{T}) where T

    out = Vector{Tuple{T,T}}(undef, 0)
    for i = 1:length(inds)
        for j = i+1:length(inds)
            push!(out, (inds[i], inds[j]))
        end
    end

    return out
end

"""
isallsame(a::Vector{T}; atol::T = convert(T, 1e-6)) where T

returns true if the entries in `a` are all within an abolute tolerance of `atol`.
"""
function isallsame(a::Vector{T}; atol::T = convert(T, 1e-6)) where T
    if length(findall(xx->isapprox(a[1], xx; atol = atol), a)) == length(a)
        return true
    end

    return false
end

"""
Returns the edges of a connected path as specified by `vertices`.
"""
function getconnectpath(vertices::Vector{Int})
    N = length(vertices)

    out = Vector{Tuple{Int,Int}}(undef, N-1)

    for i = 1:N-1
        out[i] = (vertices[i], vertices[i+1])
    end

    return out
end
