# """
# uniqueinds(a_in::Vector{T}; atol::T = convert(T, 1e-6)) where T

# Returns the unique values of `a_in`, with absolute tolerance `atol`, and the indices for each unique value.
# """
function uniqueinds(a::Vector{T}; digits::Int = 6) where T <: AbstractFloat

    atol = digits2atol(T, digits)

    b = unique(xx -> round(xx, digits = digits), a)
    idx = unique(zz -> round(a[zz], digits = digits), 1:length(a))

    inds_b = Vector{Vector{Int}}(undef, length(idx))
    for i in eachindex(idx)
        
        k = idx[i]

        inds_b[i] = findall(
            xx->isapprox(
                xx,
                round(a[k], digits = digits);
                atol = atol,
            ),
            a,
        )
    end

    return b, inds_b
end

function digits2atol(::Type{T}, digits::Int)::T where T <: AbstractFloat

    return convert(T, 10) ^ (-digits)
end

# """
# keeptargetintegers(C::Vector{Vector{Int}}, search_list::Vector{Int})

# keep an entry from C if the any of its integer array values appear in any entry of `search_list`.
# """
function keeptargetintegers(C::Vector{Vector{Int}}, search_list::Vector{Int})

    keep_flags = falses(length(C))
    for i = 1:length(C)

        if any(C[i][1] .== search_list)
            keep_flags[i] = true
        end
    end


    return C[keep_flags]
end

# """
# getpairs(inds::Vector{T})

# get all exhaustive pairwise combos without symmetry of the 1D array `inds`.
# """
function getpairs(inds::Vector{T}) where T

    out = Vector{Tuple{T,T}}(undef, 0)
    for i = 1:length(inds)
        for j = i+1:length(inds)
            push!(out, (inds[i], inds[j]))
        end
    end

    return out
end

# """
# isallsame(a::Vector{T}; atol::T = convert(T, 1e-6)) where T

# returns true if the entries in `a` are all within an abolute tolerance of `atol`.
# """
function isallsame(a::Vector{T}; atol::T = convert(T, 1e-6)) where T
    if length(findall(xx->isapprox(a[1], xx; atol = atol), a)) == length(a)
        return true
    end

    return false
end

# """
# Returns the edges of a connected path as specified by `vertices`.
# """
function getconnectpath(vertices::Vector{Int})
    N = length(vertices)

    out = Vector{Tuple{Int,Int}}(undef, N-1)

    for i = 1:N-1
        out[i] = (vertices[i], vertices[i+1])
    end

    return out
end
