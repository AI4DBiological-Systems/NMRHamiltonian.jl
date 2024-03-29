

####### generic utils.

# consider deprecating this, and use collect( Iterators.flatten() )
# won't deprecate now since combinevectors(x) is still more efficient.
# julia> @btime HAM.combinevectors($a);
#   190.265 ns (1 allocation: 176 bytes)

# julia> @btime collect(Iterators.flatten($a));
#   214.844 ns (3 allocations: 480 bytes)

# julia> @btime collect(Iterators.flatten($a));
#   219.786 ns (3 allocations: 480 bytes)

# julia> @btime HAM.testfunc($a);
#   365.442 ns (3 allocations: 480 bytes)

# julia> @btime HAM.testfunc($a);
#   353.614 ns (3 allocations: 480 bytes)

function combinevectors(x::Vector{Vector{T}})::Vector{T} where T

    if isempty(x)
        return Vector{T}(undef, 0)
    end

    N = sum(length(x[i]) for i = 1:length(x))

    y = Vector{T}(undef,N)

    st_ind = 0
    fin_ind = 0
    for i = 1:length(x)
        st_ind = fin_ind + 1
        fin_ind = st_ind + length(x[i]) - 1

        y[st_ind:fin_ind] = x[i]
    end

    return y
end

# # used as comparison for performance benchmarking combinevectors().
# function testfunc(x::Vector{Vector{T}})::Vector{T} where T
#     return collect( Iterators.flatten(x) )
# end

# """
#     convertcompactdomain(x::T, a::T, b::T, c::T, d::T)::T

# converts compact domain x ∈ [a,b] to compact domain out ∈ [c,d].
# """
function convertcompactdomain(x::T, a::T, b::T, c::T, d::T)::T where T <: AbstractFloat

    return (x-a)*(d-c)/(b-a)+c
end

function convertcompactdomain(x::Vector{T}, a::T, b::T, c::T, d::T)::Vector{T} where T <: AbstractFloat

    return collect( convertcompactdomain(x[i], a, b, c, d) for i = 1:length(x) )
end

