
import Dictionaries
import JSON3

cs_systems = Dictionaries.ArrayDictionary(Symbol.(molecule_entries), cs_sys_mixture)
cs_singlets = Dictionaries.ArrayDictionary(Symbol.(molecule_entries), cs_sys_mixture)

IterableType = Union{Vector, JSON3.Array, JSON3.Object, Dictionaries.ArrayDictionary } #Union{T, Nothing}


function ininterval(x::IterableType, position, low::T, high::T)::Bool where T

    if isempty(x)
        return false
    end

    y = x[begin+position-1]
    current_type = typeof(y)

    while current_type <: IterableType

        return ininterval(y, low, high)
    end

    return low < y < high
end

function ininterval(x::IterableType, low::T, high::T)::Bool where T

    y = collect( ininterval(x, k, low, high) for k in eachindex(x) )

    return any(y)
end

low = 7.0
high = 8.0

Q = filter(xx->ininterval(xx, low, high), cs_systems)