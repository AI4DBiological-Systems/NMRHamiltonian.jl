
# get unique values of css_sys in graph.jl.
# It is denoted here as cs.
function cstopcs(   cs::Vector{Vector{T}},
    cs_LUT::Vector{Vector{Vector{Int}}}) where T <: Real

    N_subsystems = length(cs)
    p_cs = Vector{Vector{T}}(undef, N_subsystems)

    for m = 1:N_subsystems

    p_cs[m] = Vector{T}(undef, length(cs_LUT[m]))

    for i = 1:length(cs_LUT[m])
    inds = cs_LUT[m][i]

    # take average.
    p_cs[m][i] = sum( cs[m][k] for k in inds)/length(inds)
    end

    end

    return p_cs
end


function pcstocs(   p_cs::Vector{Vector{T}},
    cs_LUT::Vector{Vector{Vector{Int}}},
    cs_len_sys::Vector{Int}) where T <: Real

    N_subsystems = length(p_cs)

    cs = Vector{Vector{T}}(undef, N_subsystems)
    pcstocs!(cs, p_cs, cs_LUT, cs_len_sys)

    return cs
end

function pcstocs!(  cs::Vector{Vector{T}},
    p_cs::Vector{Vector{T}},
    cs_LUT::Vector{Vector{Vector{Int}}},
    cs_len_sys::Vector{Int}) where T <: Real
    #
    N_subsystems = length(p_cs)

    for m = 1:N_subsystems

    cs[m] = Vector{T}(undef, cs_len_sys[m])

    # for i = 1:length(p_cs[m])
    #     inds = cs_LUT[m][i]
    #
    #     for k in inds
    #         cs[m][k] = p_cs[m][i]
    #     end
    # end
    pcstocs!(cs[m], p_cs[m], cs_LUT[m])
end

return nothing
end

function pcstocs!(  cs_m::Vector{T},
    p_cs_m::Vector{T},
    cs_LUT_m::Vector{Vector{Int}}) where T <: Real
    
    @assert length(p_cs_m) == length(cs_LUT_m) <= length(cs_m)

    for i = 1:length(p_cs_m)
    inds = cs_LUT_m[i]

    for k in inds
    cs_m[k] = p_cs_m[i]
    end
    end

    return nothing
end