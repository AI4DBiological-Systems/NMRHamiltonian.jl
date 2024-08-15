# TODO export these functions.

function get_partitions(A::SHType)
    return A.parts
end

function get_rcs(A::SHType)
    return A.cs_Δc
end

#  `cd` stands for coherence difference, and is Δc.
function get_group_cd(A::SHType)
    return A.Δc_bar
end

#  `cd` stands for coherence difference, and is Δc.
function get_component_cd(A::SHType)
    return A.Δc
end

function num_groups(A::SHType, i::Integer)
    return length(A.parts[i])
end

function get_Δc_matrix(A::SHType, i::Integer)
    num_ME_nuclei = length(A.Δc[i][begin])
    num_components = length(A.Δc[i]) # resonance comopnents in spin system i.
    return reshape(collect(Iterators.flatten(A.Δc[i])), num_ME_nuclei, num_components)
end

function get_sorted_Δc_matrix(A::SHType, i::Integer)
    #inds = sortperm(A.αs[i], rev = true)
    inds = collect(
        Iterators.flatten(
            A.parts[i]
        )
    )
    dc = A.Δc[i][inds]
    
    num_ME_nuclei = length(dc[begin])
    num_components = length(dc) # resonance comopnents in spin system i.
    U = reshape(collect(Iterators.flatten(dc)), num_ME_nuclei, num_components)

    return U, A.αs[i][inds]
end