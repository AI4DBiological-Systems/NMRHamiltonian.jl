
#### load chemical shift and J-coupling information from JSOn.
function loadcouplinginfojson(
    ::Type{T},
    save_path::String,
    )::Tuple{Vector{Int}, Vector{T}, Vector{Tuple{Int,Int}}, Vector{T}} where T <: AbstractFloat

    J = JSON3.read(read(save_path))

    H_IDs = Vector{Int}(undef, 0)
    H_css = Vector{T}(undef, 0)
    for dict in J[Symbol("chemical shift")]
        id = convert(Int, dict[Symbol("ID")])
        value = convert(Float64, dict[Symbol("value")])
        push!(H_IDs, id)
        push!(H_css, value)
    end

    J_IDs = Vector{Tuple{Int,Int}}(undef, 0)
    J_vals = Vector{T}(undef, 0)
    for dict in J[Symbol("J-coupling")]
        id1 = convert(Int, dict[Symbol("ID1")])
        id2 = convert(Int, dict[Symbol("ID2")])
        value = convert(T, dict[Symbol("value")])
        push!(J_IDs, (id1, id2))
        push!(J_vals, value)
    end

    return H_IDs, H_css, J_IDs, J_vals
end


function constructspinsystemsgraphformolecule(N_vertices::Int,
    J_inds::Vector{Tuple{Int,Int}})

    g = Graphs.SimpleGraph(N_vertices)

    for k = 1:length(J_inds)
        i, j = J_inds[k]
        Graphs.add_edge!(g, i, j)
    end

    return g
end

function convertJIDstoJinds(J_IDs::Vector{Tuple{Int,Int}},
    H_IDs::Vector{Int})::Vector{Tuple{Int,Int}}

    J_inds = Vector{Tuple{Int,Int}}(undef, 0)
    for k = 1:length(J_IDs)

        i = findfirst(xx->xx==J_IDs[k][1], H_IDs)
        j = findfirst(xx->xx==J_IDs[k][2], H_IDs)

        push!(J_inds, (i,j))
    end

    return J_inds
end

# """
# matchJlabels((label_pair_list::Vector{Tuple{Int,Int}}, search_list::Vector{Int})::Vector{Int}

# Returns the indices of entries of `label_pair_list` that has both labels in the `search_list`.
# """
function matchJlabels(label_pair_list::Vector{Tuple{Int,Int}}, search_list::Vector{Int})::Vector{Int}

    inds = Vector{Int}(undef, 0)
    for k = 1:length(label_pair_list)
        i, j = label_pair_list[k]

        if any(i .== search_list) && any(j .== search_list)
            push!(inds, k)
        end
    end

    return inds
end

function setupcsJ(
    H_IDs,
    H_css::Vector{T},
    J_IDs,
    J_vals;
    unique_cs_tol::T = convert(T, 1e-6),
    ) where T

    H_inds = collect(1:length(H_IDs))
    J_inds = convertJIDstoJinds(J_IDs, H_IDs)

    g = constructspinsystemsgraphformolecule(length(H_IDs), J_inds)
    systems_g = Graphs.connected_components(g)

    dict_H_inds_to_css = Dict(H_inds .=> H_css)

    J_inds_sys = Vector{Vector{Tuple{Int,Int}}}(undef, 0)
    J_IDs_sys = Vector{Vector{Tuple{Int,Int}}}(undef, 0)
    J_vals_sys = Vector{Vector{T}}(undef, 0)
    H_inds_sys = Vector{Vector{Int}}(undef, 0)
    cs_sys = Vector{Vector{T}}(undef, 0)

    H_inds_singlets = Vector{Vector{Int}}(undef, 0)
    cs_singlets = Vector{T}(undef, 0)

    for i in eachindex(systems_g)

        H_inds_i = systems_g[i]

        if length(H_inds_i) > 1

            cs = collect( dict_H_inds_to_css[H_inds_i[l]] for l in eachindex(H_inds_i) )

            # check if cs is a vector of the same value, up to an absolute tolerance of `unique_cs_tol`.
            if all(isapprox.(cs, cs[1], atol = unique_cs_tol))

                # add singlet.
                push!(H_inds_singlets, H_inds_i)
                push!(cs_singlets, dict_H_inds_to_css[H_inds_i[1]])
            else

                # add spin system.
                push!(H_inds_sys, H_inds_i)
                push!(cs_sys, cs)

                inds = matchJlabels(J_inds, H_inds_i)
                push!(J_inds_sys, J_inds[inds])
                push!(J_IDs_sys, J_IDs[inds])
                push!(J_vals_sys, J_vals[inds])
            end

        elseif length(H_inds_i) == 1

            # add singlet.
            push!(H_inds_singlets, H_inds_i)
            push!(cs_singlets, dict_H_inds_to_css[H_inds_i[1]])
        end
    end

    # remove singlets that have duplicate chemical shifts.
    H_inds_singlets, cs_singlets = removeredundantsinglets(H_inds_singlets, cs_singlets)

    #
    J_inds_sys_local = convertJindsglobaltolocal(H_inds_sys, J_inds_sys)

    # check.
    for i = 1:length(J_inds_sys)
        @assert length(J_inds_sys[i]) == length(J_inds_sys_local[i]) == length(J_IDs_sys[i]) == length(J_vals_sys[i])
        @assert length(H_inds_sys[i]) == length(cs_sys[i])
    end

    @assert length(H_inds_singlets) == length(cs_singlets)

    return CSJContainer(
        J_inds_sys, J_inds_sys_local, J_IDs_sys, J_vals_sys, H_inds_sys,
        cs_sys, H_inds_singlets, cs_singlets, H_inds, J_inds,
    ), g
end

function removeredundantsinglets(H_inds::Vector{Vector{Int}}, cs_singlets::Vector{T}) where T <: AbstractFloat

    cs_singlets_unique, inds_unique = uniqueinds(cs_singlets)

    H_inds_unique = Vector{Vector{Int}}(undef, 0)
    for i in eachindex(inds_unique)

        inds = inds_unique[i]
        push!(H_inds_unique, combinevectors(H_inds[inds]))
    end

    return H_inds_unique, cs_singlets_unique
end

# """
# convertJindsglobaltolocal(H_inds_sys::Vector{Vector{Int}},
#     J_inds_sys::Vector{Vector{Tuple{Int,Int}}})

# The return type is Vector{Vector{Tuple{Int,Int}}}.

# Local: Each spin system will have its own spin nucleui numbering that starts at 1.
# Global: The spins systems will together have one single spin nucleui numbering that starts at 1.
# """
function convertJindsglobaltolocal(
    H_inds_sys::Vector{Vector{Int}},
    J_inds_sys::Vector{Vector{Tuple{Int,Int}}},
    )

    @assert length(H_inds_sys) == length(J_inds_sys) # same number of spin systems.
    N_systems = length(H_inds_sys)

    J_inds_sys_local = Vector{Vector{Tuple{Int,Int}}}(undef, N_systems)

    for i = 1:N_systems

        x = H_inds_sys[i]
        J = J_inds_sys[i]

        conversion_dict = Dict(x .=> collect(1:length(x)))

        J_inds_sys_local[i] = collect(
            (conversion_dict[J[k][1]], conversion_dict[J[k][2]])
            for k in eachindex(J)
        )
    end

    return J_inds_sys_local
end

function convertJindslocaltoglobal(
    H_inds_sys::Vector{Vector{Int}},
    J_inds_sys_local::Vector{Vector{Tuple{Int,Int}}},
    )

    @assert length(H_inds_sys) == length(J_inds_sys_local) # same number of spin systems.
    N_systems = length(H_inds_sys)

    J_inds_sys = Vector{Vector{Tuple{Int,Int}}}(undef, N_systems)

    for i = 1:N_systems

        x = H_inds_sys[i]
        J = J_inds_sys_local[i]

        conversion_dict = Dict(collect(1:length(x)) .=> x)

        J_inds_sys[i] = collect(
            (conversion_dict[J[k][1]], conversion_dict[J[k][2]])
            for k in eachindex(J)
        )
    end

    return J_inds_sys
end
