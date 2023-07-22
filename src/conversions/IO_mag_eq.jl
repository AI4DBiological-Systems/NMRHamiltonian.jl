"""
output index format:
x[i][j][k][l]
i-th spin system, k-th magnetically equivalent nuclei, l-th local spin index.
"""
function getmageqmolecule(
    g,
    H_inds_sys,
    dict_ind_to_H_ID,
    dict_H_inds_to_css,
    dict_H_IDs_to_css,
    dict_J_ID_to_val,
    J_IDs;
    atol = 1e-6)

    C_g = Graphs.maximal_cliques(g)

    N_spin_systems = length(H_inds_sys)

    mag_eq_sys_inds_local = Vector{Vector{Vector{Int}}}(undef, N_spin_systems)
    mag_eq_sys_inds_global = Vector{Vector{Vector{Int}}}(undef, N_spin_systems)
    mag_eq_sys_IDs = Vector{Vector{Vector{Int}}}(undef, N_spin_systems)

    for i = 1:N_spin_systems
        C = keeptargetintegers(C_g, H_inds_sys[i])

        mag_eq_IDs0, mag_eq_inds0 = getmageqIDs(C,
            dict_ind_to_H_ID,
            dict_H_inds_to_css,
            dict_H_IDs_to_css, dict_J_ID_to_val, J_IDs;
            atol = atol)

        mag_eq_inds = combinetransitiveeqgroups(mag_eq_inds0)
        mag_eq_IDs = combinetransitiveeqgroups(mag_eq_IDs0)

        # mag_eq_inds = mag_eq_inds0
        # mag_eq_IDs = mag_eq_IDs0

        mag_eq_sys_inds_local[i] = getmageqlocalinds(H_inds_sys[i], mag_eq_inds)
        mag_eq_sys_IDs[i] = mag_eq_IDs
        mag_eq_sys_inds_global[i] = mag_eq_inds
    end

    return mag_eq_sys_inds_local, mag_eq_sys_IDs, mag_eq_sys_inds_global
end

"""
convertlabelsglobaltolocal(H_inds_sys::Vector{Vector{Int}},
    mag_eq_sys_inds::Vector{Vector{Int}})

The return type is Vector{Vector{Int}}.

Local: Each spin system will have its own spin nucleui numbering that starts at 1.
Global: The spins systems will together have one single spin nucleui numbering that starts at 1.
"""
function getmageqlocalinds(H_inds::Vector{Int},
    mag_eq_sys_inds::Vector{Vector{Int}})

    N_eqs = length(mag_eq_sys_inds)
    labels_sys_local = Vector{Vector{Int}}(undef, N_eqs)

    for i = 1:N_eqs

        x = H_inds
        conversion_dict = Dict(x .=> collect(1:length(x)))

        z = mag_eq_sys_inds[i]
        labels_sys_local[i] = collect( conversion_dict[z[k]] for k = 1:length(z) )
    end

    return labels_sys_local
end


"""
getmageqIDs(C::Vector{Vector{Int}},
    dict_ind_to_H_ID,
    dict_H_inds_to_css;
    atol::Float64 = 1e-6)

Given a list of node indices `C` of cliques of an undirected graph, with the indices of a clique being C[i],
Returns the list of node indices for each clique that are magnetically equivalent.
"""
function getmageqIDs(C::Vector{Vector{Int}},
    dict_ind_to_H_ID,
    dict_H_inds_to_css,
    dict_H_IDs_to_css,
    dict_J_ID_to_val,
    J_IDs;
    atol::Float64 = 1e-6)

    # traverse each maximally connected cliques.
    C_IDs_mag_eq = Vector{Vector{Vector{Int}}}(undef, length(C))
    C_inds_mag_eq = Vector{Vector{Vector{Int}}}(undef, length(C))
    for i = 1:length(C)

        ### partition the nucleus by unique cs.
        cs = collect( dict_H_inds_to_css[C[i][k]] for k = 1:length(C[i]) )
        cs_IDs = collect( dict_ind_to_H_ID[C[i][k]] for k = 1:length(C[i]) )
        cs_inds = C[i]

        # unique_cs = unique(round.(cs, digits = cs_round_digits))
        # unique_cs_inds = collect( inds = findall(xx->isapprox(unique_cs[k], xx; atol = atol), cs) for k = 1:length(unique_cs) )

        unique_cs, unique_cs_inds = uniqueinds(cs; atol = atol)

        # check magnetic equivalence for the chemically equivalent entries in `unique_cs_inds`.
        pass_flags = collect( checkmageq(unique_cs_inds[k], cs_IDs, dict_J_ID_to_val,
            dict_H_IDs_to_css, J_IDs; atol = atol) for k = 1:length(unique_cs_inds) )

        # store result.
        common_cs_IDs = collect( cs_IDs[unique_cs_inds[k]] for k = 1:length(unique_cs_inds) )
        C_IDs_mag_eq[i] = common_cs_IDs[pass_flags]

        common_cs_inds = collect( cs_inds[unique_cs_inds[k]] for k = 1:length(unique_cs_inds) )
        C_inds_mag_eq[i] = common_cs_inds[pass_flags]

        if length(common_cs_IDs[pass_flags]) > 2
            println("Warning for clique $(i): more than two groups of magnetically equivalent nuclei!")
        end
    end

    # remove empty entries.
    keep_flags = trues(length(C))
    for i = 1:length(C)
        if isempty(C_IDs_mag_eq[i])
            keep_flags[i] = false
        end
    end
    C_IDs_mag_eq = C_IDs_mag_eq[keep_flags]
    C_inds_mag_eq = C_inds_mag_eq[keep_flags]

    mag_eq_IDs::Vector{Vector{Int}} = combinevectors(C_IDs_mag_eq)
    mag_eq_inds::Vector{Vector{Int}} = combinevectors(C_inds_mag_eq)

    return mag_eq_IDs, mag_eq_inds
end

### routines related to checking whether the J-coupling value of A-C and B-C
# are the same, for all C's connected to either A or B. A and B are
# chemically equivalent nuclei.
"""
matchJlabels((label_pair_list::Vector{Tuple{Int,Int}}, search_list::Vector{Int})::Vector{Int}

Returns the indices of entries of `label_pair_list` that has any of the labels in the `search_list`.
"""
function matchanyJlabels(label_pair_list::Vector{Tuple{Int,Int}}, search_list::Vector{Int})::Vector{Int}

    inds = Vector{Int}(undef, 0)
    for k = 1:length(label_pair_list)
        i, j = label_pair_list[k]

        if any(i .== search_list) || any(j .== search_list)
            push!(inds, k)
        end
    end

    return inds
end

"""
getJIDstest(J_IDs::Vector{Tuple{Int,Int}},
    common_cs_IDs::Vector{Int},
    dict_H_IDs_to_css; atol = 1e-6)

Get the list of A-C, B-C, A-D, B-D, etc pairs to test.
"""
function getJIDstest(J_IDs::Vector{Tuple{Int,Int}},
    common_cs_IDs::Vector{Int},
    dict_H_IDs_to_css; atol = 1e-6)

    inds_any = matchanyJlabels(J_IDs, common_cs_IDs)
    inds_both = matchJlabels(J_IDs, common_cs_IDs)
    inds = setdiff(inds_any, inds_both)

    #inds_both =
    J_IDs_tmp = J_IDs[inds]


    target_IDs = Vector{Int}(undef, 0)
    for l = 1:length(J_IDs_tmp)
        i, j = J_IDs_tmp[l]

        # look for nuclei pairs that have the same chemical shift.
        if !isapprox(dict_H_IDs_to_css[i], dict_H_IDs_to_css[j], atol = atol)

            # find the value between i and j that isn't in common_cs_IDs.
            target_ID = i
            if any(i .== common_cs_IDs)
                target_ID = j
            end

            push!(target_IDs, target_ID)
        end
    end

    return target_IDs
end



"""
checkmageq(test_inds::Vector{Vector{Int}},
    cs_IDs::Vector{Int},
    dict_J_ID_to_val;
    atol = 1e-6)::Bool

`test_inds` is a list of indices in cs_IDs. This function checks for magnetic equivalence for the nuclei in `test_inds`.
"""
function checkmageq(test_inds::Vector{Int},
    cs_IDs::Vector{Int},
    dict_J_ID_to_val,
    dict_H_IDs_to_css,
    J_IDs;
    atol = 1e-6)::Bool

    if length(test_inds) == 1
        # need to have at least two nuclei to be magnetically equivalent.
        return false
    end


    common_cs_IDs = cs_IDs[test_inds]
    p_IDs = getpairs(common_cs_IDs)

    ## test common J within pairs.
    J_p = collect( getJfromdict(p_IDs[l][1], p_IDs[l][2], dict_J_ID_to_val) for l = 1:length(p_IDs) )
    pass_common_J_flag = isallsame(J_p; atol = atol)

    ## test common J with other IDs.

    # find all connections to first pair.
    t_IDs = getJIDstest(J_IDs, common_cs_IDs, dict_H_IDs_to_css; atol = atol)

    # check J-values.
    pass_test_J_flag = true
    for k = 1:length(t_IDs)

        J_t = collect( getJfromdict(t_IDs[k], base_ID, dict_J_ID_to_val) for base_ID in common_cs_IDs )

        pass_test_J_flag = pass_test_J_flag & isallsame(J_t; atol = atol)
    end

    return pass_common_J_flag & pass_test_J_flag
end


"""
getJfromdict(i::Int, j::Int, dict::Dict{Tuple{Int,Int},T} ) where T

query `dict` for the entry `(i,j)` and `(j,i)`. Returns the entry for `(i,j)`. Returns the entry for `(j,i)` if `(i,j)` is not found. Return zero if both entries are not found.
"""
function getJfromdict(i::Int, j::Int, dict::Dict{Tuple{Int,Int},T} ) where T
    #
    if haskey(dict, (i,j))
        return dict[(i,j)]
    end

    if haskey(dict, (j,i))
        return dict[(j,i)]
    end

    return zero(T)
end



"""
combinetransitiveeqgroups(eq_inds::Vector{Vector{Int}})

`eq_inds` contains theequivalent nuclei indices.
This function combine the groups of equivalent nuclei via transitivity.

Assumes the equivalence relation cannot be true for a group containg only one nucleus.
"""
function combinetransitiveeqgroups(eq_inds::Vector{Vector{Int}})

    if isempty(eq_inds)
        return Vector{Vector{Int}}(undef, 0)
    end

    eq_inds_flat_unique = unique(combinevectors(eq_inds))
    V = 1:length(eq_inds_flat_unique)
    dict_eq_to_V = Dict(eq_inds_flat_unique .=> V)
    dict_V_to_eq = Dict(V .=> eq_inds_flat_unique)

    eq_inds_V = collect( collect( dict_eq_to_V[eq_inds[i][k]] for k = 1:length(eq_inds[i])) for i = 1:length(eq_inds) )

    # get the edges of fully connected
    E = Vector{Tuple{Int,Int}}(undef, 0)
    for k = 1:length(eq_inds)
        E_k = getconnectpath(eq_inds_V[k])

        push!(E, E_k...)
    end

    #
    h = Graphs.SimpleGraph(V[end])
    for k = 1:length(E)
        Graphs.add_edge!(h, E[k][1], E[k][2])
    end

    H = Graphs.connected_components(h)

    # convert back to the nuclei indexing of `eq_inds`.
    H_out = collect( collect( dict_V_to_eq[H[i][k]] for k = 1:length(H[i])) for i = 1:length(H) )

    return H_out
end

### final package for a molecule.
function getmageqinfo(H_IDs, H_css::Vector{T}, J_IDs, J_vals;
    unique_cs_atol = 1e-6)::Tuple{Vector{Vector{Vector{Int}}},Vector{Vector{Vector{Int}}},Vector{Vector{Vector{Int}}}} where T <: AbstractFloat

    # parse the spin systems and singlets.
    J_inds_sys, J_inds_sys_local, J_IDs_sys, J_vals_sys, H_inds_sys,
    cs_sys, H_inds_singlets, cs_singlets, H_inds, J_inds,
    g = setupcsJ(H_IDs, H_css, J_IDs, J_vals)

    # set up look-up.
    dict_H_IDs_to_css = Dict(H_IDs .=> H_css)
    dict_H_inds_to_css = Dict(H_inds .=> H_css)
    dict_H_ID_to_ind = Dict(H_IDs .=> collect(1:length(H_IDs)))
    dict_ind_to_H_ID = Dict( collect(1:length(H_IDs)) .=> H_IDs)

    dict_J_ID_to_val = Dict(J_IDs .=> J_vals)
    dict_J_ind_to_val = Dict(J_inds .=> J_vals)

    dict_J_ID_to_ind = Dict(J_IDs .=> J_inds)
    dict_J_ind_to_ID = Dict(J_inds .=> J_IDs)

    #
    mag_eq_sys_inds_local, mag_eq_sys_IDs,
    mag_eq_sys_inds_global = getmageqmolecule(g,
    H_inds_sys, dict_ind_to_H_ID, dict_H_inds_to_css,
    dict_H_IDs_to_css, dict_J_ID_to_val, J_IDs; atol = unique_cs_atol)

    return mag_eq_sys_inds_local, mag_eq_sys_IDs,
    mag_eq_sys_inds_global
end

# unused?
function getmageqinfomixture(
    ::Type{T},
    target_entries::Vector{String},
    base_path::String,
    dict_molecule_to_filename;
    unique_cs_atol = 1e-6) where T

    N_molecules = length(target_entries)
    MEs = Vector{Vector{Vector{Vector{Int}}}}(undef, N_molecules)

    for n = 1:N_molecules

        # TODO add error-handling if name is not found in the dictionary, or filename does not exist.
        load_path = joinpath(base_path, dict_molecule_to_filename[target_entries[n]]["file name"])
        H_IDs, H_css, J_IDs, J_vals = loadcouplinginfojson(T, load_path)

        MEs[n], _ = getmageqinfo(H_IDs, H_css, J_IDs, J_vals; unique_cs_atol = unique_cs_atol)
    end

    return MEs
end
