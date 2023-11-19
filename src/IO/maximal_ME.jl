
function canmaximizeME(
    P::PhysicalParamsType{T};
    unique_cs_digits::Int = 6,
    ) where T

    #unique_cs_atol = digits2atol(T, unique_cs_digits)

    cs_sys, ME = P.cs_sys, P.ME

    for i in eachindex(cs_sys)
        cs = cs_sys[i]
        unique_cs, unique_cs_inds = uniqueinds(P.cs_sys[i]; digits = unique_cs_digits)
        
        if length(unique_cs) != length(cs)
            # there are identical chemical shifts in this spin system.
            @show unique_cs, cs

            if isempty(ME[i])
                return true
            end

            for k in eachindex(unique_cs_inds)
                k_set = Set(unique_cs_inds[k])
                ind = findfirst(xx->(Set(xx) == k_set), ME[i])

                if isnothing(ind) && length(k_set) > 1 # exclude the singleton ME, since P.ME[i] only stores non-singleton ME subsets.
                    @show k_set, ME[i]
                    return true
                end
            end
        end
    end

    return false
end

function createmaximalME(
    Phys_in::Vector{PhysicalParamsType{T}};
    unique_cs_digits::Int = 6,
    unique_J_avg_digits::Int = 6,
    ) where T <: AbstractFloat
    
    Phys = Vector{PhysicalParamsType{T}}(undef, length(Phys_in))
    for n in eachindex(Phys)

        if canmaximizeME(Phys_in[n]; unique_cs_digits = unique_cs_digits)
        
            Phys[n] = createmaximalME(
                Phys_in[n];
                unique_cs_digits = unique_cs_digits,
                unique_J_avg_digits = unique_J_avg_digits,
            )
        else
            # leave alone.
            Phys[n] = deepcopy(Phys_in[n])
        end
    end

    return Phys
end

function createmaximalME(
    P_in::PhysicalParamsType{T};
    unique_cs_digits::Int = 6,
    unique_J_avg_digits::Int = 6,
    ) where T
    
    P = deepcopy(P_in)

    N_sys = length(P.cs_sys)
    edge_weights_sys = Vector{Vector{T}}(undef, N_sys)
    edges_sys = Vector{Vector{Tuple{Int,Int}}}(undef, N_sys)

    # loop over the non-singlet spin systems.
    for i in eachindex(P.cs_sys)

        # parse, set up J-coupling connectivity graph, h, for this spin system.
        node_weights, edges, edge_weights, h = createspinsystemgraph(P, i)
        cs = node_weights
        
        # find nuclei subsets that have common chemical shifts.
        # Then make a fully connected graph between the union of these subsets.
        # Initialize new connection weights (J-coupling) to zero.
        makeidenticalcsconnected!(
            edge_weights, edges, h, cs; unique_cs_digits = unique_cs_digits,
        )

        # add edges if the average J coupling connnected to each nuclei in a subset is the same. Recall nuceli in a subset have the same chemical shift here.
        maximizeME!(
            edge_weights,
            cs,
            edges;
            unique_cs_digits = unique_cs_digits,
            unique_J_avg_digits = unique_J_avg_digits,
        )

        # remove the edges that we added that have zero edge weights (zero J-coupling)
        pruneedgeweights!(edge_weights, edges; atol = digits2atol(T, unique_J_avg_digits))
        
        edge_weights_sys[i] = edge_weights
        edges_sys[i] = edges
    end

    # overwrite P with new edges and weights.
    P = updatephy(P, edges_sys, edge_weights_sys, unique_cs_digits)

    return P
end

function updatephy(P_in::PhysicalParamsType{T}, edges_sys, edge_weights_sys, unique_cs_digits::Int) where T

    P = deepcopy(P_in)
    @assert length(P.cs_sys) == length(edges_sys) == length(edge_weights_sys)

    edges_sys_global = convertJindslocaltoglobal(P.H_inds_sys, edges_sys)

    for i in eachindex(edges_sys)

        edges, edge_weights = edges_sys[i], edge_weights_sys[i]
        edges_global = edges_sys_global[i]

        resize!(P.J_inds_sys_local[i], length(edges))
        P.J_inds_sys_local[i][:] = edges
        
        resize!(P.J_vals_sys[i], length(edge_weights))
        P.J_vals_sys[i][:] = edge_weights

        resize!(P.J_inds_sys[i], length(edges_global))
        P.J_inds_sys[i][:] = edges_global
    end

    # to get ME, just remake P.
    H_IDs, H_css, J_IDs, J_vals = extractcouplinginfo(P)
    P_new = getPhysicalParamsType(H_IDs, H_css, J_IDs, J_vals; unique_cs_digits = unique_cs_digits)

    return P_new
end


function createspinsystemgraph(P::PhysicalParamsType{T}, i::Int) where T
    @assert 1 <= i <= length(P.cs_sys)

    cs = copy(P.cs_sys[i]) # node weights. The nodes are labeled: 1:length(cs).
    edges = copy(P.J_inds_sys_local[i])
    edge_weights = copy(P.J_vals_sys[i])

    h = Graphs.SimpleGraph(length(cs))
    for k in eachindex(edges)
        Graphs.add_edge!(h, edges[k][1], edges[k][2])
    end

    return cs, edges, edge_weights, h
end

# if groups of nuclei have same cs, then add a J-coupling entry of zero to every combination between the groups of nuclei with the same cs.
function makeidenticalcsconnected!(
    edge_weights::Vector{T}, # mutates.
    edges::Vector{Tuple{Int,Int}}, # mutates.
    h, # mutates.
    cs::Vector{T};
    unique_cs_digits::Int = 6,
    ) where T

    unique_cs, unique_cs_inds = uniqueinds(cs; digits = unique_cs_digits)
    subsets = unique_cs_inds

    LUT = getpairwiseLUT(length(subsets)) # pair-wise iterators.

    for (i,j) in LUT
        A = subsets[i]
        B = subsets[j]

        makefullyconnected!(h, edges, edge_weights, A, B)
    end

    return nothing
end

function makefullyconnected!(
    h, # mutates.
    edges, # mutates.
    edge_weights::Vector{T}, # mutates.
    nodes_same_cs::Vector{Int},
    target_neighbours::Vector{Int},
    ) where T
        
    for m in nodes_same_cs

        for j in target_neighbours

            if m != j && !(Graphs.has_edge(h, m, j)) # edge direction doesn't matter since h is undirected.
                
                # use the src < dest convention.
                src = min(m,j)
                dest = max(m,j)

                Graphs.add_edge!(h, src, dest)
                push!(edges, (src, dest))
                push!(edge_weights, zero(T))
            end
        end
    end

    return nothing
end

function getedgeweight(i::Int, j::Int, edges::Vector{Tuple{Int,Int}}, ws::Vector{T}) where T
    #
    ind = findfirst(xx->(xx==(i,j) || xx==(j,i)), edges)

    if isnothing(ind)
        return zero(T)
    end

    return ws[ind]
end


function getpairwiseLUT(N::Int)

    M = div(N*(N-1), 2)
    LUT = Vector{Tuple{Int,Int}}(undef, M)

    k = 0
    for j = 2:N
        for i = 1:j-1 
            
            k += 1
            LUT[k] = (i,j)
        end
    end
    @assert k == M

    return LUT
end

function updateedgeweight!(
    ws::Vector{T}, # mutates.
    w_ij::T,
    i::Int,
    j::Int,
    edges::Vector{Tuple{Int,Int}};
    verbose = true,
    ) where T
    
    ind = findfirst(xx->(xx==(i,j) || xx==(j,i)), edges)

    if isnothing(ind)
        if verbose
            println("Warning: cannot find ($i,$j) in edges. Did not update weight.")
        end
        return nothing
    end

    ws[ind] = w_ij

    return nothing
end

function maximizeME!(
    edge_weights::Vector{T}, # mutates.
    cs::Vector{T},
    edges::Vector{Tuple{Int,Int}};
    unique_cs_digits::Int = 6,
    unique_J_avg_digits::Int = 6,
    ) where T

    unique_cs, unique_cs_inds = uniqueinds(cs; digits = unique_cs_digits)
    subsets = unique_cs_inds

    LUT = getpairwiseLUT(length(subsets)) # pair-wise iterators.

    for (i,j) in LUT
        A = subsets[i]
        B = subsets[j]
    
        Js = collect(
            collect(
                getedgeweight(A[m], B[n], edges, edge_weights)
                for n in eachindex(B)
            )
            for m in eachindex(A)
        )
    
        J_avgs = computemean.(Js)
        unique_Ja, _ = uniqueinds(J_avgs; digits = unique_J_avg_digits)

        if length(unique_Ja) == 1
            # The average agree. replace J-coupling values by their average, J_new.
            J_new = unique_Ja[begin]
    
            for m in eachindex(A)
                for n in eachindex(B)
                    #@show A[m], B[n], J_new
                    # if J_new is zero, then we might not have the edge A[m] to B[n] in edges. Therefore, we turn verbose off.
                    updateedgeweight!(edge_weights, J_new, A[m], B[n], edges; verbose = false)
                end
            end
        end
    
    end

    return nothing
end

function computemean(x::Vector{T}) where T

    if isempty(x)
        return zero(T)
    end
    return sum(x) / length(x)
end

# mutates both inputs.
function pruneedgeweights!(ws::Vector{T}, edges::Vector{Tuple{Int,Int}}; atol::T = convert(T, 1e-6)) where T
    @assert length(ws) == length(edges)

    del_inds = Vector{Int}(undef, 0)
    for i in eachindex(ws)
        if isapprox(ws[i], zero(T); atol = atol)
            push!(del_inds, i)
        end
    end
    sort!(del_inds)

    deleteat!(ws, del_inds)
    deleteat!(edges, del_inds)

    return nothing
end


# # for an undirected graph. Since we deal with small spin systems, we'll use a dense matrix.
# function createadjacencymatrix(node_weights::Vector{T}, edges::Vector{Tuple{Int,Int}}, edge_weights::Vector{T}) where T
    
#     # # # sparse matrix version.
#     # # off-diagonals (edge weights)
#     # rs = collect( edges[m][begin] for m in eachindex(edges) )
#     # cs = collect( edges[m][end] for m in eachindex(edges) )
    
#     # # make symmetric
#     # sym_rs = vcat(rs, cs)
#     # sym_cs = vcat(cs, rs)

#     # # diagonals (node weights).
#     # ws = vcat(edge_weights, edge_weights)
#     # for i in eachindex(node_weights)
        
#     #     # add the main diagonal matrix indices.
#     #     push!(sym_rs, i)
#     #     push!(sym_cs, i)

#     #     # add edge weight.
#     #     push!(ws, node_weights[i])
#     # end

#     #return sparse(sym_rc, sym_cs, ws)

#     # we deal with small spin systems (length(node_weights) < 15), so use dense matrix.
#     N = length(node_weights)
#     A = zeros(T, N, N)
#     for k in eachindex(edges)
        
#         i, j = edges[k]
#         A[i,j] = edge_weights[k]
#         A[j,i] = edge_weights[k]
#     end

#     for i in eachindex(node_weights)
#         A[i,i] = node_weights[i]
#     end

#     return A
# end