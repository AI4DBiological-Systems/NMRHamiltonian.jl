function prunecoherences0(αs_inp, Ωs_inp, Δc_m_inp, c_states_inp, ϵ)
    #
    keep_flags = collect( allabssmaller(Δc_m_inp[l], 1 + ϵ) for l = 1:length(Δc_m_inp) )
    # println("Δc_m[1]]= ", Δc_m[1]) # debug.
    # println("1+simple_coherence_atol = ", 1+simple_coherence_atol) # debug.
    # println("keep_flags[1] = ", keep_flags[1]) # debug.
    αs = αs_inp[keep_flags]
    Ωs = Ωs_inp[keep_flags]
    Δc_m = Δc_m_inp[keep_flags]
    c_states = c_states_inp[keep_flags]

    # need a bit more tolerance for the sum filter, otherwise might filter too many components.
    N_DOF = length(Δc_m_inp[1])
    keep_flags2 = collect( sumabssmaller(Δc_m[l], 1 + ϵ*N_DOF ) for l = 1:length(Δc_m) )

    return αs, Ωs, Δc_m, c_states, keep_flags2
end

function addentriestogroupsbyproximity(centers::Vector{T}, f, X::Vector{T}, radius;
    skip_flags = falses(length(X))) where T

    M = length(X)
    assignment_inds = zeros(Int, M)

    for m = 1:M
        if !skip_flags[m]
            x = X[m]
            ind = findfirst(zz->(f(x, zz) < radius), centers)

            if typeof(ind) == Int
                assignment_inds[m] = ind
            end
        end
    end

    return assignment_inds
end

function getprunedcoherences(αs_inp, Ωs_inp, Δc_m_inp, c_states_inp, keep_flags)

    αs = αs_inp[keep_flags]
    Ωs = Ωs_inp[keep_flags]
    Δc_m = Δc_m_inp[keep_flags]
    c_states = c_states_inp[keep_flags]

    return αs, Ωs, Δc_m, c_states
end

# mutates αs, Ωs, Δc_m, c_states, part_inds.
function updateresonancegroups!(αs, Ωs, Δc_m, c_states, part_inds,
    αs_inp, Ωs_inp, Δc_m_inp, c_states_inp, assignment_inds)

    @assert length(αs) == length(Ωs) == length(Δc_m) == length(c_states)

    M = length(assignment_inds)
    for m = 1:M
        ind = assignment_inds[m]

        if ind > 0
            push!(αs, αs_inp[m])
            push!(Ωs, Ωs_inp[m])
            push!(Δc_m, Δc_m_inp[m])
            push!(c_states, c_states_inp[m])

            push!(part_inds[ind], length(αs))
        end
    end

    return nothing
end

function prunecoherencemain(αs_inp::Vector{T}, Ωs_inp::Vector{T},
    Δc_m::Vector{Vector{T}}, c_states_inp, ϵ::T, option::Int,
    override_intensity::T, Δc_partition_radius::T, min_α_tol::T;
    min_threshold = -1.1) where T
    # println("min_α_tol = ", min_α_tol)
    # @assert 1==43
    αs_int, Ωs_int, Δc_m_int, c_states_int,
        keep_flags2 = prunecoherences0(αs_inp, Ωs_inp, Δc_m, c_states_inp, ϵ)

    part_inds = Vector{Vector{Int}}(undef, 0)
    Δc_centroids = Vector{Vector{T}}(undef, 0)

    αs, Ωs, Δc_m, c_states = αs_int, Ωs_int, Δc_m_int, c_states_int
    if option == 1
        αs, Ωs, Δc_m, c_states = getprunedcoherences(αs_int, Ωs_int, Δc_m_int, c_states_int, keep_flags2)

        return αs, Ωs, Δc_m, c_states, part_inds, Δc_centroids

    elseif option == 2

        keep_flags3 = αs_int .> override_intensity
        keep_flags = keep_flags2 .| keep_flags3
        αs, Ωs, Δc_m, c_states = getprunedcoherences(αs_int, Ωs_int, Δc_m_int,
            c_states_int, keep_flags)

        return αs, Ωs, Δc_m, c_states, part_inds, Δc_centroids

    elseif option == 3 || option == 4 || option == 5 || option == 6

        if !isempty(Δc_m)
            αs, Ωs, Δc_m, c_states = getprunedcoherences(αs_int, Ωs_int, Δc_m_int, c_states_int, keep_flags2)

            part_inds, Δc_centroids = partitionresonancesbyneighbors(Δc_m,
                αs, min_α_tol; radius = Δc_partition_radius)

            if option == 3 # allow Δc that are close to a Δc_bar wrt maxnorm.
                assignment_inds = addentriestogroupsbyproximity(Δc_centroids, evalmaxnorm,
                Δc_m_int, Δc_partition_radius; skip_flags = keep_flags2)

            elseif option == 4 # union of 3 and "close" to a combo coherence. wrt isapproxcombocoherence()
                assignment_inds = addentriestogroups4(Δc_centroids, evalmaxnorm,
                Δc_m_int, Δc_partition_radius; skip_flags = keep_flags2,
                atol = ϵ, cancel_tol = ϵ)

            elseif option == 5 # allow combo coherence that are similar to an existing Δc_bar.
                if min_threshold < 0
                    largest_Δc_components = collect( maximum(abs.(Δc_centroids[k])) for k = 1:length(Δc_centroids))
                    min_threshold = minimum(largest_Δc_components)/2
                    #min_threshold = 0.99 # debug.
                end

                assignment_inds = addentriestogroups5(Δc_centroids, evalmaxnorm,
                Δc_m_int, Δc_partition_radius; skip_flags = keep_flags2,
                min_threshold = min_threshold, zero_tol_per_spin = ϵ)

            elseif option == 6 # union 5 with 4.
                if min_threshold < 0
                    largest_Δc_components = collect( maximum(abs.(Δc_centroids[k])) for k = 1:length(Δc_centroids))
                    min_threshold = minimum(largest_Δc_components)/2
                    #min_threshold = 0.99 # debug.
                end

                assignment_inds = addentriestogroups6(Δc_centroids, evalmaxnorm,
                Δc_m_int, Δc_partition_radius; skip_flags = keep_flags2,
                min_threshold = min_threshold, zero_tol_per_spin = ϵ)
            end

            updateresonancegroups!(αs, Ωs, Δc_m, c_states, part_inds,
                αs_int, Ωs_int, Δc_m_int, c_states_int, assignment_inds)

            return αs, Ωs, Δc_m, c_states, part_inds, Δc_centroids
        end
    end

    return αs, Ωs, Δc_m, c_states, part_inds, Δc_centroids
end




# at least one entry is very similar.
function isapproxcombocoherence0(c_test::Vector{T},
    Δc::Vector{T}; atol = 1e-2, cancel_tol = 1e-2) where T

    # per-element relative discrepancy.
    d = Δc-c_test
    relative_abs_d = abs.(d)./abs.(Δc)
    min_val, min_ind = findmin(relative_abs_d)

    if min_val < atol

        deleteat!(d, min_ind)
        if sum(d) < cancel_tol
            return true
        end
    end

    return false
end

function isapproxcombocoherence(c_test::Vector{T},
    Δc::Vector{T}; atol = 1e-2, cancel_tol = 1e-2) where T

    # per-element relative discrepancy.
    d = Δc-c_test
    relative_abs_d = abs.(d)./abs.(Δc)
    #min_val, min_ind = findmin(relative_abs_d)

    # find elements for testing.
    keep_flags = relative_abs_d .> atol

    if count(keep_flags) < length(Δc)

        # at least one element of c_test is a good match with Δc.
        if sum(d[keep_flags]) < cancel_tol
            return true
        end
    end

    return false
end

# by proximity or is a combo coherence of an entry of centers.
function addentriestogroups4(centers::Vector{T}, f, X::Vector{T}, radius;
    skip_flags = falses(length(X)), atol = atol, cancel_tol = cancel_tol) where T

    M = length(X)
    assignment_inds = zeros(Int, M)

    for m = 1:M
        if !skip_flags[m]
            x = X[m]
            ind = findfirst(zz->(f(x, zz) < radius), centers)

            if typeof(ind) == Int
                assignment_inds[m] = ind
            else
                status = collect( isapproxcombocoherence(x, centers[j]; atol = atol, cancel_tol = cancel_tol) for j = 1:length(centers) )

                if any(status)
                    assignment_inds[m] = findfirst(status)
                    # if count(status) > 1
                    #     println("x = ", x)
                    #     println("m = ", m)
                    #     println("centers = ", centers)
                    #     println("status = ", status)
                    #     println()
                    # end
                end
            end

        end
    end

    return assignment_inds
end

function addentriestogroups6(centers::Vector{T}, f, X::Vector{T}, radius;
    skip_flags = falses(length(X)),
    min_threshold = 0.2,
    zero_tol_per_spin = 1e-2) where T

    M = length(X)
    assignment_inds = zeros(Int, M)

    zero_tol = length(X[1]) * zero_tol_per_spin

    for m = 1:M
        if !skip_flags[m]
            x = X[m]
            ind = findfirst(zz->(f(x, zz) < radius), centers)

            if typeof(ind) == Int
                assignment_inds[m] = ind
            else

                Us = getcentermatrices(centers; min_threshold = min_threshold)
                ys = collect( Us[k]*x for k = 1:length(Us) )
                min_val, min_ind = findmin(abs.(sum.(ys)))

                if min_val < zero_tol
                    assignment_inds[m] = min_ind
                end
            end

        end
    end

    return assignment_inds
end

function addentriestogroups5(centers::Vector{T}, f, X::Vector{T}, radius;
    skip_flags = falses(length(X)),
    min_threshold = 0.2,
    zero_tol_per_spin = 1e-2) where T

    M = length(X)
    assignment_inds = zeros(Int, M)

    zero_tol = length(X[1]) * zero_tol_per_spin
    Us = getcentermatrices(centers; min_threshold = min_threshold)

    for m = 1:M
        if !skip_flags[m]
            x = X[m]
            # ind = findfirst(zz->(f(x, zz) < radius), centers)
            #
            # if typeof(ind) == Int
            #     assignment_inds[m] = ind
            # else

                ys = collect( Us[k]*x for k = 1:length(Us) )
                min_val, min_ind = findmin(abs.(sum.(ys)))

                if min_val < zero_tol
                    # println("centers = ", centers)
                    # println("x = ", x)
                    # println()
                    assignment_inds[m] = min_ind
                end
            #end

        end
    end

    return assignment_inds
end

function findmaincontributors(c::Vector{T}; min_threshold::T = 0.2) where T
    inds = Vector{Int}(undef, 0)

    for k = 1:length(c)
        if abs(c[k]) > min_threshold
            push!(inds, k)
        end
    end

    return inds
end

# replace the diagonals of the identity
function getcentermatrices(cbs::Vector{Vector{T}}; min_threshold = 0.2) where T
    K = length(cbs)

    Qs = Vector{Matrix{T}}(undef, K)

    for k = 1:K
        cb = cbs[k]
        Q = collect(LinearAlgebra.I(length(cb)))

        ## single.
        # val, ind = findmax(abs.(cb))
        # Q[ind,ind] = 0

        ## multiple.
        inds = findmaincontributors(cb; min_threshold = min_threshold)
        for i in inds
            Q[i,i] = 0
        end

        Qs[k] = Q
    end

    return Qs
end

########
function allabssmaller(x, threshold)
    if maximum(abs.(x)) > threshold
        return false
    end

    return true
end

function sumabssmaller(x, threshold)
    if sum(abs.(x)) > threshold
        return false
    end

    return true
end

function prunecombocoherences!(A, relative_α_threshold, coherence_tol, Δc_partition_radius)

    for i = 1:length(A.Δc)

        keep_flags = collect( allabssmaller(A.Δc[i][l],1+coherence_tol) for l = 1:length(A.Δc[i]) )
        #keep_flags = collect( sumabssmaller(A.Δc[i][l],1+coherence_tol) for l = 1:length(A.Δc[i]) )

        A.αs[i] = A.αs[i][keep_flags]
        A.Ωs[i] = A.Ωs[i][keep_flags]
        A.Δc[i] = A.Δc[i][keep_flags]

        # remake Δc_bar.
        α_tol = relative_α_threshold*maximum(A.αs[i])
        part_inds, Δc_centroids = partitionresonancesbyneighbors(A.Δc[i],
            A.αs[i], α_tol; radius = Δc_partition_radius)

        A.Δc_bar[i] = Δc_centroids
        A.parts[i] = part_inds
    end

    return nothing
end

function prunecombocoherencesbar!(A, relative_α_threshold, coherence_tol, Δc_partition_radius)

    for i = 1:length(A.Δc)

        keep_flags = collect( allabssmaller(A.Δc[i][l],1+coherence_tol) for l = 1:length(A.Δc[i]) )
        #keep_flags = collect( sumabssmaller(A.Δc[i][l],1+coherence_tol) for l = 1:length(A.Δc[i]) )

        A.αs[i] = A.αs[i][keep_flags]
        A.Ωs[i] = A.Ωs[i][keep_flags]
        A.Δc[i] = A.Δc[i][keep_flags]

        # remake Δc_bar.
        α_tol = relative_α_threshold*maximum(A.αs[i])
        part_inds, Δc_centroids = partitionresonancesbyneighbors(A.Δc[i],
            A.αs[i], α_tol; radius = Δc_partition_radius)

        #
        keep_flags = collect( sumabssmaller(Δc_centroids[l],1+coherence_tol) for l = 1:length(Δc_centroids) )

        A.Δc_bar[i] = Δc_centroids[keep_flags]
        A.parts[i] = part_inds[keep_flags]
    end

    return nothing
end

function prunecombocoherencesbar2!(Δc,
    αs, Ωs, Δc_bar, parts;
    relative_α_threshold = 0.05, coherence_tol = 1e-2, Δc_partition_radius = 0.3)

    for i = 1:length(Δc)

        keep_flags = collect( allabssmaller(Δc[i][l],1+coherence_tol) for l = 1:length(Δc[i]) )

        αs[i] = αs[i][keep_flags]
        Ωs[i] = Ωs[i][keep_flags]
        Δc[i] = Δc[i][keep_flags]
        #println("1 count(keep_flags) = ", count(keep_flags) )

        # remake Δc_bar.
        α_tol = relative_α_threshold*maximum(αs[i])
        part_inds, Δc_centroids = partitionresonancesbyneighbors(Δc[i],
            αs[i], α_tol; radius = Δc_partition_radius)

        #
        keep_flags = collect( sumabssmaller(Δc_centroids[l],1+coherence_tol) for l = 1:length(Δc_centroids) )
        #println("2 count(keep_flags) = ", count(keep_flags) )

        Δc_bar[i] = Δc_centroids[keep_flags]
        parts[i] = part_inds[keep_flags]
    end

    return nothing
end
