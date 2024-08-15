# """
# function createorderingfromeqinds(eq_inds::Vector{Vector{Int}}, N::Int)::Tuple{Vector{Int},Int}

# Creates ordering of length N from equivalent indices.
# `eq_inds` can only contain unique values from 1:N.
# When used for a spin system, N is the number of nuclei in the system.

# Returns the indices for the N nuclei, and the number of degrees of freedom.

# Example:
# import Destruct
# tmp = collect( createorderingfromeqinds(ME[i], A.N_spins_sys[i]) for i = 1:length(A.N_spins_sys) )
# κs_β_ordering, κs_β_DOF = Destruct.destruct(tmp)
# """
function createorderingfromeqinds(eq_inds::Vector{Vector{Int}}, N::Int)::Tuple{Vector{Int},Int}

    j = 0 # degrees of freedom counter.
    out = zeros(Int, N)

    if isempty(eq_inds)
        return collect(1:N), N
    end

    # check if `eq_inds` only contains unique values from 1:N.
    if all( all( (eq_inds[k][l] in 1:N) for l = 1:length(eq_inds[k])) for k = 1:length(eq_inds) )
        for k = 1:length(eq_inds)
            j += 1

            out[eq_inds[k]] .= j
        end

        # fill the rest.
        for n = 1:N
            if out[n] == 0
                # out[n] is unassigned.
                j += 1
                out[n] = j
            end
        end

        return out, j
    end

    println("Invalid ME, using default.")
    return collect(1:N), N
end


# """
# ordering must contain integers from 1:N, where N is the number of unique entries in ordering.

# Example:
# Δc_i_l = randn(10)
# κs_β_ordering_i = [ 2; 2; 2; 3; 4; 5;  1; 1; 1; 6]
# y = condensenuclei(Δc_i_l, κs_β_ordering_i)

# """
function condensenuclei(
    x::Vector{T},
    ordering::Vector{Int},
    N::Int
    )::Vector{T} where T

    @assert length(x) == length(ordering)
    @assert norm(collect(1:N) - sort(unique(ordering))) < 1e-14

    y = zeros(T, N)
    #Ns = zeros(Int, N)
    for i in eachindex(x)

        k = ordering[i]
        y[k] += x[i] # add the coherence contributions among the equivalent resonances.
        #Ns[k] += 1
    end

    #return y ./ Ns
    return y
end

function condensenuclei(
    x::Vector{T},
    ordering::Vector{Int}
    )::Vector{T} where T
    
    @assert length(x) == length(ordering)

    N = length(unique(ordering))

    return condensenuclei(x, ordering, N)
end


# """
# returns a reduced length (in the most nested level) version of Δc.
# ME_m is equivalence indices for a spin system. Could be a spin group from an entry of MEs, the return type of getmageqinfomixture().
# Δc_m is a field in a spin group of the SHType composite type.

# Example: (work on later, see reduce/batch.jl)
# c = NMRSpectraSimulator.reduceΔc(A.Δc, ME, A.N_spins_sys)
# c2 = NMRSpectraSimulator.reduceΔc(A.Δc[1], ME[1], A.N_spins_sys[1])

# """
function reduceΔc(
    Δc_m::Vector{Vector{T}},
    ME_m::Vector{Vector{Int}},
    N_spins::Integer,
    cs::Vector{T},
    ) where T

    if isempty(ME_m)
        return Δc_m
    end

    # prepare.
    ordering, DOF = createorderingfromeqinds(ME_m, N_spins)

    cs_reduced = Vector{T}(undef, DOF)
    for i in eachindex(ordering)
        k = ordering[i]
        cs_reduced[k] = cs[i]
    end
    
    # condense.
    dc = Vector{Vector{T}}(undef, length(Δc_m))
    for l = 1:length(Δc_m) # over existing groups in Δc.

        #for k = 1:length(ME_m)
            dc[l] = condensenuclei(Δc_m[l], ordering, DOF)
        #end
    end

    return dc, cs_reduced
end


################## for testing only.


# """
# function getcoherencedropviolations(Δc::Vector{Vector{Vector{T}}}) where T

# Utility function for testing whether the `Δc_bar` and `Δc` fields from variables of data type `SHType` sum to `-1`.

# Example:
# ```
# chks_m, status_m = getcoherencedropviolations(As[1].Δc)
# chks_bar, status_bar = getcoherencedropviolations(As[1].Δc_bar)
# ```
# """
function getcoherencedropviolations(
    Δc::Vector{Vector{Vector{T}}};
    atol::T = convert(T, 1e-1),
    ) where T

    no_violations = true
    out = Vector{Tuple{Int,Int}}(undef, 0)
    for i in eachindex(Δc)

        for k in eachindex(Δc[i])
            if !isapprox(sum(Δc[i][k]), -one(T); atol = atol)
                push!(out, (i,k))
                no_violations = false
            end
        end
    end

    return out, no_violations
end

function getcoherencemagnitudeviolations(
    Δc::Vector{Vector{Vector{T}}};
    atol::T = convert(T, 1e-1),
    ) where T

    lb = -one(T) -atol
    ub = -lb

    no_violations = true

    out = Vector{Tuple{Int,Int}}(undef, 0)
    for i in eachindex(Δc)

        for k in eachindex(Δc[i])
            
            x = Δc[i][k]

            if !all(lb .< x .< ub)
                no_violations = false
                push!(out, (i,k))
            end
        end
    end

    return out, no_violations
end

function getNsys(A::SHType)::Int
    return length(A.Δc_bar)
end

# for a mixture
struct CoherenceDiagnostics

    inds_drop_Δc_m::Vector{Vector{Vector{Tuple{Int,Int}}}}
    inds_drop_Δc_bar::Vector{Vector{Vector{Tuple{Int,Int}}}}
    
    valid_drop_Δc_m::Bool
    valid_drop_Δc_bar::Bool

    inds_magnitude_Δc_m::Vector{Vector{Vector{Tuple{Int,Int}}}}
    inds_magnitude_Δc_bar::Vector{Vector{Vector{Tuple{Int,Int}}}}
    
    valid_magnitude_Δc_m::Bool
    valid_magnitude_Δc_bar::Bool
end

function extractdropstatus(A::CoherenceDiagnostics)
    return A.valid_drop_Δc_m, A.valid_drop_Δc_bar
end

function extractmagnitudestatus(A::CoherenceDiagnostics)
    return A.valid_magnitude_Δc_m, A.valid_magnitude_Δc_bar
end

# """
# function getΔcdiagnostics(
#     As::Vector{SHType{T}};
#     atol::T = convert(T, 1e-1),
#     )::CoherenceDiagnostics where T

# Container of diagnostics for whether every Δc vector sums to -1 within +/- atol, and whether the components have a magnitude that is < 1, within atol.

# Example
# ```
# coherence_diagnostics = getΔcdiagnostics(As)
# ```
# """
function getΔcdiagnostics(
    As::Vector{SHType{T}};
    atol::T = convert(T, 1e-1),
    )::CoherenceDiagnostics where T

    N = length(As)
    
    inds_drop_Δc_m = Vector{Vector{Vector{Tuple{Int,Int}}}}(undef, N)
    inds_drop_Δc_bar = Vector{Vector{Vector{Tuple{Int,Int}}}}(undef, N)
    
    valid_drop_Δc_m = true
    valid_drop_Δc_bar = true

    inds_magnitude_Δc_m = Vector{Vector{Vector{Tuple{Int,Int}}}}(undef, N)
    inds_magnitude_Δc_bar = Vector{Vector{Vector{Tuple{Int,Int}}}}(undef, N)
    
    valid_magnitude_Δc_m = true
    valid_magnitude_Δc_bar = true

    for n in eachindex(As)
        N_sys = getNsys(As[n])

        inds_drop_Δc_m[n] = Vector{Vector{Tuple{Int,Int}}}(undef, N_sys)
        inds_drop_Δc_bar[n] = Vector{Vector{Tuple{Int,Int}}}(undef, N_sys)

        inds_magnitude_Δc_m[n] = Vector{Vector{Tuple{Int,Int}}}(undef, N_sys)
        inds_magnitude_Δc_bar[n] = Vector{Vector{Tuple{Int,Int}}}(undef, N_sys)

        for i in eachindex(As[n].Δc)

            inds_drop_Δc_m[n][i], status = getcoherencedropviolations(As[n].Δc)
            valid_drop_Δc_m = valid_drop_Δc_m && status

            inds_drop_Δc_bar[n][i], status = getcoherencedropviolations(As[n].Δc_bar)
            valid_drop_Δc_bar = valid_drop_Δc_bar && status

            inds_magnitude_Δc_m[n][i], status = getcoherencemagnitudeviolations(
                As[n].Δc;
                atol = atol,
            )
            valid_magnitude_Δc_m = valid_magnitude_Δc_m && status

            inds_magnitude_Δc_bar[n][i], status = getcoherencemagnitudeviolations(
                As[n].Δc_bar;
                atol = atol,
            )
            valid_magnitude_Δc_bar = valid_magnitude_Δc_bar && status
        end
    end

    return CoherenceDiagnostics(
        inds_drop_Δc_m,
        inds_drop_Δc_bar,
        valid_drop_Δc_m,
        valid_drop_Δc_bar,
        inds_magnitude_Δc_m,
        inds_magnitude_Δc_bar,
        valid_magnitude_Δc_m,
        valid_magnitude_Δc_bar,
    )
end


function checkcoherences(
    As::Vector{SHType{T}};
    atol::T = convert(T, 1e-1),
    ) where T
    
    coherence_diagnostics = getΔcdiagnostics(As; atol = atol)
    c_drop, bar_drop = extractdropstatus(coherence_diagnostics)
    c_mag, bar_mag = extractmagnitudestatus(coherence_diagnostics)
    
    #out::Bool = c_mag && bar_mag && c_drop && bar_drop
    out::Bool = c_drop && bar_drop # the sum(Δc) = -1 is more important than the magnitude < 1.
    
    return out, coherence_diagnostics
end

function checkcoherences(
    MSPs::Vector{MoleculeSpinSystem{T}};
    coherence_sum_zero_tol::T = convert(T, 1e-14),
    ) where T

    for n in eachindex(MSPs)
        for i in eachindex(MSPs[n].spin_systems)
            sp = MSPs[n].spin_systems[i]

            LHS = norm(sum.(sp.partial_quantum_numbers)-sp.quantum_numbers)
            if LHS > coherence_sum_zero_tol
                return false
            end
        end
    end

    return true
end