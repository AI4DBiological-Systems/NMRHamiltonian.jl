"""
Creates ordering of length N from equivalent indices.
`eq_inds` can only contain unique values from 1:N.
When used for a spin system, N is the number of nuclei in the system.

Returns the indices for the N nuclei, and the number of degrees of freedom.

Example:
import Destruct
tmp = collect( createorderingfromeqinds(ME[i], A.N_spins_sys[i]) for i = 1:length(A.N_spins_sys) )
κs_β_ordering, κs_β_DOF = Destruct.destruct(tmp)
"""
function createorderingfromeqinds(eq_inds::Vector{Vector{Int}}, N::Int)

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


"""
ordering must contain integers from 1:N, where N is the number of unique entries in ordering.

Example:
Δc_m_compound_i_l = randn(10)
κs_β_ordering_i = [ 2; 2; 2; 3; 4; 5;  1; 1; 1; 6]
y = condensenuclei(Δc_m_compound_i_l, κs_β_ordering_i)

"""
function condensenuclei(x::Vector{T}, ordering::Vector{Int}, N::Int)::Vector{T} where T
    @assert length(x) == length(ordering)
    @assert norm(collect(1:N) - sort(unique(ordering))) < 1e-14

    y = zeros(T, N)
    #Ns = zeros(Int, N)
    for i = 1:length(x)

        k = ordering[i]
        y[k] += x[i]
        #Ns[k] += 1
    end

    #return y ./ Ns
    return y
end

function condensenuclei(x::Vector{T}, ordering::Vector{Int})::Vector{T} where T
    @assert length(x) == length(ordering)

    N = length(unique(ordering))

    return condensenuclei(x, ordering, N)
end

"""
unused?
returns a reduced length (in the most nested level) version of Δc_m_compound.
ME is equivalence indices. Could be an entry of MEs, the return type of getmageqinfomixture().
Δc_m_compound is a field in the SHType composite type.
"""
function reduceΔc(Δc_m_compound::Vector{Vector{Vector{T}}},
    ME::Vector{Vector{Vector{Int}}},
    N_spins_sys::Vector{Int}) where T

    N_sys = length(Δc_m_compound)
    @assert length(ME) == N_sys == length(N_spins_sys)

    # prepare.
    κs_β_ordering = Vector{Vector{Int}}(undef, N_sys)
    κs_β_DOF = Vector{Int}(undef, N_sys)
    for i = 1:N_sys
        κs_β_ordering[i], κs_β_DOF[i] = createorderingfromeqinds(ME[i], N_spins_sys[i])
    end

    # condense.
    c = Vector{Vector{Vector{T}}}(undef, N_sys)

    for i = 1:N_sys # over spin systems.

        if !isempty(Δc_m_compound[i])

            c[i] = Vector{Vector{T}}(undef, length(Δc_m_compound[i]))

            for l = 1:length(Δc_m_compound[i]) # over existing groups in Δc_m_compound.

                for k = 1:length(ME[i])
                    c[i][l] = condensenuclei(Δc_m_compound[i][l], κs_β_ordering[i], κs_β_DOF[i])
                end
            end
        else
            c[i] = Δc_m_compound[i]
        end
    end

    return c
end

"""
returns a reduced length (in the most nested level) version of Δc_m_compound.
ME_m is equivalence indices for a spin system. Could be a spin group from an entry of MEs, the return type of getmageqinfomixture().
Δc_m is a field in a spin group of the SHType composite type.

Example: (work on later, see reduce/batch.jl)
c = NMRSpectraSimulator.reduceΔc(A.Δc_m_compound, ME, A.N_spins_sys)
c2 = NMRSpectraSimulator.reduceΔc(A.Δc_m_compound[1], ME[1], A.N_spins_sys[1])

"""
function reduceΔc(Δc_m::Vector{Vector{T}},
    ME_m::Vector{Vector{Int}},
    N_spins::Int) where T

    if isempty(ME_m)
        return Δc_m
    end

    # prepare.
    ordering, DOF = createorderingfromeqinds(ME_m, N_spins)

    # condense.
    c = Vector{Vector{T}}(undef, length(Δc_m))
    for l = 1:length(Δc_m) # over existing groups in Δc_m_compound.

        for k = 1:length(ME_m)
            c[l] = condensenuclei(Δc_m[l], ordering, DOF)
        end
    end

    return c
end
