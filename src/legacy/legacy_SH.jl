##### backup.

function computeSHoldversion(
    css_sys::Vector{Vector{T}},
    J_vals_sys,
    J_inds_sys,
    intermediates_sys,
    ppm2hzfunc,
    cs_singlets,
    N_spins_singlet::Vector{Int};
    coherence_tols::Vector{T} = collect( 1e-2 for _ in eachindex(css_sys)),
    normalize_αs::Vector{Bool} = trues(length(css_sys)),
    ) where T <: AbstractFloat

    N_sys = length(css_sys)
    @assert length(J_vals_sys) == length(J_inds_sys) == N_sys == length(coherence_tols) == length(normalize_αs)
    #
    
    N_singlets = length(cs_singlets)

    N_groups = N_sys + N_singlets

    # non-singlet objects.
    states_sys = Vector{Vector{Int}}(undef, N_sys)
    coherence_state_pairs_sys = Vector{Vector{Tuple{Int,Int}}}(undef, N_sys)
    eig_vals_sys = Vector{Vector{T}}(undef, N_sys)
    Q_sys = Vector{Vector{Vector{T}}}(undef, N_sys)
    H_sys = Vector{Matrix{T}}(undef, N_sys)
    coherence_mat_sys = Vector{Matrix{T}}(undef, N_sys)
    ms_sys = Vector{Vector{Vector{T}}}(undef, N_sys)
    M_sys = Vector{Vector{T}}(undef, N_sys)

    # output buffers.
    αs = Vector{Vector{T}}(undef, N_groups)
    Ωs = Vector{Vector{T}}(undef, N_groups)

    # Spin systems.
    for i = 1:N_sys

        αs[i], Ωs[i], coherence_mat_sys[i], eig_vals_sys[i], Q_sys[i],
        coherence_state_pairs_sys[i],
        H_sys[i], H1, H2, M_sys[i] = evalSCalgorithm(
            css_sys[i],
            J_vals_sys[i],
            J_inds_sys[i],
            intermediates_sys[i],
            ppm2hzfunc;
            coherence_tol = coherence_tols[i],
        )

        # normalize αs[i] according to number of spins.
        if normalize_αs[i]
            N_spins = length(css_sys[i])
            normalizetoNspins!(αs[i], N_spins)
        end

        Id = getsingleId()
        ms_sys[i] = computequantumnumbers(Q_sys[i], Id)

        # get the unique list of all states that appear in the coherence state paris.
        c_pairs = coherence_state_pairs_sys[i]
        tmp1 = collect( c_pairs[j][1] for j in eachindex(c_pairs) )
        tmp2 = collect( c_pairs[j][2] for j in eachindex(c_pairs) )
        states_sys[i] = unique([tmp1; tmp2])
    end

    # singlets. evalRayleighproxy!() updates singlets for ΩS.
    for i = 1:N_singlets
        αs[i+N_sys] = [ N_spins_singlet[i] ]
        Ωs[i+N_sys] = [ ppm2hzfunc.(cs_singlets[i]) .* (2*π) ]
    end

    return αs, Ωs, coherence_mat_sys, eig_vals_sys, Q_sys,
    coherence_state_pairs_sys, H_sys, states_sys, ms_sys, M_sys
end


# for now, assign p is Δcs, and hard code the conversion of
#   Δcs to ΔΩ0 as unit conversion for fullsolve.
# Here, updateΩfuncs takes Δcs instead of cs.
## TODO later: generic version where p reads into p_cs, then converted to css_sys.
#   the units of p_cs needs to be defined before calling this func.
function evalfullsolveproxy!(css_sys_buffer::Vector{Vector{T}},
    cs_singlets_compact0::Vector{T},
    #N_singlets::Int,
    Ωs::Vector{Vector{T}},
    p::Vector{T},
    st_ind::Int,
    cs_LUT::Vector{Vector{Vector{Int}}},
    updateΩfuncs::Vector{Function},
    fs::T,
    SW::T,
    ppm2hzfunc) where T <: AbstractFloat
### checks.
N_sys = length(css_sys_buffer)
@assert length(updateΩfuncs) == N_sys == length(cs_LUT)

N_singlets = length(cs_singlets_compact0)

N_groups = N_sys + N_singlets
@assert length(Ωs) == N_groups

### set up.
fin_ind = st_ind -1

### multi-spin systems.
for i = 1:N_sys

## parse to css_sys.
N_cs_vars = length(cs_LUT[i])

st_ind = fin_ind +1
fin_ind = st_ind + N_cs_vars -1
p_cs = p[st_ind:fin_ind]
pcstocs!(css_sys_buffer[i], p_cs, cs_LUT[i])

## convert units.
convertΔcstoΔω0!(css_sys_buffer[i], fs, SW)

#updaterfuncs(ω0)
updateΩfuncs[i](css_sys_buffer[i])
end

### singlets.
for i = 1:N_singlets

## parse.
fin_ind += 1

# ## convert units and update.
# Ωs[i+N_sys][1] = p[fin_ind]*2*π*fs/SW

# p is Δcs, and new cs = old cs + Δcs.
tmp = cs_singlets_compact0[i] + p[fin_ind]

## update.
Ωs[i+N_sys][1] = ppm2hzfunc(tmp) * 2*π
end

return fin_ind
end

function convertΔcstoΔω0!(x::Vector{T}, fs::T, SW::T) where T
for i = 1:length(x)
x[i] = x[i]*2*π*fs/SW
end
end

function convertΔcstoΔω0(x::T, fs::T, SW::T)::T where T
return x*2*π*fs/SW
end