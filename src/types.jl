
###### l-infinite norm as a distance metric for NearestNeighbors.jl.
function evalmaxnorm(x, y)
    return maximum(abs.(x-y))
end

# Declare the metric
struct MaxNormMetric <: NearestNeighbors.Metric end
(::MaxNormMetric)(x,y) = evalmaxnorm(x,y)

######

struct SHType{T} # output of the SH simulation.

    # resonance components in non-singlet spin systems.
    αs::Vector{Vector{T}}
    Ωs::Vector{Vector{T}}

    Δc_m_compound::Vector{Vector{Vector{T}}}
    part_inds_compound::Vector{Vector{Vector{Int}}}
    Δc_bar::Vector{Vector{Vector{T}}}

    N_spins_sys::Vector{Int}

    # resonance components in singlets spin systems.
    αs_singlets::Vector{T}
    Ωs_singlets::Vector{T}

    # misc input for the SH simulation.
    fs::T
    SW::T
    ν_0ppm::T
end

struct PhysicalParamsType{T}

    H_IDs::Vector{Int}

    H_inds_sys::Vector{Vector{Int}}
    cs_sys::Vector{Vector{T}}

    H_inds_singlets::Vector{Vector{Int}}
    cs_singlets::Vector{T}

    J_inds_sys::Vector{Vector{Tuple{Int,Int}}}
    J_inds_sys_local::Vector{Vector{Tuple{Int,Int}}} # redundant?
    J_vals_sys::Vector{Vector{T}}

    ME::Vector{Vector{Vector{Int}}}
end

function getphysicalparameters(target_names::Vector{String},
    H_params_path::String,
    dict_compound_to_filename;
    unique_cs_atol = 1e-6)

    N_compounds = length(target_names)
    Phys = Vector{PhysicalParamsType{Float64}}(undef, N_compounds)

    for n = 1:N_compounds

        # TODO add error-handling if name is not found in the dictionary, or filename does not exist.
        load_path = joinpath(H_params_path, dict_compound_to_filename[target_names[n]]["file name"])
        H_IDs, H_css, J_IDs, J_vals = loadcouplinginfojson(load_path)

        J_inds_sys, J_inds_sys_local, J_IDs_sys, J_vals_sys, H_inds_sys,
            cs_sys, H_inds_singlets, cs_singlets, H_inds, J_inds,
            g = setupcsJ(H_IDs, H_css, J_IDs, J_vals)

        ME, _ = getmageqinfo(H_IDs, H_css, J_IDs, J_vals; unique_cs_atol = unique_cs_atol)

        Phys[n] = PhysicalParamsType(H_IDs, H_inds_sys, cs_sys,
            H_inds_singlets, cs_singlets, J_inds_sys, J_inds_sys_local,
            J_vals_sys, ME)
    end

    return Phys
end
