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
