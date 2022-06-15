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


# different parameterizations of the spin system FID parameters.


struct SpinSysParamsType2{T}
    κs_λ::Vector{Vector{T}} # a multiplier for each (spin group, partition element).
    κs_β::Vector{Vector{T}} # a vector coefficient for each (spin group). vector length: number of spins in the spin group.
    d::Vector{Vector{T}} # a multiplier for each (spin group, partition element).
end

function SpinSysParamsType2(x::T) where T
    return SpinSysParamsType2(Vector{Vector{T}}(undef,0), Vector{Vector{T}}(undef, 0), Vector{Vector{T}}(undef, 0))
end

struct SpinSysParamsType1{T}
    κs_λ::Vector{T} # a common multiplier for each spin group. length: number of spin groups.
    κs_β::Vector{Vector{T}} # a vector coefficient for each (spin group). vector length: number of spins in the spin group.
    d::Vector{T} # a common multiplier for each spin group. length: number of spin groups.
end

function SpinSysParamsType1(x::T) where T
    return SpinSysParamsType1(Vector{T}(undef,0), Vector{Vector{T}}(undef, 0), Vector{T}(undef, 0))
end

function constructorSSFID(x::SpinSysParamsType1{T}, y...)::SpinSysParamsType1{T} where T
    return SpinSysParamsType1(y...)
end

# function constructorSSFID(x::SpinSysFIDType2{T}, y...)::SpinSysFIDType2{T} where T
#     return SpinSysFIDType2(y...)
# end
