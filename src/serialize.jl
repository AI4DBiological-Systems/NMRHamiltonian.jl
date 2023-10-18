# serialize and deserialization to Dict(Symbol, Any) for data types that only contain native Vector{T} types, where T is a concrete type.
# This improves the chance of compatibility across Julia versions.
# The user of this library needs to write code to save the Dict{Symbol, Any} to disk. JSON3 and BSOn are used in the test bench, but not included as formal dependencies for NMRHamiltonian.jl.
# Since JSON3 is already included as dependencies, there is a routine for saving and loading the serialized object Dict{Symbol, Any} to JSOn.


function recovertype(tag::String)::DataType
    return eval(Meta.parse(tag))
end

function recovertype(tag::DataType)::DataType
    return tag
end


### mixture: As::SHType{T}.

"""
```
function serializemixture(
    As::Vector{SHType{T}},
) where T
```
Returns a `Dict{Symbol, Any}` object, which can then be used with file serialization packages for saving to disk.
For use with `saveasJSOn()`; see its docstring for more details.
"""
function serializemixture(
    As::Vector{SHType{T}},
    ) where T

    if isempty(As)
        println("input is an empty collection. Cannot serialize.")
    else

        N = length(As)

        αs_set = Vector{typeof(As[begin].αs)}(undef, N)
        Ωs_set = Vector{typeof(As[begin].Ωs)}(undef, N)

        Δc_set = Vector{typeof(As[begin].Δc)}(undef, N)
        parts_set = Vector{typeof(As[begin].parts)}(undef, N)
        Δc_bar_set = Vector{typeof(As[begin].Δc_bar)}(undef, N)
        N_spins_sys_set = Vector{typeof(As[begin].N_spins_sys)}(undef, N)

        for n in eachindex(As)
            αs_set[n] = As[n].αs
            Ωs_set[n] = As[n].Ωs

            Δc_set[n] = As[n].Δc
            parts_set[n] = As[n].parts
            Δc_bar_set[n] = As[n].Δc_bar
            N_spins_sys_set[n] = As[n].N_spins_sys
        end
        
        fs = As[begin].fs
        SW = As[begin].SW
        ν_0ppm = As[begin].ν_0ppm

        return Dict(
            # values.
            :αs_set => αs_set,
            :Ωs_set => Ωs_set,
    
            :Δc_set => Δc_set,
            :parts_set => parts_set,
            :Δc_bar_set => Δc_bar_set,
            :N_spins_sys_set => N_spins_sys_set,
    
            :fs => fs,
            :SW => SW,
            :ν_0ppm => ν_0ppm,

            # types.
            :αs_set_type => typeof(αs_set),
            :Ωs_set_type => typeof(Ωs_set),

            :Δc_set_type => typeof(Δc_set),
            :parts_set_type => typeof(parts_set),
            :Δc_bar_set_type => typeof(Δc_bar_set),
            :N_spins_sys_set_type => typeof(N_spins_sys_set),
    
            :fs_type => typeof(fs),
            :SW_type => typeof(SW),
            :ν_0ppm_type => typeof(ν_0ppm),
        )
    end

    return nothing
end


"""
```
deserializemixture(W)
```
W is of a data type that can be addressed via a key (of type `Symbol`) and returns a value. Example data types that W can be are `Dict{Symbol, Any}` or a `JSON3` Object.

Outputs a `Vector{SHType{T}}` variable, where `T` is the data type of the key `:fs` in `W`.
"""
function deserializemixture(W)

    αs_set = convert(recovertype(W[:αs_set_type]), W[:αs_set])
    Ωs_set = convert(recovertype(W[:Ωs_set_type]), W[:Ωs_set])

    Δc_set = convert(recovertype(W[:Δc_set_type]), W[:Δc_set])
    parts_set = convert(recovertype(W[:parts_set_type]), W[:parts_set])

    Δc_bar_set = convert(recovertype(W[:Δc_bar_set_type]), W[:Δc_bar_set])

    N_spins_sys_set = convert(recovertype(W[:N_spins_sys_set_type]), W[:N_spins_sys_set])
    
    fs = convert(recovertype(W[:fs_type]), W[:fs])
    SW = convert(recovertype(W[:SW_type]), W[:SW])
    ν_0ppm =  convert(recovertype(W[:ν_0ppm_type]), W[:ν_0ppm])

    # assemble output.
    T = eltype(fs)
    N = length(αs_set)

    @assert T <: AbstractFloat
    # skip error-checking on the _set fields, for convinence of mainability. TODO: later.

    As = Vector{SHType{T}}(undef, N)

    for n in eachindex(As)

        As[n] = SHType(
            αs_set[n],
            Ωs_set[n],

            Δc_set[n],
            parts_set[n],
            Δc_bar_set[n],
            N_spins_sys_set[n],

            fs,
            SW,
            ν_0ppm,
        )
    end

    return As
end

### physical parameters and corresponding entry names.

"""
```
serializephysicalparams(
    Phys::Vector{PhysicalParamsType{T}},
    molecule_entries::Vector{String},
    ) where T
```
Returns a `Dict{Symbol, Any}` object, which can then be used with file serialization packages for saving to disk.
"""
function serializephysicalparams(
    Phys::Vector{PhysicalParamsType{T}},
    molecule_entries::Vector{String},
    ) where T

    if isempty(Phys)
        println("input is an empty collection. Cannot serialize.")
    else

        N = length(Phys)

        H_IDs_set = Vector{typeof(Phys[begin].H_IDs)}(undef, N)

        H_inds_sys_set = Vector{typeof(Phys[begin].H_inds_sys)}(undef, N)
        cs_sys_set = Vector{typeof(Phys[begin].cs_sys)}(undef, N)

        H_inds_singlets_set = Vector{typeof(Phys[begin].H_inds_singlets)}(undef, N)
        cs_singlets_set = Vector{typeof(Phys[begin].cs_singlets)}(undef, N)

        J_inds_sys_set = Vector{typeof(Phys[begin].J_inds_sys)}(undef, N)
        J_inds_sys_local_set = Vector{typeof(Phys[begin].J_inds_sys_local)}(undef, N)
        J_vals_sys_set = Vector{typeof(Phys[begin].J_vals_sys)}(undef, N)

        ME_set = Vector{typeof(Phys[begin].ME)}(undef, N)

        for n in eachindex(Phys)
            H_IDs_set[n] = Phys[n].H_IDs

            H_inds_sys_set[n] = Phys[n].H_inds_sys
            cs_sys_set[n] = Phys[n].cs_sys
    
            H_inds_singlets_set[n] = Phys[n].H_inds_singlets
            cs_singlets_set[n] = Phys[n].cs_singlets
    
            J_inds_sys_set[n] = Phys[n].J_inds_sys
            J_inds_sys_local_set[n] = Phys[n].J_inds_sys_local
            J_vals_sys_set[n] = Phys[n].J_vals_sys
    
            ME_set[n] = Phys[n].ME
        end
        
        return Dict(
            :molecule_entries => molecule_entries,    
        
            # values.    
            :H_IDs_set => H_IDs_set,
            :H_inds_sys_set => H_inds_sys_set,
            :cs_sys_set => cs_sys_set,

            :H_inds_singlets_set => H_inds_singlets_set,
            :cs_singlets_set => cs_singlets_set,

            :J_inds_sys_set => J_inds_sys_set,
            :J_inds_sys_local_set => J_inds_sys_local_set,
            :J_vals_sys_set => J_vals_sys_set,

            :ME_set => ME_set,

            # types.
            :H_IDs_set_type => typeof(H_IDs_set),
            :H_inds_sys_set_type => typeof(H_inds_sys_set),
            :cs_sys_set_type => typeof(cs_sys_set),

            :H_inds_singlets_set_type => typeof(H_inds_singlets_set),
            :cs_singlets_set_type => typeof(cs_singlets_set),

            :J_inds_sys_set_type => typeof(J_inds_sys_set),
            :J_inds_sys_local_set_type => typeof(J_inds_sys_local_set),
            :J_vals_sys_set_type => typeof(J_vals_sys_set),

            :ME_set_type => typeof(ME_set),

            :T => T,
        )
    end

    return nothing
end

"""
```
deserializephysicalparams(W)
```
W is of a data type that can be addressed via a key (of type `Symbol`) and returns a value. Example data types that W can be are `Dict{Symbol, Any}` or a `JSON3` Object.

Outputs a `Vector{PhysicalParamsType{T}}` and a `Vector{String}` variable. `T` is the value stored in the `:T` key in `W`.
"""
function deserializephysicalparams(W)

    H_IDs_set = convert(recovertype(W[:H_IDs_set_type]), W[:H_IDs_set])

    H_inds_sys_set = convert(recovertype(W[:H_inds_sys_set_type]), W[:H_inds_sys_set])
    cs_sys_set = convert(recovertype(W[:cs_sys_set_type]), W[:cs_sys_set])

    H_inds_singlets_set = convert(recovertype(W[:H_inds_singlets_set_type]), W[:H_inds_singlets_set])
    cs_singlets_set = convert(recovertype(W[:cs_singlets_set_type]), W[:cs_singlets_set])

    #J_inds_sys_set = convert(recovertype(W[:J_inds_sys_set_type]), W[:J_inds_sys_set])
    A = W[:J_inds_sys_set]
    J_inds_sys_set = collect(
        collect(
            collect(
                tuple(convert(Int, A[n][i][j][begin]), convert(Int, A[n][i][j][end]))
                for j in eachindex(A[n][i])
            )
            for i in eachindex(A[n])
        )
        for n in eachindex(A)
    )
    J_inds_sys_set = convert(Vector{Vector{Vector{Tuple{Int,Int}}}}, J_inds_sys_set)

    A = W[:J_inds_sys_local_set]
    J_inds_sys_local_set = collect(
        collect(
            collect(
                tuple(convert(Int, A[n][i][j][begin]), convert(Int, A[n][i][j][end]))
                for j in eachindex(A[n][i])
            )
            for i in eachindex(A[n])
        )
        for n in eachindex(A)
    )
    J_inds_sys_local_set = convert(Vector{Vector{Vector{Tuple{Int,Int}}}}, J_inds_sys_local_set)

    #J_inds_sys_local_set = convert(recovertype(W[:J_inds_sys_local_set_type]), W[:J_inds_sys_local_set])
    J_vals_sys_set = convert(recovertype(W[:J_vals_sys_set_type]), W[:J_vals_sys_set])
    ME_set = convert(recovertype(W[:ME_set_type]), W[:ME_set])

    molecule_entries = convert(Vector{String}, W[:molecule_entries])

    # assemble output.
    T = recovertype(W[:T])
    N = length(molecule_entries)

    @assert T <: AbstractFloat
    # skip error-checking on the _set fields, for convinence of mainability. TODO: later.

    Phys = Vector{PhysicalParamsType{T}}(undef, N)

    for n in eachindex(Phys)

        Phys[n] = PhysicalParamsType(
            H_IDs_set[n],

            H_inds_sys_set[n],
            cs_sys_set[n],
        
            H_inds_singlets_set[n],
            cs_singlets_set[n],
        
            J_inds_sys_set[n],
            J_inds_sys_local_set[n],
            J_vals_sys_set[n],
            ME_set[n],
        )
    end

    return Phys, molecule_entries
end