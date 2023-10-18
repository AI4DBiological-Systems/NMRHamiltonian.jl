

# not type stable.
# args[k] much be a Vector of length(labels).
function assembledict(molecule_entries::Vector{NT},
    labels::Vector{LT},
    args...) where {NT,LT}

    N_labels = length(labels)
    @assert length(args) == N_labels # future: add error message.

    N_entries = length(molecule_entries)
    for j in eachindex(args)
        @assert length(args[j]) == N_entries # future: add error message.
    end

    dict_out = Dict()

    for n in eachindex(molecule_entries)
        dict_out[molecule_entries[n]] = Dict()

        for j in eachindex(labels)
            merge!(dict_out[molecule_entries[n]], Dict(labels[j] => args[j][n]))
        end
    end

    return dict_out
end

# makes a folder at `file_folder` if it doesn't exist.
function saveasJSON(
    output_folder::String,
    file_name::String,
    dict,
    )

    # set up save path.
    isdir(output_folder) || mkpath(output_folder)
    save_path = joinpath(output_folder, file_name)

    return saveasJSON(save_path, dict)
end

function saveasJSON(save_path::String, dict)

    #stringdata = JSON.json(dict)
    stringdata = dict

    open(save_path, "w") do f
        JSON3.pretty(f, stringdata)
        println(f)
    end

    return nothing
end

## not type safe.
# generic read with nested vectors.
# assumes:
#   - each parameter is a 1D array.
#   - every name entry has the same N_parameters_per_entry number categories of parameters.
# x[j,n] is the parameter for the j-th parameter category for the n-th name entry.
function convertnestedvectors(c::JSON3.Array, val_type::DataType)

    return collect( nestedconvert(c, i, val_type) for i in eachindex(c) )
end

# this is not type stable..
# This is for nested Vector{Any} but really have base Julia concrete types at the inner-most level.
# position takes value in {1, 2, ..., number of elements of x}.
#   - e.g., position is an element of the eachindex(x) iterator.
function nestedconvert(x::Vector{Any}, position, val_type::DataType)

    if isempty(x)
        ret = Vector{val_type}(undef, 0)
        return convert(Vector{val_type}, ret)
    end

    y = x[begin+position-1]
    current_type = typeof(y)

    while current_type <: JSON3.Array
        
        return convertnestedvectors(y, val_type) # returns a Vector{} or Vctor{Vector{T}} or higher-order nested Vector{T}'s.
    end

    return convert(val_type, y)
end

function nestedconvert(x::JSON3.Array, position, val_type::DataType)

    if isempty(x)
        ret = Vector{val_type}(undef, 0)
        return convert(Vector{val_type}, ret)
    end

    y = x[begin+position-1]
    current_type = typeof(y)

    while current_type <: JSON3.Array

        return convertnestedvectors(y, val_type) # returns a Vector{} or Vctor{Vector{T}} or higher-order nested Vector{T}'s.
    end

    return convert(val_type, y)
end

##############



# mutates `vars_vector_struct`
# `mixture_names` pair with `vars_vector_struct`.
function updatefromfile!(vars_vector_struct::Vector{ST},
    vars_target::Vector{Symbol},
    mapping_struct_to_json::Dict{Symbol,Symbol},
    J_mixture,
    target_entries::Vector{String},
    mixture_names::Vector{String}) where ST

    target_entries_sym = Symbol.(target_entries)
    
    # check.
    #J_mixture_keys = collect( i for i in eachindex(J_mixture))
    for key in target_entries_sym
        #@assert key in J_mixture_keys
        @assert key in eachindex(J_mixture)
    end

    # process.
    for n in eachindex(target_entries_sym)
        
        key = target_entries_sym[n]
        J = J_mixture[key]
        # @show key
        # @show J
        
        ind = findfirst(xx->xx==target_entries[n], mixture_names)

        if typeof(ind) != Nothing
            updatefromfile!(vars_vector_struct[ind], vars_target,
                mapping_struct_to_json, J)
        end
    end

    return nothing
end

# does not update json attribute entries that are empty.
# mutates `vars_struct`
# `mapping_json_to_struct` options are getshconfigtolabelsmapping(), 
# ST must be a mutable struct.
# the order of `all_keys` should match the order of the Vector{}-typed members in `vars_struct`. 
# J must be JSON3.Object or Dict{Symbol, .}, . means don't care.
function updatefromfile!(vars_struct::ST,
    vars_target::Vector{Symbol},
    mapping_struct_to_json::Dict{Symbol,Symbol},
    J) where ST

    ## setup.
    vars_all = fieldnames(ST)

    # get the corresponding keys in J.
    keys_target = collect( mapping_struct_to_json[v] for v in vars_target)
    
    ## for debugging.
    # @show keys_target
    # @show J_keys = collect( i for i in eachindex(J))
    ##

    ## key existence checks.
    
    # check if target fields are in the data type ST.
    for field in vars_target
        @assert field in vars_all
    end

    # check if the keys are in the dict/json3 object J
    for key in keys_target
        @assert key in eachindex(J)
    end

    ## update vars_struct.
    for i in eachindex(keys_target)
        key = keys_target[i]
        field = vars_target[i]
        
        #FT = fieldtype(ST, field)
        nested_T = nestedeltype(getfield(vars_struct, field))
        #@show getfield(vars_struct, field)
        c = J[key]
        #@show length(c)
        
        if length(c) > 0
            # conversion.
            update_value = convertnestedvectors(c, nested_T)
            setfield!(vars_struct, field, update_value)
        end
    end

    return nothing
end

function nestedeltype(x::Vector)
    y = eltype(x)
    while y <: Vector
        y = eltype(y)
    end
    return(y)
end


##### for creating cs-J-oupling JSON files.

function Jcouplingtodict(J_IDs::Vector{Tuple{Int,Int}}, J_vals::Vector{T}) where T

    return collect( Dict("ID1" => J_IDs[i][1], "ID2" => J_IDs[i][2], "value" => J_vals[i]) for i = 1:length(J_vals) )
end

function cstodict(H_IDs::Vector{Int}, H_css::Vector{T}) where T

    return collect( Dict("ID" => H_IDs[i], "value" => H_css[i]) for i = 1:length(H_css) )
end

function savecouplinginfo(
    save_path::String,
    H_IDs::Vector{Int},
    H_css::Vector{T},
    J_IDs::Vector{Tuple{Int,Int}},
    J_vals::Vector{T},
    ) where T

    @assert length(H_IDs) == length(H_css)
    @assert length(J_IDs) == length(J_vals)

    cs_dict = cstodict(H_IDs, H_css)
    J_dict = Jcouplingtodict(J_IDs, J_vals)
    dict1 = Dict("chemical shift" => cs_dict, "J-coupling" => J_dict)
    stringdata = dict1

    open(save_path, "w") do f
        JSON3.pretty(f, stringdata)
        println(f)
    end

    return nothing
end

function savecouplinginfo(
    Phys::Vector{PhysicalParamsType{T}},
    file_names::Vector{String};
    save_folder::String = pwd(),
    ) where T

    @assert length(Phys) == length(file_names)

    for n in eachindex(Phys)
        H_IDs, H_css, J_IDs, J_vals = extractcouplinginfo(Phys[n])

        savecouplinginfo(
            joinpath(save_folder, file_names[n]),
            H_IDs, H_css, J_IDs, J_vals,
        )
    end

    return nothing
end

function extractcouplinginfo(Phy::PhysicalParamsType{T}) where T

    ID_LUT = Phy.H_IDs

    H_IDs = Vector{Int}(undef, 0)
    H_css = Vector{T}(undef, 0)

    # non-singlets.
    for i in eachindex(Phy.cs_sys)
        push!(H_css, Phy.cs_sys[i]...)
        push!(H_IDs, ID_LUT[Phy.H_inds_sys[i]]...)
    end

    # singlets.
    for i in eachindex(Phy.cs_singlets)
        cs_i = Phy.cs_singlets[i]
        H_IDs_i = Phy.H_inds_singlets[i]

        cs_update = collect( cs_i for _ in eachindex(H_IDs_i))
        H_ID_update = collect( ID_LUT[H_IDs_i[j]] for j in eachindex(H_IDs_i) )

        push!(H_css, cs_update...)
        push!(H_IDs, H_ID_update...)
    end

    # J-coupling.
    J_IDs = Vector{Tuple{Int,Int}}(undef, 0)
    J_vals = Vector{T}(undef, 0)
    for i in eachindex(Phy.cs_sys)
        push!(J_vals, Phy.J_vals_sys[i]...)
        
        tmp = collect(
            (ID_LUT[Phy.J_inds_sys[i][k][1]], ID_LUT[Phy.J_inds_sys[i][k][2]])
            for k in eachindex(Phy.J_inds_sys[i])
        )
        push!(J_IDs, tmp...)
    end

    return H_IDs, H_css, J_IDs, J_vals
end

function createnamemappingJSON(
    save_path::String,
    label_names::Vector{String},
    file_names::Vector{String};
    notes::Vector{String} = collect( "" for _ in eachindex(label_names)),
    )

    @assert length(label_names) == length(file_names) == length(notes)

    dic = Dict()
    for n in eachindex(label_names)

        key, file_name, note = label_names[n], file_names[n], notes[n]

        dic[key] = Dict( "file name" => file_name, "notes" => note)
    end

    # save.
    stringdata = dic #JSON.json(dic)
    open(save_path, "w") do f
        JSON3.pretty(f, stringdata)
        println(f)
    end

    return dic
end