
######## initialialization routines.


"""
```
extractcs(Phys::Vector{PhysicalParamsType{T}}) where T
```

Assemble the J-coupling and chemical shift values from JSON files for the specified molecules.

Construction functions for molecule mixture-related configurations in NMRHamiltonian.jl depend on a provided list of the number of spin systems and singlets for each molecule. The constructor functions can be inferred this information from the output quantities of this function.

### Inputs

- `Phys`$(DOCSTRING_Phys("T"))

### Outputs

- `cs_sys_mixture:`$(DOCSTRING_cs_sys_mixture("T"))

- `cs_singlets_mixture`$(DOCSTRING_cs_singlets_mixture("T"))

"""
function extractcs(Phys::Vector{PhysicalParamsType{T}}) where T
    #
    cs_sys_mixture = collect( Phys[n].cs_sys for n in eachindex(Phys) )
    cs_singlets_mixture = collect( Phys[n].cs_singlets for n in eachindex(Phys) )

    return cs_sys_mixture, cs_singlets_mixture
end

### generate generic spin system of a mixture datastructure.

# works for sys_mixture = cs_sys_mixture for non-singlet spin systems.
# use sys_mixture = cs_singlets_mixture for singlet spin systems.
# generatorfunc is either nothing or of the form sys_mixture[n][i]::IT ↦ value::OT.
#   - set to nothing to use the provided initial_value.
function generatemixturespinsystemvector(sys_mixture::Vector{Vector{IT}},
    initial_value::OT;
    generatorfunc = nothing)::Vector{Vector{OT}} where {IT,OT}

    #
    N = length(sys_mixture)

    out = Vector{Vector{OT}}(undef, N)
    for n in eachindex(sys_mixture)
        out[n] = Vector{OT}(undef, length(sys_mixture[n]))

        for i in eachindex(out[n])
            if typeof(generatorfunc) == Nothing
                out[n][i] = deepcopy(initial_value)
            else
                out[n][i] = generatorfunc(sys_mixture[n][i])
            end
        end
    end

    return out
end


"""
```
defaultmixtureshsconfig(
    cs_sys_mixture::Vector{Vector{Vector{T}}};
    coherence_tol::T = convert(T, 0.01),
    relative_α_threshold::T = convert(T, 0.01),
    normalize_α::Bool = true
)::Vector{SHsConfigType{T}} where T
```

Assembles a default configuration data structure for the `simulate()` function for the list of molecule entries that corresponds to `cs_sys_mixture`.

### Inputs

- `cs_sys_mixture`$(DOCSTRING_cs_sys_mixture("T"))

### Optional inputs

- `coherence_tol` 
Value for each element of `coherence_tols` in `SHsConfigType`, default set to `0.01`.


- `relative_α_threshold`
Value for each element of `relative_α_thresholds` in `SHsConfigType`, default set to `0.01`.

- `normalize_α`
Value for each element of `normalize_αs` in `SHsConfigType`, default set to `true`.


### Outputs

- `out::Vector{SHsConfigType{T}}`
list of configuration settings for `simulate()`.

"""
function defaultmixtureshsconfig(
    cs_sys_mixture::Vector{Vector{Vector{T}}};
    coherence_tol::T = convert(T, 0.01),
    relative_α_threshold::T = convert(T, 0.01),
    normalize_α::Bool = true
    )::Vector{SHsConfigType{T}} where T
    
    @assert zero(T) <= relative_α_threshold <= one(T)

    x = generatemixturespinsystemvector(cs_sys_mixture, zero(T))

    out = collect( SHsConfigType(
        ones(T, length(cs_sys_mixture[n])) .* coherence_tol,
        ones(T, length(cs_sys_mixture[n])) .* relative_α_threshold,
        ones(Bool, length(cs_sys_mixture[n])) ) for n in eachindex(x) )
    
    # to ensure it remains as a Vector{Bool} instead of BitVector.
    for n in eachindex(out)
        fill!(out[n].normalize_αs, normalize_α)
    end

    return out
end

function getwriteshconfiglabels()
    return ["relative α tolerance"; "coherence tolerance"; "normalize α"] # usually singular. for one spin system.
end

# left hand side is from the field names of `SHsConfigType`.
# right hand side is from `getwriteshconfiglabels()`.
# need to manually change this Dict if SHsConfigType and getwriteshconfiglabels where changed.
function getshconfigtolabelsmapping()

    return Dict(
        Symbol("relative_α_thresholds") => Symbol("relative α tolerance"),
        Symbol("coherence_tols") => Symbol("coherence tolerance"),
        Symbol("normalize_αs") => Symbol("normalize α"))
end

##
"""
```
defaultmixturepartitionsparameters(
    cs_sys_mixture::Vector{Vector{Vector{T}}},
    default_θ::KT;
    default_γ::T = one(T),
    default_knns::Int = 0,
    default_radi::T = one(T),
    )::Vector{PartitionsParametersType{T,KT}} where {T,KT}
```
Creates a list of partition parameters configurations for the list of molecule entries corresponding to `cs_sys_mixture`. For use with the `simulate()` function. 

### Inputs
- `cs_sys_mixture`$(DOCSTRING_cs_sys_mixture("T"))

### Optional inputs
- `default_θ::KT`, `default_γ::T`, `default_knns::Int`, `default_radi::T` are the values supplied to the `PartitionsParametersType{T,KT}` constructor.

### Outputs

Single output of type `::Vector{PartitionsParametersType{T,KT}}`, the list of partition parameters configuration settings for `simulate()`.

"""
function defaultmixturepartitionsparameters(
    cs_sys_mixture::Vector{Vector{Vector{T}}},
    default_θ::KT;
    default_γ::T = one(T),
    default_knns::Int = 0,
    default_radi::T = one(T),
    )::Vector{PartitionsParametersType{T,KT}} where {T,KT}
    
    x = generatemixturespinsystemvector(cs_sys_mixture, zero(T))

    return collect( PartitionsParametersType(
        collect( deepcopy(default_θ) for _ in eachindex(cs_sys_mixture[n]) ),
        collect( default_γ for _ in eachindex(cs_sys_mixture[n]) ),
        collect( default_knns for _ in eachindex(cs_sys_mixture[n]) ),
        collect( default_radi for _ in eachindex(cs_sys_mixture[n]) ),
        ) for n in eachindex(x) )
end


function getwritepartitionconfiglabels()
    return ["θ"; "γ"; "knn"; "radius"] # usually singular. for one spin system.
end

# left hand side is from the field names of `PartitionAlgorithmType`.
# right hand side is from `getwritepartitionconfiglabels().
# LHS is usually the plural of the RHS, since LHS is multiple spin system, RHS is single spin system.
function getpartitionconfigtolabelsmapping()

    return Dict(
        Symbol("θs") => Symbol("θ"),
        Symbol("γs") => Symbol("γ"),
        Symbol("knns") => Symbol("knn"),
        Symbol("radi") => Symbol("radius"))
end

##### load routines


"""
```
loadmixtureshsconfig(
    cs_sys_mixture::Vector{Vector{Vector{T}}},
    file_path::String,
    molecule_entries::Vector{String},
    ::Type{T};
    target_entries::Vector{String} = molecule_entries,
    vars_target::Vector{Symbol} = collect(fieldnames(SHsConfigType)),
    supplied_coherence_tol::T = convert(T, 0.01),
    supplied_relative_α_threshold::T = convert(T, 0.01),
    supplied_normalize_α::Bool = true,
    )::Vector{SHsConfigType{T}} where T <: AbstractFloat
```

Loads a configuration data structure from file for the `simulate()` function. 
    
Only the targeted names and field names are updated from file. The non-targeted names and non-targeted field names are filled using the supplied values to `loadmixtureshconfig`. See `SHsConfigType` for the field names.

### Inputs

- `cs_sys_mixture`$(DOCSTRING_cs_sys_mixture("T"))

- `file_path`
Path to the configuration file.

- `molecule_entries`
list of molecule entries corresponding to `cs_sys_mixture`.

- `::Type{T}`
The floating-point data type to use. Use `Float64` if you are running a 64-bit version of Julia.


### Optional inputs

- `target_entries`
the targeted molecule entries to update from file. This should be a subset of `molecule_entries`. Default behavior is to update all molecule entries.

- `vars_target`
The targeted field names of `SHsConfigType` to update from file. The non-targeted field names will take on the values supplied. Default behavior is to update all field names.

- `supplied_coherence_tol`
If not loading from file, this value is used for each element of `relative_α_thresholds` in `SHsConfigType`, default set to `0.01`.

- `supplied_relative_α_threshold`
If not loading from file, this value is used for each element of `relative_α_thresholds` in `SHsConfigType`, default set to `0.01`.

- `supplied_normalize_α`
If not loading from file, this value is used for each element of `normalize_αs` in `SHsConfigType`, default set to `true`.

### Outputs

(WIP)
- `mixture_sh_config::Vector{SHsConfigType{T}}`
List of configuration settings for `simulate()`.

"""
function loadmixtureshsconfig(
    cs_sys_mixture::Vector{Vector{Vector{T}}},
    file_path::String,
    molecule_entries::Vector{String},
    ::Type{T};
    target_entries::Vector{String} = molecule_entries,
    vars_target::Vector{Symbol} = collect(fieldnames(SHsConfigType)), # default to update all fields.
    supplied_coherence_tol::T = convert(T, 0.01),
    supplied_relative_α_threshold::T = convert(T, 0.01),
    supplied_normalize_α::Bool = true
    )::Vector{SHsConfigType{T}} where T <: AbstractFloat

    # check.
    @assert zero(T) <= supplied_relative_α_threshold <= one(T)
    @assert isfile(file_path)
    for name in target_entries
        @assert name in molecule_entries
    end

    # create output container with default values.
    mixture_sh_config = defaultmixtureshsconfig(
        cs_sys_mixture;
        coherence_tol = supplied_coherence_tol,
        relative_α_threshold = supplied_relative_α_threshold,
        normalize_α = supplied_normalize_α
    )

    # get JSON3 object of file.
    json_string = read(file_path)
    S = JSON3.read(json_string)

    # populate container with contents of `S`.
    updatefromfile!(
        mixture_sh_config,
        vars_target, 
        getshconfigtolabelsmapping(),
        S,
        target_entries,
        molecule_entries
    )

    return mixture_sh_config
end

"""
```
loadmixturepartitionsparameters(
    cs_sys_mixture::Vector{Vector{Vector{T}}},
    file_folder::String,
    file_name::String,
    molecule_entries::Vector{String},
    default_θ::KT;
    target_entries::Vector{String} = molecule_entries,
    vars_target::Vector{Symbol} = collect(fieldnames(PartitionsParametersType)),
    default_γ::T = one(T),
    default_knns::Int = 60,
    default_radi::T = convert(T, 0.3),
    ) where {T <: AbstractFloat, KT}
```

Assembles a partition parameters configuration data structure using information from a file. The configuration is for the `simulate()` function. 
    
Only the targeted names and field names are updated from file. The non-targeted names and non-targeted field names are filled using the supplied values to `loadmixtureshconfig`. See `PartitionsParametersType` for the field names.

### Inputs

- `cs_sys_mixture`$(DOCSTRING_cs_sys_mixture("T"))

- `file_path::String`
Path to the configuration file.

- `molecule_entries::Vector{String}`
list of molecule entries corresponding to `cs_sys_mixture`.

- `default_θ::KT`
If not loading from file, each element of `θs` in `PartitionsParametersType{T,KT}` is assigned a copy of this input value is made via `deepcopy()`.

### Optional inputs

- `target_entries::Vector{String}`
the targeted molecule entries to update from file. This should be a subset of `molecule_entries`. Default behavior is to update all molecule entries.

- `vars_target::Vector{Symbol}`
The targeted field names of `SHsConfigType` to update from file. The non-targeted field names will take on the values supplied. Default behavior is to update all field names.

- `default_γ::T`
If not loading from file, this value is used for each element of `relative_α_thresholds` in `PartitionsParametersType`. Default set to `1`.

- `default_knns::Int`
If not loading from file, this value is used for each element of `knns` in `PartitionsParametersType`. Default set to `60`.

- `default_radi::T`
If not loading from file, this value is used for each element of `normalize_αs` in `PartitionsParametersType`. Default set to `0.3`.

### Outputs

The only output is of type `::Vector{SHsConfigType{T}}`, the list of configuration settings for `simulate()`.

"""
function loadmixturepartitionsparameters(
    cs_sys_mixture::Vector{Vector{Vector{T}}},
    file_path::String,
    molecule_entries::Vector{String},
    default_θ::KT;
    target_entries::Vector{String} = molecule_entries,
    vars_target::Vector{Symbol} = collect(fieldnames(PartitionsParametersType)),
    default_γ::T = one(T),
    default_knns::Int = 60,
    default_radi::T = convert(T, 0.3),
    ) where {T <: AbstractFloat, KT}

    # check.
    for name in target_entries
        @assert name in molecule_entries
    end

    # create output container with default values.
    mixture_parts_params = defaultmixturepartitionsparameters(
        cs_sys_mixture,
        default_θ;
        default_γ = default_γ,
        default_knns = default_knns,
        default_radi = default_radi
    )

    # get JSON3 object of file.
    json_string = read(file_path)
    S = JSON3.read(json_string)

    # populate container with contents of `S`.
    updatefromfile!(
        mixture_parts_params,
        vars_target, 
        getpartitionconfigtolabelsmapping(),
        S,
        target_entries,
        molecule_entries
    )

    return mixture_parts_params
end

##### save routines

function savemixtureshsconfig(file_folder::String,
    file_name::String,
    molecule_entries::Vector{String},
    α_tols::Vector{Vector{T}},
    coherence_tols::Vector{Vector{T}},
    normalize_αs::Vector{Vector{Bool}}) where T <: AbstractFloat

    sh_labels = getwriteshconfiglabels()
    dict = assembledict(molecule_entries, sh_labels, α_tols, coherence_tols, normalize_αs)
    saveasJSON(file_folder, file_name, dict)

    return nothing
end

function savemixturepartitionsparameters(file_folder::String,
    file_name::String,
    molecule_entries::Vector{String},
    θs::Vector{Vector{KT}},
    γs::Vector{Vector{T}},
    knns::Vector{Vector{Int}},
    radi::Vector{Vector{T}}) where {T<:AbstractFloat,KT}
    
    cc_labels = getwritepartitionconfiglabels()
    dict = assembledict(molecule_entries, cc_labels, θs, γs, knns, radi)
    saveasJSON(file_folder, file_name, dict)

    return nothing
end


#### load.
"""
```
getphysicalparameters(
    ::Type{T},
    target_entries::Vector{String},
    H_params_path::String,
    molecule_mapping_file_path;
    unique_cs_atol = 1e-6,
) where T <: AbstractFloat
```

Assemble the J-coupling and chemical shift values from JSON files for the specified molecules.

### Nomenclature

One set of J-coupling and one set of chemical shift values for a molecule is what NMRHamiltonian.jl calls a molecule entry.

It is common to see different J-coupling values of the same molecule reported in literature for similar experimental conditions.  NMRHamiltonian.jl requires one to record each instance as separate molecule entry.

Every molecule entry is recorded as a separate JSON file with some arbitrary file name the data collector wants to use. An additional name-mapping JSON file is needed to translate what the NMRHamiltonian.jl user wants to label the entries to what is the data collector chose to name the JSON files. The name-mapping JSON file is likely to require manual set up.

### Inputs

- `::Type{T}`       -- The `AbstractFloat` datatype to use for storing floating point data. A viable input here is `Float64`, which specifies double precision floating point numbers are to be used.
- `target_entries`  -- list of molecule entry names.
- `H_params_path`   -- path to the directory that contain the J-coupling and chemical shift JSON files. Each file corresponds to one entry of a molecule, and has a set of chemical shifts and J-coupling information in the JSON dictionary format.
    `ID1` and `ID2` are JSON dictionary keys that specify the spin nucleus label for J-coupling values, and `ID` is the nucleus label for chemical shift values. The following is an example of the JSON format for a L-Histidine entry:
```
{
       "J-coupling": [
                       {
                            "ID2": 13,
                            "ID1": 12,
                          "value": -15.350518
                       },
                       {
                            "ID2": 16,
                            "ID1": 12,
                          "value": 7.63377
                       },
                       {
                            "ID2": 16,
                            "ID1": 13,
                          "value": 5.029267
                       }
                     ],
   "chemical shift": [
                       {
                          "value": 3.18747,
                             "ID": 12
                       },
                       {
                          "value": 3.26837,
                             "ID": 13
                       },
                       {
                          "value": 7.14078,
                             "ID": 14
                       },
                       {
                          "value": 8.02487,
                             "ID": 15
                       },
                       {
                          "value": 3.9967,
                             "ID": 16
                       }
                     ]
}
```

- `molecule_mapping_file_path`   -- the files in `H_params_path` might not be named in a manner for the user to know which molecule entry corresponds to which JSON file. Therefore, the user should further define a single JSON file that encodes this mapping.
    Example: if the file `bmse000976_simulation_1.json` corresponds to the molecule entry "L-Histidine" and the file `Histidine_Govindaraju_2000.json` corresponds to the molecule entry "L-Histidine - Govindaraju", and suppose the user only wants to target these two entries, i.e. `target_entries` contain only these two entries. Then the following entry-filename-mapping JSON file should be used, and `molecule_mapping_file_path` should be the path to this entry-filename-mapping file.
```
{
    "L-Histidine": {
        "notes": "http://gissmo.bmrb.io/entry/bmse000976/simulation_1",
       "file name": "bmse000976_simulation_1.json"
    },
    "L-Histidine - Govindaraju": {
        "notes": "https://pubmed.ncbi.nlm.nih.gov/10861994/",
       "file name": "Histidine_Govindaraju_2000.json"
    }
}
```

### Optional inputs

- `unique_cs_atol`   -- Two chemical shift values in the same JSON file are assigned the same chemical shift value to both nuclei if the two values in the file are less than `unique_cs_atol`.

### Outputs

- `Phys::Vector{PhysicalParamsType{Float64}`    -- list of the data type `PhysicalParamsType` that contain chemical shift and J-coupling information for each of the molecule entries in `target_entries`.

- `dict_molecule_to_filename`                   -- the dictionary that maps molecule entries to meta information such as their corresponding JSON filenames

"""
function getphysicalparameters(
    ::Type{T},
    target_entries::Vector{String},
    H_params_path::String,
    molecule_mapping_file_path;
    unique_cs_atol::T = convert(T, 1e-6),
    ) where T
    
    load_paths, dict_molecule_to_filename = getloadpaths(target_entries, H_params_path, molecule_mapping_file_path)

    return getphysicalparameters(
        T,
        load_paths;
        unique_cs_atol = unique_cs_atol,
        )
end

function getphysicalparameters(
    ::Type{T},
    load_paths::Vector{String};
    unique_cs_atol::T = convert(T, 1e-6),
    ) where T <: AbstractFloat

    Phys = Vector{PhysicalParamsType{T}}(undef, length(load_paths))

    for n in eachindex(load_paths)
        H_IDs, H_css, J_IDs, J_vals = loadcouplinginfojson(T, load_paths[n])

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

function getloadpaths(
    target_entries::Vector{String},
    H_params_path::String,
    molecule_mapping_file_path::String
    )

    dict_molecule_to_filename = JSON3.read(read(molecule_mapping_file_path)) # map molecule entries to coupling information file names.

    N_molecules = length(target_entries)
    load_paths = Vector{String}(undef, N_molecules)

    file_name_symbol = Symbol("file name")

    for n in eachindex(load_paths)

        ## if all keys to dictionary were type `Symbol`.
        name_key = Symbol(target_entries[n])
        load_paths[n] = joinpath(H_params_path, dict_molecule_to_filename[name_key][file_name_symbol])
    end

    return load_paths, dict_molecule_to_filename
end
