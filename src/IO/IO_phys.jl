
############# IO ζ to cs-sys to Phy::PhysicalParamsType.

# modified from condensenuclei(), to not accumulate the equivalence variables.
"""
ordering must contain integers from 1:N, where N is the number of unique entries in ordering.

Example:
using Linearalgebra
cs = randn(10)
ordering = [ 2; 2; 2; 3; 4; 5;  1; 1; 1; 6]
z = nuclei2vars(cs, ordering)

"""
function nuclei2vars(
    x::Vector{T},
    ordering::Vector{Int},
    )::Vector{T} where T

    @assert length(x) == length(ordering)

    N = length(unique(ordering))
    @assert norm(collect(1:N) - sort(unique(ordering))) < 1e-14

    y = zeros(T, N)
    for i in eachindex(x)

        k = ordering[i]
        y[k] = x[i] # replace when indicated by ordering.
    end

    return y
end

function vars2nuclei(
    y::Vector{T},
    ordering::Vector{Int},
    )::Vector{T} where T

    @assert length(y) == maximum(ordering)
    for i in eachindex(y)
        @assert 1 <= i
    end

    N = length(unique(ordering))
    @assert norm(collect(1:N) - sort(unique(ordering))) < 1e-14

    x = zeros(T, length(ordering))
    for i in eachindex(x)

        k = ordering[i]
        x[i] = y[k] # replace when indicated by ordering.
    end

    return x
end
# 

### reading and writing chemical shift values that correspond to Δc_bar.
# account for any magnetic equivalence relationships here.
# creates copy.
function readbasechemshifts(
    cs_sys::Vector{Vector{T}},
    cs_singlets::Vector{T},
    ME::Vector{Vector{Vector{Int}}},
    )::Vector{Vector{T}} where T
    
    N_sys = length(ME)
    @assert N_sys == length(cs_sys)

    N = N_sys + length(cs_singlets)
    cs_shifts = Vector{Vector{T}}(undef, N)
    
    N_spins_sys = collect( length(cs_sys[i]) for i in eachindex(cs_sys))

    for i in eachindex(ME)
        ordering, DOF = createorderingfromeqinds(ME[i], N_spins_sys[i])
        cs_shifts[i] = nuclei2vars(cs_sys[i], ordering)
        # if length(ME[i]) > 0

        #     ordering, DOF = createorderingfromeqinds(ME[i], N_spins_sys[i])
        #     cs_shifts[i] = nuclei2vars(cs_sys[i], ordering)
        # else
        #     #
        #     cs_shifts[i] = copy(cs_sys[i])
        # end
    end

    for i in eachindex(cs_singlets)
        cs_shifts[N_sys+i] = [cs_singlets[i];]
    end

    return cs_shifts
end

"""
function readbasechemshifts(P::PhysicalParamsType{T})::Vector{Vector{T}} where T

Outputs `shifts`, with index structure [spin system][cs index]. cs is in ppm. Singlet spin systems are at the end of `shifts`.
"""
function readbasechemshifts(P::PhysicalParamsType{T})::Vector{Vector{T}} where T
    return readbasechemshifts(P.cs_sys, P.cs_singlets, P.ME)
end

"""
function writebasechemshifts!(P::PhysicalParamsType{T}, shifts::Vector{Vector{T}}) where T

`shifts`: index structure [spin system][cs index]. cs has to be in ppm. Singlet spin systems are at the end of `shifts`.
"""
function writebasechemshifts!(P::PhysicalParamsType{T}, shifts::Vector{Vector{T}}) where T
    return writebasechemshifts!(P.cs_sys, P.cs_singlets, shifts, P.ME)
end

function writebasechemshifts!(
    cs_sys::Vector{Vector{T}}, # mutates. [spin system][var ind]
    cs_singlets::Vector{T}, # mutates. [spin system]
    cs_shifts::Vector{Vector{T}}, # in ppm. [spin system][var ind]. singlet spin systems are at the end.
    ME::Vector{Vector{Vector{Int}}},
    ) where T
    
    N_sys = length(ME)
    @assert N_sys == length(cs_sys)

    N = N_sys + length(cs_singlets)
    @assert length(cs_shifts) == N

    N_spins_sys = collect( length(cs_sys[i]) for i in eachindex(cs_sys))

    for i in eachindex(ME)
        ordering, DOF = createorderingfromeqinds(ME[i], N_spins_sys[i])
        cs_sys[i] = vars2nuclei(cs_shifts[i], ordering)    

        # if length(ME[i]) > 0

        #     N_shift_vars = length(ME[i])
        #     cs_shifts[i] = Vector{T}(undef, N_shift_vars)
        
        #     for j in eachindex(ME[i])
        #         k = ME[i][j][begin] # any nucleus ID work. we take the first 1.
        #         cs_sys[i][k] = cs_shifts[i][j]
        #     end
        # else
            
        #     cs_sys[i][:] = cs_shifts[i]
        #end
    end

    for i in eachindex(cs_singlets)
        cs_singlets[i] = cs_shifts[N_sys+i][begin]
    end

    return nothing

end



#### put into table.

function locatenucleisys(
    inds::Vector{Vector{Int}},
    values::Vector{Vector{T}},
    target_ind::Int,
    ) where T

    @assert length(inds) == length(values)

    for i in eachindex(inds)

        #@assert length(inds[i]) == length(values[i]) # work work for singlets.

        j = findfirst(xx->xx==target_ind, inds[i])
        
        if typeof(j) <: Integer
            return values[i][j]
        end
    end

    return convert(T, NaN)
end

function locatenucleisinglets(
    inds::Vector{Vector{Int}},
    values::Vector{T},
    target_ind::Int,
    ) where T

    @assert length(inds) == length(values)

    for i in eachindex(inds)

        #@assert length(inds[i]) == length(values[i]) # work work for singlets.

        j = findfirst(xx->xx==target_ind, inds[i])
        
        if typeof(j) <: Integer
            return values[i]
        end
    end

    return convert(T, NaN)
end

# need to document this and explain the difference between this and extractcs()
function getcs(Phy::PhysicalParamsType{T}, ID::Int) where T # hydrogen ID.

    H_ind = findfirst(xx->xx==ID, Phy.H_IDs)
    if !(typeof(H_ind) <: Integer)
        return convert(T, NaN)
    end

    for i in eachindex(Phy.H_inds_sys)
        
        out = locatenucleisys(Phy.H_inds_sys, Phy.cs_sys, H_ind)
        #@show i, out

        if isfinite(out)
            # success.
            return out
        end
    end

    for i in eachindex(Phy.H_inds_singlets)
        
        out = locatenucleisinglets(Phy.H_inds_singlets, Phy.cs_singlets, H_ind)
        #@show i, out

        if isfinite(out)
            # success.
            return out
        end
    end

    # fail.
    return convert(T, NaN)
end

"""
extractMEnuclei(
    Phys::Vector{PhysicalParamsType{T}},
) where T

Outputs: IDs, cs, entry_IDs
- IDs::Vector{Vector{Int}}
Each entry contains the set of nuclei labels that is associated with the corresponding entry in `cs`.

- cs::Vector{T}
Each entry contains the chemical shift entry of a magnetically equivalent group of nuclei. The magnetic equivalent is based on the data stored in the `ME` field of the entries in Phy`.

- entry_IDs::Vector{Tuple{Int,Int}}
Each entry contains a pair of integers that is associated with the corresponding entry in `cs`. The first integer is the molecule entry number. The second is the spin system number that the corresponding nuclei in `IDs` belongs to.
"""
function extractMEnuclei(Phys::Vector{PhysicalParamsType{T}}) where T

    IDs = Vector{Vector{Int}}(undef, 0)
    cs = Vector{T}(undef, 0)
    entry_IDs = Vector{Tuple{Int,Int}}(undef, 0) # (molecule, spin system). non-singlet spin systems appear before singlet spin systems.

    for n in eachindex(Phys)

        # non-singlet spin systems.
        for i in eachindex(Phys[n].H_inds_sys)

            all_inds = deepcopy(Phys[n].H_inds_sys[i])
            #@show n, i, all_inds
            
            # ME nuclei
            for j in eachindex(Phys[n].ME[i])

                nuclei_inds = Phys[n].ME[i][j]
                all_inds = setdiff(all_inds, nuclei_inds)

                nuclei_IDs = Phys[n].H_IDs[nuclei_inds]
                push!(IDs, nuclei_IDs)

                nuclei_cs = getcs(Phys[n], nuclei_IDs[begin])
                push!(cs, nuclei_cs)

                push!(entry_IDs, (n,i))
            end

            # non-ME nuclei
            for j in eachindex(all_inds)

                ind = all_inds[j]

                nuclei_ID = Phys[n].H_IDs[ind]
                push!(IDs, [nuclei_ID;])

                nuclei_cs = getcs(Phys[n], nuclei_ID)
                push!(cs, nuclei_cs)

                push!(entry_IDs, (n,i))
            end
        end

        # singlet spin systems.

        N_non_singlet_sys = length(Phys[n].H_inds_sys)

        for i in eachindex(Phys[n].H_inds_singlets)

            inds = Phys[n].H_inds_singlets[i]

            nuclei_IDs = Phys[n].H_IDs[inds]
            push!(IDs, nuclei_IDs)

            #@show n, i, nuclei_IDs
            nuclei_cs = getcs(Phys[n], nuclei_IDs[begin])
            push!(cs, nuclei_cs)

            push!(entry_IDs, (n, N_non_singlet_sys+i))
        end
    end

    return IDs, cs, entry_IDs
end

########################## Phys IO front_end


## initialialization routines.


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




#### load.
"""
```
getphysicalparameters(
    ::Type{T},
    target_entries::Vector{String},
    H_params_path::String,
    molecule_mapping_file_path;
    unique_cs_atol::T = convert(T, 1e-6),
) where T <: AbstractFloat
```

Assemble the J-coupling and chemical shift values from JSON files for the specified molecules.

### Nomenclature

One set of J-coupling and one set of chemical shift values for a molecule is what NMRHamiltonian.jl calls a molecule entry.

It is common to see different J-coupling values of the same molecule reported in literature for similar experimental conditions.  NMRHamiltonian.jl requires one to record each instance as separate molecule entry.

Every molecule entry is recorded as a separate JSON file with some arbitrary file name the data collector wants to use. An additional name-mapping JSON file is needed to translate what the NMRHamiltonian.jl user wants to label the entries to what is the data collector chose to name the JSON files. The name-mapping JSON file is likely to require manual set up.

### Inputs

- `::Type{T}`       -- The `AbstractFloat` datatype to use for storing floating-point data. A viable input here is `Float64`, which specifies double precision floating-point numbers are to be used.
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
    ) where T <: AbstractFloat
    
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
