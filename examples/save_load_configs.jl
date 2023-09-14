# include("../src/NMRHamiltonian.jl")
# import .NMRHamiltonian
# # remove later. dependencies.
# import Kronecker, Graphs, ConvexClustering, Distances
# # end remove.

import NMRHamiltonian

using DataDeps, CodecZlib, Tar

using LinearAlgebra
import PythonPlot
import JSON3

PythonPlot.close("all")
fig_num = 1

# PythonPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

#### user inputs.

molecule_entries = ["D-(+)-Glucose"; "Ethanol"; "L-Methionine"; 
    "L-Phenylalanine"; "L-Glutathione reduced"; "L-Glutathione oxidized";
    "HEPES"; "L-Glutamine"; "DSS"; "D2O"]

# machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
fs = 14005.602240896402
SW = 20.0041938620844
ν_0ppm = 10656.011933076665

# # machine values for the BMRB 500 MHz glucose experiment.
# ν_0ppm = 6752.490995937095
# SW = 16.0196917451925
# fs = 9615.38461538462

## pull the sample coupling values into dictionary structures.
include("./helpers/data.jl") # use DataDeps.jl and Tar.jl to download and extract the sample coupling values.

root_data_path = getdatapath() # coupling values data repository root path

H_params_path = joinpath(root_data_path, "coupling_info") # folder of coupling values. # replace with your own values in actual usage.

molecule_mapping_root_path = joinpath(root_data_path, "molecule_name_mapping")
molecule_mapping_file_path = joinpath(molecule_mapping_root_path, "select_molecules.json")
#molecule_mapping_file_path = joinpath(molecule_mapping_root_path, "GISSMO_names.json")

dict_molecule_to_filename = JSON3.read(read(molecule_mapping_file_path)) # map molecule entries to coupling information file names.

### end inputs.


#### set up

# get the chemical shift values from file.
println("Timing: getphysicalparameters")
@time Phys = NMRHamiltonian.getphysicalparameters(
    Float64,
    molecule_entries,
    H_params_path,
    molecule_mapping_file_path;
    unique_cs_atol = 1e-6)

# extract chemical shifts for spin systems and singlets.
cs_sys_mixture, cs_singlets_mixture = NMRHamiltonian.extractcs(Phys)



#### generate default values. Use to assemble configs later.

# one file.
α_tols = NMRHamiltonian.generatemixturespinsystemvector(cs_sys_mixture, 0.05)
coherence_tols = NMRHamiltonian.generatemixturespinsystemvector(cs_sys_mixture, 0.01)
normalize_αs = NMRHamiltonian.generatemixturespinsystemvector(cs_sys_mixture, true)

# one file.
manual_radi = NMRHamiltonian.generatemixturespinsystemvector(cs_sys_mixture, 0.1)
manual_knns = NMRHamiltonian.generatemixturespinsystemvector(cs_sys_mixture, 64)
manual_θs = NMRHamiltonian.generatemixturespinsystemvector(cs_sys_mixture, 0.1) # for SqExp
manual_γs = NMRHamiltonian.generatemixturespinsystemvector(cs_sys_mixture, 0.4)

####
# insert code to modify the default values.
####

#### save configs.
file_folder = "./configs"

NMRHamiltonian.savemixtureshsconfig(file_folder, "mixture_SH.json", molecule_entries, α_tols, coherence_tols, normalize_αs)
NMRHamiltonian.savemixturepartitionsparameters(file_folder, "mixture_partition.json", molecule_entries, manual_θs, manual_γs, manual_knns, manual_radi)

#### load configs

# create container that stores the loaded results. Populate with default values.
mixture_sh_config = NMRHamiltonian.defaultmixtureshsconfig(
    cs_sys_mixture;
    normalize_α = false
)
mixture_parts_params = NMRHamiltonian.defaultmixturepartitionsparameters(
    cs_sys_mixture,
    1.0
)


import JSON3

file_name = "mixture_SH.json"
json_string = read(joinpath(file_folder, file_name))
S = JSON3.read(json_string)


vars_target = [:coherence_tols; :relative_α_thresholds] # should be a subset of propertynames(sh_config)

# leave out default_normalize_α on purpose to see if it updates. the file has it set to true.
NMRHamiltonian.updatefromfile!( mixture_sh_config,
    vars_target, 
    NMRHamiltonian.getshconfigtolabelsmapping(), S,
    molecule_entries[1:3],
    molecule_entries)

println("mixture_sh_config")
display(mixture_sh_config)
#@assert 1==2

##

file_name = "mixture_partition.json"
json_string = read(joinpath(file_folder, file_name))
S = JSON3.read(json_string)

# leave out default_normalize_α on purpose to see if it updates. the file has it set to true.
vars_target = [:θs; :radi]
NMRHamiltonian.updatefromfile!( mixture_parts_params,
    vars_target,
    NMRHamiltonian.getpartitionconfigtolabelsmapping(), S,
    molecule_entries[4:end],
    molecule_entries)
println("mixture_parts_params")
display(mixture_parts_params)


mixture_sh_config2 = NMRHamiltonian.loadmixtureshsconfig(
    cs_sys_mixture,
    joinpath(file_folder, "mixture_SH.json"),
    molecule_entries,
    Float64,
)

mixture_parts_params2 = NMRHamiltonian.loadmixturepartitionsparameters(
    cs_sys_mixture,
    joinpath(file_folder, "mixture_partition.json"),
    molecule_entries,
    2.0
)

# try mixture_cs_sys and mixture_cs_sys.
# consolidate this back to mixture_type1.jl

# variables = load_data()                   # Load data from function load_data()

# name_list =[:A, :C, :F, :No]                # typically 3 element string vector although could be 4 or 5 (user defined from a read in file - different every time)
# master_list = fieldnames(MyStruct)
# base_case = Vector{Any}(undef, 3)
# for (i,n) in enumerate(name_list)
#     if n in master_list
#         base_case[i] = getfield(variables, n)
#     end
# end