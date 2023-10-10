

# continues from the end of simulate_SH.jl

### test serialization.
using Test
import JSON3
import BSON
include("../test/helpers/roundtrip.jl")
roundtripJSON(As)
roundtripBSON(As)



include("../test/helpers/roundtrip.jl")
roundtripJSON(Phys, molecule_entries)
roundtripBSON(Phys, molecule_entries)

# # Table readout
n_select = 1 # 3.
h = xx->NMRHamiltonian.getcs(Phys[n_select], xx)
tmp = h.(Phys[n_select].H_IDs)

col_nuclei_IDs, col_cs, col_entry_IDs = NMRHamiltonian.extractMEnuclei(Phys)
tab1 = [col_entry_IDs col_nuclei_IDs col_cs]


# # Phys write to and read from JSON files.
# This allows us to re-do the SH simulation from updated chemical shift or J-coupling values from a new file.

# We'll reassign some chemical shifts on purpose, and see if roundtrip serialization works.
Phys2 = deepcopy(Phys)

# Valine, 2nd ME nuclei group (which is nuclei 1, 2, 3)
molecule_select = 2
spin_system_select = 1
subsystem_select = 2
cs_shifts = NMRHamiltonian.readbasechemshifts(Phys[molecule_select])
cs_shifts[spin_system_select][subsystem_select] = 0.14

# save.
NMRHamiltonian.writebasechemshifts!(Phys2[molecule_select], cs_shifts)

# Valine, 2nd ME nuclei group.
molecule_select = 4
spin_system_select = 2
subsystem_select = 1
cs_shifts = NMRHamiltonian.readbasechemshifts(Phys[molecule_select])
cs_shifts[spin_system_select][subsystem_select] = -0.1

# save.
NMRHamiltonian.writebasechemshifts!(Phys2[molecule_select], cs_shifts)


# extract the complete list of chemical shifts and J-coupling for all nuclei. We can inspect these to see what will be saved to file.
H_IDs_set = collect( NMRHamiltonian.extractcouplinginfo(Phys2[n])[1] for n in eachindex(Phys2) )
H_css_set = collect( NMRHamiltonian.extractcouplinginfo(Phys2[n])[2] for n in eachindex(Phys2) )
J_IDs_set = collect( NMRHamiltonian.extractcouplinginfo(Phys2[n])[3] for n in eachindex(Phys2) )
J_vals_set = collect( NMRHamiltonian.extractcouplinginfo(Phys2[n])[4] for n in eachindex(Phys2) )



# file names, not file paths. These files should be in the same folder.
coupling_JSON_file_names = [
    "glutamine_new.json";
    "isoleucine_new.json";
    "phenylalanine_new.json";
    "DSS_new.json";
]

NMRHamiltonian.savecouplinginfo(
    Phys2,
    coupling_JSON_file_names;
    save_folder = joinpath(pwd(), "files/new_coupling_info"),
)

# create compound name mapping dictionary JSON:

load_paths, dict_molecule_to_filename = NMRHamiltonian.getloadpaths(
    molecule_entries, H_params_path, molecule_mapping_file_path,
)

molecule_mapping_file_path_new = joinpath(pwd(),"files/new_molecule_names_mapping.json")
H_params_path_new = joinpath(pwd(), "files/new_coupling_info")

# save to JSON.
label_names = [
    "L-Glutamine new";
    "L-Isoleucine new";
    "L-Phenylalanine new";
    "DSS new";
]

notes = collect( "fitted result" for _ in eachindex(label_names) )

dic_name_mapping_new = NMRHamiltonian.createnamemappingJSON(
    molecule_mapping_file_path_new,
    label_names,
    coupling_JSON_file_names;
    notes = notes,
)

# load from JSON.
@time Phys_new = NMRHamiltonian.getphysicalparameters(
    Float64,
    label_names,
    H_params_path_new,
    molecule_mapping_file_path_new;
    unique_cs_atol = 1e-6,
)

@show Phys[2].cs_sys # original.
@show Phys2[2].cs_sys # what should be saved.
@show Phys_new[2].cs_sys # what is loaded. should match what is saved.
println()

@show Phys[4].cs_singlets # original.
@show Phys2[4].cs_singlets # what should be saved.
@show Phys_new[4].cs_singlets # what is loaded. should match what is saved.
println()
