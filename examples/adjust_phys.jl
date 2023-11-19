include("a.jl")


import Random
Random.seed!(25)


PLT.close("all")
fig_num = 1

#T = Float32
T = Float64

### user inputs.
fs, SW, ν_0ppm = HAM.getpresetspectrometer(T, "700")
#fs, SW, ν_0ppm = 9615.38461538462, 16.0196917451925, 6752.905945857455

λ0 = convert(T, 3.4)

## pull the sample coupling values into dictionary structures.

root_data_path = DS.getdatapath(DS.NMR2023()) # coupling values data repository root path

H_params_path = joinpath(root_data_path, "coupling_info") # folder of coupling values. # replace with your own values in actual usage.

molecule_mapping_root_path = joinpath(
    root_data_path,
    "molecule_name_mapping",
)
molecule_mapping_file_path = joinpath(
    molecule_mapping_root_path,
    "select_molecules.json",
)

config = HAM.SHConfig{T}(
    coherence_tol = convert(T, 0.01),
    relative_α_threshold = convert(T, 0.005),
    tol_radius_1D = convert(T, 0.01), # strictly between 0 and 1. The lower, the better the approximation, but would a larger partition (i.e. more resonance groups).
    nuc_factor = convert(T, 2.0),
)
unique_cs_atol = convert(T, 1e-6)
unique_J_avg_atol = convert(T, 1e-6)

molecule_entries = [
    "L-Phenylalanine";
    "L-Valine";
    "L-Leucine";
    "L-Isoleucine";
    "L-Glutamine";
    "L-Proline";
    "L-"
    "alpha-D-Glucose";
    "beta-D-Glucose";
    "DSS";
]

println("Timing: getphysicalparameters")
@time Phys = NMRHamiltonian.getphysicalparameters(
    T,
    molecule_entries,
    H_params_path,
    molecule_mapping_file_path;
    unique_cs_atol = convert(T, 1e-6),
)

# DSS and Phenylalanine might have magnetic equivalence, but the coupling constants from the sources of PublicationDatasets do not reflect this.

Phys2 = HAM.createmaximalME(
    Phys;
    unique_cs_atol = unique_cs_atol,
    unique_J_avg_atol = unique_J_avg_atol,
)

# the first and last molecules, Phenylalanine and DSS, should have ME, but doesn't.
println("Non-adjusted:")
for i in eachindex(Phys)
    println("$(molecule_entries[i])")
    @show Phys[i].cs_sys
    @show Phys[i].ME
    println()
end
println()

# after adjustment, they do.
println("Adjusted:")
for i in eachindex(Phys2)
    println("$(molecule_entries[i])")
    @show Phys2[i].cs_sys
    @show Phys2[i].ME
    println()
end

# save the adjusted/checked physical chemistry parameters to JSON.
save_file_names = collect( "$(molecule_entries[n])_checked.json" for n in eachindex(molecule_entries) )
save_folder = "./files/"
HAM.savecouplinginfo(Phys2, save_file_names; save_folder = save_folder)

nothing