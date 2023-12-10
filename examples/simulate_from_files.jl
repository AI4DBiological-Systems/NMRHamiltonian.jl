# run a.jl first.

include("./helpers/utils.jl")

PLT.close("all")
fig_num = 1

#PLT.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

T = Float64
#T = Float32
#T = Float16 # will result in numerical precision-related errors. Please use Float32 or above.

### user inputs.

molecule_entries = [
    "L-Isoleucine";
    "L-Leucine";
    "L-Valine";
    "DSS";
    "Singlet - 4.9 ppm";
]

# machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
fs, SW, ν_0ppm = HAM.getpresetspectrometer(T, "700")

config = HAM.SHConfig{T}(
    coherence_tol = convert(T, 0.01),
    relative_α_threshold = convert(T, 0.005),
    max_deviation_from_mean = convert(T, 0.2), # positive number. Larger means less resonance groups, but worse approximation for each group.
    acceptance_factor = convert(T, 0.99), # keep this close to 1. Takes values from (0,1).
    total_α_threshold = convert(T, 0.01), # final intensity pruning.
)

# ## Identify where to get the reference J-coupling and chemical shift files.

# we should how to use individual files in this example.
molecule_mapping_file_path = "./files/molecule_name_mapping/demo_compounds.json"
H_params_path = "./files/coupling_info/"

### end inputs.



println("Timing: getphysicalparameters")
@time Phys = NMRHamiltonian.getphysicalparameters(
    T,
    molecule_entries,
    H_params_path,
    molecule_mapping_file_path;
    unique_cs_digits = 6,
)

println("Timing: simulate()")
@time As, MSPs = HAM.simulate(
    Phys,
    molecule_entries,
    fs,
    SW,
    ν_0ppm,
    config,
)
println()

println()

cs_values = collect( HAM.readbasechemshifts(Phys[n]) for n in eachindex(Phys) )
println("The reference chemical shift values (in ppm) for the Δc features: [compoune][spin system][ME nuclei]")
display(cs_values)
println()

del_c_bars = collect( As[n].Δc_bar for n in eachindex(As) )
N_res_groups = collect( length.(As[n].Δc_bar) for n in eachindex(As) )
N_components = collect( length.(As[n].αs) for n in eachindex(As) )
N_nucs = collect( 
    collect(
        length(As[n].Δc[i][begin])
        for i in eachindex(As[n].Δc)
    ) 
    for n in eachindex(As)
)
println("[molecule_entries N_nucs N_res_groups N_components]:")
display([molecule_entries N_nucs N_res_groups N_components])
println()

nothing