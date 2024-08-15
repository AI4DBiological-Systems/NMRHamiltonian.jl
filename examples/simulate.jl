# run a.jl first.

include("./helpers/utils.jl")

PLT.close("all")
fig_num = 1

#PLT.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

T = Float64
#T = Float32

molecule_mapping_file_path = joinpath(pwd(), "ref_files", "molecule_name_mapping", "demo_compounds.json")
H_params_path = joinpath(pwd(), "ref_files", "coupling_info")

### user inputs.

molecule_entries = [
    "L-Serine";
    "alpha-D-Glucose";
    "beta-D-Glucose";  
    "L-Isoleucine";
    "L-Glutamine";
    "L-Valine";
    "L-Leucine";
    "DSS";
    "Singlet - 4.9 ppm";
]

# machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
fs, SW, ν_0ppm = HAM.getpresetspectrometer(T, "700")

config = HAM.SHConfig{T}(
    coherence_tol = convert(T, 0.01),
    relative_α_threshold = convert(T, 0.001),
    max_deviation_from_mean = convert(T, 0.05),
    acceptance_factor = convert(T, 0.99),
    total_α_threshold = zero(T),
)
unique_cs_digits = 6

println("Timing: HAM.loadandsimulate")
@time Phys, As, MSPs = HAM.loadandsimulate(
    fs, SW, ν_0ppm,
    molecule_entries,
    H_params_path,
    molecule_mapping_file_path,
    config;
    unique_cs_digits = unique_cs_digits,
);

rcs = collect(HAM.get_rcs(A) for A in As)

nothing