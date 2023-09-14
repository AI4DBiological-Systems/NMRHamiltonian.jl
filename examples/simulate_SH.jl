# include("../src/NMRHamiltonian.jl")
# import .NMRHamiltonian
# # remove later. dependencies.
# import Kronecker, Graphs, ConvexClustering, Distances
# # end remove.

import NMRHamiltonian

using DataDeps, CodecZlib, Tar

using LinearAlgebra
import PythonPlot # do Pkg.add("PythonPlot") if this is missing.
#import JSON3

#import Optim


include("./helpers/utils.jl")

PythonPlot.close("all")
fig_num = 1

#PythonPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

T = Float64

### user inputs.

# molecule_entries = ["D-(+)-Glucose"; "Ethanol"; "L-Methionine";
#     "L-Phenylalanine"; "L-Glutathione reduced"; "L-Glutathione oxidized";
#     "HEPES"; "L-Glutamine"; "DSS"; "D2O"]

#molecule_entries = ["L-Methionine"; "L-Phenylalanine"; "DSS"; "Ethanol"; "L-Isoleucine"]
#molecule_entries = ["D-(+)-Glucose"; "Ethanol"; "DSS"; "D2O"; ]

#molecule_entries = ["D2O, 4.8ppm"]
#molecule_entries = ["DSS"]
#molecule_entries = ["L-Histidine - 2000 study";]
#molecule_entries = ["L-Histidine, 500 MHz";]
#molecule_entries = ["L-Histidine";]
#molecule_entries = ["DL-Serine"]

#molecule_entries = ["3-Hydroxybutyrate"]
#molecule_entries = ["L-Glutathione reduced"]
#molecule_entries = ["L-Glutathione oxidized"]
#molecule_entries = ["Uridine"]
#molecule_entries = ["HEPES";]

#molecule_entries = ["Ethanol";]
#molecule_entries = ["L-Phenylalanine";]

molecule_entries = ["L-Glutamine"; "L-Valine"; "L-Phenylalanine"; "DSS"]
#γ_base = 0.01

#molecule_entries = ["L-Serine";]

# molecule_entries = ["Agmatine";]
# γ_base = 1.0

#molecule_entries = ["Gamma-Aminobutyric acid"]
#molecule_entries = ["L-Valine";]

#molecule_entries = ["L-Isoleucine";]
#γ_base = 1e-6

# molecule_entries = ["L-Leucine";]
# γ_base = 0.1

# molecule_entries = [
#     "D-(+)-Glucose";
#     "Ethanol";
#     "L-Methionine";     
#     "L-Phenylalanine";
#     "L-Glutathione reduced";
#     "L-Glutathione oxidized";       
#     "HEPES";
#     "L-Glutamine";
#     "DSS";
#     "D2O";
# ]
# γ_base = 0.01
# γ_rate = 1.1
# max_iter_γ = 5

# molecule_entries = ["DSS";]
# γ_base = 0.1
# γ_rate = 1.05
# max_iter_γ = 100

# machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
fs, SW, ν_0ppm = fetchsamplemachinesettings("700")

## pull the sample coupling values into dictionary structures.
include("./helpers/data.jl") # use DataDeps.jl and Tar.jl to download and extract the sample coupling values.

root_data_path = getdatapath() # coupling values data repository root path

H_params_path = joinpath(root_data_path, "coupling_info") # folder of coupling values. # replace with your own values in actual usage.

molecule_mapping_root_path = joinpath(root_data_path, "molecule_name_mapping")
molecule_mapping_file_path = joinpath(molecule_mapping_root_path, "select_molecules.json")
#molecule_mapping_file_path = joinpath(molecule_mapping_root_path, "GISSMO_names.json")

## doesn't produce partition tree due to numerical accuracy. 
# cc_gap_tol = 1e-6
# cc_max_iters = 300
# assignment_zero_tol = 1e-6

# These produce a aprtition tree over the searched γ values for many molecule entries.

config_SH = HAM.SHConfig(
    T;
    max_partition_size_offset = 0,
    partition_rate = 2.0,
    γ_base = 0.1,
    γ_rate = 1.05,
    max_iter_γ = 100,
    fully_connected_convex_clustering = false, #  overides all knn-related optional inputs
    length_scale_base = 10.0,
    length_scale_rate = 0.7,
    min_dynamic_range = 0.95,
    cc_gap_tol = 1e-8,
    cc_max_iters = 300,
    assignment_zero_tol = 1e-3,
    report_γ_cost = true,
    verbose_kernel = true,
)

### end inputs.



println("Timing: getphysicalparameters")
@time Phys = NMRHamiltonian.getphysicalparameters(
    Float64,
    molecule_entries,
    H_params_path,
    molecule_mapping_file_path;
    unique_cs_atol = 1e-6,
)

println("Timing: simulate()")
@time As, Rs, MSPs = HAM.simulate(
    Phys,
    molecule_entries,
    fs,
    SW,
    ν_0ppm,
    config_SH,
)

# view the knn or radius actually used for molecule 1, spin system 1.
@show Rs[1][1].connectivity_label
connectivity_parameter = eval(Meta.parse(string(Rs[1][1].connectivity_parameter)))
@show connectivity_parameter

println("Are the returned partitions nested successively? ", all(NMRHamiltonian.isnestedsuccessively(Rs[1][1].Gs)))

### visualize a target molecule and spin group, given a T-2 decay parameter.
λ0 = 3.4 # T-2 decay parameter.
molecule_select = 1
spin_system_select = 1
ppm_offset = 0.2 # for display padding.
N_viz = 50000

a = As[molecule_select].αs[spin_system_select]
F = As[molecule_select].Ωs[spin_system_select]

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

ΩS_ppm = hz2ppmfunc.( F ./ (2*π) )
ΩS_ppm_sorted = sort(ΩS_ppm)

u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - ppm_offset)
u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + ppm_offset)

P_min = hz2ppmfunc(u_min)
P_max = hz2ppmfunc(u_max)

P = LinRange(P_min, P_max, N_viz)
U = ppm2hzfunc.(P)
U_rad = U .* (2*π)

# absorption Lorentzian.

function evalzerophasecl1Darray(u_rad, αs::Vector{T}, Ωs::Vector{T}, λ::T)::Complex{T} where T <: AbstractFloat

    out = zero(Complex{T})
    for l = 1:length(αs)
        out += evalzerophaseclpartitionelement(u_rad, αs[l], Ωs[l], λ)
    end

    return out
end

function evalzerophaseclpartitionelement(r,
    α::T, Ω::T, λ::T)::Complex{T} where T <: AbstractFloat

    out = α/(λ+im*(r-Ω))

    return out
end

q = uu->evalzerophasecl1Darray(uu, a, F, λ0)
q_U = q.(U_rad)

# plot.
PythonPlot.figure(fig_num)
fig_num += 1

PythonPlot.plot(P, real.(q_U))
PythonPlot.gca().invert_xaxis()

PythonPlot.ylabel("real part")
#PythonPlot.legend()
PythonPlot.title("spectrum of $(molecule_entries[molecule_select]), spin system $(spin_system_select)")



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

nothing