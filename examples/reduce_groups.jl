
include("../src/NMRHamiltonian.jl")
import .NMRHamiltonian
#import NMRHamiltonian
import NMRSignalSimulator

using LinearAlgebra
using FFTW
import PyPlot
import JSON

import Statistics

#using PyPlot, PyCall
import PyPlot

include("./helpers/utils.jl")
#include("./helpers/coherence.jl")

import Random
Random.seed!(25)

PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])


### user inputs.

save_folder = "/home/roy/MEGAsync/outputs/NMR/groups/publication/static"
save_fig_flag = true
#save_fig_flag = false

λ0 = 3.4

SH_config_path = "/home/roy/Documents/repo/NMRData/input/SH_configs/select_compounds_SH_configs_reduce.json"
#SH_config_path = "/home/roy/Documents/repo/NMRData/input/SH_configs/select_compounds_SH_configs_low_intensity_threshold.json"
surrogate_config_path = "/home/roy/Documents/repo/NMRData/input/surrogate_configs/select_compounds_SH_configs.json"

#molecule_names = ["L-Serine"; ]
#molecule_names = ["Agmatine"; ]
#molecule_names = ["L-Glutamine"; ]
molecule_names = ["alpha-D-Glucose"; ]
#molecule_names = ["beta-D-Glucose"; ]; # undertest.
#molecule_names = ["Ethanol"; ];
#molecule_names = ["L-Leucine"; ]
#molecule_names = ["L-Phenylalanine"; ]

#molecule_names = ["alpha-D-Glucose-brain"; ]
#molecule_names = ["beta-D-Glucose-brain"; ];
#molecule_names = ["alpha-D-Glucose-brain"; "beta-D-Glucose-brain";]

P_min = Inf
P_max = Inf
# P_min = 3.6
# P_max = 4

#fs, SW, ν_0ppm = fetchsamplemachinesettings("400")
#fs, SW, ν_0ppm = fetchsamplemachinesettings("500")
#fs, SW, ν_0ppm = fetchsamplemachinesettings("600")
fs, SW, ν_0ppm = fetchsamplemachinesettings("700")
#fs, SW, ν_0ppm = fetchsamplemachinesettings("900")

# path to the json file that provides the mapping from a compound name to its spin system info file name.
H_params_path = "/home/roy/Documents/repo/NMRData/input/coupling_info"
dict_compound_to_filename = JSON.parsefile("/home/roy/Documents/repo/NMRData/input/compound_mapping/select_compounds.json")

### end inputs.

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)
spectrometer_freq = round(fs/SW, digits = 2)

println("Timing: getphysicalparameters")
@time Phys = NMRHamiltonian.getphysicalparameters(molecule_names,
    H_params_path,
    dict_compound_to_filename;
    unique_cs_atol = 1e-6)

#
println("Timing: setupmixtureproxies()")
@time As = NMRHamiltonian.setupmixtureSH(molecule_names,
    fs, SW, ν_0ppm,
    Phys;
    config_path = SH_config_path,
    prune_Δc_option = 0,
    normalize_α_for_spin_sys = false)
#
@time Bs = NMRHamiltonian.setupmixtureSH(molecule_names,
    fs, SW, ν_0ppm,
    Phys;
    config_path = SH_config_path,
    prune_Δc_option = 1,
    normalize_α_for_spin_sys = false)
#
@time Cs = NMRHamiltonian.setupmixtureSH(molecule_names,
    fs, SW, ν_0ppm,
    Phys;
    config_path = SH_config_path,
    prune_Δc_option = 2,
    normalize_α_for_spin_sys = false)
#
@time Ds = NMRHamiltonian.setupmixtureSH(molecule_names,
    fs, SW, ν_0ppm,
    Phys;
    config_path = SH_config_path,
    prune_Δc_option = 3,
    normalize_α_for_spin_sys = false)
#
@time Es = NMRHamiltonian.setupmixtureSH(molecule_names,
    fs, SW, ν_0ppm,
    Phys;
    config_path = SH_config_path,
    prune_Δc_option = 4,
    normalize_α_for_spin_sys = false)

#
i = 1
println("A: No combo reducing: Spin group $(i):")
println("Number of resonance groups: ", length(As[1].Δc_bar[i]))
println("Number of components: ", length(As[1].αs[i]))
println()

println("B: Combo reduced: Spin group $(i):")
println("Number of resonance groups: ", length(Bs[1].Δc_bar[i]))
println("Number of components: ", length(Bs[1].αs[i]))
println()

println("C: Combo reduced: Spin group $(i):")
println("Number of resonance groups: ", length(Cs[1].Δc_bar[i]))
println("Number of components: ", length(Cs[1].αs[i]))
println()

println("D: Combo reduced: Spin group $(i):")
println("Number of resonance groups: ", length(Ds[1].Δc_bar[i]))
println("Number of components: ", length(Ds[1].αs[i]))
println()

println("E: Combo reduced: Spin group $(i):")
println("Number of resonance groups: ", length(Es[1].Δc_bar[i]))
println("Number of components: ", length(Es[1].αs[i]))
println()

fig_num = plotabsorptionlorentzians(Bs[1].αs[i], Bs[1].Ωs[i], λ0, fs, SW, ν_0ppm, fig_num;
    a_ref = As[1].αs[i],
    F_ref = As[1].Ωs[i],
    title_string = "$(spectrometer_freq) MHz: B vs reference",
    P_min = P_min, P_max = P_max)
#
fig_num = plotabsorptionlorentzians(Cs[1].αs[i], Cs[1].Ωs[i], λ0, fs, SW, ν_0ppm, fig_num;
    a_ref = As[1].αs[i],
    F_ref = As[1].Ωs[i],
    title_string = "$(spectrometer_freq) MHz: C vs reference",
    P_min = P_min, P_max = P_max)
#
fig_num = plotabsorptionlorentzians(Ds[1].αs[i], Ds[1].Ωs[i], λ0, fs, SW, ν_0ppm, fig_num;
    a_ref = As[1].αs[i],
    F_ref = As[1].Ωs[i],
    title_string = "$(spectrometer_freq) MHz: D vs reference",
    P_min = P_min, P_max = P_max)
#
#
fig_num = plotabsorptionlorentzians(Es[1].αs[i], Es[1].Ωs[i], λ0, fs, SW, ν_0ppm, fig_num;
    a_ref = As[1].αs[i],
    F_ref = As[1].Ωs[i],
    title_string = "$(spectrometer_freq) MHz: E vs reference",
    P_min = P_min, P_max = P_max)
#
fig_num = plotabsorptionlorentzians(Ds[1].αs[i], Ds[1].Ωs[i], λ0, fs, SW, ν_0ppm, fig_num;
    a_ref = Bs[1].αs[i],
    F_ref = Bs[1].Ωs[i],
    title_string = "$(spectrometer_freq) MHz: D vs B",
    P_min = P_min, P_max = P_max)
#
fig_num = plotabsorptionlorentzians(Ds[1].αs[i], Ds[1].Ωs[i], λ0, fs, SW, ν_0ppm, fig_num;
    a_ref = Es[1].αs[i],
    F_ref = Es[1].Ωs[i],
    title_string = "$(spectrometer_freq) MHz: D vs E",
    P_min = P_min, P_max = P_max)

##### some more info for option 1 (As) and option 2 (Bs).

F_A = As[1].Ωs[i] ./ (2*π)
P_A = hz2ppmfunc.(F_A)
inds_A = As[1].part_inds_compound[i]

F_B = Bs[1].Ωs[i] ./ (2*π)
P_B = hz2ppmfunc.(F_B)
inds_B = Bs[1].part_inds_compound[i]

cb_B = Bs[1].Δc_bar[i]
cb_A = As[1].Δc_bar[i]

P_A_groups = collect( P_A[inds_A[k]] for k = 1:length(inds_A) )
P_A_avgs = collect( Statistics.mean(P_A[inds_A[k]]) for k = 1:length(inds_A) )

P_B_groups = collect( P_B[inds_B[k]] for k = 1:length(inds_B) )
P_B_avgs = collect( Statistics.mean(P_B[inds_B[k]]) for k = 1:length(inds_B) )
sort_inds_B = sortperm(P_B_avgs)

MB = [array2matrix(cb_B); P_B_avgs']
MA = [array2matrix(cb_A); P_A_avgs']
