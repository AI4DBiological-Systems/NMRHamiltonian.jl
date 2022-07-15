
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
molecule_names = ["Agmatine"; ]
#molecule_names = ["L-Glutamine"; ]
#molecule_names = ["alpha-D-Glucose"; ]
#molecule_names = ["beta-D-Glucose"; ];
#molecule_names = ["Ethanol"; ];
#molecule_names = ["L-Leucine"; ]

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
    prune_combo_Δc_flag = false)
#
@time Bs = NMRHamiltonian.setupmixtureSH(molecule_names,
    fs, SW, ν_0ppm,
    Phys;
    config_path = SH_config_path,
    prune_combo_Δc_flag = true)

#
i = 1
println("No combo reducing: Spin group $(i):")
println("Number of resonance groups: ", length(As[1].Δc_bar[i]))
println("Number of components: ", length(As[1].αs[i]))
println()

println("Combo reduced: Spin group $(i):")
println("Number of resonance groups: ", length(Bs[1].Δc_bar[i]))
println("Number of components: ", length(Bs[1].αs[i]))
println()

fig_num = plotabsorptionlorentzians(Bs[1].αs[i], Bs[1].Ωs[i], λ0, fs, SW, ν_0ppm, fig_num;
    a_ref = As[1].αs[i],
    F_ref = As[1].Ωs[i],
    title_string = "$(spectrometer_freq) MHz: Pruned vs reference",
    P_min = P_min, P_max = P_max)

#cb = Δc_bar[i]

# same_flags = collect( isapprox(cb[k], -cb[15]; atol = 1e-1) for k = 1:length(cb) )
# count(same_flags)

#cb_mat = array2matrix(cb)

# currently trying to figure out if no-pruning leads to too many components.
# might be true for alpha-Glucose.
# try brai metabolite study for J-coupling.
