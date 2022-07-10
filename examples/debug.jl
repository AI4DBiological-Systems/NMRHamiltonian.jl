
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

# include("./helpers/group.jl")
include("./helpers/utils.jl")

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

#molecule_names = ["alpha-D-Glucose"; ]
#molecule_names = ["beta-D-Glucose"; ];
molecule_names = ["Ethanol"; ];

prune_combo_Δc_flag = true
#prune_combo_Δc_flag = false

simple_coherence_atol = 1e-2

#fs, SW, ν_0ppm = fetchsamplemachinesettings("400")
#fs, SW, ν_0ppm = fetchsamplemachinesettings("500")
#fs, SW, ν_0ppm = fetchsamplemachinesettings("600")
#fs, SW, ν_0ppm = fetchsamplemachinesettings("700")
fs, SW, ν_0ppm = fetchsamplemachinesettings("900")

# path to the json file that provides the mapping from a compound name to its spin system info file name.
H_params_path = "/home/roy/Documents/repo/NMRData/input/coupling_info"
dict_compound_to_filename = JSON.parsefile("/home/roy/Documents/repo/NMRData/input/compound_mapping/select_compounds.json")

### end inputs.

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

println("Timing: getphysicalparameters")
@time Phys = NMRHamiltonian.getphysicalparameters(molecule_names,
    H_params_path,
    dict_compound_to_filename;
    unique_cs_atol = 1e-6)

#
println("Timing: setupmixtureSH()")
@time mixture_params = NMRHamiltonian.setupmixtureSH(molecule_names,
    fs, SW, ν_0ppm,
    Phys;
    config_path = SH_config_path,
    # tol_coherence = tol_coherence_default,
    # α_relative_threshold = α_relative_threshold_default,
    # Δc_partition_radius = Δc_partition_radius_default,
    prune_combo_Δc_flag = prune_combo_Δc_flag,
    simple_coherence_atol = simple_coherence_atol)
As = mixture_params


# allabssmaller(cb[4],1+1e-2)
# sumabssmaller(cb[4],1+(1e-2)*length(cb[4]))


### plot.
# This is the frequency range that we shall work with.
P = LinRange(-0.5, 4, 50000)
U = ppm2hzfunc.(P)
U_rad = U .* (2*π)

# absorption Lorentzian.
q = uu->evalzerophaseclmixture(uu, As, λ0)
q_singlets = uu->evalzerophaseclsinglets(uu, A, λ0)



A = As[1]
qs = getqs(A, λ0)

# function plotallgroups(q, qs, q_singlets, U_rad, fig_num)
#     PyPlot.figure(fig_num)
#     fig_num += 1
#
#     PyPlot.plot(P, real.(q.(U_rad)), label = "sum", "--")
#
#     for i = 1:length(qs)
#         for k = 1:length(qs[i])
#             PyPlot.plot(P, real.(qs[i][k].(U_rad)), label = "group ($(i),$(k))")
#         end
#     end
#
#     PyPlot.plot(P, real.(q_singlets.(U_rad)), label = "singlets")
#
#     PyPlot.gca().invert_xaxis()
#
#     PyPlot.legend()
#     PyPlot.title("groups")
#
#     return fig_num
# end
# fig_num = plotallgroups(q, qs, q_singlets, U_rad, fig_num)

spectrometer_freq = round(fs/SW, digits = 2)

using Plots; plotly()
plots_save_path = "/home/roy/Documents/workspaces/WIP/debug/$(molecule_names[1])_$(spectrometer_freq)-$(prune_combo_Δc_flag).html"


title_string = "$(spectrometer_freq)-$(prune_combo_Δc_flag)"
canvas_size = (1500, 600)
plot_obj, q_U, qs_U, q_singlets_U = plotgroups(title_string, P, U, q, qs, q_singlets, real, 1.0;
canvas_size = canvas_size)

if save_fig_flag
    Plots.savefig(plot_obj, plots_save_path)
end
display(plot_obj)

S = collect( sum(A.Δc_bar[1][k]) for k = 1:length(A.Δc_bar[1]) )
F = collect( NMRHamiltonian.allabssmaller(A.Δc_bar[1][k], 1+1e-2) for k = 1:length(A.Δc_bar[1]) )

c = deepcopy(A.Δc_m_compound)
cb = deepcopy(A.Δc_bar)

i = 1
M = array2matrix(cb[i])
M = [collect(1:size(M,2))'; M]
println("$(spectrometer_freq) MHz, Δc_bar[$(i)] (rows 2 to end)")
display(M)
