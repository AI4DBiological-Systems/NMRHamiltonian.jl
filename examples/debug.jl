
include("/home/roy/Documents/repo/NMRHamiltonian.jl/src/NMRHamiltonian.jl")
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

# tol_coherence_default = 1e-2 # resonances are pairs of eigenvalues of the Hamiltonian that have quantum numbers that differ by -1. This is the tolerance away from -1 that is allowed.
# #α_relative_threshold_default = 0.05 # resonances with relative amplitude less than this factor compared to the maximum resonance in the spin group will be removed. Set to 0.0 to see every single resonance component.
# α_relative_threshold_default = 0.0
# Δc_partition_radius_default = 0.3 # determines how many resonances get grouped together. Larger number means less groups and thus more resonances per group.
λ0 = 3.4

# Δr_default = 1.0 # the samples used to build the surrogate is taken every `Δr` radian on the frequency axis. Decrease for improved accuracy at the expense of computation resources.
# Δκ_λ_default = 0.05 # the samples used to build thes urrogate for κ_λ are taken at this sampling spacing. Decrease for improved accuracy at the expense of computation resources.
# Δcs_max_scalar_default = 0.2 # In units of ppm. interpolation border that is added to the lowest and highest resonance frequency component of the mixture being simulated.
# κ_λ_lb_default = 0.5 # interpolation lower limit for κ_λ.
# κ_λ_ub_default = 2.5 # interpolation upper limit for κ_λ.

SH_config_path = "/home/roy/Documents/repo/NMRData/input/SH_configs/select_compounds_SH_configs_reduce.json"
#SH_config_path = "/home/roy/Documents/repo/NMRData/input/SH_configs/select_compounds_SH_configs_low_intensity_threshold.json"
surrogate_config_path = "/home/roy/Documents/repo/NMRData/input/surrogate_configs/select_compounds_SH_configs.json"

#molecule_names = ["alpha-D-Glucose"; ]; skip_sum_flag = false;
molecule_names = ["beta-D-Glucose"; ]; skip_sum_flag = false;
#molecule_names = ["L-Serine"; ]; skip_sum_flag = true;

#molecule_names = ["L-Isoleucine"; ]; skip_sum_flag = true;
#molecule_names = ["Epinephrine"; ]; skip_sum_flag = true;
#molecule_names = ["L-Leucine"; ]; skip_sum_flag = true;
#molecule_names = ["L-Lysine"; ]; skip_sum_flag = false;
#molecule_names = ["L-Phenylalanine"; ]; skip_sum_flag = true;
#molecule_names = ["L-Glutamine"; ]; skip_sum_flag = false;


#molecule_names = ["L-Valine"; ]; skip_sum_flag = true;
#molecule_names = ["Ethanol"; ]; skip_sum_flag = true;
#molecule_names = ["DSS"; ]; skip_sum_flag = true;
#molecule_names = ["ATP"; ] ; skip_sum_flag = true;
#molecule_names = ["L-Glutathione oxidized"; ]; skip_sum_flag = false;
#molecule_names = ["L-Glutathione reduced"; ]; skip_sum_flag = false;
#molecule_names = ["beta-Alanine"; ] ; skip_sum_flag = true;
#molecule_names = ["L-Alanine"; ] ; skip_sum_flag = true;

### serine.

#skip_sum_flag = false

#fs, SW, ν_0ppm = fetchsamplemachinesettings("400")
# fs, SW, ν_0ppm = fetchsamplemachinesettings("500")
 fs, SW, ν_0ppm = fetchsamplemachinesettings("600")
# fs, SW, ν_0ppm = fetchsamplemachinesettings("700")
# fs, SW, ν_0ppm = fetchsamplemachinesettings("900")

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
    prune_combo_Δc_bar_flag = true)
As = mixture_params

##prunecombocoherences!(As[1], α_relative_threshold_default, tol_coherence_default, Δc_partition_radius_default)
#prunecombocoherencesbar!(As[1], α_relative_threshold_default, tol_coherence_default, Δc_partition_radius_default)

### plot.

function evalzerophaseclpartitionelement(r,
    α::T, Ω::T, λ::T)::Complex{T} where T <: Real

    out = α/(λ+im*(r-Ω))

    return out
end

function evalzerophaseclresonancegroup(u_rad, A,
    i::Int, k::Int, λ::T)::Complex{T} where T <: Real

    out = zero(Complex{T})
    for l in A.part_inds_compound[i][k]
        out += evalzerophaseclpartitionelement(u_rad, A.αs[i][l], A.Ωs[i][l], λ)
    end

    return out
end

function evalzerophaseclsinglets(u_rad, A, λ::T)::Complex{T} where T <: Real

    out = zero(Complex{T})
    for i = 1:length(A.Ωs_singlets)
        out += evalzerophaseclpartitionelement(u_rad, A.αs_singlets[i], A.Ωs_singlets[i], λ)
    end

    return out
end

function evalzerophaseclmixture(u_rad, As, λ::T)::Complex{T} where T <: Real

    out = zero(Complex{T})
    for n = 1:length(As)
        for i = 1:length(As[n].part_inds_compound)
            for k = 1:length(As[n].part_inds_compound[i])

                out += evalzerophaseclresonancegroup(u_rad, As[n], i, k, λ)
            end
        end

        out += evalzerophaseclsinglets(u_rad, As[n], λ)
    end

    return out
end


# This is the frequency range that we shall work with.
P = LinRange(3, 4, 50000)
U = ppm2hzfunc.(P)
U_rad = U .* (2*π)

# absorption Lorentzian.
q = uu->evalzerophaseclmixture(uu, As, λ0)
q_singlets = uu->evalzerophaseclsinglets(uu, A, λ0)

function getqs(A, λ::T) where T <: Real
    #
    N_sys = length(A.part_inds_compound)

    qs = Vector{Vector{Function}}(undef, N_sys)
    for i = 1:N_sys

        N_groups = length(A.part_inds_compound[i])
        qs[i] = collect( uu->evalzerophaseclresonancegroup(uu, A, i, k, λ) for k = 1:N_groups)
    end

    return qs
end

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

save_fig_flag = true

using Plots; plotly()
plots_save_path = "/home/roy/Documents/workspaces/WIP/debug"

title_string = "groups"
canvas_size = (1500, 600)
plot_obj, q_U, qs_U, q_singlets_U = plotgroups(title_string, P, U, q, qs, q_singlets, real, 1.0;
canvas_size = canvas_size)

if save_fig_flag
    Plots.savefig(plot_obj, plots_save_path)
end
display(plot_obj)

@assert 1==2


## parameters that affect qs.
# A.d, A.κs_λ, A.κs_β
# A.d_singlets, A.αs_singlets, A.Ωs_singlets, A.β_singlets, A.λ0, A.κs_λ_singlets
q = uu->NMRSignalSimulator.evalclproxymixture(uu, mixture_params, Bs)

# create the functions for each resonance group.
B = Bs[1]
A = As[1]
qs = collect( collect( ωω->B.qs[i][k](ωω-B.ss_params.d[i], B.ss_params.κs_λ[i]) for k = 1:length(B.qs[i]) ) for i = 1:length(B.qs) )
q_singlets = ωω->NMRSignalSimulator.evalclsinglets(ωω, B.d_singlets, A.αs_singlets, A.Ωs_singlets, B.β_singlets, B.λ0, B.κs_λ_singlets)

# create the function for the entire compound.
q = uu->NMRSignalSimulator.evalclproxymixture(uu, As[1:1], Bs[1:1])

# evaluate at the plotting positions.
q_U = q.(U_rad)

qs_U = collect( collect( qs[i][k].(U_rad) for k = 1:length(qs[i]) ) for i = 1:length(qs) )
q_singlets_U = q_singlets.(U_rad)



function plotgroups(P, q_U, qs_U, q_singlets_U, fig_num)

    PyPlot.figure(fig_num)
    fig_num += 1

    PyPlot.plot(P, q_U, label = "sum", "--")

    for i = 1:length(qs_U)
        for k = 1:length(qs_U[i])

            PyPlot.plot(P, qs_U[i][k], label = "group ($(i),$(k))")
        end
    end

    PyPlot.plot(P, q_singlets_U, label = "singlets")

    return fig_num
end


fig_num = plotgroups(P, q_U, qs_U, q_singlets_U, fig_num)

#### sanity check.
q_check_U = q_singlets_U
if !isempty(qs) # some compounds only have singlets.
    q_check_U += sum( sum( qs[i][k].(U_rad) for k = 1:length(qs[i]) ) for i = 1:length(qs) )
end

discrepancy = norm(q_check_U- q_U)
println("sanity check. This should be numerically zero: ", discrepancy)
