# Demo: code walk-through
For those who prefer working with source files instead of the Julia REPL and notebooks, this demo follows the script `examples/demo.jl`, which is derived from `examples/simulate_SH.jl` in the NMRHamiltonian repository. In this demo, we load the physical chemistry parameters and alias mapping files, then simulate the resonance frequency and intensity, followed by an application of the convex clustering to partition the resonance components into resonance groups.

First, load the libraries. You can use `using` (loads NMRHamiltonian into current namespace) or `import` (loads NMRHamiltonian into its own namespace `NMRHamiltonian`). We use the latter in this demo.

```julia
import NMRHamiltonian
HAM = NMRHamiltonian # so we don't need to type NMRHamiltonian. everytime we reference something from this library.


using DataDeps, Tar, CodecZlib

using LinearAlgebra
import PyPlot

# actual values from Bruker machines. Curated from BMRB experiments.
function fetchsamplemachinesettings(tag)

    fs = 14005.602240896402
    SW = 20.0041938620844
    ν_0ppm = 10656.011933076665

    if tag == "700"
        # machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
        fs = 14005.602240896402 # in Hz.
        SW = 20.0041938620844 # spectral window, in ppm.
        ν_0ppm = 10656.011933076665 # in Hz, where the zero ppm is detected in the experiment. Not used for this demo.

    elseif tag == "600"
        ## machine values from a 600 MHz experiment: bmse000915, methionine.
        fs = 9615.38461538462
        SW = 16.022093454391
        ν_0ppm = 6685.791496181313

    elseif tag == "900"
        ## machine values from a 900 MHz experiment: GISSMO, leucine entry.
        fs = 14423.0769230769
        SW = 16.0300195009073
        ν_0ppm = 10160.027322585376

    elseif tag == "500"
        fs = 6493.50649350649
        SW = 12.9911090156122
        ν_0ppm = 4035.6644246816795

    elseif tag == "400"

        ### 400 MHz, bmse000297, ethanol.
        fs = 4807.69230769231
        SW = 12.0152693165838
        ν_0ppm = 2884.905244600881
    end

    return fs, SW, ν_0ppm
end

# the same function for ensuring the dataset at 10.5281/zenodo.8174261 is downloaded.
function getdatapath()::String

    dataset_alias = "AI4DBiological-Systems_NMR_data" # don't use spaces or 'strange' symbols like commas, colons, etc.
    archive_file_name = "nmr_physical_parameters_dataset.tar.gz" # the filename on the data repository that we download.
    url = "https://zenodo.org/record/8174261/files/nmr_physical_parameters_dataset.tar.gz?download=1"

    register(DataDep("$dataset_alias",
        """
        Dataset: Selected BMRB 1D 1H NMR data and physical chemistry values compiled from literature
        Author: Roy Chih Chung Wang
        License: [Creative Commons Attribution Non Commercial Share Alike 4.0 International](https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode)
        DOI: 10.5281/zenodo.8174261

        Description:
        The experiments were downloaded from BMRB,
        The physical chemistry parameters data and sample configuration-related files were collected, organized, converted or manually entered in as JSON files by Roy.
        The physical chemistry parameters, i.e., chemical shift and J-coupling values, are from two sources:
        - GISSMO: (https://gissmo.bmrb.io/ acccessed Mar. 2023). See Dashti, et. al. DOI: 10.1021/acs.analchem.8b02660 and DOI: 10.1021/acs.analchem.7b02884 for more details.
        - From Govindaraju, et. al.'s work: DOI: 10.1002/1099-1492(200005)13:3<129::AID-NBM619>3.0.CO;2-V

        Please cite the data sources and this data repository if you find the contents helpful for your work. See the Zenodo DOI entry for more description.
        """,
        url
    ));

    #readdir(datadep"AI4DBiological-Systems NMR data") # have to manually type out the alias. Does not allow string variable substitution.
    local_dataset_archive_path = @datadep_str("$dataset_alias") # call the actual macro to allow string variable substitution.

    # extract archive, then delete. Do this only if archive file still exists.

    root_data_path = joinpath(local_dataset_archive_path, "contents")
    
    if isfile(joinpath(local_dataset_archive_path, archive_file_name))
        t = @task begin; ispath(root_data_path) || mkpath(root_data_path); end
        schedule(t); wait(t)
    
        #t = @task begin; Tar.extract(joinpath(local_dataset_archive_path, archive_file_name), root_data_path); end
        t = @task begin; extractuncompress(joinpath(local_dataset_archive_path, archive_file_name), root_data_path); end
        schedule(t); wait(t)
        rm(joinpath(local_dataset_archive_path, archive_file_name)) # delete the archive file.
    end

    return root_data_path

    # # return root_data_path. however, this unpacks in the current working directory!
    # archive_file_path = joinpath(local_dataset_archive_path, archive_file_name)
    # if isfile(archive_file_path)
    #     DataDeps.unpack(archive_file_path)
    # end
    #return local_dataset_archive_path
end

function extractuncompress(src_path, dest_path)
    tar_gz = open(src_path)
    tar = GzipDecompressorStream(tar_gz)
    dir = Tar.extract(tar, dest_path)
    close(tar)
    
    return dir
end

# This was tested in a REPL environment, not a notebook environment. You might need to modify this for notebooks to get it to display the plots.
PyPlot.close("all")
fig_num = 1
PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

T = Float64
```

The variable `molecule_entries` contain compound aliases. For this demo, any alias in the `molecule_name_mapping/select_molecules.json` file from [10.5281/zenodo.8174261](https://zenodo.org/record/8174261) would work.
```julia
### user inputs.

molecule_entries = ["L-Glutamine"; "L-Valine"; "L-Phenylalanine"; "DSS"]

# machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
fs, SW, ν_0ppm = fetchsamplemachinesettings("700")

## pull the sample coupling values into dictionary structures.
# use DataDeps.jl and Tar.jl to download and extract the sample coupling values.
root_data_path = getdatapath() # coupling values data repository root path

H_params_path = joinpath(root_data_path, "coupling_info") # folder of coupling values. # replace with your own values in actual usage.

molecule_mapping_root_path = joinpath(root_data_path, "molecule_name_mapping")
molecule_mapping_file_path = joinpath(molecule_mapping_root_path, "select_molecules.json")
#molecule_mapping_file_path = joinpath(molecule_mapping_root_path, "GISSMO_names.json")


# # These produce a aprtition tree over the searched γ values for many molecule entries.
config_SH = HAM.SHConfig(
    T;
    coherence_tol = 0.01, # set to a small number. must be larger than zero. 0.01 seemed to work across our test simulation configurations.
    relative_α_threshold = 0.01, # the larger it is, the more resonance components would be dropped for computation speed. However, setting it to 0 would keep every resonance component even if it had a very small intensity. 0.01 seemed to work across our test simulation configurations.
    max_partition_size_offset = 0,
    partition_rate = 2.0,
    γ_base = 0.01,
    γ_rate = 1.05,
    max_iter_γ = 100,
    fully_connected_convex_clustering = false, #  overides all knn-related optional inputs
    length_scale_base = 10.0,
    length_scale_rate = 0.7,
    min_dynamic_range = 0.95,
    cc_gap_tol = 1e-8, # set this much smaller than assignment_zero_tol to get partition tree, this is the stopping gradient residual of the convex clustering algorithm.
    cc_max_iters = 300,
    assignment_zero_tol = 1e-3,
    report_γ_cost = true,
    verbose_kernel = true,
)
```

Now, load the parameters, simulate the resonance components, and partition into resonance groups.
```julia
# read the physical parameters from file, from via the DataDeps routine from "helpers/data.jl".
Phys = HAM.getphysicalparameters(
    Float64,
    molecule_entries,
    H_params_path,
    molecule_mapping_file_path;
    unique_cs_atol = 1e-6,
)
# You can tweak the J-coupling or chemical shift values of Phys manually if desired, before calling simulate(). We won't tweak anything in this demo.

# Assemble and eigendecompose the Hamiltonian, get frequencies and intensities for each resonance component, and partition the comonents into resonance groups.
As, Rs, MSPs = HAM.simulate(
    Phys,
    molecule_entries,
    fs,
    SW,
    ν_0ppm,
    config_SH,
)
```

Now, we'll plot the absorption Lorentzian spectrum of the simulated resonance frequencies and intensities. We'll have to make up a T2 relaxation value, which corresponds to the component width.
```julia
### visualize a target molecule and spin group, given a T-2 decay parameter.
λ0 = 3.4 # make up a T-2 decay parameter for this plotting example.
molecule_select = 1
spin_system_select = 1
ppm_offset = 0.2 # for display padding.
N_viz = 50000 # this many points to plot

a = As[molecule_select].αs[spin_system_select]
F = As[molecule_select].Ωs[spin_system_select]

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs # Hz to ppm conversion. This is useful to remember.
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW) # ppm to Hz conversion. This is useful to remember.

ΩS_ppm = hz2ppmfunc.( F ./ (2*π) ) # convert to radial frequency.
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
PyPlot.figure(fig_num)
fig_num += 1

PyPlot.plot(P, real.(q_U))
PyPlot.gca().invert_xaxis()

PyPlot.ylabel("real part")
#PyPlot.legend()
PyPlot.title("spectrum of $(molecule_entries[molecule_select]), spin system $(spin_system_select)")
```
