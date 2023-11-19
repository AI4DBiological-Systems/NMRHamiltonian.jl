# Demo: code walk-through
For those who prefer working with source files instead of the Julia REPL and notebooks, this demo follows the script `examples/simulate.jl`, which is derived from `examples/simulate_SH.jl` in the NMRHamiltonian repository. In this demo, we load the physical chemistry parameters and alias mapping files, then simulate the resonance frequency and intensity, followed by an application of the convex clustering to partition the resonance components into resonance groups.

First, load the libraries. You can use `using` (loads NMRHamiltonian into current namespace) or `import` (loads NMRHamiltonian into its own namespace `NMRHamiltonian`). We use the latter in this demo.

The code involving DataDeps, Tar, CodecZlib were developed on `DataDeps v0.7.11`, `Tar v1.10.0`, `CodecZlib v0.7.3`. You might need to modify it if you use a different version of `DataDeps`.

```julia
import NMRHamiltonian
HAM = NMRHamiltonian # so we don't need to type NMRHamiltonian. everytime we reference something from this library.


using DataDeps, Tar, CodecZlib

using LinearAlgebra
import PythonPlot as PLT # do Pkg.add("PythonPlot") if you're missing this library.

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
end

function extractuncompress(src_path, dest_path)
    tar_gz = open(src_path)
    tar = GzipDecompressorStream(tar_gz)
    dir = Tar.extract(tar, dest_path)
    close(tar)
    
    return dir
end

# This was tested in a REPL environment, not a notebook environment. You might need to modify this for notebooks to get it to display the plots.
PLT.close("all")
fig_num = 1

# select the AbstractFloat data type you want to use. Anything lower than Float32 will result in numerical precision-related errors in our tests.
#T = Float64
T = Float32
```

The variable `molecule_entries` contain compound aliases. For this demo, any alias in the `molecule_name_mapping/select_molecules.json` file from [10.5281/zenodo.8174261](https://zenodo.org/record/8174261) would work.
```julia
### user inputs.

molecule_entries = [
    "L-Serine";
    "alpha-D-Glucose";
    "beta-D-Glucose";
    "Ethanol";
    "L-Methionine";     
    "L-Phenylalanine";
    "L-Glutathione reduced";
    "L-Glutathione oxidized";       
    "Uridine";
    "L-Glutamine";
    "L-Valine";
    "DSS";
]

# machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
fs, SW, ν_0ppm = HAM.getpresetspectrometer(T, "700")

## pull the sample coupling values into dictionary structures.
# use DataDeps.jl and Tar.jl to download and extract the sample coupling values.
root_data_path = getdatapath() # coupling values data repository root path

H_params_path = joinpath(root_data_path, "coupling_info") # folder of coupling values. # replace with your own values in actual usage.


molecule_mapping_root_path = joinpath(
    root_data_path,
    "molecule_name_mapping",
)
molecule_mapping_file_path = joinpath(
    molecule_mapping_root_path,
    "select_molecules.json",
)
#molecule_mapping_file_path = joinpath(molecule_mapping_root_path, "GISSMO_names.json")


# Type ?help HAM.SHConfig for details.
config_SH = HAM.SHConfig{T}(
    coherence_tol = convert(T, 0.01),
    relative_α_threshold = convert(T, 0.005),
    tol_radius_1D = convert(T, 0.1), # strictly between 0 and 1. The lower, the better the approximation, but would a larger partition (i.e. more resonance groups).
    nuc_factor = convert(T, 1.5),
)

```

Now, load the parameters, simulate the resonance components, and partition into resonance groups.
```julia
# read the physical parameters from file, from via the DataDeps routine from "helpers/data.jl".
Phys = NMRHamiltonian.getphysicalparameters(
    T,
    molecule_entries,
    H_params_path,
    molecule_mapping_file_path;
    unique_cs_atol = convert(T, 1e-6),
)
# You can tweak the J-coupling or chemical shift values of Phys manually if desired, before calling simulate(). We won't tweak anything in this demo.

# Assemble and eigendecompose the Hamiltonian, get frequencies and intensities for each resonance component, and partition the comonents into resonance groups.
As, MSPs = HAM.simulate(
    Phys,
    molecule_entries,
    fs,
    SW,
    ν_0ppm,
    config_SH,
)
```

Now, we'll plot the absorption Lorentzian spectrum of the simulated resonance frequencies and intensities. We'll have to make up a T2 relaxation value, which corresponds to the component width.

We will focus in on the first molecule (change `molecule_select = 1` to select a different molecule, the molecule names are in `molecule_entries`) and its first spin system (change `spin_system_select = 1` to select a different spin system).
If `T` was set to a floating-point data type other than the default one from your operating system's architecture, then you need to explicitly type-cast each floating-point number to type `T`, via `convert`. To find out what the default floating-point data type is, run `typeof(1.0)`.

```julia
### visualize a target molecule and spin group, given a T-2 decay parameter.
λ0 = convert(T, 3.4) # make up a T-2 decay parameter for this plotting example.
molecule_select = 1
spin_system_select = 1
ppm_offset = convert(T, 0.2) # for display padding.
N_viz = 50000 # this many points to plot



println("Visualizing $(molecule_entries[molecule_select]) and spin system $(spin_system_select).")

a = As[molecule_select].αs[spin_system_select]
F = As[molecule_select].Ωs[spin_system_select]

hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs # Hz to ppm conversion. This is useful to remember.
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW) # ppm to Hz conversion. This is useful to remember.

ΩS_ppm = hz2ppmfunc.( F ./ HAM.twopi(T) ) # convert to radial frequency.
ΩS_ppm_sorted = sort(ΩS_ppm)

u_min = ppm2hzfunc(ΩS_ppm_sorted[1] - ppm_offset)
u_max = ppm2hzfunc(ΩS_ppm_sorted[end] + ppm_offset)

P_min = hz2ppmfunc(u_min)
P_max = hz2ppmfunc(u_max)

P = LinRange(P_min, P_max, N_viz)
U = ppm2hzfunc.(P)
U_rad = U .* HAM.twopi(T)

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
PLT.figure(fig_num)
fig_num += 1

PLT.plot(P, real.(q_U))

PLT.gca().invert_xaxis()
PLT.ylabel("real part")
PLT.title("spectrum of $(molecule_entries[molecule_select]), spin system $(spin_system_select)")
```
