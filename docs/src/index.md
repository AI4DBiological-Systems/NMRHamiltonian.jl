# NMRHamiltonian.jl
NMRHamiltonian is a library for simulating the resonance frequencies and relative intensities of compounds in isotropic liquid. It is intended for 1D 1H NMR, which is a spin-1/2 NMR. 

## Table of contents
```@contents
Pages = ["index.md"]
Depth = 5
```

## Install
NMRHamiltonian is hosted on a custom Julia registry. We need to add that before installing NMRHamiltonian: 
``` julia
using Pkg
Pkg.Registry.add(RegistrySpec(url = "https://github.com/AI4DBiological-Systems/PublicJuliaRegistry"))
Pkg.add("NMRHamiltonian")
```

To update this package once it is installed, do
``` julia
using Pkg
Pkg.update("NMRHamiltonian")
```

# Important exported functions
The following functions are the focus of NMRHamiltonian:

* `getphysicalparameters()` given a list of compound aliases, parse the physical chemistry parameters from file to data structure.

* `SHConfig()` assembles the simulation configuration container variable. This is a constructor to a composite data type (i.e., this data type is like the C programming language's `struct`).

* `simulate()` the outputs of `getphysicalparameters()` and `SHConfig()`, and returns a data structure that contain the simulated frequencies and relative intensities for each resonance comonents, and their resonance groups.

This is not a full quantum simulation of the NMR spectrum; only the frequencies, intensities, and this partition of the resonance components called a set of *resonance groups* is returned.

See the demo pages for an example walk-through of how to use this library.

## Julia Basics

- Read the [Julia documentation](https://docs.julialang.org/en/v1/stdlib/Pkg/) for installing Julia packages.

- We use 1-indexing for addressing arrays here, but you can use `[begin + n]` instead of `[n]` if you wish to use a 0-indexing index `n`. 

# A compiled list of physical chemistry parameters
The `getphysicalparameters()` function in NMRHamiltonian requires the user to supply it with a list of compound aliases, each of which is a string called an `entry`, as well as the corresponding file path of each entry to a JSON file that stores its physical chemistry parameters. These parameters are the chemical shift and J-coupling values of an compound, as we are only concerned with spin-1/2 NMR in isotropic liquid in NMRHamiltonian. As some of these values change as a function of ambient temperature, pH, spectrometer frequency, one can have multiple entries that describe the same compound for different experimental or spectrometer configurations.

The demo and test scripts for NMRHamiltonian.jl requires the compiled list of parameters at [10.5281/zenodo.8174261](https://zenodo.org/record/8174261), but you are encouraged to devise your own that follows its JSON format. You'll need a alias-to-file-path mapping file, which is the `molecule_name_mapping/selected_molecules.json` file from this dataset. You'll also need the individual parameter files for each compound at a certain experimental or spectrometer configurations, which are the files in the `coupling_info` folder from this dataset. See these files to get a sense of how the format is like.

For the parameter files, the JSON file structure goes like this:
```json
{
       "J-coupling": [
                       {
                            "ID2": 9,
                            "ID1": 8,
                          "value": -12.531894
                       },
                       {
                            "ID2": 10,
                            "ID1": 8,
                          "value": 3.511369
                       },
                       {
                            "ID2": 10,
                            "ID1": 9,
                          "value": 5.979
                       }
                     ],
   "chemical shift": [
                       {
                          "value": 3.9744,
                             "ID": 8
                       },
                       {
                          "value": 3.9352,
                             "ID": 9
                       },
                       {
                          "value": 3.83016,
                             "ID": 10
                       }
                     ]
}
```
This is the `L-Serine_mod.json` file in `coupling_info` from the dataset. This indicates that there are J-couplings (nucleus 1 ID, nucleus 2 ID, J-coupling value in Hz) = (8, 9, -12.531894), (8, 10, 3.511369), (9, 10, 5.979). There are chemical shifts for each nuclei: (8, 3.9744 ppm), (9, 3.9352 ppm), (10, 3.83016 ppm).

You can construct ficticious compound entries and simulate it, just to see what the resonance frequencies and intensities might look like.

For the mapping file, it follows the format:
```json
{
    "Dopamine - 500 MHz": {
        "notes": "http://gissmo.bmrb.io/entry/bmse000909/simulation_1",
        "file name": "bmse000909_simulation_1.json"
    },
    "L-Cystine": {
        "notes": "http://gissmo.bmrb.io/entry/bmse000035/simulation_1",
        "file name": "bmse000035_simulation_1.json"
    },
}
```
It means that the parameter file for compound alias "Dopamine - 500 MHz" is `bmse000909_simulation_1.json` in the *nuclei parameter file path folder*, which together with the file path for the mapping file, are inputs to the `getphysicalparameters()` function in this library. We will go over this function in the demo.

## Citation
Please cite [10.5281/zenodo.8174261](https://zenodo.org/record/8174261) and the appropriate references within if you use any physical chemistry values from it. The references to the original data sources from which these values were curated from are also in the instructions.

## Download a NMR parameters dataset via DataDeps.jl

Run the following commands to load the `getdatapath()` function into Julia REPL, which calls [DataDeps.jl](https://github.com/oxinabox/DataDeps.jl) to download and uses [Tar.jl](https://github.com/JuliaIO/Tar.jl) and [CodecZlib](https://github.com/JuliaIO/CodecZlib.jl) to extract the files on your local machine. If you usually use a source file workflow, these commands are stored in `examples/helpers/data.jl` in the NMRHamiltonian repository.

```julia
using DataDeps, Tar, CodecZlib


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
```

Now, run this function to download and extract the data.
```julia
getdatapath()
```

# Terminology

- A *ME-modulo nuclei index* corresponds to a maximal set of nuclei where each nucleus is magnetically equivalent (ME) with all other nucleus in that maximal set. For example: if compound X has 10 nuclei contributing to its 1D 1H NMR spectrum, but nucleus #2 and #5 are ME, the other nucleus do not form ME with each other, then we can re-label the nuclei indices `i = 1, 2, ..., 10` to the ME-modulo nuclei index `j = 1, 2, ..., 9`, with `j = 2` denoting the set of nuclei `{i == 2, i == 5}`. The order of coherences-related quantities in NMRHamiltonian use ME-modulo nuclei indices as oppose to the nuclei indices.

- A *radial* frequency in units of radians is a frequency in Hz multiplied by 2π.

- A *resonance group* is a set of resonance components that have similar degrees-of-freedom. The set of resonance groups of a spin system is a partition of the resonance components for that spin system.

- A *part* is an element of a partition of a set `X`. This means a part is a subset of `X`.

- Throughout this documentation, we use `x == y` to mean the value in variable `x` equals to the value in variable `y`. We use `x = y` to mean we assign the variable `x` the value  `y`.


# License
NMRHamiltonian has the Mozilla Public License Version 2.0.

# Authors
Code author:
- Roy Chih Chung Wang

Supervisors:
- Dave Campbell (Carleton University, Bank of Canada)
- Miroslava Čuperlović-Culf (National Research Council of Canada)

# Funding
This projected was funded by the AI-for-Design Challenge Program from the National Research Council of Canada.