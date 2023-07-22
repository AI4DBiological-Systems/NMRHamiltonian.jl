# include("../src/NMRHamiltonian.jl")
# import .NMRHamiltonian
# # remove later. dependencies.
# import Kronecker, Graphs, ConvexClustering, Distances
# # end remove.

import NMRHamiltonian

using DataDeps, CodecZlib, Tar

using LinearAlgebra
import PyPlot
import JSON3

import Optim


PyPlot.close("all")
fig_num = 1

PyPlot.matplotlib["rcParams"][:update](["font.size" => 22, "font.family" => "serif"])

### user inputs.

molecule_entries = ["D-(+)-Glucose"; "Ethanol"; "L-Methionine";
    "L-Phenylalanine"; "L-Glutathione reduced"; "L-Glutathione oxidized";
    "L-Leucine"; "L-Isoleucine"; "L-Valine";
    "HEPES"; "L-Glutamine"; "DSS"; "D2O"]

#molecule_entries = ["L-Methionine"; "L-Phenylalanine"; "DSS"; "Ethanol"; "L-Isoleucine"]
#molecule_entries = ["D-(+)-Glucose"; "Ethanol"; "DSS"; "D2O"; ]

#molecule_entries = ["D2O, 4.8ppm"]
#molecule_entries = ["DSS"]
#molecule_entries = ["L-Histidine - 2000 study";]
#molecule_entries = ["L-Histidine, 500 MHz";]
#molecule_entries = ["L-Histidine";]
#molecule_entries = ["DL-Serine"]
#molecule_entries = ["Gamma-Aminobutyric acid"]
#molecule_entries = ["3-Hydroxybutyrate"]
#molecule_entries = ["L-Glutathione reduced"]
#molecule_entries = ["L-Glutathione oxidized"]
#molecule_entries = ["Uridine"]
#molecule_entries = ["HEPES";]
#molecule_entries = ["L-Leucine";]
#molecule_entries = ["L-Isoleucine";]
#molecule_entries = ["Ethanol";]
#molecule_entries = ["L-Valine";]
#molecule_entries = ["L-Glutamine";]
#molecule_entries = ["L-Phenylalanine";]

# machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
fs = 14005.602240896402
SW = 20.0041938620844
ν_0ppm = 10656.011933076665


# # machine values for the BMRB 500 MHz glucose experiment.
# ν_0ppm = 6752.490995937095
# SW = 16.0196917451925
# fs = 9615.38461538462

## pull the sample coupling values into dictionary structures.
include("./helpers/data.jl") # use DataDeps.jl and Tar.jl to download and extract the sample coupling values.

root_data_path = getdatapath() # coupling values data repository root path

H_params_path = joinpath(root_data_path, "coupling_info") # folder of coupling values. # replace with your own values in actual usage.

molecule_mapping_root_path = joinpath(root_data_path, "molecule_name_mapping")
molecule_mapping_file_path = joinpath(molecule_mapping_root_path, "select_molecules.json")
#molecule_mapping_file_path = joinpath(molecule_mapping_root_path, "GISSMO_names.json")


### end inputs.

println("Timing: getphysicalparameters")
@time Phys = NMRHamiltonian.getphysicalparameters(
    Float64,
    molecule_entries,
    H_params_path,
    molecule_mapping_file_path;
    unique_cs_atol = 1e-6)

# see if any molecules have ME in multiple spin systems.
MEs = collect( Phys[n].ME for n in eachindex(Phys) )

# ethanol.
n = 2
i = 1

# Leucine
n = 7
i = 1

ME_m = MEs[n][i]
N_spins = length(Phys[n].cs_sys[i])

ordering, DOF = NMRHamiltonian.createorderingfromeqinds(ME_m, N_spins)
@assert length(ordering) == N_spins

c = randn(length(ordering))
println("before, c: ", c)

cp = NMRHamiltonian.condensenuclei(c, ordering, DOF)
@assert length(cp) == length(unique(ordering))


# reverse.
function applyreverseME(
    y::Vector{T},
    ordering::Vector{Int},
    M::Int,
    )::Vector{T} where T

    N = length(y)
    #@assert length(y) == length(unique(ordering))
    @assert norm(collect(1:N) - sort(unique(ordering))) < 1e-14

    x = zeros(T, M)
    #fill!(x, Inf)
    for i in eachindex(x)

        k = ordering[i]
        x[i] = y[k]
    end

    return x
end

# not really used in practice.
function applyforwardME(
    x::Vector{T},
    ordering::Vector{Int},
    N::Int
    )::Vector{T} where T

    @assert length(x) == length(ordering)
    @assert norm(collect(1:N) - sort(unique(ordering))) < 1e-14

    y = zeros(T, N)
    #fill!(y, Inf)
    for i in eachindex(x)

        k = ordering[i]

        ## if no averaging, just overwrite.
        y[k] = x[i] # add the coherence contributions among the equivalent resonances.

        ## if averaging.
        #occurances = length(findall(xx->xx==k, ordering))
        #y[k] = y[k]/occurances
    end

    return y
end

x = applyreverseME(cp, ordering, N_spins)

xf = applyforwardME(x, ordering, DOF)
xr = applyreverseME(xf, ordering, N_spins)

@assert norm(xr - x) < 1e-14

#@assert 1==2


###### magnetic equivalence in multiple spin systems in the same molecule.
molecule_mapping_file_path = joinpath(molecule_mapping_root_path, "GISSMO_names.json")
dict_molecule_to_filename = JSON3.read(read(molecule_mapping_file_path)) # map molecule entries to coupling 

db_names = String.(keys(dict_molecule_to_filename))
molecule_entries = db_names

println("Timing: getphysicalparameters")
@time Phys = NMRHamiltonian.getphysicalparameters(
    Float64,
    molecule_entries,
    H_params_path,
    molecule_mapping_file_path;
    unique_cs_atol = 1e-6)

# see if any molecules have ME in multiple spin systems.
MEs = collect( Phys[n].ME for n in eachindex(Phys) )
lens = collect( collect( !isempty(MEs[n][i]) for i in eachindex(MEs[n]) ) for n in eachindex(MEs) )
Q = count.(lens)

# get entries X that have magnetic equivalence in more than one spin system.
inds = findall(xx->xx>1, Q)
X = molecule_entries[inds]

p0 = Phys[inds[17]]
me0 = MEs[inds[17]]
p0.cs_sys

load_paths = collect( 
    joinpath(
        H_params_path,
        dict_molecule_to_filename[Symbol(X[n])][Symbol("file name")],
     ) for n in eachindex(X)
)

@time Phys = NMRHamiltonian.getphysicalparameters(
    Float64,
    load_paths;
    unique_cs_atol = 1e-6,
)

p = Phys[1]
me = MEs[1]
p.cs_sys