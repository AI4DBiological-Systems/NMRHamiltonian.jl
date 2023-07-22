module NMRHamiltonian


using LinearAlgebra
import Kronecker, Graphs
import JSON3

# using Pkg
# Pkg.add(url="https://github.com/RoyCCWang/ConvexClustering.jl")
import ConvexClustering
import Distances

const CC = ConvexClustering

include("./doc_strings/doc_types.jl")
include("./doc_strings/doc_configs.jl")


include("./types.jl")
include("./utils.jl")
include("./serialize.jl")

include("./settings/config_IO.jl")
include("./settings/part_configs.jl")
include("./settings/front_end.jl")
include("./settings/cvxc_config.jl")

include("./conversions/conversion_helpers.jl")
include("./conversions/IO_mag_eq.jl")
include("./conversions/IO_csJ.jl")
include("./conversions/cs_parse.jl")
include("./conversions/IO_phys.jl")

#include("../src/SH/SH_front_end.jl")
#include("../src/SH/molecule.jl")
include("./SH/SH.jl")
include("./SH/Hamiltonian.jl")
include("./SH/operators.jl")
include("./SH/singlets.jl")
include("./SH/front_end.jl") # front end for scripts in the SH folder.
include("./SH/prune_coherence.jl")
include("./SH/reduce_nuclei.jl")


include("./partition/partition.jl")
include("./partition/cvxc.jl")

# front end for auto-generating configs given hyperparameters.
include("./autogen.jl")

export 
getphysicalparameters,
PhysicalParamsType,
extractcs,

defaultmixtureshsconfig,
loadmixtureshsconfig,

defaultmixturepartitionsparameters,
loadmixturepartitionsparameters,

generatemixturepartitionalgorithm,
defaultknnconfig,
defaultknnsearchconfig,
defaultradiusconfig,
defaultradiussearchconfig,
createsearchθconfigs,
createsearchγconfigs,

setupconstantparameteroptions,

serializemixture,
simulate,
SHConfig,
SHType,

readbasechemshifts,
extractMEnuclei

end
