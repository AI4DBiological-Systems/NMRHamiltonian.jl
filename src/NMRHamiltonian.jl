module NMRHamiltonian


using LinearAlgebra
import Kronecker, Graphs
import JSON3

import SingleLinkagePartitions
const SL = SingleLinkagePartitions

include("./doc_strings/doc_types.jl")
include("./doc_strings/doc_configs.jl")


include("./types.jl")
include("./utils.jl")
include("./serialize.jl")

# data conversion, disk input/output.
include("./IO/config_IO.jl")
include("./IO/conversion_helpers.jl")
include("./IO/IO_mag_eq.jl")
include("./IO/IO_csJ.jl")
include("./IO/cs_parse.jl")
include("./IO/IO_phys.jl")

# spin Hamiltonian.
include("./SH/SH.jl")
include("./SH/Hamiltonian.jl")
include("./SH/operators.jl")
include("./SH/singlets.jl")


# overall simulator.
include("./engine/reduce_nuclei.jl")
include("./engine/partition.jl")
include("./engine/front_end.jl") # front end for scripts in the SH folder.



export 
getphysicalparameters,
PhysicalParamsType,
extractcs,

serializephysicalparams,
deserializephysicalparams,
serializemixture,
deserializemixture,
saveasJSON,
readJSON,

simulate,
SHConfig,
SHType,

readbasechemshifts,
writebasechemshifts!,
extractMEnuclei

end
