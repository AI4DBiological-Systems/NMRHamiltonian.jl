module NMRHamiltonian

using LinearAlgebra
import Kronecker, Graphs
import JSON3

import SingleLinkagePartitions
const SL = SingleLinkagePartitions

using Serialization
using Statistics

# constant values.
function twopi(::Type{Float32})::Float32
    return 6.2831855f0 #convert(T, 2*π)
end

function twopi(::Type{Float64})::Float64
    return 6.283185307179586 #convert(T, 2*π)
end

function twopi(::Type{T})::T where T <: AbstractFloat
    return convert(T, 2*π)
end


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


include("./IO/maximal_ME.jl") # generate Phys from an existing Phys, such that ME is maximized.

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
