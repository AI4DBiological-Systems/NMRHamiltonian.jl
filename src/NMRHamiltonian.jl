module NMRHamiltonian

# Write your package code here.
using Distributed, LinearAlgebra, FFTW, SharedArrays, OffsetArrays, Statistics

import Kronecker, Interpolations, NearestNeighbors, JSON, Graphs

include("../src/types.jl")
include("../src/utils.jl")

include("./conversions/conversion_helpers.jl")
include("./conversions/IO_mag_eq.jl")
include("./conversions/IO_csJ.jl")
include("./conversions/cs_parse.jl")

#include("../src/SH/SH_front_end.jl")
#include("../src/SH/molecule.jl")
include("./SH/SH.jl")
include("./SH/Hamiltonian.jl")
include("./SH/operators.jl")
include("./SH/resonance_partitioning.jl")
include("./SH/prune_coherence.jl")
include("./SH/reduce_nuclei.jl")

end
