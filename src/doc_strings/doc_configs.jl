


DOCSTRING_cs_sys_mixture(T, show_type = true) = """
$(show_type ? "`::Vector{Vector{Vector{$(T)}}}`" : "")


The set of chemical shift values for non-singlet spin systems for multiple molecule entries. See `PhysicalParamsType`.
"""

DOCSTRING_cs_singlets_mixture(T, show_type = true) = """
$(show_type ? "`::Vector{Vector{$(T)}}`" : "")


Singlet chemical shift values for multiple molecule entries. See `PhysicalParamsType`.
"""

DOCSTRING_Phys(T, show_type = true) = """
$(show_type ? "`::Vector{PhysicalParamsType{$(T)}}`" : "")


Physical parameters for multiple molecule entries.
"""



############## functions.

DOCSTRING_cc_intro = """
`NMRHamiltonian` poses the resonance group partitioning for each non-singlet spin system as a weighted convex clustering problem. The `ConvexClustering` package is used. Each resonance component of a spin system is a point, and we want to partition the set of all resonance components in a spin system. The resultant partition is the set of resonance groups.
"""

DOCSTRING_cc_connectivity = """
To reduce computational burden, it is customary to only allow a certain amount of non-zero weights in the set up of a weighted convex clustering problem, which is imposing a certain fixed connectivity between the points. The two connectivity construction methods `ConvexClustering` implements is construction by k-nearest neigbors, or construction by ϵ-radius ball between pairs of points.
"""


DOCSTRING_getgraphconfigfunc(show_type = true) = """
$(show_type ? "`::Function`" : "")

A function that specifies the connectivity construction configuration for the resonance group partition algorithm. This function has the form:
    
```
(n::String, i::Int, Δc_m::Vector{Vector{AbstractFloat}}, αs::Vector{AbstractFloat})->graph_config::NMRHamiltonian.CC.WeightedGraphConfigType
```
This function is provided so that one can create their own custom algorithm that utilize the molecule entry (`n`) spin system index (`i`), set of resonance components in the spin system (`Δc_m`), and set of resonance intensities for the components (`αs`) to generate a different connectivity configuration.

### Background
$(DOCSTRING_cc_intro)
$(DOCSTRING_cc_connectivity)

### Pre-made examples

See the source code of the following for an example of pre-made options.
-- For an iterative search over a range of graph connectivity parameters, see `defaultknnsearchconfig` or `defaultradiussearchconfig`. These start with values that return a large partition size, gradually decreases as the search goes on, and stops after the returned partition size is below a specified threshold size. You can copy and modify these functions with your own algorithm for the threshold size that uses the information of the manditory inputs of these functions.

-- For examples that are non-iterative searching algorithm of a graph connectivity parameter, see `defaultknnconfig` or `defaultradiusconfig`.
"""

####

DOCSTRING_getsearchθconfigfunc(show_type = true) = """
$(show_type ? "`::Function`" : "")

A function that specifies the search sequence and stopping condition for searching for a suitable hyperparameter θ of the kernel weight function at the beginning of the resonance group partitioning algorithm. No convex clustering algorithm is solved in this iterative process.
    
If searching is to be enabled, pass a function of the following form
```
(n::String, i::Int, Δc_m::Vector{Vector{AbstractFloat}}, αs::Vector{AbstractFloat})->graph_config::NMRHamiltonian.ConvexClustering.SearchθConfigType
```
A pre-made example is `createsearchθconfigs`, where the search terminantes after a minimum dynamic range is reached among the weights between the connected points in the set to be partitioned.

To disable iterative searching for θ, this function should be set to `disablesearch()`.

"""

DOCSTRING_getsearchγconfigfunc(show_type = true) = """
$(show_type ? "`::Function`" : "")

A function that specifies the search sequence and stopping condition for searching for a suitable sparsity-inducing regularization parameter γ by iteratively setting up and solving a sequence weighted convex clustering problems. This could be time consuming.
    
If searching is to be enabled, pass a function of the following form
```
(n::String, i::Int, Δc_m::Vector{Vector{AbstractFloat}}, αs::Vector{AbstractFloat})->graph_config::NMRHamiltonian.ConvexClustering.SearchγConfigType
```
A pre-made example is `createsearchγconfigs`, where the search terminantes after a solved partition has a size (the number of resonance groups) is lower than an adaptive threshold value `max_partition_size`.

To disable iterative searching for γ, this function should be set to `disablesearch()`.

"""