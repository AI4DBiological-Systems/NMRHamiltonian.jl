### search parameter configuration types and usage methods.

const isnestedsuccessively = ConvexClustering.isnestedsuccessively
const isnestedin = ConvexClustering.isnestedin

# would nee connectivity config to instantiate the cluster graph.
struct PartitionSearchRecord{T}
    
    # γ-traversal of nested partition tree.
    Gs::Vector{Vector{Vector{Int}}}
    γs::Vector{T}
    chosen_ind_γ::Int # the index selected for the surrogate model, out of Gs and γs.

    # affects the constructed nested partition tree.
    θs::Vector{T}
    #connectivity::AbstractSearchConnectivityType{T} # contains knn or radius info.
    connectivity_label::Symbol
    connectivity_parameter::Symbol

    # the convex clustering problem, given θ.
    #edge_pairs::Vector{Tuple{Int64, Int64}}
    #w::Vector{Float64}
end


# function extractGs(Rs::Vector{Vector{PartitionSearchRecord{T}}}) where T

#     return collect(
#         collect(
#             Rs[n][i].Gs for i in eachindex(Rs[n])
#         ) for n in eachindex(Rs)
#     )
# end

# function extractallγs(Rs::Vector{Vector{PartitionSearchRecord{T}}}) where T

#     return collect(
#         collect(
#             Rs[n][i].γs for i in eachindex(Rs[n])
#         ) for n in eachindex(Rs)
#     )
# end

function extractθγ(Rs::Vector{Vector{PartitionSearchRecord{T}}}) where T

    γs = collect(
        collect(
            Rs[n][i].γs[Rs[n][i].chosen_ind_γ] for i in eachindex(Rs[n])
        ) for n in eachindex(Rs)
    )

    θs = collect(
        collect(
            Rs[n][i].θs[end] for i in eachindex(Rs[n])
        ) for n in eachindex(Rs)
    )

    return θs, γs
end

# function extractγs(Rs_set::PartitionSearchRecord{T})
    

# end
# for a single spin system, unlike PartitionsParametersType is for multiple spin systems.
# it is intended that end users don't save this to file (unless using BSON.jl). They write source code to modify this.
# it is intended that end users do save PartitionsParametersType (hyperparameters for multiple spin systems) to file.
struct PartitionAlgorithmType

    getgraphconfigfunc::Function # (name,i,cm,as)-> configs for graph, optim, assignment, θ, γ.
    getoptimconfigfunc::Function
    getassignmentconfigfunc::Function
    getsearchθconfigfunc::Function
    getsearchγconfigfunc::Function

    store_trace::Bool
    report_cost::Bool
    verbose_kernel::Bool

    getdefaultγfunc::Function # a function: (name::String, i::Int, Δc_m::Vector{Vector{T}}, αs::Vector{T}) ↦ γ::T
    getdefaultθfunc::Function # a function: (name::String, i::Int, Δc_m::Vector{Vector{T}}, αs::Vector{T}) ↦ θ::KT
end

function returnonefunc(
    name::String,
    i::Int,
    Δc_m::Vector{Vector{T}},
    αs_i_prune::Vector{T},
    )::T where T

    return one(T)
end

function PartitionAlgorithmType(
    getgraphconfigfunc::Function, # pre-made options are: defaultknnsearchconfig, other options: defaultknnconfig, defaultradiusconfig, defaultradiussearchconfig
    getoptimconfigfunc::Function,
    getsearchθconfigfunc::Function, # options are  createsearchθconfigs or disablesearch()
    getsearchγconfigfunc::Function; # options are createsearchγconfigs or disablesearch()

    getassignmentconfigfunc::Function = defeaultassignmentconfig,

    store_trace::Bool = false,
    report_cost::Bool = false,
    verbose_kernel::Bool = false,

    getdefaultγfunc::Function = returnonefunc,
    getdefaultθfunc::Function = returnonefunc)

    #
    return PartitionAlgorithmType(
        getgraphconfigfunc,
        getoptimconfigfunc,
        getassignmentconfigfunc,
        getsearchθconfigfunc,
        getsearchγconfigfunc,

        store_trace,
        report_cost,
        verbose_kernel,

        getdefaultγfunc,
        getdefaultθfunc
        )
end

# universal disable function that returns nothing instead of the intended default return type.
function disablesearch(args...)
    return nothing
end

### front end

"""
```
generatemixturepartitionalgorithm(
    names::Vector{String},
    θs::Vector{Vector{KT}},
    γs::Vector{Vector{T}},
    Phys::Vector{PhysicalParamsType{T}},
    getgraphconfigfunc::Function,
    getsearchθconfigfunc::Function,
    getsearchγconfigfunc::Function,
    runcgsolverlib::Function;
    report_cost::Bool = false,
    verbose_kernel::Bool = false
    ) where {T <: AbstractFloat, KT}
```

Assembles a partition parameters configuration data structure using information from a file. The configuration is for the `simulate()` function.

Only the targeted names and field names are updated from file. The non-targeted names and non-targeted field names are filled using the supplied values to `loadmixtureshconfig`.

See `PartitionsParametersType` for the field names.

### Inputs

- `names::Vector{String}`
List of molecule entries.

- `θs::Vector{Vector{KT}}`
With an abuse of notation, each element has the same function as the `θs` field from the `PartitionsParametersType` data type.

- `γs::Vector{Vector{T}}`
With an abuse of notation, each element has the same function as the `γs` field from the `PartitionsParametersType` data type.


### Optional inputs

- `Phys`$(DOCSTRING_Phys("T"))

- `getgraphconfigfunc`$(DOCSTRING_getgraphconfigfunc())

- `getsearchθconfigfunc`$(DOCSTRING_getsearchθconfigfunc())

- `getsearchγconfigfunc`$(DOCSTRING_getsearchγconfigfunc())

- `report_cost::Bool`
Output to REPL the optimization cost and the sparsity-inducing regularization value γ that was used in the partition problem set up.

- `verbose_kernel::Bool`
Output to REPL diagnostics about the kernel parameter for the weight function in the partition problem set up.

### Outputs

The only output is of type `::Vector{Vector{PartitionAlgorithmType{T,KT}}}`, the list of partition configuration settings for `simulate()`.

"""
function generatemixturepartitionalgorithm(
    names::Vector{String},
    θs::Vector{Vector{KT}},
    γs::Vector{Vector{T}},
    Phys::Vector{PhysicalParamsType{T}};
    getgraphconfigfunc::Function = defaultknnconfig, # options: defaultknnconfig, defaultradiusconfig, defaultknnsearchconfig, defaultradiussearchconfig

    getsearchθconfigfunc::Function = disablesearch, # options:  createsearchθconfigs, disablesearch
    getsearchγconfigfunc::Function = disablesearch, # options: createsearchγconfigs, disablesearch
    getassignmentfunc::Function = defeaultassignmentconfig,

    #store_trace::Bool = false,
    report_cost::Bool = false,
    verbose_kernel::Bool = false,
    
    # for optimization.
    gap_tol::T = convert(T, 1e-8),
    max_iters::Int = 300,
    ) where {T <: AbstractFloat, KT}

    @assert length(γs) == length(θs) == length(Phys) == length(names)

    # getoptimconfigfunc = (nn,ii,cc,aa)->defeaultoptimconfig(nn,ii,cc,aa,runcgsolverlib)
    getoptimconfigfunc = (nn,ii,cc,aa)->defeaultoptimconfig(
        nn, ii, cc, aa;
        gap_tol = gap_tol,
        max_iters = max_iters,
    )

    out = Vector{Vector{PartitionAlgorithmType}}(undef, length(names))

    for n in eachindex(Phys)
        p = Phys[n]

        out[n] = Vector{PartitionAlgorithmType}(undef, length(p.cs_sys))

        for i in eachindex(p.cs_sys)

            θ_default = θs[n][i]
            γ_default = γs[n][i]

            # won't get used if getsearchγconfigfunc = disablesearch()
            getdefaultγfunc = (nn,ii,cc,aa)->γ_default

            # won't get used if getsearchθconfigfunc = disablesearch()
            getdefaultθfunc = (nn,ii,cc,aa)->θ_default

            # use this to search both θ then γ, ignoring both γ_default and θ_default.
            out[n][i] = PartitionAlgorithmType(
                getgraphconfigfunc,
                getoptimconfigfunc,
                getsearchθconfigfunc,
                getsearchγconfigfunc;

                getdefaultγfunc = getdefaultγfunc,
                getdefaultθfunc = getdefaultθfunc,
                getassignmentconfigfunc = getassignmentfunc,

                #store_trace = store_trace,
                report_cost = report_cost,
                verbose_kernel = verbose_kernel)
        end
    end

    return out
end

"""
```
setupconstantparameteroptions(
    molecule_entries::Vector{String},
    Ps::Vector{PartitionsParametersType{T,KT}};
    metric::Distances.Metric = Distances.Euclidean(),
    kernelfunc::Function = evalSqExpkernel,
    ) where {T,KT}
```

Assembles constant partition parameters configuration generating functions using information given a list of `PartitionsParametersType`. These configuration generating functions can be used by the `simulate()` function's `getgraphconfigfunc` optional input argument.

The `knns` and `radi` fields of the list of `PartitionsParametersType` Julia variables are used to create configuration generating functions that force the partition algorithm in `simulate()` to use those values.

### Inputs

- `molecule_entries::Vector{String}`
List of molecule entries.

- `Ps::Vector{PartitionsParametersType{T,KT}}`
List of partition parameters data structure. Aside from manually constructing this, one can use `defaultmixturepartitionsparameters()` to construct them with default values, or `loadmixturepartitionsparameters()` to construct them with default values and values from a file.


### Optional inputs

- `report_cost::Bool`
Output to REPL the optimization cost and the sparsity-inducing regularization value γ that was used in the partition problem set up.

- `verbose_kernel::Bool`
Output to REPL diagnostics about the kernel parameter for the weight function in the partition problem set up.

### Outputs

- `constantknnconfig::Function`
Use this with the  `getgraphconfigfunc` optional input argument of `generatemixturepartitionalgorithm()` to tell the partition algorithm to use the `knns` values in `Ps` to set up the convex clustering problem in the partition algorithm.

- `constantradiusconfig::Function`
Use this with the  `getgraphconfigfunc` optional input argument of `generatemixturepartitionalgorithm()` to tell the partition algorithm to use the `radi` values in `Ps` to set up the convex clustering problem in the partition algorithm.

- `mixture_θs::Vector{Vector{KT}}`
The kernel hyperparameters in `Ps`, in a nested `Vector` form. This can be used as one of the inputs of `generatemixturepartitionalgorithm()`.

- `mixture_γs::Vector{Vector{T}}`
The regularization parameters in `Ps`, in a nested `Vector` form. This can be used as one of the inputs of `generatemixturepartitionalgorithm()`.

"""
function setupconstantparameteroptions(
    molecule_entries::Vector{String},
    Ps::Vector{PartitionsParametersType{T,KT}}
    ) where {T,KT}

    N = length(Ps)
    #cpknngraphconfigfuncs = Vector{Vector{Function}}(undef, N)
    #cpradiusgraphconfigfuncs = Vector{Vector{Function}}(undef, N)
    mixture_θs = Vector{Vector{KT}}(undef, N)
    mixture_γs = Vector{Vector{T}}(undef, N)

    for n in eachindex(Ps)

        N_sys = length(Ps[n].knns)
        @assert length(Ps[n].radi) == length(Ps[n].θs) == length(Ps[n].γs) == N_sys

        #cpknngraphconfigfuncs[n] = Vector{Function}(undef, N_sys)
        #cpradiusgraphconfigfuncs[n] = Vector{Function}(undef, N_sys)
        mixture_θs[n] = Vector{KT}(undef, N_sys)
        mixture_γs[n] = Vector{T}(undef, N_sys)

        # if N_sys < 1
        #     # if the n-th molecule entry does not have any non-singlet spin systems, allocate a dummy entry anyways for code compatibility.
        #     cpknngraphconfigfuncs[n] = Vector{Function}(undef, 1)
        #     cpradiusgraphconfigfuncs[n] = Vector{Function}(undef, 1)
        #     mixture_θs[n] = Vector{KT}(undef, 1)
        #     mixture_γs[n] = Vector{T}(undef, 1)
        #
        #     # this initalization should never be used, so set to values that will cause an error.
        #     cpknngraphconfigfuncs[n][begin] = xx->nothing
        #     cpradiusgraphconfigfuncs[n][begin] = xx->nothing
        #     mixture_θs[n][begin] = -1.0
        #     mixture_γs[n][begin] = -1.0
        # end

        # populate container.
        for i in eachindex(Ps[n].knns)
            # cpknngraphconfigfuncs[n][i] = (nn,ii,cc,aa)->getgraphonfigsupplyknn(
            #     nn,ii,cc,aa;
            #     knn = Ps[n].knns[i],
            #     metric = metric,
            #     kernelfunc = kernelfunc)
            #
            # cpradiusgraphconfigfuncs[n][i] = (nn,ii,cc,aa)->getgraphconfigsupplyradius(
            #     nn,ii,cc,aa;
            #     manual_radius = Ps[n].radi[i],
            #     metric = metric,
            #     kernelfunc = kernelfunc)

            mixture_θs[n][i] = Ps[n].θs[i]
            mixture_γs[n][i] = Ps[n].γs[i]
        end

    end

    constantknnconfig = (nn,ii,cc,aa)->getgraphonfigsupplyknn(nn,ii,cc,aa, Ps, molecule_entries)
    constantradiusconfig = (nn,ii,cc,aa)->getgraphconfigsupplyradius(nn,ii,cc,aa, Ps, molecule_entries)

    #return cpknngraphconfigfuncs, cpradiusgraphconfigfuncs, mixture_θs, mixture_γs
    return constantknnconfig, constantradiusconfig, mixture_θs, mixture_γs
end


### algorithm

function defeaultassignmentconfig(name::String, i::Int, c::Vector{Vector{T}},
    a::Vector{T};
    metric::Distances.Metric = Distances.Euclidean(),
    assignment_zero_tol::T = convert(T, 1e-3)) where T <: AbstractFloat

    return ConvexClustering.AssignmentConfigType(metric, assignment_zero_tol)
end

function defeaultoptimconfig(
    name::String,
    i::Int,
    c::Vector{Vector{T}},
    a::Vector{T};
    # runcgsolverlib::Function = (xx,ff,dff,gg_tol)->ConvexClustering.runOptimjl(xx, ff, dff, gg_tol),
    gap_tol::T = convert(T, 1e-8),
    σ_base::T = convert(T, 0.4),
    σ_rate::T = convert(T, 1.05),
    updateσfunc::Function = nn::Int->convert(T, σ_base*σ_rate^(nn-1)),
    updateϵfunc::Function = nn::Int->convert(T, 1/(nn-1)^2),
    #runoptimfunc::Function = (xx,ff,dff,gg_tol)->ConvexClustering.runOptimjl(xx, ff, dff, gg_tol; verbose = verbose_subproblem),
    max_iters::Int = 300,
    ) where T <: AbstractFloat

    ## anonmymous functions do not take keyword inputs. cannot control verbose option of runcgsolver from here.
    #runoptimfunc = (xx,ff,dff,gg_tol)->runcgsolverlib(xx, ff, dff, gg_tol; verbose = verbose_subproblem),
    #runoptimfunc = runcgsolverlib

    return ConvexClustering.ALMConfigType(
        gap_tol;
        #runoptimfunc;
        max_iters = max_iters,
        updateσfunc = updateσfunc,
        updateϵfunc = updateϵfunc,
        )
end


### hyperparameters

# weight function hyperparameter search.

"""
```
function createsearchθconfigs(
    name,
    i,
    Δc_m::Vector{Vector{T}},
    αs::Vector{T};
    length_scale_base::T = convert(T, 10.0),
    length_scale_rate::T = convert(T, 0.7),
    min_dynamic_range::T = convert(T, 0.95),
    length_scale_max_iters::Int = 1000,
    getθfunc::Function = nn::Int->convert(T, lengthscale2θ(evalgeometricsequence(nn-1, length_scale_base, length_scale_rate))),
    ) where T <: AbstractFloat
```
Returns `NMRHamiltonian.ConvexClustering.SearchθConfigType(length_scale_max_iters, min_dynamic_range, getθfunc)`.

The default input parameters are for use with the `evalSqExpkernel()` kernel function.

`min_dynamic_range` should take a value in the range of the kernel function, which in the default case of `evalSqExpkernel()` is between 0 and 1.

"""
function createsearchθconfigs(
    name,
    i,
    Δc_m::Vector{Vector{T}},
    αs::Vector{T};
    length_scale_base::T = convert(T, 10.0),
    length_scale_rate::T = convert(T, 0.7),
    min_dynamic_range::T = convert(T, 0.95),
    length_scale_max_iters::Int = 1000,
    getθfunc::Function = nn::Int->convert(T, lengthscale2θ(evalgeometricsequence(nn-1, length_scale_base, length_scale_rate))),
    ) where T <: AbstractFloat

    return ConvexClustering.SearchθConfigType(length_scale_max_iters, min_dynamic_range, getθfunc)
end

"""
```
createsearchγconfigs(
    name,
    i,
    Δc_m::Vector{Vector{T}},
    αs::Vector{T};
    γ_base::T = convert(T, 0.1),
    γ_rate::T = convert(T, 1.05),
    max_iter_γ::Int = 100,
    getγfunc::Function = nn::Int->convert(T, evalgeometricsequence(nn-1, γ_base, γ_rate)),
    max_partition_size_offset::Int = 2, # set to large number to force max_partition_size to be length(Δc_m), i.e., every part is a singleton, i.e. every resonance component is in its own resonance group.
    partition_rate = 2.0, # only used when max_partition_size_offset is 0. Above 1 means more than the dimension of a Δc vector, below 1 means less than the dimension.
    verbose = false,
    )where T <: AbstractFloat
```
- Returns `NMRHamiltonian.ConvexClustering.SearchγConfigType(max_iter_γ, max_partition_size, getγfunc)`.

- The default input parameters are for use with the `evalSqExpkernel()` kernel function.

- `γ_base` is the starting sparsity-inducing regularization parameter in the search sequence. The overall trend is that larger γ produces a smaller partition size, but this trend might not be true for γ values that close to each other.

- `max_partition_size` is computed as:
```
D = length(first(Δc_m))
max_partition_size = clamp( D + max_partition_size_offset, 1, length(Δc_m) )
```

"""
function createsearchγconfigs(
    name,
    i,
    Δc_m::Vector{Vector{T}},
    αs::Vector{T};
    γ_base::T = convert(T, 0.1),
    γ_rate::T = convert(T, 1.05),
    max_iter_γ::Int = 100,
    getγfunc::Function = nn::Int->convert(T, evalgeometricsequence(nn-1, γ_base, γ_rate)),
    max_partition_size_offset::Integer = 3, # set to large number to force max_partition_size to be length(Δc_m), i.e., every part is a singleton, i.e. every resonance component is in its own resonance group.
    partition_rate::Real = 2.0, # only used when max_partition_size_offset is 0. Above 1 means more than the dimension of a Δc vector, below 1 means less than the dimension.
    # verbose::Bool = false,
    )where T <: AbstractFloat

    # stop searching once the size of the returned partition is less than `max_partition_size`.
    @assert !isempty(Δc_m)
    D = length(first(Δc_m))

    # assume manual as default.
    max_partition_size::Int = clamp( D + max_partition_size_offset, 1, length(Δc_m) ) # TODO make sure this guard is set for all max_partition_size variables in the package.

    # it not using manual.
    if max_partition_size_offset == 0

        proposed_size = round(Int, D*partition_rate)
        max_partition_size = clamp( proposed_size, 1, length(Δc_m) )
    end

    return ConvexClustering.SearchγConfigType(max_iter_γ, max_partition_size, getγfunc)
end

# build graph using knn, ignore radius.
# """
# usage: a way to create an anonmymous function for the `getgraphconfigfunc` input argument for `generatemixturepartitionalgorithm()`.

# Given `knn::Int`, the anomymous function is created by:
# getgraphconfigfunc = (nn,ii,cc,aa)->getgraphonfigsupplyknn(nn,ii,cc,aa; knn = knn).
# """
# function getgraphonfigsupplyknn(
#     name,
#     i,
#     Δc_m::Vector{Vector{T}},
#     αs::Vector{T};
#     knn::Int = 64,
#     metric::Distances.Metric = Distances.Euclidean(),
#     kernelfunc::Function = evalSqExpkernel) where T

#     connectivity = ConvexClustering.KNNType(knn)

#     return CC.WeightedGraphConfigType(connectivity, metric, kernelfunc)
# end
function getgraphonfigsupplyknn(
    name,
    i,
    Δc_m::Vector{Vector{T}},
    αs::Vector{T},
    Ps::Vector{PartitionsParametersType{T,KT}},
    molecule_entries::Vector{String};
    fallback_knn::Int = 64,
    metric::Distances.Metric = Distances.Euclidean(),
    kernelfunc::Function = evalSqExpkernel,
    ) where {T,KT}

    n = findfirst(xx->xx==name, molecule_entries)

    # use fallback if problem with input.
    knn = fallback_knn
    if !isnothing(n)
        if isassigned(Ps[n].knns, i)
            knn = Ps[n].knns[i]
        else
            println("There is an issue with the partition configuration. Cannot find spin index $(i) for entry $(name). Using Fallback knn: $(fallback_knn).")
        end
    else
        println("There is an issue with the partition configuration. Cannot find $(name) in provided molecule entries. Using Fallback knn: $(fallback_knn).")
    end

    # make sure knn is valid.
    N = length(Δc_m)
    knn = clamp(knn, 1, N-1)

    connectivity = ConvexClustering.KNNType(knn)

    return CC.WeightedGraphConfigType(connectivity, metric, kernelfunc)
end

# """
# ```
# defaultknnconfig(
#     name,
#     i,
#     Δc_m::Vector{Vector{T}},
#     αs::Vector{T};
#     metric::Distances.Metric = Distances.Euclidean(),
#     kernelfunc::Function = evalSqExpkernel
#     ) where T <: AbstractFloat
# ```

# This is an example of a valid custom function

# Usage:
# Then pass this function to `generatemixturepartitionalgorithm()`:
# getgraphconfigfunc = defaultknnconfig

# """
# function getgraphconfigsupplyradius(name, i,
#     Δc_m::Vector{Vector{T}}, αs::Vector{T};
#     manual_radius::T = one(T),
#     metric::Distances.Metric = Distances.Euclidean(),
#     kernelfunc::Function = evalSqExpkernel) where T <: AbstractFloat

#     connectivity = ConvexClustering.RadiusType(manual_radius)

#     return CC.WeightedGraphConfigType(connectivity, metric, kernelfunc)
# end
function getgraphconfigsupplyradius(
    name,
    i,
    Δc_m::Vector{Vector{T}},
    αs::Vector{T},
    Ps::Vector{PartitionsParametersType{T,KT}},
    molecule_entries::Vector{String};
    fallback_radius::T = one(T),
    metric::Distances.Metric = Distances.Euclidean(),
    kernelfunc::Function = evalSqExpkernel) where {T <: AbstractFloat, KT}

    n = findfirst(xx->xx==name, molecule_entries)

    # use fallback if problem with input.
    radius = fallback_radius
    if !isnothing(n)
        if isassigned(Ps[n].radi, i)
            radius = Ps[n].radi[i]
        else
            println("There is an issue with the partition configuration. Cannot find spin index $(i) for entry $(name). Using Fallback radius: $(fallbacK_radius).")
        end
    else
        println("There is an issue with the partition configuration. Cannot find $(name) in provided molecule entries. Using Fallback radius: $(fallbacK_radius).")
    end

    # make sure the radius is non-negative.
    radius = max(radius, zero(T))

    connectivity = ConvexClustering.RadiusType(radius)

    return CC.WeightedGraphConfigType(connectivity, metric, kernelfunc)
end


### utility methods

"""
'''
evalgeometricsequence(n::Int, a0::T, r::T)::T where T <: AbstractFloat
'''
returns `a0*r^n`
"""
function evalgeometricsequence(n::Int, a0::T, r::T)::T where T <: AbstractFloat
    return a0*r^n
end

# non-negative-valued kernels, used as weight function.
"""
```
function evalSqExpkernel(x::Vector{T}, z::Vector{T}, θ::T)::T where T <: AbstractFloat
```
returns `exp(-θ*norm(x-z)^2)`
"""
function evalSqExpkernel(x::Vector{T}, z::Vector{T}, θ::T)::T where T <: AbstractFloat
    return exp(-θ*norm(x-z)^2)
end

"""
```
function lengthscale2θ(l::T)::T where T <: AbstractFloat
```
returns `1/(2*l^2)`
"""
function lengthscale2θ(l::T)::T where T <: AbstractFloat
    return 1/(2*l^2)
end

"""
```
function θ2lengthscale(θ::T)::T where T<: AbstractFloat
```
returns `1/sqrt(2*θ)`
"""
function θ2lengthscale(θ::T)::T where T<: AbstractFloat
    return 1/sqrt(2*θ)
end
