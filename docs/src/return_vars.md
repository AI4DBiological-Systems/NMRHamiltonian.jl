# Demo: return variables
Run the code on the `Demo: code walk-through` first.

We'll explore the various field names of the return variables for the simulation. We'll set look at the first compound (i.e., `n = 1`) and the first spin system (i.e., `i = 1`) for this demo. In Julia, set these variables:

```julia
n = 1
i = 1
```

### Resonance components result for the n-th compound entry, `As[n]`
First, we find out what data type the `As` variable is. 
```julia
typeof(As)
```
This should show up in the REPL.
```
Vector{SHType{Float64}}
```
Next, we display the list of field names for the variable `As[n]`, which corresponds to the 1-st compound entry.
```julia
ns = propertynames(As[n])
```
This should show up in the REPL.
```
(:αs, :Ωs, :Δc, :parts, :Δc_bar, :N_spins_sys, :fs, :SW, :ν_0ppm)
```
Let's see their respective data types:
```julia
collect( typeof( getfield(As[n], field_name) ) for field_name in ns )
```
This should show up in the REPL, and we added additional comments to show what the nested arrays index.
```
9-element Vector{DataType}:
 Vector{Vector{Float64}} # αs
 Vector{Vector{Float64}} # Ωs
 Vector{Vector{Vector{Float64}}} # Δc
 Vector{Vector{Vector{Int64}}} # part
 Vector{Vector{Vector{Float64}}} # Δc_bar
 Vector{Int64} # N_spins_sys
 Float64 # fs, sampling frequency (in Hz) of the spectrometer used for this simulation (i.e., it is an input to the simulation, returned for book keeping).
 Float64 # SW, the spectral window (in ppm) of the spectrometer used for this simulation.
 Float64 # ν_0ppm, the 0 ppm resonance frequency used for this simulation.
```

Some of the field names are described below:

- `αs[i][l]`: resonance intensities for the i-th spin system, l-th resonance component. The intensities are relative to the number of nuclei for the n-th compound entry.
- `Ωs[i][l]`: resonance frequencies for the i-th spin system, l-th resonance component. This is a radial frequency, which has units of radians. Divide by 2π to convert to Hz.

NMRHamiltonian scales `αs` for a given n-th entry such that we should recover approximately the number of nuclei if we sum over the relative intensities. For example: 
```julia
sum(As[1].αs[1])
```
gives
```
4.994790368138611
```
This number approaches 5 as `relative_α_threshold` is set to approach zero, the number of nuclei for the `n == 1` compound entry of this demo, which is `L-Glutamine`. However, setting `relative_α_threshold = 0.0` incurs significant computational burden, as every resonance component is kept, even if it has a tiny intensity. In practice, we recommend leaving `relative_α_threshold == 0.01` or some other small number, but keep in mind that the sum of the relative intensities would not be exactly equal to the number of nuclei used to assemble the Hamiltonian.

- `Δc[i][l][j]`: order of coherence for the i-th spin system, l-th resonance component, and j-th ME-modulo nuclei index. The set of 1-D arrays `Δc[i]` is the input set we pass to convex clustering to construct a partition. You can use your own clustering/partition construction algorithm; just overwrite `Δc_bar[i]` and `part[i]` accordingly with your results. 


- `part`: a part is an element of a partition in combinatorics literature terminology. `part[i][l]` contains the resonance component indices (the `l` index from `Δc`) that belong in the part for the i-th spin system, j-th part. In NMRHamiltonian, a part of a spin system is called a resonance group, as it is a set of resonance components that share a similar degrees of freedom.

- `Δc_bar[i][k][j]`: representative order of coherence for the i-th spin system, k-th resonance group, and j-th ME-modulo nuclei index. The number of `k` index values is the number of parts that convex clustering returned.

- `N_spins_sys[i]`: the number of nuclei in the i-th spin system. These are not ME_modulo nuclei. This is used for normalization purposes for finalizing `αs`.


### Spin system information for the n-th compound entry, `MSPs[n]`

First, we find out what data type the `MSPs` variable is. 
```julia
typeof(MSPs)
```
This should show up in the REPL.
```
Vector{MoleculeSpinSystem{Float64}}
```
Next, we display the list of field names for the variable `Rs[n][i]`, which corresponds to the convex clustering results for the first compound entry, first spin system.
```julia
ns = propertynames(MSPs[n])
```
This should show up in the REPL.
```
(:spin_systems, :singlet_intensities, :singlet_frequencies)
```

- `spin_systems`: Discussed shortly.

- `singlet_intensities[m]`: The m-th singlet intensity of the n-th compound entry.

- `singlet_frequencies[m]`: The m-th singlet (radial) frequency (in radians) of the n-th compound entry.

### Hamiltonian for the n-th compound entry, i-th spin system, `MSPs[n].spin_system[i]`

First, we find out what data type the `MSPs[n].spin_systems[i]` variable is. 
```julia
typeof(MSPs[n].spin_systems[i])
```
This should show up in the REPL.
```
NMRHamiltonian.SpinSystem{Float64}
```
Next, we display the list of field names for the variable `MSPs[n].spin_systems[i]`, which corresponds to the convex clustering results for the first compound entry, first spin system.
```julia
ns = propertynames(MSPs[n].spin_systems[i])
```
This should show up in the REPL.
```
(:intensities, :frequencies, :H, :coherence_mat, :coherence_state_pairs, :states, :partial_quantum_numbers, :quantum_numbers, :coherence_tol)
```
Let's see their respective data types:
```julia
collect( typeof( getfield(MSPs[n].spin_systems[i], field_name) ) for field_name in ns )
```
This should show up in the REPL, and we added additional comments to show what the nested arrays index.
```
9-element Vector{DataType}:
 Vector{Float64} (alias for Array{Float64, 1}) # intensities
 Vector{Float64} (alias for Array{Float64, 1}) # frequencies
 NMRHamiltonian.Hamiltonian{Float64} # H
 Matrix{Float64} (alias for Array{Float64, 2}) # coherence_mat
 Vector{Tuple{Int64, Int64}} (alias for Array{Tuple{Int64, Int64}, 1}) # coherence_state_pairs
 Vector{Int64} (alias for Array{Int64, 1}) # states
 Vector{Vector{Float64}} (alias for Array{Array{Float64, 1}, 1}) # partial_quantum_numbers
 Vector{Float64} (alias for Array{Float64, 1}) # quantum_numbers
 Float64 # coherence_tol
```

- `intensities`: These are all the resonance intensities that passed the (-1)-quantum coherence condition (within a numerical tolerance of `coherence_tol`). The `αs` intensities were created from intensities.

- `frequencies`: Similar to `intensities`, but for resonance (radial) frequencies, in radians.

- `H`: information for the Hamiltonian matrix. The `H.matrix` is the Hamiltonian matrix. `H.contributions` is a 1-D array of two matrices that sum to `H.matrix`, each of the matrices corresponds to a term in the Hamiltonian equation in our manuscript. The `H.eigenvalues` and `H.eigenvectors` are the eigenpairs of `H.matrix`.

- `coherence_mat`: this is the order of coherence as computed using the z-axis spin angular momentum.

- `coherence_state_pairs`: this is the eigenstate pairs (sometimes denoted by `{(r,s)}` in our manuscript) that corresponds to each resonance component in `intensities` and `frequencies`.

- `states`: the unique states in `coherence_state_pairs`.

- `partial_quantum_numbers`: the contributions to the quantum number of each state in `states`. One contribution from each quantum sub-system, where each sub-system is a ME-modulo nuclei index. These generate the order of coherence features `Δc`.

- `quantum_numbers`: the quantum numbers of an eigenstate, used to determine if a state-pair from `coherence_state_pairs` satisfies the (-1)-quantum coherence condition (within a numerical tolerance of `coherence_tol`). If it does, the state-pair creates a valid resonance component, and its intensity and frequency is computed and stored as an element of `intensities` and `frequencies`, respectively.

### Physical chemistry parameters for the n-th compound entry, `Phys[n]`

First, we find out what data type the `Phys[n]` variable is. 
```julia
typeof(Phys[n])
```
This should show up in the REPL.
```
NMRHamiltonian.PhysicalParamsType{Float64}
```
Next, we display the list of field names for the variable `Phys[n]`, which corresponds to the convex clustering results for the first compound entry, first spin system.
```julia
ns = propertynames(Phys[n])
```
This should show up in the REPL.
```
(:H_IDs, :H_inds_sys, :cs_sys, :H_inds_singlets, :cs_singlets, :J_inds_sys, :J_inds_sys_local, :J_vals_sys, :ME)
```
Let's see their respective data types:
```julia
collect( typeof( getfield(Phys[n], field_name) ) for field_name in ns )
```
This should show up in the REPL, and we added additional comments to show what the nested arrays index.
```
9-element Vector{DataType}:
 Vector{Int64} # H_IDs
 Vector{Vector{Int64}} # H_inds_sys
 Vector{Vector{Float64}} # cs_sys
 Vector{Vector{Int64}} # H_inds_singlets
 Vector{Float64} # cs_singlets
 Vector{Vector{Tuple{Int64, Int64}}} # J_inds_sys
 Vector{Vector{Tuple{Int64, Int64}}} # J_inds_sys_local
 Vector{Vector{Float64}} # J_vals_sys
 Vector{Vector{Vector{Int64}}} # ME
```

- `H_IDs`: The nuclei labels parsed from the physical chemistry file for the n-th compound entry.

- `H_inds_sys[i]`: The nuclei from a re-labelled `H_IDs` that are in the i-th spin system. The re-label starts from 1, and increments by 1.

- `cs_sys[i]`: The chemical shifts (in units ppm) that corresponds to the nuclei in `H_inds_sys[i]`, for the i-th spin system.

- `H_inds_singlets`: The nuclei from a re-labelled `H_IDs` that generate singlet spin systems.

- `cs_singlets`: The chemical shifts (in units ppm) that corresponds to the nuclei in `H_inds_singlets`.

- `J_inds_sys[i]`: Pairs of nuclei from `H_inds_sys[i]` such that each pair has a J-coupling value for the i-th spin system.

- `J_inds_sys_local[i]`: This is a version of `J_inds_sys[i]` but the nuclei indices are re-labelled to start from 1, i.e. this is a re-labelling of the nuclei for the i-th spin system.

- `J_vals_sys`: J-coupling values in Hz, corresponds to the nuclei in `J_inds_sys[i]` and `J_inds_sys_local[i]`.

- `ME[i][m]`: The nuclei indices from a re-labelled `H_IDs` (so that it starts from 1, increments by 1) for the m-th non-singleton ME-modulo nuclei index in the i-th spin system. This means `ME[i][m]` contains some nuclei labels that together form the m-th set of nuclei such that it contains only nuclei that are magnetically equivalent to each other, and this m-th set of nuclei is the largest possible set. Non-singleton means that `ME[i][m]` does not include nucleus that are not magnetically equivalent to anyone else; any nucleus is magnetically equivalent to itself, and we exclude it to be in `ME[i][m]`. For example, if `Phys[n].ME[i]` contains 
```
2-element Vector{Vector{Int64}}:
 [7, 5, 6] # assign to j == 1
 [1, 2, 3] # assign to j == 2.
 # The other values of j's are assigned according to ascending numerical value of whatever nuclei labels are left in this spin system.
```
then it means the n-th compound entry, i-th spin system has two non-singleton magnetically equivalent set of nuclei. The first set of nuclei is nuclei `{7,5,6}` and the second set is nuclei `{1,2,3}`. If there are 8 nuclei for this compound entry, then nuclei 4 and 8 are both not magnetically equivalent to any other nuclei that isn't itself. The ME-module index `j` for this spin system would be `{1,2,3,4}`, with `j == 1` representing nuclei `{7,5,6}`, `j == 2` representing nuclei `{1,2,3}`, `j == 3` representing nucleus `4`, `j == 4` representing nucleus `8`.

### Store simulated variables

JSON: in this demo, we use the [JSON3.jl](https://github.com/quinnj/JSON3.jl) library for loading JSON. Please make sure it is installed before running the following code.
```
# store.
file_name = "test.json"
S = NMRHamiltonian.serializemixture(As)
NMRHamiltonian.saveasJSON(file_name, S)

# load.
import JSON 3
file_path = file_name
json_string = read(file_path)
W = JSON3.read(json_string)
```

You can also use [BSON.jl](https://github.com/JuliaIO/BSON.jl).
```
# store.
file_name = "test.bson"
S = NMRHamiltonian.serializemixture(As)
BSON.bson(file_name, S)

# load.
W = BSON.load(file_name)
```

A round-trip test routine is in the `test` folder of the NMRHamiltonian repository. See the functions `roundtripJSON()` and `roundtripBSON()` in that folder for more details on the test.

An alternative method to serialize data is with the Julia Base `Serialization` library. My understanding is that compatibility is not guaranteed between Julia versions, so it should be used for short-term storage.