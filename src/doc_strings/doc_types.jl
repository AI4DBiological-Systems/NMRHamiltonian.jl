################# SHType
DOCSTRING_SHType_αs(T, show_type = true) = """
$(show_type ? "`::Vector{Vector{$(T)}}`" : "")

Resonance intensities for each spin systems. Index `[i][l]` refers to the i-th spin system, l-th resonnace component.
"""

DOCSTRING_rad2hz = """
To convert to Hz, divide by `2π`.
"""

DOCSTRING_SHType_Ωs(T, show_type = true) = """
$(show_type ? "`::Vector{Vector{$(T)}}`" : "")

Resonance frequencies for each  spin systems. In units radians. $DOCSTRING_rad2hz Index `[i][l]` refers to the i-th spin system, l-th resonnace component.
"""

DOCSTRING_Δc = """
The Δc feature of a resonance component can be interpreted as the partial contributions to the quantum number of the resonance component from each subsystem of the the composite spin system. Its eleemnts take on values betwee `-1` and `1` for spin-1/2 NMR.
"""

DOCSTRING_SHType_Δc(T, show_type = true) = """
$(show_type ? "`::Vector{Vector{Vector{$(T)}}}`" : "")

Index `[i][l]` is a 1-D array of the Δc feature of the l-th resonance component in the i-th spin system. $DOCSTRING_Δc
"""

DOCSTRING_SHType_parts(show_type = true) = """
$(show_type ? "`::Vector{Vector{Vector{Int}}}`" : "")

For a fixed spin system `i`, `parts[i][k]` is a 1-D array of indices that specify which of the resonance components (e.g., which elements of αs[i] and Ωs[i]) belong in the resonance group (i,k). The following creates `α_ik`, the resonance components for the (i,k)-th resonance group of the molecule entry:
```
α_ik = αs[parts[i][k]]
```
The indices in `parts` follow the 1-indexing scheme.
"""


DOCSTRING_SHType_Δc_bar(T, show_type = true) = """
$(show_type ? "`::Vector{Vector{Vector{$(T)}}}`" : "")

Index `[i][k]` is the weighted average of the features of the (i,k)-th resonance group. The features are 
```
Δc_ik = Δc[parts[i][k]]
```
and the weights are 
```
αs[parts[i][k]]
```
The length of each `Δc_bar[i]` in a spin system should be the same, and should be less than or equal to the number of spins in the i-th spin system. It won't be equal if there are exhibits magnetic equivalence between the spins.
"""

DOCSTRING_SHType_N_spins_sys(show_type = true) = """
$(show_type ? "`::VectorInt}`" : "")

The i-th element is the number of spins in the i-th spin system.
"""


DOCSTRING_SHType_αs_singlets(T, show_type = true) = """
$(show_type ? "`::Vector{$(T)}`" : "")

The i-th element is the resonance intensity for the i-th singlet resonance component.
"""

DOCSTRING_SHType_Ωs_singlets(T, show_type = true) = """
$(show_type ? "`::Vector{$(T)}`" : "")

The i-th element is the resonance frequency for the i-th singlet resonance component. $DOCSTRING_rad2hz
"""

DOCSTRING_hz2ppm_ppm2hz = """
To convert from Hz to ppm frequency units, and the inverse, use the following anonymous functions:
```
hz2ppmfunc = uu->(uu - ν_0ppm)*SW/fs
ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)
```
"""

DOCSTRING_spectrometer_freq = """
Spectrometer frequency (in MHz) is calculated as SW/fs, where `SW` is in ppm, and `fs` is in Hz. An example value is 700.14 ≈ 14005.6/20.004.
"""

DOCSTRING_fs = """
Sampling frequency of the spectrometer, in Hz. This is the rate at which the complex-valued free-induction decay signal is sampled at. An example value is around 14005.6 Hz for a 700 MHz spectrometer.
"""

DOCSTRING_SW = """
Spectral window in ppm for the spectrometer. An example value is around 20.004 ppm for a 700 MHz spectrometer.
"""

DOCSTRING_ν_0ppm = """
The 0 ppm reference frequency, in Hz. An example value is around 10656.01 Hz for a spectrometer with `SW` ≈ 20.004 ppm and `fs` ≈ 14005.6 Hz.
"""

################# SHsConfigType

# DOCSTRING_coherence_tol(T, show_type = true) = """
# $(show_type ? "`::$(T)`" : "")

# Must be a positive number, as it is a proximity tolerance. Resonance components that doesn't numerically match the -1 quantum coherence condition but falls within this tolerance are kept in the simulation. Larger tolerances produce more resonance components. A value of `0.01` worked well for internal testing of `NMRHamiltonian`.
# """

# DOCSTRING_relative_α_threshold(T, show_type = true) = """
# $(show_type ? "`::$(T)`" : "")

# This should be between 0 and 1. A value of `0.01` means the resonance components that have an intensity less than 0.01 times the maximum intensity in the spin system would be discarded. Larger values means less resonance components are kept.
# """

# DOCSTRING_normalize_α(show_type = true) = """
# $(show_type ? "`::Bool`" : "")

# whether to normalize the intensities of each spin system such that their sum (before discarding via `relative_α_threshold`) adds up to the number of nuclei in the spin system.
# """


################## PhysicalParamsType

DOCSTRING_PhysicalParamsType_H_IDs = """
`::Vector{Int}`

Numerical (positive integer) numbering of nuclei for a molecule entry. This numbering is taken from a file, so it is not guaranteed to be consecutive nor start at 1.
"""

DOCSTRING_PhysicalParamsType_H_inds_sys = """
`::Vector{Vector{Int}}`

Nested list of nuclei labels in consecutive positive integers, which `NMRHamiltonian` refers to as nuclei indices. The element at `[i]` is a list of nuclei indices for the i-th non-singlet spin system.
"""

DOCSTRING_PhysicalParamsType_cs_sys(T, show_type = true) = """
$(show_type ? "`::Vector{Vector{$(T)}}`" : "")

Nested list of non-singlet chemical shifts. The element at `[i]` is a list of chemical shift values for the i-th non-singlet spin system.
"""

DOCSTRING_PhysicalParamsType_H_inds_singlets = """
`::Vector{Vector{Int}}`

Nested list of nuclei labels in consecutive positive integers, which `NMRHamiltonian` refers to as nuclei indices. The element at `[i]` is a list of nuclei indices that contribute to the same i-th singlet.
"""

DOCSTRING_PhysicalParamsType_cs_singlet(T, show_type = true) = """
$(show_type ? "`::Vector{$(T)}`" : "")

List of chemical shift values for singlets.
"""


DOCSTRING_PhysicalParamsType_J_inds_sys = """
`::Vector{Vector{Tuple{Int,Int}}}`

Nested list of pairs of nuclei indices. The element at `[i]` is a list of J-coupling indices for the i-th non-singlet spin system. J-coupling index is a pair of nuclei indices.
"""


DOCSTRING_PhysicalParamsType_J_inds_sys_local = """
`::Vector{Vector{Tuple{Int,Int}}}`

Nested list of pairs of nuclei indices. The element at `[i]` is a list of J-coupling indices for the i-th non-singlet spin system. J-coupling index is a pair of nuclei indices. The nuclei indices are re-numbered to start from 1 for each spin system, which `NMRHamiltonian` refers to as local nuclei indices.

Internally, spin Hamiltonian matrices for each spin system use this particular indexing scheme in their construction.
"""

DOCSTRING_PhysicalParamsType_J_vals_sys(T, show_type = true) = """
$(show_type ? "`::Vector{Vector{$(T)}}`" : "")

Nested list of J-coupling values. The element at `[i]` is a list J-coupling values for the i-th non-singlet spin system.
"""

DOCSTRING_PhysicalParamsType_ME = """
`::Vector{Vector{Vector{Int}}}`

The element at `[i][k]` is a list of nuclei indices that are magnetically equivalent. The indices in this list collectively make up the k-th magnetically equivalent group in the i-th spin system. The nuclei indices are re-numbered to start from 1 for each spin system, which `NMRHamiltonian` refers to as local nuclei indices.

Note that `NMRHamiltonian` defines magnetically equivalent groups must have two or more nuclei, and these equivalence are determined only based on the J-coupling and chemical shift values provided in the corresponding file for the molecule entry at hand.
"""