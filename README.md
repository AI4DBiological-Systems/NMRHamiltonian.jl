# NMRHamiltonian.jl
Given chemical shifts and J-coupling values, produce simulated amplitude and frequencies, and partitions of them.

# Documentation
[https://AI4DBiological-Systems.github.io/NMRHamiltonian.jl/](https://AI4DBiological-Systems.github.io/NMRHamiltonian.jl/)

See `/examples/adjust_phys.jl` for how to modify and save J-coupling values read from file, such that the magnetic equivalence of the molecule is maximized.

The tutorial in the documentation is based on `/examples/simulate.jl`.

# Install
Add the custom registries, then the package.
```
using Pkg

pkg"registry add https://github.com/RoyCCWang/RWPublicJuliaRegistry"

pkg"registry add https://github.com/AI4DBiological-Systems/PublicJuliaRegistry"

pkg"add NMRHamiltonian"
```

# Citation
Our work is undergoing peer review. Please cite our ChemRxiv version if you use this software.
[https://doi.org/10.26434/chemrxiv-2023-0s196](https://doi.org/10.26434/chemrxiv-2023-0s196)

# License
See LICENSE.md for the licence.
