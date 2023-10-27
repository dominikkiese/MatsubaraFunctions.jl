[![CI](https://github.com/dominikkiese/MatsubaraFunctions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/dominikkiese/MatsubaraFunctions.jl/actions/workflows/CI.yml)
[![Documentation](https://github.com/dominikkiese/MatsubaraFunctions.jl/actions/workflows/Documentation.yml/badge.svg)](https://dominikkiese.github.io/MatsubaraFunctions.jl/dev/)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![DOI](https://zenodo.org/badge/585688328.svg)](https://zenodo.org/doi/10.5281/zenodo.10048247)

# MatsubaraFunctions.jl

This package aims at providing a convenient interface to rapidly prototype algorithms for multivariable Green's functions of the form $G_{i_1 ... i_n}(i\omega_1, ..., i\omega_m)$,
where $i_k$ could denote lattice or orbital indices and $\omega_l$ are fermionic/bosonic Matsubara frequencies.

# Installation

The package is not yet registered, but available from

```julia
pkg> add https://github.com/dominikkiese/MatsubaraFunctions.jl
```
