# Advanced Usage: Automated Symmetry Reduction

In many cases, symmetries can be used to drastically reduce the numerical effort for self-consistent calculations. For one-particle Green's functions $G_{i_1 i_2}(i\omega)$, for example, hermicity of the Hamiltonian dictates that $G_{i_1 i_2}(i\omega) = G^{\star}_{i_2 i_1}(-i\omega)$, relating positive and negative Matsubara frequencies. We offer an automated way to compute the set of irreducible (i.e. unrelatable by symmetries) `MeshFunction` components, as is illustrated in the following example

```julia
using MatsubaraFunctions

ξ  = 0.5
T  = 1.0
N_ = 128
g  = MatsubaraMesh(T, N_, Fermion)
f  = MeshFunction(g)

for v in g
    f[v] = 1.0 / (im * plain_value(v) - ξ)
end 

# complex conjugation acting on Green's function
function conj(
    w :: Tuple{MatsubaraFrequency}
    ) :: Tuple{Tuple{MatsubaraFrequency}, Operation}

    return (-w[1],), Operation(false, true)
end 

# compute the symmetry group 
SG = SymmetryGroup([Symmetry{1}(conj)], f)

# obtain another Green's function by symmetrization
function init(
    w :: Tuple{MatsubaraFrequency}
    ) :: ComplexF64

    return f(w...)
end 

InitFunc = InitFunction{1, ComplexF64}(init)
h        = MeshFunction(g)
SG(h, InitFunc)
@assert h == f
```

## Functions and Functors for symmetry reduction

```@docs
Operation
```

```@docs
Symmetry
```

```@docs
SymmetryGroup
```

```@docs
InitFunction
```

```@docs
get_reduced
```

```@docs
init_from_reduced!
```

```@docs
sgn
```

```@docs
con
```

## I/O to HDF5 files

TBD