# Advanced Usage: Automated Symmetry Reduction

In many cases, the numerical effort of computing functions in the Matsubara domain can be drastically reduced by the use of symmetries. For one-particle Green's functions $G_{i_1 i_2}(i\omega)$, for example, hermicity of the Hamiltonian dictates that $G_{i_1 i_2}(i\omega) = G^{\star}_{i_2 i_1}(-i\omega)$, relating positive and negative Matsubara frequencies. We offer an automated way to compute the set of irreducible (i.e. unrelatable by symmetries) `MatsubaraFunction` components, as is illustrated in the following example

```julia
ξ = 0.5
T = 1.0
N_= 128
g = MatsubaraGrid(T, N_, Fermion)
f = MatsubaraFunction(g, 1)

for v in g
    f[v, 1] = 1.0 / (im * value(v) - ξ)
end 

# complex conjugation acting on Green's function
function conj(
    w :: Tuple{MatsubaraFrequency},
    x :: Tuple{Int64}
    ) :: Tuple{Tuple{MatsubaraFrequency}, Tuple{Int64}, MatsubaraOperation}

    return (-w[1],), (x[1],), MatsubaraOperation(false, true)
end 

# compute the symmetry group 
SG = MatsubaraSymmetryGroup([MatsubaraSymmetry{1, 1}(conj)], f)

# obtain another Green's function by symmetrization
function init(
    w :: Tuple{MatsubaraFrequency},
    x :: Tuple{Int64}
    ) :: ComplexF64

    return f[w, x...]
end 

InitFunc = MatsubaraInitFunction{1, 1, ComplexF64}(init)
h        = MatsubaraFunction(g, 1)
SG(h, InitFunc)
@assert h == f
```

## Functions and Functors for symmetry reduction


```@docs
MatsubaraOperation
```

```@docs
MatsubaraSymmetry
```

```@docs
MatsubaraSymmetryGroup
```

```@docs
MatsubaraInitFunction
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



# I/O to HDF5 files

```@docs
save_matsubara_symmetry_group!
```    


