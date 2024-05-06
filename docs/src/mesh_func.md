# MeshFunction

```@meta
CurrentModule = MatsubaraFunctions
```

# Constructor

A `MeshFunction` is a collection of $D$ `Mesh` instances together with an associated data array $G_{i_1...i_D}$ for each point $(\omega_1, ..., \omega_D)$ in the cartesian product of the meshes. 
For more info on meshes see [General info on meshes](@ref).

```@docs
MeshFunction
```

Some possible constructors are

```julia
T = 1.0
N = 128
m1 = MatsubaraMesh(T, N, Fermion)   # constructs Mesh of Matsubara frequencies
m2 = IndexMesh(5)                   # constructs Mesh of indices

# single-frequency Matsubara function (scalar):
f1_complex = MeshFunction(m1) 
f1_real    = MeshFunction(m1; data_t=Float64) 

# single-frequency Matsubara function with additional index structure:
## Matsubara frequencies in the 1. dimension and index structure in the 2.
f2_complex = MeshFunction(m1, m2) 
f2_real    = MeshFunction(m1, m2; data_t=Float64) 

# Matsubara frequencies in the 1. dimension and index structure in the 2. and 3.
f3_complex = MeshFunction(m1, m2, m2) 
f3_real    = MeshFunction(m1, m2, m2; data_t=Float64) 

# two-frequency Matsubara function with additional index structure:
f4_complex = MeshFunction(m1, m1, m2, m2) 
f4_real    = MeshFunction(m1, m1, m2, m2; data_t=Float64) 
```




# Indexing and assignment

There are two possible ways to access the data of a `MeshFunction`, using either 
* the bracket `[]` for indexing and assignment: It will return (or assign) the value of the function precisely for the given arguments (indices or `MeshPoint` objects).
* or the parenthesis `()` operator: In addition, `()` is well defined even for out of bounds access. `()` also allows to insert `Float64` for the frequency arguments, in which case multilinear interpolation is performed.


```julia
ξ = 0.5
T = 1.0
N = 128
g = MatsubaraMesh(T, N, Fermion)    
f = MeshFunction(g)                 

# assign a value at each Matsubara frequency in g:
for v in g
    f[v] = 1.0 / (im * plain_value(v) - ξ)
end 

# access MeshFunction data
v = g[rand(eachindex(g))]
println(f[v])               # fast data access, throws error if v is out of bounds
println(f(v; lim=1.))       # fast data access, defined even if v is out of bounds in which case lim is returned (standard value for lim is 0)
println(f(plain_value(v)))  # slow data access, uses interpolation 

```


#### Internals about Indexing 

```@docs
LinearIndex
```


```@docs
LinearIndex_bc
```

```@docs
to_meshes
```

#### Internals about Linear interpolation

```@docs
InterpolationParam
```

```@docs
indices(:: InterpolationParam{N}) where N
```

```@docs
weights
```


# I/O to HDF5 files

`MeshFunction` objects can be saved in HDF5 file format as

```julia
using MeshFunctions 
using HDF5

file = h5open("test.h5", "w")
T    = 1.0
N    = 128
g    = MatsubaraMesh(T, N, Fermion)
f    = MeshFunction(g, 1)

save_matsubara_function!(file, "func", f)
fp = load_matsubara_function(file, "func")
close(file)
```

```@docs
save!(:: HDF5.File, :: String, :: MeshFunction)
```    

```@docs
load_mesh_function
```    


# Interfacing with Julia's solvers

```@docs
flatten
```

```@docs
flatten!
```

```@docs
unflatten!
```



# Getter and setter functions

```@docs
meshes
```

```@docs
set!
```

# Arithmetic operations

```@docs
add
```

```@docs
add!
```

```@docs
subtract
```

```@docs
subtract!
```

```@docs
mult
```

```@docs
mult!
```


```@docs
mult_add!
```




# Miscellaneous utilities


```@docs
absmax
```

```@docs
arg_absmax
```





