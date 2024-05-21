# MeshFunction

```@meta
CurrentModule = MatsubaraFunctions
```

## Constructor

A `MeshFunction` is a collection of $D$ `Mesh` instances together with an associated data array $G_{i_1...i_D}$ for each point 
$(x_1, ..., x_D)$ in the cartesian product of the meshes. For more detail, see [General mesh information](@ref).

```@docs
MeshFunction
```

Examples of possible constructors are:

```julia
T = 1.0
N = 128
m1 = MatsubaraMesh(T, N, Fermion) # a mesh of Matsubara frequencies
m2 = IndexMesh(5)                 # a mesh of indices

# scalar-valued mesh function with a single frequency argument:
f1_complex = MeshFunction(m1) 
f1_real    = MeshFunction(m1; data_t = Float64) 

# vector-valued mesh function with a single frequency argument:
f2_complex = MeshFunction(m1, m2) 
f2_real    = MeshFunction(m1, m2; data_t = Float64) 

# matrix-valued mesh function with a single frequency argument:
f3_complex = MeshFunction(m1, m2, m2) 
f3_real    = MeshFunction(m1, m2, m2; data_t = Float64) 

# matrix-valued mesh function with two frequency arguments:
f4_complex = MeshFunction(m1, m1, m2, m2) 
f4_real    = MeshFunction(m1, m1, m2, m2; data_t=Float64) 
```

## Indexing and assignment

There are two possible ways to access the data of a `MeshFunction`: 
* the bracket `[]` for indexing and assignment. It will return (or assign) the value of the function for the given arguments 
  (indices or `MeshPoint` objects).
* the parenthesis `()`. As opposed to `[]`, `()` is defined even for out-of-bounds access and can be called with plain value types 
  (e.g. Float64 for Matsubara frequencies), in which case a multilinear interpolation is performed.

```julia
ξ = 0.5
T = 1.0
N = 128
g = MatsubaraMesh(T, N, Fermion)    
f = MeshFunction(g)                 

# assign a value to each Matsubara frequency in g:
for v in g
    f[v] = 1.0 / (im * plain_value(v) - ξ)
end 

# access MeshFunction data
v = g[rand(eachindex(g))]
println(f[v])              # fast data access, throws error if v is out of bounds
println(f(v; lim = 1.0))   # fast data access, returns lim if v is out of bounds (lim defaults to 0)
println(f(plain_value(v))) # slow data access, uses interpolation 
```

### Indexing details

```@docs
LinearIndex
```

```@docs
LinearIndex_bc
```

```@docs
to_meshes
```

### Linear interpolation details

```@docs
InterpolationParam
```

```@docs
indices(:: InterpolationParam{N}) where N
```

```@docs
weights
```

## I/O

`MeshFunction` objects can be saved in HDF5 file format as

```julia
using MeshFunctions 
using HDF5

file = h5open("test.h5", "w")
T    = 1.0
N    = 128
g    = MatsubaraMesh(T, N, Fermion)
f    = MeshFunction(g, 1)

save!(file, "func", f)
fp = load_mesh_function(file, "func")
close(file)
```

```@docs
save!(:: HDF5.File, :: String, :: MeshFunction)
```    

```@docs
load_mesh_function
```    

## Interfacing with other Julia packages

```@docs
flatten
```

```@docs
flatten!
```

```@docs
unflatten!
```

## Getter and setter functions

```@docs
meshes
```

```@docs
set!
```

## Arithmetic operations

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

## Utility

```@docs
absmax
```

```@docs
arg_absmax
```