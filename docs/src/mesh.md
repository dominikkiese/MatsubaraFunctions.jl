# General mesh information

`Mesh` instances (aka sorted sets of `MeshPoint`s) provide the axes for `MeshFunction`s. 
Indexing a `MeshFunction` with a `MeshPoint` is as fast as with regular indices. 
To ensure that `MeshPoint`s are not used for the wrong axis a hash comparison is performed if 
`MatsubaraFunctions.DEBUG() = true` is set.

A `MeshPoint` $m_i$ may, for example, represent:
 * a Matsubara frequency, see [Matsubara mesh](@ref),
 * a momentum point for a periodic lattice, see [Brillouin zone mesh](@ref),
 * an ordinary index, see [Index mesh](@ref),
 * other user-defined types, see [Custom mesh implementation](@ref).

## Usage

There are specialized constructors for each mesh type. For example:

```julia
T  = 1.0
N  = 128
m1 = MatsubaraMesh(T, N, Fermion) # range of 2N Matsubara frequencies (symmetric around zero)

M  = 16
m2 = IndexMesh(M) # indices from 1 to M
```

`Mesh` instances are iterable

```julia
T = 1.0
N = 128
m = MatsubaraMesh(T, N, Fermion)

for v in m
  println(index(v))       # prints linear mesh index
  println(value(v))       # prints value associated with the mesh point (here: a MatsubaraFrequency struct)
  println(plain_value(v)) # prints bare value of the mesh point (here: a Float64)
end
```

`Mesh` objects can be saved in HDF5 file format as

```julia
using MatsubaraFunctions 
using HDF5

file = h5open("test.h5", "w")
T    = 1.0
N    = 128
g    = MatsubaraMesh(T, N, Fermion)

save!(file, "grid", g) 
gp = load_mesh(file, "grid")
close(file)
```

## API

TBD

## Types

```@docs
Mesh
```

```@docs
AbstractMesh
```

```@docs
MeshPoint
```

```@docs
AbstractMeshPoint
```

```@docs
MatsubaraFunctions.AbstractDomain
```

```@docs
AbstractValue
```

## Functions

```@docs
points
```

```@docs
domain
```

```@docs
index(:: T) where {T <: AbstractMeshPoint}
```

```@docs
value( :: MeshPoint{T}) where {T <: AbstractValue}
```

```@docs
plain_value
```

```@docs
load_mesh(:: HDF5.File, :: String)
```    