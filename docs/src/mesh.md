# General info on meshes
`Mesh` instances describe the axes of MeshFunctions.
A `Mesh` $\vec{m}$ is a sorted set of $N$ `MeshPoint`s $m_i$ corresponding to the indices $i\in 1,...,N$.
Indexing a `MeshFunction` with a `MeshPoint` is as fast as with regular indices (e.g. `Int`).
To make sure that a `MeshPoint` is not used for the wrong axis both `MeshPoint` and `Mesh` contain a hash which are compared if `MatsubaraFunctions.DEBUG() = true` is set.

A `MeshPoint` $m_i$ may represent
 * Matsubara frequencies $\omega$, see [Matsubara frequency mesh](@ref),
 * a point in momentum space $\vec{k}$, see [Brillouin zone mesh](@ref),
 * or ordinary indices, see [Index mesh](@ref).
For convenience the package offers corresponding constructors, arithmetric operations, I/O and many more functionalities. 



## Usage
There are specialized constructors for each type of `Mesh`.
For example:

```julia
T  = 1.0
N  = 128
m1 = MatsubaraMesh(T, N, Fermion) # represents range of 2N Matsubara frequencies (symmetric around zero)

M = 16
m2 = IndexMesh(M)                 # represents indices 1,...,M
```

`Mesh` instances are iterable

```julia
T = 1.0
N = 128
m = MatsubaraMesh(T, N, Fermion)

for v in m
  println(index(v))         # prints index 1,...,length(m)
  println(value(v))         # prints object of type MatsubaraFrequency
  println(plain_value(v))   # prints Float64 corresponding to the MatsubaraFrequency; equivalent to value(value(v))
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


# API


### Types

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





### Functions
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




### To be implemented for a new mesh type
Imagine we develop a new type of mesh called `NewMeshT` (subtype of `AbstractMesh`).
The mesh points have a new value type called `NewValueT` (subtype of `AbstractValue`).
We also need an abstract domain which encodes all the details of a mesh.

As a reference, consider our most simple mesh, `IndexMesh`, for which we have

| Abstract type            | Concrete type                            |
|--------------------------|------------------------------------------|
| AbstractMesh             | Mesh{MeshPoint{Index}, IndexDomain}      |
| AbstractValue            | Index                                    |
| AbstractDomain           | IndexDomain                              |

Additional to the new structs, one needs to implement:
* specialized constructor(s)
* `value(::NewValueT <: AbstractValue)`: needed for calling `plain_value(mp::MeshPoint{NewValueT})`.
* `is_inbounds(::NewValueT <: AbstractValue, ::NewMeshT <: AbstractMesh)`: Checks if the value is within the bounds of the mesh
* `save!()`, `load_mesh()`: for I/O to HDF5 files
* `mesh_index(::NewValueT <: AbstractValue, ::NewMeshT <: AbstractMesh)`: maps `NewValueT` instance to `Int` for indexing (without considering boundary conditions)
* `mesh_index_bc(::NewValueT <: AbstractValue, ::NewMeshT <: AbstractMesh)`: same as `mesh_index(...)` __with__ boundary conditions.
