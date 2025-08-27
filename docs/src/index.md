# MatsubaraFunctions.jl

This package aims at providing a convenient interface to rapidly prototype algorithms for multivariable Green's functions of the form $G_{i_1 ... i_n}(i\omega_1, ..., i\omega_m)$, where $i_k$ denote lattice or orbital indices and $\omega_l$ are fermionic/bosonic Matsubara frequencies.

# Installation

The package is not yet registered, but available from

```julia
pkg> add https://github.com/dominikkiese/MatsubaraFunctions.jl
```

# Related
For an advanced example for the usage of MatsubaraFunctions, see `https://github.com/dominikkiese/MBEsolver.jl`

# Quick start User Guide

## Constructing Meshes
Meshes define the domains (frequency, momentum, index) for your functions.

Matsubara Frequency Mesh
```julia
using MatsubaraFunctions

# Fermionic mesh: temperature = 0.1, N = 10
fermion_mesh = MatsubaraMesh(0.1, 10, Fermion)

# Bosonic mesh: temperature = 0.1, N = 10
boson_mesh = MatsubaraMesh(0.1, 10, Boson)
```

Brillouin Zone Mesh
```julia
using StaticArrays

# 2D square lattice, L = 8
basis = SMatrix{2,2,Float64}([1.0 0.0; 0.0 1.0])
bz = BrillouinZone(8, basis)
bz_mesh = BrillouinZoneMesh(bz)
```

Index Mesh
```julia
index_mesh = IndexMesh(5)  # indices 1:5
```

## Creating a MeshFunction
A MeshFunction stores data on one or more meshes.
```julia
# Create a function on a Matsubara mesh and an index mesh
f = MeshFunction(fermion_mesh, index_mesh; data_t=ComplexF64)

# Assign values
f[1, 2] = 1.0 + 0.0im
```

## Evaluating and Interpolating Functions
```julia
# Direct evaluation at mesh indices
val = f[1, 2]

# Evaluation at mesh points
mp = points(fermion_mesh, 1)
val2 = f(mp, points(index_mesh, 2))

# Interpolation (e.g., at a non-mesh frequency)
val3 = f(0.15, 3)
```

## Using Symmetries
Symmetries can reduce computation and enforce physical constraints.
```julia
# Define a symmetry operation (e.g., sign change)
sym = Symmetry{2}((w) -> (w, Operation{ComplexF64}(sgn=true)))

# Apply symmetry to a tuple of mesh values
result = sym((Index(1), Index(2)))
```

## Saving and Loading Functions
```julia
using HDF5

# Save to HDF5
h = h5open("func_data.h5", "w")
save!(h, "myfunc", f)
close(h)

# Load from HDF5
h = h5open("func_data.h5", "r")
f2 = load_mesh_function(h, "myfunc")
close(h)
```

## Parallelization (MPI)
For large-scale computations, use MPI helpers.
```julia
using MatsubaraFunctions

# Split a range among MPI ranks
my_range = mpi_split(1:100)

# Reduce a MeshFunction across ranks
mpi_allreduce!(f)
```

## Pade Approximation (Analytic Continuation)
```julia
xdat = [1.0, 2.0, 3.0]
ydat = [0.5, 0.2, -0.1]
pade = PadeApprox(xdat, ydat)

# Evaluate the Pade approximant at a new point
val = pade(2.5)
```

## Tips

- Use info(f) to print details about any MeshFunction.
- Use points(mesh) to access mesh points, and value(mp) to get their coordinates.
- Combine multiple meshes for multidimensional functions.