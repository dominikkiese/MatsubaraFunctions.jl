# Go through the publication of MatsubaraFunctions.jl a write a version of the code examples compatible with the unstable branch.

using MatsubaraFunctions
using StaticArrays

# Section 3.1
T = 1.0
k = 10
v = MatsubaraFrequency(T, k, Fermion)
w = MatsubaraFrequency(T, k, Boson)


### OUTPUT FROM CHATGPT

# 1. Basic types and meshes


## Fermionic and bosonic Matsubara meshes:

T = 0.1     # temperature
N = 32      # number of positive frequencies

mf_fermi = MatsubaraMesh(T, N, Fermion)
mf_boson = MatsubaraMesh(T, N, Boson)

length(mf_fermi), length(mf_boson)
collect(eachindex(mf_fermi))  # linear indices

## The i-th mesh point and its numeric value
i = 33
ωi = points(mf_fermi, i)   # MeshPoint
ω = value(ωi)              # MatsubaraFrequency
ω.value                    # numeric value, Float64
value(ω)                   # same as above
ωi == mf_fermi[i]          # indexing


## Index meshes (for discrete orbital/spin/lattice indices):
idx_mesh = IndexMesh(5)          # indices 1:5
length(idx_mesh)
points(idx_mesh, 2)  

## Brillouin-zone meshes:

### 2D square lattice reciprocal basis (identity)
L = 16  # mesh size along each direction
basis = @SMatrix [1.0 0.0; 0.0 1.0]  # 2x2 static matrix

bz = BrillouinZone(L, basis)
bzmesh = BrillouinZoneMesh(bz)

length(bzmesh)

# Convert CartesianIndex to linear index
cart_idx = CartesianIndex(2, 3)
lin_idx = LinearIndices((L, L))[cart_idx]

kpt = points(bzmesh, lin_idx)
k = value(kpt)
k.value[1], k.value[2]  # k-vector components


# 2. MeshFunction basics

Gω = MeshFunction(mf_fermi; data_t=ComplexF64)

# Assign via cartesian indices
Gω[3] = 1.0 + im
Gω[2] = 0.5 - im

# Read back and function-style evaluation
Gω[3]                          # at linear index      
Gω(points(mf_fermi, 3))        # at mesh point
Gω(value(points(mf_fermi, 3))) # off-mesh numeric evaluation (interpolated)
Gω(-19.1)                       # off-mesh numeric evaluation (interpolated)
Gω(1000.0)                      # off-mesh numeric evaluation works also out of bounds


# Vector-valued function using an index mesh:
vec_idx = IndexMesh(5)
Σωi = MeshFunction(mf_fermi, vec_idx; data_t=ComplexF64)
Σωi[3, 2] = 0.5 - 0.1im
Σωi[3, 2]
Σωi(points(mf_fermi, 3), points(vec_idx, 2))

# Structure and quick info:
size(Gω), length(Gω), collect(eachindex(Gω))
info(Gω)

# Basic arithmetic and comparisons:
Hω = MeshFunction(Gω.meshes... ; data_t=eltype(Gω.data))

Hω[3] = 0.25im

Sω = Gω + Hω
Tω = 2 * Gω - Hω

# alternative: manipulate the data arrays directly
Tω_data = 2 .* Gω.data .- Hω.data
Tω = MeshFunction(Gω.meshes...; data_t=eltype(Gω.data))
Tω.data .= Tω_data

Sω == Tω                    #   equality of MeshFunction objects
all(Sω.data .== Tω.data)    # elementwise equality
# should be false in this example :)


# 3. Evaluation, interpolation, and off-mesh access

## Function-call syntax supports MeshPoints or numeric values:

ω_test = value(points(mf_fermi, 5)).value + 0.1
val = Gω(ω_test)   # off-mesh; uses interpolation along Matsubara axis

