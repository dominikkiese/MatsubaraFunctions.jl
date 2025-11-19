# Practical Examples (Updated to current API)

This page updates the examples from the original paper to the current MatsubaraFunctions.jl API. It covers mesh construction, MeshFunction usage, evaluation/interpolation, I/O, symmetries, MPI helpers, Pade approximation, and two physics examples (HF and GW) in the atom limit.

Prerequisites:
```julia
using MatsubaraFunctions
using HDF5
```

## 1. Basic types and meshes

Fermionic and bosonic Matsubara meshes:
```julia
β = 10.0    # inverse temperature
N = 32      # number of positive frequencies

mf_fermi = MatsubaraMesh(β, N, Fermion)
mf_boson = MatsubaraMesh(β, N, Boson)

length(mf_fermi), length(mf_boson)
collect(eachindex(mf_fermi))  # linear indices

# The i-th mesh point and its numeric value
i = 3
ωi = points(mf_fermi, i)   # MeshPoint
ω = value(ωi)              # Float64
```

Index meshes (for discrete orbital/spin/lattice indices):
```julia
im = IndexMesh(5)          # indices 1:5
length(im)
points(im, 2)              # MeshPoint representing index 2
```

Brillouin-zone meshes:
```julia
# 2D square lattice reciprocal basis (identity)
basis = [1.0 0.0; 0.0 1.0]
bz = BrillouinZone(basis)
bzmesh = BrillouinZoneMesh(bz; size=(16, 16))

length(bzmesh)
kpt = points(bzmesh, CartesianIndex(2, 3))
k = value(kpt)             # SVector-like or vector of Float64
```

## 2. MeshFunction basics

Create MeshFunction objects (on one or more meshes), assign values, and inspect.

Scalar 1D function on fermionic frequencies:
```julia
Gω = MeshFunction(mf_fermi; data_t=ComplexF64)

# Assign via cartesian indices
Gω[3] = 1.0 + 0.0im

# Read back and function-style evaluation
Gω[3]
Gω(3)                          # at mesh index
Gω(points(mf_fermi, 3))        # at mesh point
Gω(value(points(mf_fermi, 3))) # off-mesh numeric evaluation (interpolated)
```

Vector-valued function using an index mesh:
```julia
vec_idx = IndexMesh(5)
Σωi = MeshFunction(mf_fermi, vec_idx; data_t=ComplexF64)
Σωi[3, 2] = 0.5 - 0.1im
Σωi(3, 2)
```

Structure and quick info:
```julia
size(Gω), length(Gω), collect(eachindex(Gω))
info(Gω)
```

Basic arithmetic and comparisons:
```julia
Hω = similar(Gω)
Hω[3] = 0.25im

Sω = Gω + Hω
Tω = 2 .* Gω .- Hω

Sω == Tω  # elementwise equality
```

## 3. Evaluation, interpolation, and off-mesh access

Function-call syntax supports MeshPoints or numeric values:
```julia
ω_test = value(points(mf_fermi, 5)) + 0.1
val = Gω(ω_test)   # off-mesh; uses interpolation along Matsubara axis
```

Note: For fine control over interpolation, see the docstrings in func_itp.jl.

## 4. I/O with HDF5

Save and load meshes and MeshFunction objects:
```julia
# Save
h5open("example.h5", "w") do h
    save!(h, "mesh/fermi", mf_fermi)
    save!(h, "func/Gω", Gω)
    save!(h, "func/Σωi", Σωi)
end

# Load
mf_fermi2 = h5open("example.h5", "r") do h
    load_mesh(h, "mesh/fermi")
end

Gω2 = h5open("example.h5", "r") do h
    load_mesh_function(h, "func/Gω")
end
```

## 5. Symmetries

Define and use symmetry operations with SymmetryGroup to reduce computation. Example: complex conjugation symmetry G*(ν) = G(-ν).

```julia
using MatsubaraFunctions: Symmetry, SymmetryGroup, Operation

# A one-arg symmetry mapping (frequency -> -frequency), with conjugation
conj_sym = Symmetry{1}() do args
    ω = args[1]                   # MeshPoint
    return ( -ω, ), Operation(sgn=false, con=true)
end

# Build a MeshFunction and a group
G = MeshFunction(mf_fermi; data_t=ComplexF64)
SG = SymmetryGroup([conj_sym], G)

# Provide an initializer for irreducible elements
init = x -> begin
    ω = x[1]                      # MeshPoint
    1.0 / (im * value(ω) + 0.5)
end

# Fill G over irreducible set, propagate via symmetries
SG(G, init)
```

You can query equivalence classes or irreducible indices:
```julia
classes(SG, length(mf_fermi))
irreducible(SG, length(mf_fermi))
```

## 6. Parallelization helpers (MPI)

Use MPI to split work and reduce results:
```julia
using MPI
MPI.Init()

chunk = mpi_split(1:100)          # rank-local subrange
local_sum = sum(chunk)
global_sum = mpi_allreduce(local_sum)

# Reduce MeshFunction data in-place
mpi_allreduce!(Gω)

MPI.Finalize()
```

Run with:
```sh
mpiexec -n 4 julia --project -e 'using Pkg; Pkg.instantiate(); include("your_script.jl")'
```

## 7. Pade approximation (analytic continuation)

Fit a Pade approximant and evaluate:
```julia
x = collect(range(-1, 1; length=21))
y = @. 1 / (1 + x^2) + 0.01randn()
P = PadeApprox(x, y)
P(0.3)
```

## 8. Hartree-Fock (HF) in the atomic limit (updated)

Self-consistent solution for density n in the HF approximation. We avoid deprecated APIs and implement the Matsubara sum directly.

```julia
using MatsubaraFunctions
using NLsolve

β = 1/0.3
U = 0.9
N = 1000

mf = MatsubaraMesh(β, N, Fermion)
G  = MeshFunction(mf; data_t=ComplexF64)

# Helper: Matsubara sum of a MeshFunction f(ν)
function matsubara_sum(f::MeshFunction)
    T = 1/β
    s = zero(eltype(f))
    for i in eachindex(mf)
        s += f(i)                 # f(i) == f(points(mf, i))
    end
    return T*s
end

# Initialize G0(ν) = 1/(iν)
for i in eachindex(mf)
    ω = value(points(mf, i))
    G[i] = 1.0 / (im*ω)
end

# Fixed-point: n = T*Σν G(ν) e^{iν0+}; and G(ν) = 1/(iν - U n)
function fixed_point!(F, nvec)
    n = nvec[1]
    for i in eachindex(mf)
        ω = value(points(mf, i))
        G[i] = 1.0 / (im*ω - U*n)
    end
    F[1] = real(matsubara_sum(G)) - n
    return nothing
end

n0 = [real(matsubara_sum(G))]
res = nlsolve((F, n) -> fixed_point!(F, n), n0, method=:anderson)
n_HF = res.zero[1]
```

Notes:
- We compute the Matsubara sum explicitly from the mesh. If you have a convenience function in your environment for density, you can swap it in.
- This example uses 1D MeshFunction; for multi-index objects, use index meshes.

## 9. GW in the atomic limit (updated)

A compact GW sketch leveraging the new Mesh/MeshFunction API and symmetry helpers.

```julia
β = 1/0.3
U = 0.9
N = 1000

mfF = MatsubaraMesh(β, N, Fermion)
mfB = MatsubaraMesh(β, N, Boson)

G     = MeshFunction(mfF; data_t=ComplexF64)
Σ     = MeshFunction(mfF; data_t=ComplexF64)
Pω    = MeshFunction(mfB; data_t=ComplexF64)
ηD, ηM = MeshFunction(mfB; data_t=ComplexF64), MeshFunction(mfB; data_t=ComplexF64)

# Symmetry: conjugation G*(ν) = G(-ν) => reuse for Σ and P
using MatsubaraFunctions: Symmetry, SymmetryGroup, Operation

conjF = Symmetry{1}() do args
    ( -args[1], ), Operation(sgn=false, con=true)
end
conjB = Symmetry{1}() do args
    ( -args[1], ), Operation(sgn=false, con=true)
end

SGF = SymmetryGroup([conjF], G)
SGB = SymmetryGroup([conjB], Pω)

# Initialize G0
for i in eachindex(mfF)
    ω = value(points(mfF, i))
    G[i] = 1.0 / (im*ω)
end

# A larger summation mesh can help tails
sum_mesh = MatsubaraMesh(β, 4N, Fermion)

# Compute P(Ω) = T Σν G(Ω+ν) G(ν), using symmetries
function init_P(args)
    Ω = value(args[1])
    T = 1/β
    s = 0.0 + 0.0im
    for j in eachindex(sum_mesh)
        ν = value(points(sum_mesh, j))
        s += G(ν + Ω) * G(ν)
    end
    return T*s
end
SGB(Pω, init_P)

# Screened interactions
for i in eachindex(mfB)
    valP = Pω(i)
    ηD[i] = +U / (1 - U*valP)
    ηM[i] = -U / (1 + U*valP)
end

# Self-energy with simple symmetry-accelerated initializer
function init_Σ(args)
    ν = value(args[1])
    T = 1/β
    # Un = U * n; here approximate n from G
    n = real(T*sum(G(j) for j in eachindex(mfF)))
    s = U*n
    # Convolution with screened interactions
    acc = 0.0 + 0.0im
    for j in eachindex(sum_mesh)
        νp = value(points(sum_mesh, j))
        acc += G(νp) * ( 0.25*ηD(ν - νp) + 0.75*ηM(ν - νp) + 0.5*U )
    end
    return s - T*acc
end
SGF(Σ, init_Σ)

# Dyson update for G
for i in eachindex(mfF)
    ω = value(points(mfF, i))
    G[i] = 1.0 / (im*ω - Σ[i])
end
```

Notes:
- This is a compact single-pass sketch to illustrate the updated interfaces. In practice, wrap it in a non-linear fixed-point loop (e.g., Anderson) updating G, P, η, Σ.

## 10. Multiboson exchange (MBE/SBE) building blocks

The current API supports efficient index computation and direct data access. Below are updated skeletons of the performance-critical parts.

Efficient index computation and constant extrapolation along Matsubara axes:
```julia
# Advanced/internal: grid_index_extrp is internal; call via module
# It finds the nearest in-bounds linear index for extrapolation.
# Use with care and profile in your application.

w1 = ω + ν + νp
η_idx = MatsubaraFunctions.grid_index_extrp(w1, grids(ηD, 1))  # 1st mesh of ηD
val = ηD[η_idx]
```

Batched precomputation of irreducible vertex tensors T^D, T^M:
```julia
function calc_T_ph!(
    T_D::MeshFunction, T_M::MeshFunction,
    η_S::MeshFunction, λ_S::MeshFunction,
    η_D::MeshFunction, λ_D::MeshFunction,
    η_M::MeshFunction, λ_M::MeshFunction,
    M_S::MeshFunction, M_T::MeshFunction, M_D::MeshFunction, M_M::MeshFunction,
    U::Float64,
)
    Threads.@threads for ip in eachindex(grids(T_D, 3))   # ν'
        vp = value(points(grids(T_D, 3), ip))
        λ1_idx2 = MatsubaraFunctions.grid_index_extrp(vp, grids(λ_D, 2))
        vp_idx  = MatsubaraFunctions.grid_index_extrp(vp, grids(M_S, 2))

        for iv in eachindex(grids(T_D, 2))                 # ν
            v  = value(points(grids(T_D, 2), iv))
            w2 = vp - v

            λ1_idx3 = MatsubaraFunctions.grid_index_extrp(v,  grids(λ_D, 2))
            η2_idx  = MatsubaraFunctions.grid_index_extrp(w2, grids(η_D, 1))
            λ2_idx1 = MatsubaraFunctions.grid_index_extrp(w2, grids(λ_D, 1))
            v_idx   = MatsubaraFunctions.grid_index_extrp(v,  grids(M_S, 2))
            w2_idx  = MatsubaraFunctions.grid_index_extrp(w2, grids(M_S, 1))

            for iw in eachindex(grids(T_D, 1))             # Ω
                w  = value(points(grids(T_D, 1), iw))
                w1 = w + v + vp
                v2 = w + v

                η1_idx  = MatsubaraFunctions.grid_index_extrp(w1, grids(η_D, 1))
                λ1_idx1 = MatsubaraFunctions.grid_index_extrp(w1, grids(λ_D, 1))
                λ2_idx2 = MatsubaraFunctions.grid_index_extrp(v2, grids(λ_D, 2))
                w_idx   = MatsubaraFunctions.grid_index_extrp(w,  grids(M_S, 1))
                w1_idx  = MatsubaraFunctions.grid_index_extrp(w1, grids(M_S, 1))
                v2_idx  = MatsubaraFunctions.grid_index_extrp(v2, grids(M_S, 2))

                # SBE parts
                p1 = λ_S[λ1_idx1, λ1_idx2, 1] * η_S[η1_idx, 1] * λ_S[λ1_idx1, λ1_idx3, 1]
                p2 = λ_D[λ2_idx1, λ1_idx3, 1] * η_D[η2_idx, 1] * λ_D[λ2_idx1, λ2_idx2, 1]
                p3 = λ_M[λ2_idx1, λ1_idx3, 1] * η_M[η2_idx, 1] * λ_M[λ2_idx1, λ2_idx2, 1]

                # MBE parts
                m1 = M_S[w1_idx, v_idx, vp_idx, 1]
                m2 = M_T[w1_idx, v_idx, vp_idx, 1]
                m3 = M_D[w2_idx, v_idx, v2_idx, 1]
                m4 = M_M[w2_idx, v_idx, v2_idx, 1]

                T_D[iw, iv, ip, 1] = -2U + M_D[w_idx, v_idx, vp_idx, 1] + 0.5*(p1 + m1 - p2 - m3) + 1.5*(m2 - p3 - m4)
                T_M[iw, iv, ip, 1] = +2U + M_M[w_idx, v_idx, vp_idx, 1] - 0.5*(p1 + m1 + p2 + m3) + 0.5*(m2 + p3 + m4)
            end
        end
    end
    return nothing
end
```

Bethe-Salpeter-like MBE step with symmetry acceleration and views:
```julia
function calc_M!(
    M::MeshFunction,    # target M^D or M^M
    Π::MeshFunction,    # bubble slice
    T::MeshFunction,    # precomputed irreducible
    M_D::MeshFunction,  # for channel coupling
    SG::SymmetryGroup   # symmetry group of M
)
    Tmesh1, Tmesh2, Tmesh3 = grids(T, 1), grids(T, 2), grids(T, 3)
    Πmesh2 = grids(Π, 2)

    init = function(args)
        Ω, ν, νp = value(args[1]), value(args[2]), value(args[3])

        # Construct views for fast contractions along a fermionic axis
        v1, v2 = Πmesh2(points(Tmesh3, 1)), Πmesh2(points(Tmesh3, length(Tmesh3)))
        Π_slice  = view(Π, points(grids(Π, 1), Π(grids(Π, 1)(Ω))), v1:v2)
        M_slice  = view(M_D, points(Tmesh1, Tmesh1(Ω)),
                             points(Tmesh2, Tmesh2(ν)),
                             :)
        TL_slice = view(T, points(Tmesh1, Tmesh1(Ω)),
                           points(Tmesh2, Tmesh2(ν)),
                           :)
        TR_slice = view(T, points(Tmesh1, Tmesh1(Ω)),
                           points(Tmesh2, Tmesh2(νp)),
                           :)

        # Handle offsets consistently (illustrative)
        acc = zero(eltype(Π))
        @inbounds for i in 1:length(TL_slice)
            acc -= (TL_slice[i] - M_slice[min(i, size(M_slice, 1))]) * Π_slice[i] * TR_slice[i]
        end
        return (1/β) * acc
    end

    SG(M, init; mode=:hybrid)
    return nothing
end
```

Notes:
- The MBE code above is a template illustrating the updated APIs and common optimization patterns (preindexing, views, symmetry).
- grid_index_extrp is an internal helper; profile and validate carefully before relying on it.

---

Tips:
- Use info(obj) to inspect Mesh and MeshFunction objects.
- Use points(mesh, i) to retrieve MeshPoints and value(::MeshPoint) to get coordinates.
- Save and load with save!(h5, path, obj) and load_mesh/load_mesh_function.
- Prefer broadcasting and direct buffer reuse for performance-critical sections.