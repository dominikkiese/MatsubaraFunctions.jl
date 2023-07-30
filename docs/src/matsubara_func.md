# MatsubaraFunction

A `MatsubaraFunction` is a collection of `MatsubaraGrid` instances together with an associated tensor structure $G_{i_1...i_n}$ for each point $(\omega_1, ..., \omega_m)$ in the cartesian product of the grids. Some possible constructors are

```julia
T = 1.0
N = 128
g = MatsubaraGrid(T, N, Fermion)

# 1D grid, rank 1 tensor with index dimension 1 (scalar valued)
f1_complex = MatsubaraFunction(g, 1) 
f1_real    = MatsubaraFunction(g, 1, Float64) 

# 1D grid, rank 1 tensor with index dimension 5 (vector valued)
f2_complex = MatsubaraFunction(g, 5) 
f2_real    = MatsubaraFunction(g, 5, Float64) 

# 1D grid, rank 2 tensor with index dimension 5 (matrix valued)
f3_complex = MatsubaraFunction(g, (5, 5)) 
f3_real    = MatsubaraFunction(g, (5, 5), Float64) 

# 2D grid, rank 2 tensor with index dimension 5 (matrix valued)
f4_complex = MatsubaraFunction((g, g), (5, 5)) 
f4_real    = MatsubaraFunction((g, g), (5, 5), Float64) 
```

There are two possible ways to access the data of a `MatsubaraFunction`, using either the bracket `[]` or the parenthesis `()` operator. The former can be used together with a set of linear indices or with a combination of `MatsubaraFrequency` objects and linear indices (for the tensor structure). It will return the value of the function precisely for the given arguments. `()` allows to substitute `Float64` for the frequency arguments, in which case a multilinear interpolation is performed. In addition, `()` is well defined even for out of bounds access, using either a custom boundary condition or, for 1D grids, polynomial extrapolation.

```julia
ξ = 0.5
T = 1.0
N = 128
g = MatsubaraGrid(T, N, Fermion)
f = MatsubaraFunction(g, 1)

for v in g
    # if there is only one index of dimension 1, it does not need to be specified, i.e. 
    # f[v] can be used instead of f[v, 1] (also works for the '()' operator)
    f[v] = 1.0 / (im * value(v) - ξ)
end 

# access MatsubaraFunction data
v = g[rand(eachindex(g))]
println(f[v])        # fast data access, throws error if v is out of bounds
println(f(v))        # fast data access, defined even if v is out of bounds
println(f(value(v))) # slow data access, uses interpolation 

# fallback methods for out of bounds access
vp = MatsubaraFrequency(T, 256, Fermion)
println(f(vp))                                  # default x -> 0.0
println(f(vp; bc = x -> 1.0))                   # custom boundary condition x -> 1.0
println(f(vp; bc = x -> 1.0 / im / value(x)))   # custom boundary condition x -> 1.0 / im / value(x)
println(f(value(vp); bc = x -> 1.0 / im / x))   # custom boundary condition x -> 1.0 / im / x
println(f(vp; extrp = (true, ComplexF64(0.0)))) # polynomial extrapolation in 1D, constant term set to 0.0
```

`MatsubaraFunction` objects can be saved in HDF5 file format as

```julia
using MatsubaraFunctions 
using HDF5

file = h5open("test.h5", "w")
T    = 1.0
N    = 128
g    = MatsubaraGrid(T, N, Fermion)
f    = MatsubaraFunction(g, 1)

save_matsubara_function!(file, "func", f)
fp = load_matsubara_function(file, "func")
close(file)
```

# Advanced Usage: Matsubara Sums 

For `MatsubaraFunction` objects $G_{i_1 ... i_n}(i\omega)$ defined on 1D grids, we export the function `sum_me`, which computes the series $\Sigma_m G_{i_1 ... i_n}(i\omega_{m}) e^{i\omega_m 0^+}$ for $m \in \mathbb{Z}$ using tail fits of $G$ together with analytic formulas for summations of the form $\Sigma_m \frac{1}{(i\omega_m)^\alpha}e^{i\omega_m 0^+}$ with $\alpha \in \mathbb{N}$. This, however, requires $G$ to be representable by a Laurent series in an elongated annulus about the imaginary axis.

```julia
ξ = 0.5
T = 1.0
N = 128
g = MatsubaraGrid(T, N, Fermion)
f = MatsubaraFunction(g, 1)

for v in g
    f[v] = 1.0 / (im * value(v) - ξ)
end 

# evaluate the series and compare to analytic result
ρ(x, T) = 1.0 / (exp(x / T) + 1.0)
println(abs(sum_me(f, ComplexF64(0.0)) - (ρ(+ξ, T) - 1.0)))
```

# Advanced Usage: Automated Symmetry Reduction

In many cases, the numerical effort of computing functions in the Matsubara domain can be drastically reduced by the use of symmetries. For one-particle Green's functions $G_{i_1 i_2}(i\omega)$, for example, hermicity of the Hamiltonian dictates that $G_{i_1 i_2}(i\omega) = G^{\star}_{i_2 i_1}(-i\omega)$, relating positive and negative Matsubara frequencies. We offer an automated way to compute the set of irreducible (i.e. unrelatable by symmetries) `MatsubaraFunction` components, as is illustrated in the following example

```julia
ξ = 0.5
T = 1.0
N = 128
g = MatsubaraGrid(T, N, Fermion)
f = MatsubaraFunction(g, 1)

for v in g
    f[v] = 1.0 / (im * value(v) - ξ)
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

# Advanced Usage: Running in parallel

To simplify the parallelization of algorithms when using the package, we export some useful methods based on the MPI.jl wrapper. For further information on how to set up MPI with Julia see [https://github.com/JuliaParallel/MPI.jl](https://github.com/JuliaParallel/MPI.jl).

```julia
using MatsubaraFunctions 
using MPI 

MPI.Init()
mpi_info()
mpi_println("I print on main.")
ismain = mpi_ismain() # ismain = true if rank is 0

T = 1.0
N = 128
g = MatsubaraGrid(T, N, Fermion)
f = MatsubaraFunction(g, 1)

# simple loop parallelization for UnitRange
for vidx in mpi_split(1 : length(g))
    println("My rank is $(mpi_rank()): $(vidx)")
end

# simple (+) allreduce
mpi_allreduce!(f)
```

In addition, calls of `MatsubaraSymmetryGroup` with an initialization function have an opt-in switch (`mode`) to enable parallel evaluation of the `MatsubaraInitFunction` (by default `mode = :serial`). If `mode = :polyester`, shared memory multithreading via the `Polyester` ([https://github.com/JuliaSIMD/Polyester.jl](https://github.com/JuliaSIMD/Polyester.jl)) Julia package is used. This mode is recommended for initialization functions that are cheap to evaluate and are unlikely to benefit from `Threads.@threads` due to the overhead from invoking the Julia scheduler. For more expensive functions, users can choose between `mode = :threads`, which simply uses `Threads.@threads`, and `mode = :hybrid`. The latter combines both MPI and native Julia threads and can therefore be used to run calculations on multiple compute nodes.

# Types

```@docs
MatsubaraFunction
```

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
PadeApprox
```

# Functions

```@docs
grids
```

```@docs
grids_shape
```

```@docs
shape
```

```@docs
data_shape
```

```@docs
absmax
```

```@docs
argmax
```

```@docs
mpi_split
```

```@docs
mpi_allreduce!
```

```@docs
mpi_ismain
```

```@docs
mpi_println
```

```@docs
mpi_info
```

```@docs
mpi_barrier
```

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
set!
```

```@docs
flatten
```

```@docs
flatten!
```

```@docs
unflatten!
```

```@docs
LinearIndex
```

```@docs
to_Matsubara
```

```@docs
upper_tail_moments
```

```@docs
lower_tail_moments
```

```@docs
sum_me
```

```@docs
density
```

```@docs
sgn
```

```@docs
con
```

```@docs
save_matsubara_function!
```    

```@docs
load_matsubara_function
```    

```@docs
save_matsubara_symmetry_group!
```    

```@docs
load_matsubara_symmetry_group
```    

```@docs
coeffs
```    

```@docs
xdat
```    