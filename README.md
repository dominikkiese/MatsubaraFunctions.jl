# MatsubaraFunctions.jl

This package aims at providing a convenient interface to rapidly prototype algorithms for multivariable Green's functions of the form $G_{i_1 ... i_n}(i\omega_1, ..., i\omega_m)$ where $i_k$ denotes a lattice (orbital) index and $\omega_l$ a Matsubara frequency. The latter can either be of fermionic or bosonic type, that is $\omega^{n}_l = \frac{\pi}{\beta}(2n + 1)$ or $\omega^{n}_l = \frac{2\pi n}{\beta}$ where $n \in \mathbb{Z}$. Here, $\beta = 1/ T$ is the inverse temperature.


# Installation

```julia
pkg> add https://github.com/dominikkiese/MatsubaraFunctions.jl
```


# Basic Usage 

MatsubaraFunctions.jl provides three elementary classes: `MatsubaraFrequency`, `MatsubaraGrid` and `MatsubaraFunction`, each of them shipping with its own set of convenience constructors, getter functions and allowed operations which are partially detailed below. For further examples, see our tests. <br /> 


A `MatsubaraFrequency` at temperature $T$ and with index $n$ can be generated using

```julia
T = 1.0
n = 5
v = MatsubaraFrequency(T, n, Fermion)
w = MatsubaraFrequency(T, n, Boson) 
```

Matsubara frequencies can be added, subtracted and their sign can be reversed, producing a new instance of `MatsubaraFrequency`.

```julia
w1 = v + v # type(v1) = :Boson
v2 = w - v # type(v2) = :Fermion
v3 = -v    # type(v3) = :Fermion 
```


A `MatsubaraGrid` is a sorted (symmetric) set of `MatsubaraFrequency` objects and can be constructed by

```julia
T  = 1.0
N  = 128
g1 = MatsubaraGrid(T, N, Fermion) # no. frequencies is 2N
g2 = MatsubaraGrid(T, N, Boson)   # no. frequencies is 2N - 1
info(g1)
println()
info(g2)
```

`MatsubaraGrid` instances are iterable

```julia
T = 1.0
N = 128
g = MatsubaraGrid(T, N, Fermion)

for v in g
  println(value(v)) 
  println(index(v))
end
```

and can be evaluated using either a `MatsubaraFrequency` or `Float64`. As long as the input argument is in bounds, this will return the corresponding linear index of the grid in the former case and the linear index of the closest frequency in the latter case 

```julia
T   = 1.0
N   = 128
g   = MatsubaraGrid(T, N, Fermion)
idx = rand(1 : length(g))
@assert g(g[idx]) == idx 
@assert g(value(g[idx])) == idx 
```


A `MatsubaraFunction` is a collection of `MatsubaraGrid` instances together with an associated tensor structure $G_{i_1...i_n}$ for each point $(\omega_1, ..., \omega_m)$ in the cartesian product of the grids. Some possible constructors are

```julia
T = 1.0
N = 128
g = MatsubaraGrid(T, N, Fermion)

# 1D grid, rank 1 tensor with index dimension 1 (scalar valued)
f1_complex = MatsubaraFunction(g, 1) 
f1_real    = MatsubaraFunction(g, 1, Float64) 
info(f1_complex)
println()
info(f1_real)
println()

# 1D grid, rank 1 tensor with index dimension 5 (vector valued)
f2_complex = MatsubaraFunction(g, 5) 
f2_real    = MatsubaraFunction(g, 5, Float64) 
info(f2_complex)
println()
info(f2_real)
println()

# 1D grid, rank 2 tensor with index dimension 5 (matrix valued)
f3_complex = MatsubaraFunction(g, (5, 5)) 
f3_real    = MatsubaraFunction(g, (5, 5), Float64) 
info(f3_complex)
println()
info(f3_real)
println()

# 2D grid, rank 2 tensor with index dimension 5 (matrix valued)
f4_complex = MatsubaraFunction((g, g), (5, 5)) 
f4_real    = MatsubaraFunction((g, g), (5, 5), Float64) 
info(f4_complex)
println()
info(f4_real)
println()
```

There are two possible ways to access the data of a `MatsubaraFunction`, using either the bracket `[]` or the parenthesis `()` operator. The former can be used together with a set of linear indices or with a combination of `MatsubaraFrequency` objects and linear indices. It will return the value of the function exactly for the given arguments. The latter allows to substitute `Float64` for the frequency arguments, in which case a multilinear interpolation is performed. In addition, `()` is well defined even for out of bounds access, providing the user with the option to provide a boundary condition or (for 1D grids) use polynomial extrapolation.

```julia
ξ = 0.5
T = 1.0
N = 128
g = MatsubaraGrid(T, N, Fermion)
f = MatsubaraFunction(g, 1)

for v in g
    f[v, 1] = 1.0 / (im * value(v) - ξ)
end 

# access MatsubaraFunction data
println(f[v, 1])        # fast data access, throws error if v is out of bounds
println(f(v, 1))        # fast data access, defined even if v is out of bounds
println(f(value(v), 1)) # slow data access, uses interpolation 

# fallback methods for out of bounds access
vp = MatsubaraFrequency(T, 256, Fermion)
println(f(vp, 1))                                # default x -> 0.0
println(f(vp, 1; bc = x -> 1.0))                 # custom boundary condition x -> 1.0
println(f(vp, 1; bc = x -> 1.0 / im / value(x))) # custom boundary condition x -> 1.0 / im / value(x)
println(f(value(vp), 1; bc = x -> 1.0 / im / x)) # custom boundary condition x -> 1.0 / im / x
println(f(vp, 1; extrp = true))                  # polynomial extrapolation in 1D 
```

`MatsubaraGrid` and `MatsubaraFunction` objects can be saved in HDF5 file format as

```julia
using MatsubaraFunctions 
using HDF5

file = h5open("test.h5", "w")
T    = 1.0
N    = 128
g    = MatsubaraGrid(T, N, Fermion)
f    = MatsubaraFunction(g, 1)

save_matsubara_grid!(file, "grid", g) 
save_matsubara_function!(file, "func", f)
gp = load_matsubara_grid(file, "grid")
fp = load_matsubara_function(file, "func")
```

# Advanced Usage: Matsubara Sums 

For `MatsubaraFunction` objects $G_{i_1 ... i_n}(i\omega)$ defined on 1D grids, we export the function `sum_me`, which computes the series $\Sigma_m G_{i_1 ... i_n}(i\omega^{m}) e^{i\omega^m 0^+}$ for $m \in \mathbb{Z}$ using tail fits of $G$ together with analytic formulas for summations of the form $\Sigma_m \frac{1}{(i\omega^m)^\alpha}e^{i\omega^m 0^+}$ with $\alpha \in \mathbb{N}$. This, however, requires $G$ to be representable by a Laurent series in an elongated annulus about the imaginary axis.

```julia
ξ = 0.5
T = 1.0
N = 128
g = MatsubaraGrid(T, N, Fermion)
f = MatsubaraFunction(g, 1)

for v in g
    f[v, 1] = 1.0 / (im * value(v) - ξ)
end 

# evaluate the series and compare to analytic result
ρ(x, T) = 1.0 / (exp(x / T) + 1.0)
println(abs(sum_me(f1_complex, 1) - (ρ(+ξ, T) - 1.0)))
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
    f[v, 1] = 1.0 / (im * value(v) - ξ)
end 

# complex conjugation acting on Green's function
function conj(
    w :: Tuple{MatsubaraFrequency},
    x :: Tuple{Int64}
    ) :: Tuple{Tuple{MatsubaraFrequency}, Tuple{Int64}, Operation}

    return (-w[1],), (x[1],), Operation(false, true)
end 

# compute the symmetry group 
SG = SymmetryGroup([Symmetry{1, 1}(conj)], f)

# symmetrize and compare to f
ftest = deepcopy(f)

for class in SG.classes 
    ftest[class[1][1], class[1][2]...] = f[class[1][1], class[1][2]...]
end 

SG(ftest)
println(maximum(abs.(ftest.data .- f.data)))
```

# Advanced Usage: MPI Helpers

To simplify the parallelization of algorithms involving `MatsubaraFunction` instances, we export some useful methods based on the MPI.jl wrapper. For further information on how 
to set up MPI with Julia see https://github.com/JuliaParallel/MPI.jl.

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
  comm = MPI.COMM_WORLD 
  println("My rank is $(MPI.Comm_rank(comm)): $(vidx)")
end

# simple (+) allreduce
mpi_allreduce!(f)
```

















