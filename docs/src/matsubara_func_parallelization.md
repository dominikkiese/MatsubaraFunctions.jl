

# Advanced Usage: Running in parallel

To simplify the parallelization of algorithms when using the package, we export some preliminary methods based on the MPI.jl wrapper. For further information on how to set up MPI with Julia see [https://github.com/JuliaParallel/MPI.jl](https://github.com/JuliaParallel/MPI.jl).

```julia
using MatsubaraFunctions 
using MPI 

MPI.Init()
mpi_info()
mpi_println("I print on main.")
ismain = mpi_ismain() # ismain = true if rank is 0

T = 1.0
N_= 128
g = MatsubaraGrid(T, N_, Fermion)
f = MatsubaraFunction(g, 1)

# simple loop parallelization for UnitRange
for vidx in mpi_split(1 : length(g))
    println("My rank is $(mpi_rank()): $(vidx)")
end

# simple (+) allreduce
mpi_allreduce!(f)
```

In addition, calls of `MatsubaraSymmetryGroup` with an initialization function have an opt-in switch (`mode`) to enable parallel evaluation of the `MatsubaraInitFunction` (by default `mode = :serial`). If `mode = :polyester`, shared memory multithreading via the `Polyester` ([https://github.com/JuliaSIMD/Polyester.jl](https://github.com/JuliaSIMD/Polyester.jl)) Julia package is used. This mode is recommended for initialization functions that are cheap to evaluate and are unlikely to benefit from `Threads.@threads` due to the overhead from invoking the Julia scheduler. For more expensive functions, users can choose between `mode = :threads`, which simply uses `Threads.@threads`, and `mode = :hybrid`. The latter combines both MPI and native Julia threads and can therefore be used to run calculations on multiple compute nodes.


# Functions

```@docs
mpi_comm
```

```@docs
mpi_rank
```

```@docs
mpi_size
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
