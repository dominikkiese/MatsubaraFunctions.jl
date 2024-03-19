# MatsubaraGrid

A `MatsubaraGrid{PT}` is a sorted (symmetric) set of `MatsubaraFrequency{PT}` objects and can be constructed by

```julia
T  = 1.0
N  = 128
g1 = MatsubaraGrid(T, N, Fermion) # total no. frequencies is 2N 
g2 = MatsubaraGrid(T, N, Boson)   # total no. frequencies is 2N - 1
```

where N is the number of non-negative frequencies which are defined as follows:

|Particle type:             | Fermion | Boson    |
|---------------------------|---------|----------|
|total no. of frequencies   | 2N      | 2N-1     |
|Range of Matsubara index n | -N: N-1 | -N+1:N-1 |
|def. of frequencies v_n    | (2n+1)πT| 2nπT     |

`MatsubaraGrid{PT}` instances are iterable

```julia
T = 1.0
N = 128
g = MatsubaraGrid(T, N, Fermion)

for v in g
  println(value(v)) 
  println(index(v))
end
```

and can be evaluated using either a `MatsubaraFrequency{PT}`, `MatsubaraIndex{PT}` or `Float64`. As long as the input argument is in bounds, this will return the corresponding linear index of the grid in the two former cases and the linear index of the closest frequency in the latter case 

```julia
T   = 1.0
N   = 128
g   = MatsubaraGrid(T, N, Fermion)
idx = rand(eachindex(g))
@assert g(g[idx]) == idx 
@assert g(value(g[idx])) == idx 
```

`MatsubaraGrid{PT}` objects can be saved in HDF5 file format as

```julia
using MatsubaraFunctions 
using HDF5

file = h5open("test.h5", "w")
T    = 1.0
N    = 128
g    = MatsubaraGrid(T, N, Fermion)

save_matsubara_grid!(file, "grid", g) 
gp = load_matsubara_grid(file, "grid")
close(file)
```

# Types

```@docs
AbstractMatsubaraGrid
```

```@docs
MatsubaraGrid
```

# Functions

```@docs
firstindex
```

```@docs
lastindex
```

```@docs
axes
```

```@docs
is_inbounds
```

```@docs
N
```

```@docs
firstvalue
```

```@docs
lastvalue
```

```@docs
indices
```

```@docs
values
```

```@docs
info
```

```@docs
save_matsubara_grid!
```       

```@docs
load_matsubara_grid
```    