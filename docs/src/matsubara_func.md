# MatsubaraFunction

## Constructor

A `MatsubaraFunction` is a collection of `MatsubaraGrid` instances together with an associated tensor structure $G_{i_1...i_n}$ for each point $(\omega_1, ..., \omega_m)$ in the cartesian product of the grids. 

```@docs
MatsubaraFunction
```

Some possible constructors are

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



# Indexing and assignment

There are two possible ways to access the data of a `MatsubaraFunction`, using either the bracket `[]` or the parenthesis `()` operator. The former can be used together with a set of linear indices or with a combination of `MatsubaraFrequency` objects and linear indices (for the tensor structure). It will return the value of the function precisely for the given arguments. `()` also allows to substitute `Float64` for the frequency arguments, in which case a multilinear interpolation is performed. In addition, `()` is well defined even for out of bounds access, since it makes use of either polynomial or constant extrapolation in this case.

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

# polynomial extrapolation in 1D, constant term set to 1 (default is 0)
vp = MatsubaraFrequency(T, 256, Fermion)
println(f(vp; extrp = ComplexF64(1.0))) 
```



```@docs
set!
```


```@docs
grid_axes
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

## Getter Functions

```@docs
grids
```

```@docs
shape
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



# Miscellaneous utilities
```@docs
Base.length
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
absmax
```

```@docs
argmax
```





# I/O to HDF5 files


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



```@docs
save_matsubara_function!
```    

```@docs
load_matsubara_function
```    

