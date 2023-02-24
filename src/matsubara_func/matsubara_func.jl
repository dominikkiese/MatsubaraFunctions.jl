"""
    struct MatsubaraFunction{GD, SD, DD, Q <: Number}

MatsubaraFunction type with fields:
* `grids :: NTuple{GD, MatsubaraGrid}` : collection of MatsubaraGrid
* `shape :: NTuple{SD, Int64}`         : shape of the tensor structure on every grid point
* `data  :: Array{Q, DD}`              : data array

Examples:
```julia 
# construction
T = 1.0
N = 128
g = MatsubaraGrid(T, N, Fermion)

# 1D grid, rank 1 tensor with index dimension 1 (scalar valued)
f1_complex = MatsubaraFunction(g, 1)                    # complex valued (default)
f1_real    = MatsubaraFunction(g, 1, Float64)           # other data type

# 1D grid, rank 1 tensor with index dimension 5 (vector valued)
f2_complex = MatsubaraFunction(g, 5)                    # complex valued (default)
f2_real    = MatsubaraFunction(g, 5, Float64)           # other data type

# 1D grid, rank 2 tensor with index dimension 5 (matrix valued)
f3_complex = MatsubaraFunction(g, (5, 5))               # complex valued (default)
f3_real    = MatsubaraFunction(g, (5, 5), Float64)      # other data type

# 2D grid, rank 2 tensor with index dimension 5 (matrix valued)
f4_complex = MatsubaraFunction((g, g), (5, 5))          # complex valued (default)
f4_real    = MatsubaraFunction((g, g), (5, 5), Float64) # other data type

# usage 
ξ = 0.5
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
println(f(vp, 1; extrp = true))                  # polynomial extrapolation (only for 1D grids)
```
"""
struct MatsubaraFunction{GD, SD, DD, Q <: Number}
    grids :: NTuple{GD, MatsubaraGrid}
    shape :: NTuple{SD, Int64}          
    data  :: Array{Q, DD}

    # basic constructor
    function MatsubaraFunction(
        grids  :: NTuple{GD, MatsubaraGrid}, 
        shape  :: NTuple{SD, Int64}, 
        data   :: Array{Q, DD}
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

        # error for integer data type
        if Q <: Integer || Q <: Complex{Int} error("Integer data type not supported") end
        
        if checks
            # check dimensions
            @assert GD + SD == DD "Dimensions do not match"
            
            # check grids
            for g in grids
                @assert temperature(g) ≈ temperature(grids[1]) "Grids must have same temperature"
                @assert issorted(Float64[value(w) for w in g]) "Grids must be sorted"
            end
        end

        return new{GD, SD, DD, Q}(grids, shape, data)
    end

    # convenience constructors
    function MatsubaraFunction(
        grids  :: NTuple{GD, MatsubaraGrid},
        shape  :: NTuple{SD, Int64},
               :: Type{Q}
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{GD, SD, GD + SD, Q} where {GD, SD, Q <: Number}

        dims = length.(grids)..., shape...
        data = Array{Q, GD + SD}(undef, dims)

        return MatsubaraFunction(grids, shape, data; checks)
    end

    function MatsubaraFunction(
        grid   :: MatsubaraGrid,
        shape  :: NTuple{SD, Int64},
               :: Type{Q}
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{1, SD, 1 + SD, Q} where {SD, Q <: Number}

        return MatsubaraFunction((grid,), shape, Q; checks)
    end

    function MatsubaraFunction(
        grids  :: NTuple{GD, MatsubaraGrid},
        shape  :: Int64,
               :: Type{Q}
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{GD, 1, GD + 1, Q} where {GD, Q <: Number}

        return MatsubaraFunction(grids, (shape,), Q; checks)
    end

    function MatsubaraFunction(
        grid   :: MatsubaraGrid,
        shape  :: Int64,
               :: Type{Q}
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{1, 1, 2, Q} where {Q <: Number}

        return MatsubaraFunction((grid,), (shape,), Q; checks)
    end

    # fallback methods if data type Q is not provided
    function MatsubaraFunction(
        grids  :: NTuple{GD, MatsubaraGrid},
        shape  :: NTuple{SD, Int64}
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{GD, SD, GD + SD, ComplexF64} where {GD, SD}

        return MatsubaraFunction(grids, shape, ComplexF64; checks)
    end

    function MatsubaraFunction(
        grid   :: MatsubaraGrid,
        shape  :: NTuple{SD, Int64},
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{1, SD, 1 + SD, ComplexF64} where {SD}

        return MatsubaraFunction(grid, shape, ComplexF64; checks)
    end

    function MatsubaraFunction(
        grids  :: NTuple{GD, MatsubaraGrid},
        shape  :: Int64,
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{GD, 1, GD + 1, ComplexF64} where {GD}

        return MatsubaraFunction(grids, shape, ComplexF64; checks)
    end

    function MatsubaraFunction(
        grid   :: MatsubaraGrid,
        shape  :: Int64,
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{1, 1, 2, ComplexF64}

        return MatsubaraFunction(grid, shape, ComplexF64; checks)
    end
end



"""
    function grids_shape(
        f :: MatsubaraFunction{GD, SD, DD, Q}
        ) :: NTuple{GD, Int64} where {GD, SD, DD, Q <: Number}

Returns length of grids
"""
function grids_shape(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: NTuple{GD, Int64} where {GD, SD, DD, Q <: Number}

    return length.(f.grids)
end

"""
    function grids_shape(
        f   :: MatsubaraFunction{GD, SD, DD, Q},
        idx :: Int64
        )   :: Int64 where {GD, SD, DD, Q <: Number}

Returns length of f.grids[idx]
"""
function grids_shape(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    idx :: Int64
    )   :: Int64 where {GD, SD, DD, Q <: Number}

    return length(f.grids[idx])
end

"""
    function shape(
        f :: MatsubaraFunction{GD, SD, DD, Q}
        ) :: NTuple{SD, Int64} where {GD, SD, DD, Q <: Number}

Returns f.shape
"""
function shape(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: NTuple{SD, Int64} where {GD, SD, DD, Q <: Number}

    return f.shape 
end 

"""
    function shape(
        f   :: MatsubaraFunction{GD, SD, DD, Q},
        idx :: Int64
        )   :: Int64 where {GD, SD, DD, Q <: Number}

Returns f.shape[idx]
"""
function shape(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    idx :: Int64
    )   :: Int64 where {GD, SD, DD, Q <: Number}

    return f.shape[idx]
end 

"""
    function data_shape(
        f :: MatsubaraFunction{GD, SD, DD, Q}
        ) :: NTuple{DD, Int64} where {GD, SD, DD, Q <: Number}

Returns shape of f.data
"""
function data_shape(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: NTuple{DD, Int64} where {GD, SD, DD, Q <: Number}

    return size(f.data)
end

"""
    function data_shape(
        f   :: MatsubaraFunction{GD, SD, DD, Q},
        idx :: Int64
        )   :: Int64 where {GD, SD, DD, Q <: Number}

Returns length of dimension idx of f.data
"""
function data_shape(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    idx :: Int64
    )   :: Int64 where {GD, SD, DD, Q <: Number}

    return size(f.data, idx)
end

"""
    function absmax(
        f :: MatsubaraFunction{GD, SD, DD, Q}
        ) :: Float64 where {GD, SD, DD, Q <: Number}

Returns largest element of f.data (in absolute terms)
"""
function absmax(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: Float64 where {GD, SD, DD, Q <: Number}

    return maximum(abs.(f.data))
end

"""
    function Base.:argmax(
        f :: MatsubaraFunction{GD, SD, DD, Q}
        ) :: CartesianIndex{DD} where {GD, SD, DD, Q <: Number}

Returns position of largest element of f.data (in absolute terms)
"""
function Base.:argmax(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: CartesianIndex{DD} where {GD, SD, DD, Q <: Number}

    return argmax(abs.(f.data))
end

"""
    function info(
        f :: MatsubaraFunction{GD, SD, DD, Q}
        ) :: Nothing where {GD, SD, DD, Q <: Number}

Prints some properties of f
"""
function info(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: Nothing where {GD, SD, DD, Q <: Number}

    println("MatsubaraFunction properties")
    println("----------------------------")
    println("Temperature     : $(temperature(f.grids[1]))")
    println("Grid dimension  : $(GD)")
    println("Shape dimension : $(SD)")
    println("Data dimension  : $(DD)")
    println("Data type       : $(Q)")

    return nothing
end



# load methods 
include("mpi_helpers.jl")
include("matsubara_func_ops.jl")
include("matsubara_func_index.jl")
include("matsubara_func_eval.jl")
include("matsubara_func_symmetries.jl")