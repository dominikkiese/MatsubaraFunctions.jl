"""
    struct MatsubaraFunction{GD, SD, DD, Q <: Number}

MatsubaraFunction type with fields:
* `grids :: NTuple{GD, MatsubaraGrid}` : collection of MatsubaraGrid
* `shape :: NTuple{SD, Int64}`         : shape of the tensor structure on every grid point
* `data  :: Array{Q, DD}`              : data array
"""
struct MatsubaraFunction{GD, SD, DD, Q <: Number}
    grids :: NTuple{GD, MatsubaraGrid}
    shape :: NTuple{SD, Int64}          
    data  :: Array{Q, DD}

    # default constructor
    function MatsubaraFunction(
        grids :: NTuple{GD, MatsubaraGrid}, 
        shape :: NTuple{SD, Int64}, 
        data  :: Array{Q, DD}
        )     :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

        # error for integer data type
        if Q <: Integer || Q <: Complex{Int} error("Integer data type not supported") end
        
        # check dimensions
        @check GD > 0 "Grid dimension cannot be zero"
        @check SD > 0 "Shape dimension cannot be zero"
        @check GD + SD == DD "Dimensions do not match"
        
        # check grids
        for g in grids
            @check temperature(g) ≈ temperature(grids[1]) "Grids must have same temperature"
            @check issorted(values(g)) "Grids must be sorted"
        end

        return new{GD, SD, DD, Q}(grids, shape, data)
    end

    # convenience constructors
    function MatsubaraFunction(
        grids :: NTuple{GD, MatsubaraGrid},
        shape :: NTuple{SD, Int64},
              :: Type{Q}
        )     :: MatsubaraFunction{GD, SD, GD + SD, Q} where {GD, SD, Q <: Number}

        dims = length.(grids)..., shape...
        data = Array{Q, GD + SD}(undef, dims)

        return matsubara_function_grid_cp(grids, shape, data)
    end

    function MatsubaraFunction(
        grid  :: MatsubaraGrid,
        shape :: NTuple{SD, Int64},
              :: Type{Q}
        )     :: MatsubaraFunction{1, SD, 1 + SD, Q} where {SD, Q <: Number}

        return MatsubaraFunction((grid,), shape, Q)
    end

    function MatsubaraFunction(
        grids :: NTuple{GD, MatsubaraGrid},
        shape :: Int64,
              :: Type{Q}
        )     :: MatsubaraFunction{GD, 1, GD + 1, Q} where {GD, Q <: Number}

        return MatsubaraFunction(grids, (shape,), Q)
    end

    function MatsubaraFunction(
        grid  :: MatsubaraGrid,
        shape :: Int64,
              :: Type{Q}
        )     :: MatsubaraFunction{1, 1, 2, Q} where {Q <: Number}

        return MatsubaraFunction((grid,), (shape,), Q)
    end

    # fallback methods if data type Q is not provided
    function MatsubaraFunction(
        grids :: NTuple{GD, MatsubaraGrid},
        shape :: NTuple{SD, Int64}
        )     :: MatsubaraFunction{GD, SD, GD + SD, ComplexF64} where {GD, SD}

        return MatsubaraFunction(grids, shape, ComplexF64)
    end

    function MatsubaraFunction(
        grid  :: MatsubaraGrid,
        shape :: NTuple{SD, Int64},
        )     :: MatsubaraFunction{1, SD, 1 + SD, ComplexF64} where {SD}

        return MatsubaraFunction(grid, shape, ComplexF64)
    end

    function MatsubaraFunction(
        grids :: NTuple{GD, MatsubaraGrid},
        shape :: Int64,
        )     :: MatsubaraFunction{GD, 1, GD + 1, ComplexF64} where {GD}

        return MatsubaraFunction(grids, shape, ComplexF64)
    end

    function MatsubaraFunction(
        grid  :: MatsubaraGrid,
        shape :: Int64,
        )     :: MatsubaraFunction{1, 1, 2, ComplexF64}

        return MatsubaraFunction(grid, shape, ComplexF64)
    end

    # copy constructor 
    function MatsubaraFunction(
        f :: MatsubaraFunction{GD, SD, DD, Q}
        ) :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

        return matsubara_function_grid_cp(grids(f), shape(f), copy(f.data))
    end
end

# safer version of MatsubaraFunction default constructor using grid copies
function matsubara_function_grid_cp(
    grids :: NTuple{GD, MatsubaraGrid}, 
    shape :: NTuple{SD, Int64}, 
    data  :: Array{Q, DD}
    )     :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

    return MatsubaraFunction(copy.(grids), shape, data)
end

"""
    function grids(
        f :: MatsubaraFunction{GD, SD, DD, Q}
        ) :: NTuple{GD, MatsubaraGrid} where {GD, SD, DD, Q <: Number}

Returns `f.grids`
"""
function grids(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: NTuple{GD, MatsubaraGrid} where {GD, SD, DD, Q <: Number}

    return f.grids
end

"""
    function grids(
        f   :: MatsubaraFunction{GD, SD, DD, Q},
        idx :: Int64
        )   :: MatsubaraGrid where {GD, SD, DD, Q <: Number}

Returns `f.grids[idx]`
"""
function grids(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    idx :: Int64
    )   :: MatsubaraGrid where {GD, SD, DD, Q <: Number}

    return f.grids[idx]
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

    return length.(grids(f))
end

"""
    function grids_shape(
        f   :: MatsubaraFunction{GD, SD, DD, Q},
        idx :: Int64
        )   :: Int64 where {GD, SD, DD, Q <: Number}

Returns length of `f.grids[idx]`
"""
function grids_shape(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    idx :: Int64
    )   :: Int64 where {GD, SD, DD, Q <: Number}

    return length(grids(f, idx))
end

"""
    function shape(
        f :: MatsubaraFunction{GD, SD, DD, Q}
        ) :: NTuple{SD, Int64} where {GD, SD, DD, Q <: Number}

Returns `f.shape`
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

Returns `f.shape[idx]`
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

Returns shape of `f.data`
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

Returns length of dimension `idx` of `f.data`
"""
function data_shape(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    idx :: Int64
    )   :: Int64 where {GD, SD, DD, Q <: Number}

    return size(f.data, idx)
end

"""
    function length(
        f :: MatsubaraFunction{GD, SD, DD, Q}
        ) :: Int64 where {GD, SD, DD, Q <: Number}

Returns length of `f.data`
"""
function Base.:length(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: Int64 where {GD, SD, DD, Q <: Number}

    return length(f.data)
end

"""
    function temperature(
        f :: MatsubaraFunction{GD, SD, DD, Q}
        ) :: Float64 where {GD, SD, DD, Q <: Number}

Returns temperature for which `f.grids` are defined
"""
function temperature(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: Float64 where {GD, SD, DD, Q <: Number}

    return temperature(grids(f, 1))
end

function Base.:copy(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

    return MatsubaraFunction(f)
end

"""
    function absmax(
        f :: MatsubaraFunction{GD, SD, DD, Q}
        ) :: Float64 where {GD, SD, DD, Q <: Number}

Returns largest element of `f.data` (in absolute terms)
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

Returns position of largest element of `f.data` (in absolute terms)
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

Prints summary of function properties
"""
function info(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: Nothing where {GD, SD, DD, Q <: Number}

    println("MatsubaraFunction properties")
    println("----------------------------")
    println("Temperature     : $(temperature(f))")
    println("Grid dimension  : $(GD)")
    println("Shape dimension : $(SD)")
    println("Data dimension  : $(DD)")
    println("Data type       : $(Q)")

    return nothing
end

# load methods 
include("matsubara_func_ops.jl")
include("matsubara_func_index.jl")
include("matsubara_func_eval.jl")
include("matsubara_func_symmetries.jl")
include("matsubara_func_io.jl")