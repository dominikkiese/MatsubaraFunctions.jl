"""
    struct MatsubaraFunction{GD, SD, DD, Q <: Number}

MatsubaraFunction type with fields:
* `grids :: NTuple{GD, AbstractMatsubaraGrid}` : collection of MatsubaraGrid
* `shape :: NTuple{SD, Int64}`                 : shape of the tensor structure on every grid point
* `data  :: Array{Q, DD}`                      : multidimensional data array
"""
struct MatsubaraFunction{GD, SD, DD, Q <: Number}
    grids :: NTuple{GD, AbstractMatsubaraGrid}
    shape :: NTuple{SD, Int64}          
    data  :: Array{Q, DD}

    function MatsubaraFunction(
        grids :: NTuple{GD, AbstractMatsubaraGrid}, 
        shape :: NTuple{SD, Int64}, 
        data  :: Array{Q, DD}
        )     :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

        if Q <: Integer || Q <: Complex{Int} error("Integer data type not supported") end
        
        # check dimensions
        @check GD > 0 "Grid dimension cannot be zero"
        @check GD + SD == DD "Dimensions do not match"
        
        # check grids
        T = temperature(grids[1])

        for g in grids
            @check temperature(g) ≈ T "Grids must have same temperature"
            @check issorted(values(g)) "Grids must be sorted"
        end

        return new{GD, SD, DD, Q}(grids, shape, data)
    end

    function MatsubaraFunction(
        grids  :: NTuple{GD, AbstractMatsubaraGrid},
        shape  :: Vararg{Int64, SD}
        ;
        data_t :: DataType = ComplexF64
        )      :: MatsubaraFunction{GD, SD, GD + SD, data_t} where {GD, SD}
        
        data = Array{data_t, GD + SD}(undef, length.(grids)..., shape...)
        return MatsubaraFunction(copy.(grids), ntuple(i -> shape[i], SD), data)
    end

    function MatsubaraFunction(
        grid   :: AbstractMatsubaraGrid,
        shape  :: Vararg{Int64, SD}
        ;
        data_t :: DataType = ComplexF64
        )      :: MatsubaraFunction{1, SD, 1 + SD, Type{data_t}} where {SD}

        return MatsubaraFunction((grid,), shape...; data_t)
    end

    function MatsubaraFunction(
        grids  :: Vararg{AbstractMatsubaraGrid, GD},
        ;
        data_t :: DataType = ComplexF64
        )      :: MatsubaraFunction{GD, 0, GD, Type{data_t}} where {GD}

        return MatsubaraFunction(ntuple(i -> grids[i], GD); data_t)
    end

    function MatsubaraFunction(
        f :: MatsubaraFunction
        ) :: MatsubaraFunction

        return MatsubaraFunction(grids(f), shape(f), copy(f.data))
    end

    # specialization
    function MatsubaraFunction(
        grid   :: AbstractMatsubaraGrid
        ;
        data_t :: DataType = ComplexF64
        )      :: MatsubaraFunction{1, 0, 1, Type{data_t}}

        return MatsubaraFunction((grid,); data_t)
    end
end

"""
    function grids(
        f :: MatsubaraFunction{GD, SD, DD, Q}
        ) :: NTuple{GD, AbstractMatsubaraGrid} where {GD, SD, DD, Q <: Number}

Returns `f.grids`
"""
function grids(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: NTuple{GD, AbstractMatsubaraGrid} where {GD, SD, DD, Q <: Number}

    return f.grids
end

"""
    function grids(
        f   :: MatsubaraFunction,
        idx :: Int64
        )   :: AbstractMatsubaraGrid

Returns `f.grids[idx]`
"""
function grids(
    f   :: MatsubaraFunction,
    idx :: Int64
    )   :: AbstractMatsubaraGrid

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
        f   :: MatsubaraFunction,
        idx :: Int64
        )   :: Int64

Returns length of `f.grids[idx]`
"""
function grids_shape(
    f   :: MatsubaraFunction,
    idx :: Int64
    )   :: Int64

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
        f   :: MatsubaraFunction,
        idx :: Int64
        )   :: Int64

Returns `f.shape[idx]`
"""
function shape(
    f   :: MatsubaraFunction,
    idx :: Int64
    )   :: Int64

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
        f   :: MatsubaraFunction,
        idx :: Int64
        )   :: Int64

Returns length of dimension `idx` of `f.data`
"""
function data_shape(
    f   :: MatsubaraFunction,
    idx :: Int64
    )   :: Int64

    return size(f.data, idx)
end

"""
    function length(
        f :: MatsubaraFunction
        ) :: Int64

Returns length of `f.data`
"""
function Base.:length(
    f :: MatsubaraFunction
    ) :: Int64

    return length(f.data)
end

"""
    function temperature(
        f :: MatsubaraFunction
        ) :: Float64

Returns temperature for which `f.grids` are defined
"""
function temperature(
    f :: MatsubaraFunction
    ) :: Float64

    return temperature(grids(f, 1))
end

function Base.:copy(
    f :: MatsubaraFunction
    ) :: MatsubaraFunction

    return MatsubaraFunction(f)
end

"""
    function absmax(
        f :: MatsubaraFunction
        ) :: Float64

Returns largest element of `f.data` (in absolute terms)
"""
function absmax(
    f :: MatsubaraFunction
    ) :: Float64

    return maximum(abs.(f.data))
end

"""
    function argmax(
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

#----------------------------------------------------------------------------------------------#

include("matsubara_func_ops.jl")
include("matsubara_func_index.jl")
include("matsubara_func_eval.jl")
include("matsubara_func_symmetries.jl")
include("matsubara_func_io.jl")

#----------------------------------------------------------------------------------------------#

export 
    MatsubaraFunction,
    grids, 
    grids_shape,
    shape,
    data_shape,
    temperature,
    absmax,
    info