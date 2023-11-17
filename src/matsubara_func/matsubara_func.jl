"""
    struct MatsubaraFunction{GD, SD, DD, Q <: Number}

MatsubaraFunction type with fields:
* `grids  :: NTuple{GD, AbstractMatsubaraGrid}` : collection of Matsubara grids
* `shape  :: NTuple{SD, Int64}`                 : shape of tensor structure on every grid point
* `offset :: NTuple{DD, Int64}`                 : offsets for data access
* `data   :: OffsetArray{Q, DD, Array{Q, DD}}`  : multidimensional data array
"""
struct MatsubaraFunction{GD, SD, DD, Q <: Number}
    grids  :: NTuple{GD, AbstractMatsubaraGrid}
    shape  :: NTuple{SD, Int64}          
    offset :: NTuple{DD, Int64}          
    data   :: OffsetArray{Q, DD, Array{Q, DD}}

    function MatsubaraFunction(
        grids :: NTuple{GD, AbstractMatsubaraGrid}, 
        shape :: NTuple{SD, Int64},
        data  :: OffsetArray{Q, DD, Array{Q, DD}}
        )     :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

        if Q <: Integer || Q <: Complex{Int} error("Integer data type not supported") end
        
        # check dimensions
        @DEBUG GD > 0 "Grid dimension cannot be zero"
        @DEBUG GD + SD == DD "Dimensions do not match"
        
        # check grids
        T = temperature(grids[1])

        for i in eachindex(grids)
            @DEBUG temperature(grids[i]) ≈ T "Grids must have same temperature"
            @DEBUG issorted(values(grids[i])) "Grids must be sorted"
            @DEBUG first_index(grids[i]) == firstindex(data, i) "First index must agree between grids and data"
        end
        
        return new{GD, SD, DD, Q}(grids, shape, ntuple(i -> firstindex(data, i) - 1, DD), data)
    end

    function MatsubaraFunction(
        grids        :: NTuple{GD, AbstractMatsubaraGrid},
        shape        :: Vararg{Int64, SD}
        ;
        data_t       :: DataType         = ComplexF64,
        shape_offset :: NTuple{SD,Int64} = ntuple(i -> 0, SD)
        )            :: MatsubaraFunction{GD, SD, GD + SD, data_t} where {GD, SD}
        
        data = OffsetArray(Array{data_t, GD + SD}(undef, length.(grids)..., shape...), .-div.(length.(grids) .+ 2, 2)..., shape_offset...)
        return MatsubaraFunction(copy.(grids), ntuple(i -> shape[i], SD), data)
    end

    function MatsubaraFunction(
        grid         :: AbstractMatsubaraGrid,
        shape        :: Vararg{Int64, SD}
        ;
        data_t       :: DataType         = ComplexF64,
        shape_offset :: NTuple{SD,Int64} = ntuple(i -> 0, SD)
        )            :: MatsubaraFunction{1, SD, 1 + SD, data_t} where {SD}

        return MatsubaraFunction((grid,), shape...; data_t, shape_offset)
    end

    function MatsubaraFunction(
        grids  :: Vararg{AbstractMatsubaraGrid, GD},
        ;
        data_t :: DataType = ComplexF64
        )      :: MatsubaraFunction{GD, 0, GD, data_t} where {GD}

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
        )      :: MatsubaraFunction{1, 0, 1, data_t}

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
    function data_shape(f :: MatsubaraFunction{GD, SD, DD, Q}) where {GD, SD, DD, Q <: Number}

Returns shape of `f.data`
"""
function data_shape(f :: MatsubaraFunction{GD, SD, DD, Q}) where {GD, SD, DD, Q <: Number}
    return axes(f.data)
end

"""
    function data_shape(f :: MatsubaraFunction, idx :: Int64)

Returns length of dimension `idx` of `f.data`
"""
function data_shape(f :: MatsubaraFunction, idx :: Int64)
    return axes(f.data, idx)
end

function Base.size(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: NTuple{DD, Int64} where {GD, SD, DD, Q <: Number}
    
    return size(f.data)
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