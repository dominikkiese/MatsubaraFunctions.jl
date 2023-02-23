#== 
    Indexing for MatsubaraFunctions:
        -> CartesianIndex
        -> LinearIndex
        -> getindex 
        -> setindex!
==#

# CartesianIndex from MatsubaraFrequency and sites
function Base.CartesianIndex(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    w :: NTuple{GD, MatsubaraFrequency},
    x :: Vararg{Int64, SD} 
    ) :: CartesianIndex{DD} where {GD, SD, DD, Q <: Number}

    idxs = ntuple(i -> f.grids[i](w[i]), GD)
    return CartesianIndex(idxs..., x...)
end

# CartesianIndex from LinearIndex
function Base.CartesianIndex(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    idx :: Int64
    )   :: CartesianIndex{DD} where {GD, SD, DD, Q <: Number}

    # bounds check performed by Base
    return CartesianIndices(data_shape(f))[idx]
end



# LinearIndex from MatsubaraFrequency and sites
"""
    function LinearIndex(
        f :: MatsubaraFunction{GD, SD, DD, Q},
        w :: NTuple{GD, MatsubaraFrequency},
        x :: Vararg{Int64, SD} 
        ) :: Int64 where {GD, SD, DD, Q <: Number}

Returns linear index for access to f.data
"""
function LinearIndex(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    w :: NTuple{GD, MatsubaraFrequency},
    x :: Vararg{Int64, SD} 
    ) :: Int64 where {GD, SD, DD, Q <: Number}

    idxs = ntuple(i -> f.grids[i](w[i]), GD)

    # bounds check for x performed by Base
    return LinearIndices(data_shape(f))[idxs..., x...]
end

# LinearIndex from CartesianIndex
"""
    function LinearIndex(
        f    :: MatsubaraFunction{GD, SD, DD, Q},
        cidx :: CartesianIndex{DD}
        )    :: Int64 where {GD, SD, DD, Q <: Number}

Returns linear index for access to f.data
"""
function LinearIndex(
    f    :: MatsubaraFunction{GD, SD, DD, Q},
    cidx :: CartesianIndex{DD}
    )    :: Int64 where {GD, SD, DD, Q <: Number}

    # bounds check performed by Base
    return LinearIndices(data_shape(f))[cidx]
end

# LinearIndex from Vararg
"""
    function LinearIndex(
        f :: MatsubaraFunction{GD, SD, DD, Q},
        x :: Vararg{Int64, DD}
        ) :: Int64 where {GD, SD, DD, Q <: Number}

Returns linear index for access to f.data    
"""
function LinearIndex(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    x :: Vararg{Int64, DD}
    ) :: Int64 where {GD, SD, DD, Q <: Number}

    # bounds check performed by Base
    return LinearIndices(data_shape(f))[x...]
end



# Conversion from CartesianIndex back to MatsubaraFrequency + sites
"""
    function to_Matsubara(
        f    :: MatsubaraFunction{GD, SD, DD, Q},
        cidx :: CartesianIndex{DD}
        )    :: Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}} where {GD, SD, DD, Q <: Number}

Returns coordinates in grids and index of tensor structure
"""
function to_Matsubara(
    f    :: MatsubaraFunction{GD, SD, DD, Q},
    cidx :: CartesianIndex{DD}
    )    :: Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}} where {GD, SD, DD, Q <: Number}

    # bounds check performed by Base
    return ntuple(i -> f.grids[i][cidx[i]], GD), ntuple(i -> cidx[GD + i], SD)
end

# Conversion from LinearIndex back to MatsubaraFrequency + sites
"""
    function to_Matsubara(
        f   :: MatsubaraFunction{GD, SD, DD, Q},
        idx :: Int64 
        )   :: Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}} where {GD, SD, DD, Q <: Number}

Returns coordinates in grids and index of tensor structure
"""
function to_Matsubara(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    idx :: Int64 
    )   :: Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}} where {GD, SD, DD, Q <: Number}

    cidx = CartesianIndex(f, idx)
    return to_Matsubara(f, cidx)
end



# getindex from MatsubaraFrequency + sites
function Base.:getindex(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    w :: NTuple{GD, MatsubaraFrequency},
    x :: Vararg{Int64, SD} 
    ) :: Q where {GD, SD, DD, Q <: Number}

    # bounds check perfomed by Base
    return f.data[CartesianIndex(f, w, x...)]
end

function Base.:getindex(
    f :: MatsubaraFunction{1, SD, DD, Q},
    w :: MatsubaraFrequency,
    x :: Vararg{Int64, SD} 
    ) :: Q where {SD, DD, Q <: Number}

    # bounds check perfomed by Base
    return f.data[CartesianIndex(f, (w,), x...)]
end

# getindex from CartesianIndex
function Base.:getindex(
    f    :: MatsubaraFunction{GD, SD, DD, Q},
    cidx :: CartesianIndex{DD},
    )    :: Q where {GD, SD, DD, Q <: Number}

    # bounds check perfomed by Base
    return f.data[cidx]
end

# getindex from LinearIndex
function Base.:getindex(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    idx :: Int64,
    )   :: Q where {GD, SD, DD, Q <: Number}

    # performed by Base
    return f.data[idx]
end

# getindex from Vararg
function Base.:getindex(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    x :: Vararg{Int64, DD}
    ) :: Q where {GD, SD, DD, Q <: Number}

    # bounds check performed by Base
    return f.data[x...]
end



# setindex! from MatsubaraFrequency + sites
function Base.:setindex!(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    val :: Qp,
    w   :: NTuple{GD, MatsubaraFrequency},
    x   :: Vararg{Int64, SD} 
    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    # bounds check performed by Base
    f.data[CartesianIndex(f, w, x...)] = val

    return nothing
end

function Base.:setindex!(
    f   :: MatsubaraFunction{1, SD, DD, Q},
    val :: Qp,
    w   :: MatsubaraFrequency,
    x   :: Vararg{Int64, SD} 
    )   :: Nothing where {SD, DD, Q <: Number, Qp <: Number}

    # bounds check performed by Base
    f.data[CartesianIndex(f, (w,), x...)] = val

    return nothing
end

# setindex! from CartesianIndex
function Base.:setindex!(
    f    :: MatsubaraFunction{GD, SD, DD, Q},
    val  :: Qp,
    cidx :: CartesianIndex{DD},
    )    :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    # bounds check performed by Base
    f.data[cidx] = val

    return nothing
end

# setindex! from LinearIndex
function Base.:setindex!(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    val :: Qp,
    idx :: Int64,
    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    # bounds check performed by Base
    f.data[idx] = val

    return nothing
end

# setindex! from Vararg
function Base.:setindex!(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    val :: Qp,
    x   :: Vararg{Int64, DD}
    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    # bounds check performed by Base
    f.data[x...] = val

    return nothing
end