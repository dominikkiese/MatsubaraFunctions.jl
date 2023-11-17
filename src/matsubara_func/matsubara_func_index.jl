function Base.:CartesianIndex(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    w :: NTuple{GD, AbstractMatsubaraFrequency},
    x :: Vararg{Int64, SD} 
    ) :: CartesianIndex{DD} where {GD, SD, DD, Q <: Number}

    idxs = ntuple(i -> grids(f, i)(w[i]), GD)
    @DEBUG all(ntuple(i -> 1 <= x[i] <= shape(f, i), SD)) "Tensor indices invalid, shape is $(shape(f))"
    return CartesianIndex(idxs..., x...)
end

function CartesianIndex_extrp(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    w :: NTuple{GD, AbstractMatsubaraFrequency},
    x :: Vararg{Int64, SD} 
    ) :: CartesianIndex{DD} where {GD, SD, DD, Q <: Number}

    idxs = ntuple(i -> grid_index_extrp(w[i], grids(f, i)), GD)
    @DEBUG all(ntuple(i -> firstindex(f.data, GD + i) <= x[i] <= lastindex(f.data, GD + i), SD)) "Tensor indices invalid, shape is $(shape(f))"

    return CartesianIndex(idxs..., x...)
end

function Base.:CartesianIndex(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    idx :: Int64
    )   :: CartesianIndex{DD} where {GD, SD, DD, Q <: Number}

    return CartesianIndices(data_shape(f))[idx]
end

#----------------------------------------------------------------------------------------------#

"""
    function LinearIndex(
        f :: MatsubaraFunction{GD, SD, DD, Q},
        w :: NTuple{GD, AbstractMatsubaraFrequency},
        x :: Vararg{Int64, SD} 
        ) :: Int64 where {GD, SD, DD, Q <: Number}

Returns linear index for access to `f.data`
"""
function LinearIndex(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    w :: NTuple{GD, AbstractMatsubaraFrequency},
    x :: Vararg{Int64, SD} 
    ) :: Int64 where {GD, SD, DD, Q <: Number}

    idxs = ntuple(i -> grids(f, i)(w[i]), GD)
    return LinearIndices(data_shape(f))[idxs..., x...]
end

"""
    function LinearIndex(
        f    :: MatsubaraFunction{GD, SD, DD, Q},
        cidx :: CartesianIndex{DD}
        )    :: Int64 where {GD, SD, DD, Q <: Number}

Returns linear index for access to `f.data`
"""
function LinearIndex(
    f    :: MatsubaraFunction{GD, SD, DD, Q},
    cidx :: CartesianIndex{DD}
    )    :: Int64 where {GD, SD, DD, Q <: Number}

    return LinearIndices(data_shape(f))[cidx]
end

"""
    function LinearIndex(
        f :: MatsubaraFunction{GD, SD, DD, Q},
        x :: Vararg{Int64, DD}
        ) :: Int64 where {GD, SD, DD, Q <: Number}

Returns linear index for access to `f.data`
"""
function LinearIndex(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    x :: Vararg{Int64, DD}
    ) :: Int64 where {GD, SD, DD, Q <: Number}

    return LinearIndices(data_shape(f))[x...]
end

#----------------------------------------------------------------------------------------------#

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

    return ntuple(i -> grids(f, i)[cidx[i]], GD), ntuple(i -> cidx[GD + i], SD)
end

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

#----------------------------------------------------------------------------------------------#

function grid_index(
    w :: Union{UnitRange, Colon, Base.IdentityUnitRange}
    ) :: Union{UnitRange, Colon, Base.IdentityUnitRange}

    return w
end

function Base.:getindex(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    w :: NTuple{GD, Union{AbstractMatsubaraFrequency, UnitRange, Colon}},
    x :: Vararg{Union{Int64, UnitRange, Colon}, SD} 
    ) :: Union{Q, AbstractArray{Q}} where {GD, SD, DD, Q <: Number}

    return f.data[ntuple(i -> grid_index(w[i]), GD)..., x...]
end

function Base.:getindex(
    f :: MatsubaraFunction{1, SD, DD, Q},
    w :: Union{AbstractMatsubaraFrequency, UnitRange, Colon},
    x :: Vararg{Union{Int64, UnitRange, Colon}, SD} 
    ) :: Union{Q, AbstractArray{Q}} where {SD, DD, Q <: Number}

    return f.data[grid_index(w), x...]
end

function Base.:getindex(
    f :: MatsubaraFunction{GD, 0, DD, Q},
    w :: Vararg{Union{AbstractMatsubaraFrequency, UnitRange, Colon}, GD} 
    ) :: Union{Q, AbstractArray{Q}} where {GD, DD, Q <: Number}

    return f.data[ntuple(i -> grid_index(w[i]), GD)...]
end

function Base.:getindex(
    f    :: MatsubaraFunction{GD, SD, DD, Q},
    cidx :: CartesianIndex{DD},
    )    :: Q where {GD, SD, DD, Q <: Number}

    return f.data[cidx]
end

function Base.:getindex(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    idx :: Int64,
    )   :: Q where {GD, SD, DD, Q <: Number}

    return f.data[idx]
end

function Base.:getindex(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    x :: Vararg{Union{Int64, UnitRange, Colon}, DD}
    ) :: Union{Q, AbstractArray{Q}} where {GD, SD, DD, Q <: Number}

    return f.data[x...]
end

# specialization
function Base.:getindex(
    f   :: MatsubaraFunction{GD, SD, 1, Q},
    idx :: Int64,
    )   :: Q where {GD, SD, Q <: Number}

    return f.data[idx]
end

function Base.:getindex(
    f :: MatsubaraFunction{1, 0, DD, Q},
    w :: Union{AbstractMatsubaraFrequency, UnitRange, Colon}
    ) :: Union{Q, AbstractArray{Q}} where {DD, Q <: Number}

    return f.data[grid_index(w)]
end

function Base.:getindex(
    f :: MatsubaraFunction{1, 0, 1, Q},
    w :: Union{UnitRange, Colon}
    ) :: Union{Q, AbstractArray{Q}} where {Q <: Number}

    return f.data[grid_index(w)]
end

#----------------------------------------------------------------------------------------------#

function Base.:view(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    w :: NTuple{GD, Union{AbstractMatsubaraFrequency, UnitRange, Colon, Base.IdentityUnitRange}},
    x :: Vararg{Union{Int64, UnitRange, Colon}, SD} 
    ) :: SubArray{Q} where {GD, SD, DD, Q <: Number}

    return view(f.data, ntuple(i -> grid_index(w[i]), GD)..., x...)
end

function Base.:view(
    f :: MatsubaraFunction{1, SD, DD, Q},
    w :: Union{AbstractMatsubaraFrequency, UnitRange, Colon, Base.IdentityUnitRange},
    x :: Vararg{Union{Int64, UnitRange, Colon}, SD} 
    ) :: SubArray{Q} where {SD, DD, Q <: Number}

    return view(f.data, grid_index(w), x...)
end

function Base.:view(
    f :: MatsubaraFunction{GD, 0, DD, Q},
    w :: Vararg{Union{AbstractMatsubaraFrequency, UnitRange, Colon, Base.IdentityUnitRange}, GD} 
    ) :: SubArray{Q} where {GD, DD, Q <: Number}

    return view(f.data, ntuple(i -> grid_index(w[i]), GD)...)
end

function Base.:view(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    x :: Vararg{Union{Int64, UnitRange, Colon, Base.IdentityUnitRange}, DD}
    ) :: SubArray{Q} where {GD, SD, DD, Q <: Number}

    return view(f.data, x...)
end

# specialization
function Base.:view(
    f :: MatsubaraFunction{1, 0, DD, Q},
    w :: Union{AbstractMatsubaraFrequency, UnitRange, Colon, Base.IdentityUnitRange},
    ) :: SubArray{Q} where {DD, Q <: Number}

    return view(f.data, grid_index(w))
end

function Base.:view(
    f :: MatsubaraFunction{1, 0, 1, Q},
    w :: Union{UnitRange, Colon, Base.IdentityUnitRange},
    ) :: SubArray{Q} where {Q <: Number}

    return view(f.data, grid_index(w))
end

#----------------------------------------------------------------------------------------------#

function Base.:setindex!(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    val :: Qp,
    w   :: NTuple{GD, Union{AbstractMatsubaraFrequency}},
    x   :: Vararg{Int64, SD} 
    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    f.data[ntuple(i -> grid_index(w[i]), GD)..., x...] = val
    return nothing
end

function Base.:setindex!(
    f   :: MatsubaraFunction{1, SD, DD, Q},
    val :: Qp,
    w   :: Union{AbstractMatsubaraFrequency},
    x   :: Vararg{Int64, SD} 
    )   :: Nothing where {SD, DD, Q <: Number, Qp <: Number}

    f.data[grid_index(w), x...] = val
    return nothing
end

function Base.:setindex!(
    f   :: MatsubaraFunction{GD, 0, DD, Q},
    val :: Qp,
    w   :: Vararg{Union{AbstractMatsubaraFrequency}, GD}
    )   :: Nothing where {GD, DD, Q <: Number, Qp <: Number}

    f.data[ntuple(i -> grid_index(w[i]), GD)...] = val
    return nothing
end

function Base.:setindex!(
    f    :: MatsubaraFunction{GD, SD, DD, Q},
    val  :: Qp,
    cidx :: CartesianIndex{DD},
    )    :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    f.data[cidx] = val
    return nothing
end

function Base.:setindex!(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    val :: Qp,
    idx :: Int64,
    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    f.data[idx] = val
    return nothing
end

function Base.:setindex!(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    val :: Qp,
    x   :: Vararg{Int64, DD}
    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    f.data[x...] = val
    return nothing
end

# specialization
function Base.:setindex!(
    f   :: MatsubaraFunction{1, 0, DD, Q},
    val :: Qp,
    w   :: Union{AbstractMatsubaraFrequency}
    )   :: Nothing where {DD, Q <: Number, Qp <: Number}

    f.data[grid_index(w)] = val
    return nothing
end

#----------------------------------------------------------------------------------------------#

export 
    LinearIndex,
    to_Matsubara