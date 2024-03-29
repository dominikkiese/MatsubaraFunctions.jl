# CartesianIndex from MatsubaraFrequency or MatsubaraIndex and sites
function Base.:CartesianIndex(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    w :: NTuple{GD, Union{MatsubaraFrequency, MatsubaraIndex}},
    x :: Vararg{Int64, SD} 
    ) :: CartesianIndex{DD} where {GD, SD, DD, Q <: Number}

    # calling grid with frequency or index performs inbounds check
    idxs = ntuple(i -> grids(f, i)(w[i]), GD)

    # check if tensor indices are inbounds 
    @check any(ntuple(i -> !(1 <= x[i] <= shape(f, i)), SD)) == false "Tensor indices invalid, shape is $(shape(f))"

    return CartesianIndex(idxs..., x...)
end

# extrapolated CartesianIndex from MatsubaraFrequency or MatsubaraIndex and sites
# (i.e. frequencies which are out of bounds will be reset to mesh boundaries)
function CartesianIndex_extrp(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    w :: NTuple{GD, Union{MatsubaraFrequency, MatsubaraIndex}},
    x :: Vararg{Int64, SD} 
    ) :: CartesianIndex{DD} where {GD, SD, DD, Q <: Number}

    # determine grid indices
    idxs = ntuple(i -> grid_index_extrp(w[i], grids(f, i)), GD)

    # check if tensor indices are inbounds 
    @check any(ntuple(i -> !(1 <= x[i] <= shape(f, i)), SD)) == false "Tensor indices invalid, shape is $(shape(f))"

    return CartesianIndex(idxs..., x...)
end

# CartesianIndex from LinearIndex
function Base.:CartesianIndex(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    idx :: Int64
    )   :: CartesianIndex{DD} where {GD, SD, DD, Q <: Number}

    # bounds check performed by Base
    return CartesianIndices(data_shape(f))[idx]
end

"""
    function LinearIndex(
        f :: MatsubaraFunction{GD, SD, DD, Q},
        w :: NTuple{GD, Union{MatsubaraFrequency, MatsubaraIndex}},
        x :: Vararg{Int64, SD} 
        ) :: Int64 where {GD, SD, DD, Q <: Number}

Returns linear index for access to `f.data`
"""
function LinearIndex(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    w :: NTuple{GD, Union{MatsubaraFrequency, MatsubaraIndex}},
    x :: Vararg{Int64, SD} 
    ) :: Int64 where {GD, SD, DD, Q <: Number}

    # calling grid with frequency performs inbounds check
    idxs = ntuple(i -> grids(f, i)(w[i]), GD)

    # bounds check for x performed by Base
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

    # bounds check performed by Base
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

    # bounds check performed by Base
    return LinearIndices(data_shape(f))[x...]
end

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

# overload to generalize getindex methods
function grid_index(
    w    :: Union{UnitRange, Colon},
    grid :: MatsubaraGrid
    )    :: Union{UnitRange, Colon}

    return w
end

# getindex methods
function Base.:getindex(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    w :: NTuple{GD, Union{MatsubaraFrequency, MatsubaraIndex, UnitRange, Colon}},
    x :: Vararg{Union{Int64, UnitRange, Colon}, SD} 
    ) :: Union{Q, AbstractArray{Q}} where {GD, SD, DD, Q <: Number}

    # bounds check performed by Base
    return f.data[ntuple(i -> grid_index(w[i], grids(f, i)), GD)..., x...]
end

function Base.:getindex(
    f :: MatsubaraFunction{1, SD, DD, Q},
    w :: Union{MatsubaraFrequency, MatsubaraIndex, UnitRange, Colon},
    x :: Vararg{Union{Int64, UnitRange, Colon}, SD} 
    ) :: Union{Q, AbstractArray{Q}} where {SD, DD, Q <: Number}

    # bounds check performed by Base
    return f.data[grid_index(w, grids(f, 1)), x...]
end

function Base.:getindex(
    f :: MatsubaraFunction{GD, 1, DD, Q},
    w :: Vararg{Union{MatsubaraFrequency, MatsubaraIndex, UnitRange, Colon}, GD} 
    ) :: Union{Q, AbstractArray{Q}} where {GD, DD, Q <: Number}

    # bounds check performed by Base
    @check shape(f, 1) == 1 "MatsubaraFunction is not scalar but vector valued!"
    return f.data[ntuple(i -> grid_index(w[i], grids(f, i)), GD)..., 1]
end

function Base.:getindex(
    f    :: MatsubaraFunction{GD, SD, DD, Q},
    cidx :: CartesianIndex{DD},
    )    :: Q where {GD, SD, DD, Q <: Number}

    # bounds check perfomed by Base
    return f.data[cidx]
end

function Base.:getindex(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    idx :: Int64,
    )   :: Q where {GD, SD, DD, Q <: Number}

    # performed by Base
    return f.data[idx]
end

function Base.:getindex(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    x :: Vararg{Union{Int64, UnitRange, Colon}, DD}
    ) :: Union{Q, AbstractArray{Q}} where {GD, SD, DD, Q <: Number}

    # bounds check performed by Base
    return f.data[x...]
end

# specialization to avoid method ambiguity 
function Base.:getindex(
    :: MatsubaraFunction{GD, SD, 1, Q}, 
    :: Int64
    ) where {GD, SD, Q <: Number}

    error("The type MatsubaraFunction{GD, SD, 1, Q} is not supported, data dimension cannot be smaller than 2") 
end

# views
function Base.:view(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    w :: NTuple{GD, Union{MatsubaraFrequency, MatsubaraIndex, UnitRange, Colon}},
    x :: Vararg{Union{Int64, UnitRange, Colon}, SD} 
    ) :: SubArray{Q} where {GD, SD, DD, Q <: Number}

    # bounds check performed by Base
    return view(f.data, ntuple(i -> grid_index(w[i], grids(f, i)), GD)..., x...)
end

function Base.:view(
    f :: MatsubaraFunction{1, SD, DD, Q},
    w :: Union{MatsubaraFrequency, MatsubaraIndex, UnitRange, Colon},
    x :: Vararg{Union{Int64, UnitRange, Colon}, SD} 
    ) :: SubArray{Q} where {SD, DD, Q <: Number}

    # bounds check performed by Base
    return view(f.data, grid_index(w, grids(f, 1)), x...)
end

function Base.:view(
    f :: MatsubaraFunction{GD, 1, DD, Q},
    w :: Vararg{Union{MatsubaraFrequency, MatsubaraIndex, UnitRange, Colon}, GD} 
    ) :: SubArray{Q} where {GD, DD, Q <: Number}

    # bounds check performed by Base
    @check shape(f, 1) == 1 "MatsubaraFunction is not scalar but vector valued!"
    return view(f.data, ntuple(i -> grid_index(w[i], grids(f, i)), GD)..., 1)
end

function Base.:view(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    x :: Vararg{Union{Int64, UnitRange, Colon}, DD}
    ) :: SubArray{Q} where {GD, SD, DD, Q <: Number}

    # bounds check performed by Base
    return view(f.data, x...)
end

# setindex! methods
function Base.:setindex!(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    val :: Qp,
    w   :: NTuple{GD, Union{MatsubaraFrequency, MatsubaraIndex}},
    x   :: Vararg{Int64, SD} 
    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    # bounds check performed by Base
    f.data[ntuple(i -> grid_index(w[i], grids(f, i)), GD)..., x...] = val

    return nothing
end

function Base.:setindex!(
    f   :: MatsubaraFunction{1, SD, DD, Q},
    val :: Qp,
    w   :: Union{MatsubaraFrequency, MatsubaraIndex},
    x   :: Vararg{Int64, SD} 
    )   :: Nothing where {SD, DD, Q <: Number, Qp <: Number}

    # bounds check performed by Base
    f.data[grid_index(w, grids(f, 1)), x...] = val

    return nothing
end

function Base.:setindex!(
    f   :: MatsubaraFunction{GD, 1, DD, Q},
    val :: Qp,
    w   :: Vararg{Union{MatsubaraFrequency, MatsubaraIndex}, GD}
    )   :: Nothing where {GD, DD, Q <: Number, Qp <: Number}

    # bounds check performed by Base
    @check shape(f, 1) == 1 "MatsubaraFunction is not scalar but vector valued"
    f.data[ntuple(i -> grid_index(w[i], grids(f, i)), GD)..., 1] = val

    return nothing
end

function Base.:setindex!(
    f    :: MatsubaraFunction{GD, SD, DD, Q},
    val  :: Qp,
    cidx :: CartesianIndex{DD},
    )    :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    # bounds check performed by Base
    f.data[cidx] = val

    return nothing
end

function Base.:setindex!(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    val :: Qp,
    idx :: Int64,
    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    # bounds check performed by Base
    f.data[idx] = val

    return nothing
end

function Base.:setindex!(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    val :: Qp,
    x   :: Vararg{Int64, DD}
    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    # bounds check performed by Base
    f.data[x...] = val

    return nothing
end