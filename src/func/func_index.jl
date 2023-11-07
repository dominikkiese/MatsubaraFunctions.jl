# cartesian index
#----------------------------------------------------------------------------------------------#

function Base.:CartesianIndex(
    f :: MeshFunction{MD, SD, DD, Q},
    p :: NTuple{MD, AbstractMeshPoint},
    x :: Vararg{Int64, SD} 
    ) :: CartesianIndex{DD} where {MD, SD, DD, Q <: Number}

    @DEBUG all(ntuple(i -> 1 <= x[i] <= shape(f, i), SD)) "Indices invalid"
    return CartesianIndex(ntuple(i -> index_valid(meshes(f, i), p[i]), MD)..., x...)
end

function Base.:CartesianIndex(
    f   :: MeshFunction{MD, SD, DD, Q},
    idx :: Int64
    )   :: CartesianIndex{DD} where {MD, SD, DD, Q <: Number}

    return CartesianIndices(size(f.data))[idx]
end

# linear index
#----------------------------------------------------------------------------------------------#

"""
    function LinearIndex(
        f :: MeshFunction{MD, SD, DD, Q},
        p :: NTuple{MD, AbstractMeshPoint},
        x :: Vararg{Int64, SD} 
        ) :: Int64 where {MD, SD, DD, Q <: Number}

Returns linear index for access to `f.data`
"""
function LinearIndex(
    f :: MeshFunction{MD, SD, DD, Q},
    p :: NTuple{MD, AbstractMeshPoint},
    x :: Vararg{Int64, SD} 
    ) :: Int64 where {MD, SD, DD, Q <: Number}

    return LinearIndices(size(f.data))[ntuple(i -> index_valid(meshes(f, i), p[i]), MD)..., x...]
end

"""
    function LinearIndex(
        f    :: MeshFunction{MD, SD, DD, Q},
        cidx :: CartesianIndex{DD}
        )    :: Int64 where {MD, SD, DD, Q <: Number}

Returns linear index for access to `f.data`
"""
function LinearIndex(
    f    :: MeshFunction{MD, SD, DD, Q},
    cidx :: CartesianIndex{DD}
    )    :: Int64 where {MD, SD, DD, Q <: Number}

    return LinearIndices(size(f.data))[cidx]
end

"""
    function LinearIndex(
        f :: MeshFunction{MD, SD, DD, Q},
        x :: Vararg{Int64, DD}
        ) :: Int64 where {MD, SD, DD, Q <: Number}

Returns linear index for access to `f.data`
"""
function LinearIndex(
    f :: MeshFunction{MD, SD, DD, Q},
    x :: Vararg{Int64, DD}
    ) :: Int64 where {MD, SD, DD, Q <: Number}

    return LinearIndices(size(f.data))[x...]
end

# conversion to meshes
#----------------------------------------------------------------------------------------------#

"""
    function to_meshes(
        f    :: MeshFunction{MD, SD, DD, Q},
        cidx :: CartesianIndex{DD}
        )    :: Tuple{NTuple{MD, AbstractMeshPoint}, NTuple{SD, Int64}} where {MD, SD, DD, Q <: Number}

Returns mesh points and indices
"""
function to_meshes(
    f    :: MeshFunction{MD, SD, DD, Q},
    cidx :: CartesianIndex{DD}
    )    :: Tuple{NTuple{MD, AbstractMeshPoint}, NTuple{SD, Int64}} where {MD, SD, DD, Q <: Number}

    return ntuple(i -> meshes(f, i)[cidx[i]], MD), ntuple(i -> cidx[MD + i], SD)
end

"""
    function to_meshes(
        f   :: MeshFunction{MD, SD, DD, Q},
        idx :: Int64 
        )   :: Tuple{NTuple{MD, AbstractMeshPoint}, NTuple{SD, Int64}} where {MD, SD, DD, Q <: Number}

Returns mesh points and indices
"""
function to_meshes(
    f   :: MeshFunction{MD, SD, DD, Q},
    idx :: Int64 
    )   :: Tuple{NTuple{MD, AbstractMeshPoint}, NTuple{SD, Int64}} where {MD, SD, DD, Q <: Number}

    cidx = CartesianIndex(f, idx)
    return to_meshes(f, cidx)
end

# getindex
#----------------------------------------------------------------------------------------------#

function index_valid(
    m :: AbstractMesh,
    p :: Union{UnitRange, Colon}
    ) :: Union{UnitRange, Colon}

    return p
end

function Base.:getindex(
    f :: MeshFunction{MD, SD, DD, Q},
    p :: NTuple{MD, Union{AbstractMeshPoint, UnitRange, Colon}},
    x :: Vararg{Union{Int64, UnitRange, Colon}, SD} 
    ) :: Union{Q, AbstractArray{Q}} where {MD, SD, DD, Q <: Number}

    return f.data[ntuple(i -> index_valid(meshes(f, i), p[i]), MD)..., x...]
end

function Base.:getindex(
    f :: MeshFunction{1, SD, DD, Q},
    p :: Union{AbstractMeshPoint, UnitRange, Colon},
    x :: Vararg{Union{Int64, UnitRange, Colon}, SD} 
    ) :: Union{Q, AbstractArray{Q}} where {SD, DD, Q <: Number}

    return f.data[index_valid(meshes(f, 1), p), x...]
end

function Base.:getindex(
    f :: MeshFunction{MD, 0, DD, Q},
    p :: Vararg{Union{AbstractMeshPoint, UnitRange, Colon}, MD} 
    ) :: Union{Q, AbstractArray{Q}} where {MD, DD, Q <: Number}

    return f.data[ntuple(i -> index_valid(meshes(f, i), p[i]), MD)...]
end

function Base.:getindex(
    f    :: MeshFunction{MD, SD, DD, Q},
    cidx :: CartesianIndex{DD},
    )    :: Q where {MD, SD, DD, Q <: Number}

    return f.data[cidx]
end

function Base.:getindex(
    f   :: MeshFunction{MD, SD, DD, Q},
    idx :: Int64,
    )   :: Q where {MD, SD, DD, Q <: Number}

    return f.data[idx]
end

function Base.:getindex(
    f :: MeshFunction{MD, SD, DD, Q},
    x :: Vararg{Union{Int64, UnitRange, Colon}, DD}
    ) :: Union{Q, AbstractArray{Q}} where {MD, SD, DD, Q <: Number}

    return f.data[x...]
end

# specialization
function Base.:getindex(
    f   :: MeshFunction{MD, SD, 1, Q},
    idx :: Int64,
    )   :: Q where {MD, SD, Q <: Number}

    return f.data[idx]
end

function Base.:getindex(
    f :: MeshFunction{1, 0, DD, Q},
    p :: Union{AbstractMeshPoint, UnitRange, Colon}
    ) :: Union{Q, AbstractArray{Q}} where {DD, Q <: Number}

    return f.data[index_valid(meshes(f, 1), p)]
end

function Base.:getindex(
    f :: MeshFunction{1, 0, 1, Q},
    p :: Union{UnitRange, Colon}
    ) :: Union{Q, AbstractArray{Q}} where {Q <: Number}

    return f.data[p]
end

# views
#----------------------------------------------------------------------------------------------#

function Base.:view(
    f :: MeshFunction{MD, SD, DD, Q},
    p :: NTuple{MD, Union{AbstractMeshPoint, UnitRange, Colon}},
    x :: Vararg{Union{Int64, UnitRange, Colon}, SD} 
    ) :: SubArray{Q} where {MD, SD, DD, Q <: Number}

    return view(f.data, ntuple(i -> index_valid(meshes(f, i), p[i]), MD)..., x...)
end

function Base.:view(
    f :: MeshFunction{1, SD, DD, Q},
    p :: Union{AbstractMeshPoint, UnitRange, Colon},
    x :: Vararg{Union{Int64, UnitRange, Colon}, SD} 
    ) :: SubArray{Q} where {SD, DD, Q <: Number}

    return view(f.data, index_valid(meshes(f, 1), p), x...)
end

function Base.:view(
    f :: MeshFunction{MD, 0, DD, Q},
    p :: Vararg{Union{AbstractMeshPoint, UnitRange, Colon}, MD} 
    ) :: SubArray{Q} where {MD, DD, Q <: Number}

    return view(f.data, ntuple(i -> index_valid(meshes(f, i), p[i]), MD)...)
end

function Base.:view(
    f :: MeshFunction{MD, SD, DD, Q},
    x :: Vararg{Union{Int64, UnitRange, Colon}, DD}
    ) :: SubArray{Q} where {MD, SD, DD, Q <: Number}

    return view(f.data, x...)
end

# specialization
function Base.:view(
    f :: MeshFunction{1, 0, DD, Q},
    p :: Union{AbstractMeshPoint, UnitRange, Colon},
    ) :: SubArray{Q} where {DD, Q <: Number}

    return view(f.data, index_valid(meshes(f, 1), p))
end

function Base.:view(
    f :: MeshFunction{1, 0, 1, Q},
    p :: Union{UnitRange, Colon},
    ) :: SubArray{Q} where {Q <: Number}

    return view(f.data, p)
end

# setindex
#----------------------------------------------------------------------------------------------#

function Base.:setindex!(
    f   :: MeshFunction{MD, SD, DD, Q},
    val :: Qp,
    p   :: NTuple{MD, Union{AbstractMeshPoint}},
    x   :: Vararg{Int64, SD} 
    )   :: Nothing where {MD, SD, DD, Q <: Number, Qp <: Number}

    f.data[ntuple(i -> index_valid(meshes(f, i), p[i]), MD)..., x...] = val
    return nothing
end

function Base.:setindex!(
    f   :: MeshFunction{1, SD, DD, Q},
    val :: Qp,
    p   :: Union{AbstractMeshPoint},
    x   :: Vararg{Int64, SD} 
    )   :: Nothing where {SD, DD, Q <: Number, Qp <: Number}

    f.data[index_valid(meshes(f, 1), p), x...] = val
    return nothing
end

function Base.:setindex!(
    f   :: MeshFunction{MD, 0, DD, Q},
    val :: Qp,
    p   :: Vararg{Union{AbstractMeshPoint}, MD}
    )   :: Nothing where {MD, DD, Q <: Number, Qp <: Number}

    f.data[ntuple(i -> index_valid(meshes(f, i), p[i]), MD)...] = val
    return nothing
end

function Base.:setindex!(
    f    :: MeshFunction{MD, SD, DD, Q},
    val  :: Qp,
    cidx :: CartesianIndex{DD},
    )    :: Nothing where {MD, SD, DD, Q <: Number, Qp <: Number}

    f.data[cidx] = val
    return nothing
end

function Base.:setindex!(
    f   :: MeshFunction{MD, SD, DD, Q},
    val :: Qp,
    idx :: Int64,
    )   :: Nothing where {MD, SD, DD, Q <: Number, Qp <: Number}

    f.data[idx] = val
    return nothing
end

function Base.:setindex!(
    f   :: MeshFunction{MD, SD, DD, Q},
    val :: Qp,
    x   :: Vararg{Int64, DD}
    )   :: Nothing where {MD, SD, DD, Q <: Number, Qp <: Number}

    f.data[x...] = val
    return nothing
end

# specialization
function Base.:setindex!(
    f   :: MeshFunction{1, 0, DD, Q},
    val :: Qp,
    p   :: Union{AbstractMeshPoint}
    )   :: Nothing where {DD, Q <: Number, Qp <: Number}

    f.data[index_valid(meshes(f, 1), p)] = val
    return nothing
end

# export
#----------------------------------------------------------------------------------------------#

export 
    LinearIndex,
    to_meshes