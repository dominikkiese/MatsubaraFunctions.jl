# cartesian index
#----------------------------------------------------------------------------------------------#

function Base.:CartesianIndex(
    f :: MeshFunction{MD, SD, DD, Q, AT},
    p :: NTuple{MD, Union{<: AbstractValue, <: AbstractMeshPoint}},
    x :: Vararg{Int64, SD} 
    ) where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    @DEBUG all(ntuple(i -> 1 <= x[i] <= shape(f, i), SD)) "Indices invalid"
    return CartesianIndex(map((y, m) -> mesh_index(y, m), p, meshes(f))..., x...)
end

function CartesianIndex_bc(
    f :: MeshFunction{MD, SD, DD, Q, AT},
    p :: NTuple{MD, Union{<: AbstractValue, <: AbstractMeshPoint}},
    x :: Vararg{Int64, SD} 
    ) where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    @DEBUG all(ntuple(i -> 1 <= x[i] <= shape(f, i), SD)) "Indices invalid"
    return CartesianIndex(map((y, m) -> mesh_index_bc(y, m), p, meshes(f))..., x...)
end

function Base.:CartesianIndex(f :: MeshFunction{MD, SD, DD, Q, AT}, idx :: Int64) where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}
    return CartesianIndices(size(f.data))[idx]
end

# linear index
#----------------------------------------------------------------------------------------------#

"""
    function LinearIndex(
        f :: MeshFunction{MD, SD, DD, Q, AT},
        p :: NTuple{MD, Union{<: AbstractValue, <: AbstractMeshPoint}},
        x :: Vararg{Int64, SD} 
        ) :: Int64 where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

Returns linear index for access to `f.data`
"""
function LinearIndex(
    f :: MeshFunction{MD, SD, DD, Q, AT},
    p :: NTuple{MD, Union{<: AbstractValue, <: AbstractMeshPoint}},
    x :: Vararg{Int64, SD} 
    ) :: Int64 where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return LinearIndices(size(f.data))[map((y, m) -> mesh_index(y, m), p, meshes(f))..., x...]
end

"""
    function LinearIndex_bc(
        f :: MeshFunction{MD, SD, DD, Q, AT},
        p :: NTuple{MD, Union{<: AbstractValue, <: AbstractMeshPoint}},
        x :: Vararg{Int64, SD} 
        ) :: Int64 where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

Returns linear index for access to `f.data` under boundary conditions
"""
function LinearIndex_bc(
    f :: MeshFunction{MD, SD, DD, Q, AT},
    p :: NTuple{MD, Union{<: AbstractValue, <: AbstractMeshPoint}},
    x :: Vararg{Int64, SD} 
    ) :: Int64 where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return LinearIndices(size(f.data))[map((y, m) -> mesh_index_bc(y, m), p, meshes(f))..., x...]
end

"""
    function LinearIndex(f :: MeshFunction{MD, SD, DD, Q, AT}, cidx :: CartesianIndex{DD}
        ) :: Int64 where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

Returns linear index for access to `f.data`
"""
function LinearIndex(f :: MeshFunction{MD, SD, DD, Q, AT}, cidx :: CartesianIndex{DD}
    ) :: Int64 where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return LinearIndices(size(f.data))[cidx]
end

"""
    function LinearIndex(f :: MeshFunction{MD, SD, DD, Q, AT}, x :: Vararg{Int64, DD}
        ) :: Int64 where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

Returns linear index for access to `f.data`
"""
function LinearIndex(f :: MeshFunction{MD, SD, DD, Q, AT}, x :: Vararg{Int64, DD}
    ) :: Int64 where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return LinearIndices(size(f.data))[x...]
end

# conversion to meshes
#----------------------------------------------------------------------------------------------#

"""
    function to_meshes(f :: MeshFunction{MD, SD, DD, Q, AT}, cidx :: CartesianIndex{DD}
        ) :: Tuple{NTuple{MD, Union{<: AbstractMeshPoint}}, NTuple{SD, Int64}} where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

Returns mesh points and indices
"""
function to_meshes(f :: MeshFunction{MD, SD, DD, Q, AT}, cidx :: CartesianIndex{DD}
    ) :: Tuple{NTuple{MD, Union{<: AbstractMeshPoint}}, NTuple{SD, Int64}} where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return ntuple(i -> meshes(f, i)[cidx[i]], MD), ntuple(i -> cidx[MD + i], SD)
end

"""
    function to_meshes(f :: MeshFunction{MD, SD, DD, Q, AT}, idx :: Int64
        ) :: Tuple{NTuple{MD, Union{<: AbstractMeshPoint}}, NTuple{SD, Int64}} where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

Returns mesh points and indices
"""
function to_meshes(f :: MeshFunction{MD, SD, DD, Q, AT}, idx :: Int64
    ) :: Tuple{NTuple{MD, Union{<: AbstractMeshPoint}}, NTuple{SD, Int64}} where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    cidx = CartesianIndex(f, idx)
    return to_meshes(f, cidx)
end

# getindex
#----------------------------------------------------------------------------------------------#

function Base.:getindex(
    f :: MeshFunction{MD, SD, DD, Q, AT},
    p :: NTuple{MD, Union{<: AbstractValue, <: AbstractMeshPoint, UnitRange, Colon}},
    x :: Vararg{Union{Int64, UnitRange, Colon}, SD} 
    ) where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return f.data[map((y, m) -> mesh_index(y, m), p, meshes(f))..., x...]
end

function Base.:getindex(
    f :: MeshFunction{1, SD, DD, Q, AT},
    p :: Union{<: AbstractValue, <: AbstractMeshPoint, UnitRange, Colon},
    x :: Vararg{Union{Int64, UnitRange, Colon}, SD} 
    ) where {SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return f.data[mesh_index(p, meshes(f, 1)), x...]
end

function Base.:getindex(f :: MeshFunction{MD, 0, DD, Q, AT}, p :: Vararg{Union{<: AbstractValue, <: AbstractMeshPoint, UnitRange, Colon}, MD}
    ) where {MD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return f.data[map((y, m) -> mesh_index(y, m), p, meshes(f))...]
end

function Base.:getindex(f :: MeshFunction{MD, SD, DD, Q, AT}, cidx :: CartesianIndex{DD}
    ) where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return f.data[cidx]
end

function Base.:getindex(f :: MeshFunction{MD, SD, DD, Q, AT}, idx :: Int64
    ) where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return f.data[idx]
end

function Base.:getindex(f :: MeshFunction{MD, SD, DD, Q, AT}, x :: Vararg{Union{Int64, UnitRange, Colon}, DD}
    ) where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return f.data[x...]
end

# specialization
function Base.:getindex(f :: MeshFunction{MD, SD, 1, Q, AT}, idx :: Int64) where {MD, SD, Q <: Number, AT <: AbstractArray{Q, 1}}
    return f.data[idx]
end

function Base.:getindex(f :: MeshFunction{1, 0, DD, Q, AT}, p :: Union{<: AbstractValue, <: AbstractMeshPoint, UnitRange, Colon}
    ) where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return f.data[mesh_index(p, meshes(f, 1))]
end

function Base.:getindex(f :: MeshFunction{1, 0, 1, Q, AT}, p :: Union{UnitRange, Colon}
    ) where {Q <: Number, AT <: AbstractArray{Q, 1}}

    return f.data[p]
end

# views
#----------------------------------------------------------------------------------------------#

function Base.:view(
    f :: MeshFunction{MD, SD, DD, Q, AT},
    p :: NTuple{MD, Union{<: AbstractValue, <: AbstractMeshPoint, UnitRange, Colon}},
    x :: Vararg{Union{Int64, UnitRange, Colon}, SD} 
    ) where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return view(f.data, map((y, m) -> mesh_index(y, m), p, meshes(f))..., x...)
end

function Base.:view(
    f :: MeshFunction{1, SD, DD, Q, AT},
    p :: Union{<: AbstractValue, <: AbstractMeshPoint, UnitRange, Colon},
    x :: Vararg{Union{Int64, UnitRange, Colon}, SD} 
    ) where {SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return view(f.data, mesh_index(p, meshes(f, 1)), x...)
end

function Base.:view(f :: MeshFunction{MD, 0, DD, Q, AT}, p :: Vararg{Union{<: AbstractValue, <: AbstractMeshPoint, UnitRange, Colon}, MD}
    ) where {MD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return view(f.data, map((y, m) -> mesh_index(y, m), p, meshes(f))...)
end

function Base.:view(f :: MeshFunction{MD, SD, DD, Q, AT}, x :: Vararg{Union{Int64, UnitRange, Colon}, DD}
    ) where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return view(f.data, x...)
end

# specialization
function Base.:view(f :: MeshFunction{1, 0, DD, Q, AT}, p :: Union{<: AbstractValue, <: AbstractMeshPoint, UnitRange, Colon}
    ) where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return view(f.data, mesh_index(p, meshes(f, 1)))
end

function Base.:view(f :: MeshFunction{1, 0, 1, Q, AT}, p :: Union{UnitRange, Colon}) where {Q <: Number, AT <: AbstractArray{Q, 1}}
    return view(f.data, p)
end

# setindex
#----------------------------------------------------------------------------------------------#

function Base.:setindex!(
    f   :: MeshFunction{MD, SD, DD, Q, AT},
    val :: Qp,
    p   :: NTuple{MD, Union{<: AbstractValue, <: AbstractMeshPoint}},
    x   :: Vararg{Int64, SD} 
    ) where {MD, SD, DD, Q <: Number, Qp <: Number, AT <: AbstractArray{Q, DD}}

    f.data[map((y, m) -> mesh_index(y, m), p, meshes(f))..., x...] = val
    return nothing
end

function Base.:setindex!(
    f   :: MeshFunction{1, SD, DD, Q, AT},
    val :: Qp,
    p   :: Union{<: AbstractValue, <: AbstractMeshPoint},
    x   :: Vararg{Int64, SD} 
    ) where {SD, DD, Q <: Number, Qp <: Number, AT <: AbstractArray{Q, DD}}

    f.data[mesh_index(p, meshes(f, 1)), x...] = val
    return nothing
end

function Base.:setindex!(
    f   :: MeshFunction{MD, 0, DD, Q, AT},
    val :: Qp,
    p   :: Vararg{Union{<: AbstractValue, <: AbstractMeshPoint}, MD}
    ) where {MD, DD, Q <: Number, Qp <: Number, AT <: AbstractArray{Q, DD}}

    f.data[map((y, m) -> mesh_index(y, m), p, meshes(f))...] = val
    return nothing
end

function Base.:setindex!(
    f    :: MeshFunction{MD, SD, DD, Q, AT},
    val  :: Qp,
    cidx :: CartesianIndex{DD},
    ) where {MD, SD, DD, Q <: Number, Qp <: Number, AT <: AbstractArray{Q, DD}}

    f.data[cidx] = val
    return nothing
end

function Base.:setindex!(
    f   :: MeshFunction{MD, SD, DD, Q, AT},
    val :: Qp,
    idx :: Int64,
    ) where {MD, SD, DD, Q <: Number, Qp <: Number, AT <: AbstractArray{Q, DD}}

    f.data[idx] = val
    return nothing
end

function Base.:setindex!(
    f   :: MeshFunction{MD, SD, DD, Q, AT},
    val :: Qp,
    x   :: Vararg{Int64, DD}
    ) where {MD, SD, DD, Q <: Number, Qp <: Number, AT <: AbstractArray{Q, DD}}

    f.data[x...] = val
    return nothing
end

# specialization
function Base.:setindex!(
    f   :: MeshFunction{1, 0, DD, Q, AT},
    val :: Qp,
    p   :: Union{<: AbstractValue, <: AbstractMeshPoint}
    ) where {DD, Q <: Number, Qp <: Number, AT <: AbstractArray{Q, DD}}

    f.data[mesh_index(p, meshes(f, 1))] = val
    return nothing
end

# export
#----------------------------------------------------------------------------------------------#

export 
    LinearIndex,
    LinearIndex_bc,
    to_meshes