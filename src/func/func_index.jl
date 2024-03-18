# cartesian index
#----------------------------------------------------------------------------------------------#

function Base.:CartesianIndex(
    f :: MeshFunction{DD, Q, AT},
    p :: Vararg{Union{<: AbstractValue, <: AbstractMeshPoint}, DD}
    ) where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return CartesianIndex(map((y, m) -> mesh_index(y, m), p, meshes(f))...)
end

function CartesianIndex_bc(
    f :: MeshFunction{DD, Q, AT},
    p :: Vararg{Union{<: AbstractValue, <: AbstractMeshPoint}, DD}
    ) where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return CartesianIndex(map((y, m) -> mesh_index_bc(y, m), p, meshes(f))...)
end

function Base.:CartesianIndex(f :: MeshFunction{DD, Q, AT}, idx :: Int64) where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}
    return CartesianIndices(size(f.data))[idx]
end

# linear index
#----------------------------------------------------------------------------------------------#

"""
function LinearIndex(
    f :: MeshFunction{DD, Q, AT},
    p :: Vararg{Union{<: AbstractValue, <: AbstractMeshPoint}, DD}
    ) :: Int64 where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

Returns linear index for access to `f.data`
"""
function LinearIndex(
    f :: MeshFunction{DD, Q, AT},
    p :: Vararg{Union{<: AbstractValue, <: AbstractMeshPoint}, DD}
    ) :: Int64 where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return LinearIndices(size(f.data))[map((y, m) -> mesh_index(y, m), p, meshes(f))...]
end

"""
    function LinearIndex_bc(
        f :: MeshFunction{DD, Q, AT},
        p :: Vararg{Union{<: AbstractValue, <: AbstractMeshPoint}, DD}
        ) :: Int64 where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

Returns linear index for access to `f.data` under boundary conditions
"""
function LinearIndex_bc(
    f :: MeshFunction{DD, Q, AT},
    p :: Vararg{Union{<: AbstractValue, <: AbstractMeshPoint}, DD}
    ) :: Int64 where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return LinearIndices(size(f.data))[map((y, m) -> mesh_index_bc(y, m), p, meshes(f))...]
end

"""
    function LinearIndex(f :: MeshFunction{DD, Q, AT}, cidx :: CartesianIndex{DD}
        ) :: Int64 where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

Returns linear index for access to `f.data`
"""
function LinearIndex(f :: MeshFunction{DD, Q, AT}, cidx :: CartesianIndex{DD}
    ) :: Int64 where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return LinearIndices(size(f.data))[cidx]
end

"""
    function LinearIndex(f :: MeshFunction{DD, Q, AT}, x :: Vararg{Int64, DD}
        ) :: Int64 where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

Returns linear index for access to `f.data`
"""
function LinearIndex(f :: MeshFunction{DD, Q, AT}, x :: Vararg{Int64, DD}
    ) :: Int64 where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return LinearIndices(size(f.data))[x...]
end

# to avoid ambiguities:
function LinearIndex(::MatsubaraFunctions.MeshFunction{0, Q, AT}) where {Q<:Number, AT<:AbstractArray{Q, 0}}
    return nothing
end

# conversion to meshes
#----------------------------------------------------------------------------------------------#

"""
    function to_meshes(f :: MeshFunction{DD, Q, AT}, cidx :: CartesianIndex{DD}
        ) :: NTuple{DD, Union{<: AbstractMeshPoint}} where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

Returns mesh points and indices
"""
function to_meshes(f :: MeshFunction{DD, Q, AT}, cidx :: CartesianIndex{DD}
    ) :: NTuple{DD, Union{<: AbstractMeshPoint}} where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return ntuple(i -> meshes(f, i)[cidx[i]], DD)
end

"""
    function to_meshes(f :: MeshFunction{DD, Q, AT}, idx :: Int64
        ) :: NTuple{DD, Union{<: AbstractMeshPoint}} where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

Returns mesh points and indices
"""
function to_meshes(f :: MeshFunction{DD, Q, AT}, idx :: Int64
    ) :: NTuple{DD, Union{<: AbstractMeshPoint}} where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    cidx = CartesianIndex(f, idx)
    return to_meshes(f, cidx)
end

# getindex
#----------------------------------------------------------------------------------------------#

function Base.:getindex(
    f :: MeshFunction{DD, Q, AT},
    p :: Vararg{Union{<: AbstractValue, <: AbstractMeshPoint, UnitRange, Colon}, DD}
    ) where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return f.data[map((y, m) -> mesh_index(y, m), p, meshes(f))...]
end


function Base.:getindex(f :: MeshFunction{DD, Q, AT}, cidx :: CartesianIndex{DD}
    ) where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return f.data[cidx]
end

function Base.:getindex(f :: MeshFunction, idx :: Int64
    )

    return f.data[idx]
end

function Base.:getindex(f :: MeshFunction{DD, Q, AT}, x :: Vararg{Union{Int64, UnitRange, Colon}, DD}
    ) where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return f.data[x...]
end

function Base.:getindex(f :: MeshFunction{1, Q, AT}, x :: Int64
    ) where {Q <: Number, AT <: AbstractArray{Q, 1}}

    return f.data[x]
end

# views
#----------------------------------------------------------------------------------------------#

function Base.:view(
    f :: MeshFunction{DD, Q, AT},
    p :: Vararg{Union{<: AbstractValue, <: AbstractMeshPoint, UnitRange, Colon}, DD}
    ) where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return view(f.data, map((y, m) -> mesh_index(y, m), p, meshes(f))...)
end

#function Base.:view(f :: MeshFunction{MD, SD, DD, Q, AT}, x :: Vararg{Union{Int64, UnitRange, Colon}, DD}
#    ) where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}
#
#    return view(f.data, x...)
#end


# setindex
#----------------------------------------------------------------------------------------------#

function Base.:setindex!(
    f   :: MeshFunction{DD, Q, AT},
    val :: Qp,
    p   :: Vararg{Union{Int, <: AbstractValue, <: AbstractMeshPoint}, DD}
    ) where {DD, Q <: Number, Qp <: Number, AT <: AbstractArray{Q, DD}}

    f.data[map((y, m) -> mesh_index(y, m), p, meshes(f))...] = val
    return nothing
end

function Base.:setindex!(
    f    :: MeshFunction{DD, Q, AT},
    val  :: Qp,
    cidx :: CartesianIndex{DD},
    ) where {DD, Q <: Number, Qp <: Number, AT <: AbstractArray{Q, DD}}

    f.data[cidx] = val
    return nothing
end

function Base.:setindex!(
    f   :: MeshFunction{DD, Q, AT},
    val :: Qp,
    idx :: Int64,
    ) where {DD, Q <: Number, Qp <: Number, AT <: AbstractArray{Q, DD}}

    f.data[idx] = val
    return nothing
end

# to avoid ambiguities:
function Base.:setindex!(f::MatsubaraFunctions.MeshFunction{1, Q, AT}, val::Qp, idx::Int64) where {Q<:Number, Qp<:Number, AT<:AbstractVector{Q}}
    return f.data[idx] = val
end

#function Base.:setindex!(
#    f   :: MeshFunction{DD, Q, AT},
#    val :: Qp,
#    x   :: Vararg{Int64, DD}
#    ) where {DD, Q <: Number, Qp <: Number, AT <: AbstractArray{Q, DD}}
#
#    f.data[x...] = val
#    return nothing
#end

## to avoid ambiguities:
#function setindex!(::MatsubaraFunctions.MeshFunction{0, Q, AT}, ::Qp) where {Q<:Number, Qp<:Number, AT<:AbstractArray{Q, 0}}
#    return nothing
#end

# export
#----------------------------------------------------------------------------------------------#

export 
    LinearIndex,
    LinearIndex_bc,
    to_meshes