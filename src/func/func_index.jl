# eachindex method
#----------------------------------------------------------------------------------------------#

function Base.:eachindex(f :: MeshFunction{DD, Q, MT, AT}
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}
    
    return eachindex(f.data)
end

# cartesian index
#----------------------------------------------------------------------------------------------#

function Base.:CartesianIndex(f :: MeshFunction{DD, Q, MT, AT}, x :: Vararg{Union{MeshPoint, <: AbstractValue}, DD}
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return CartesianIndex(_mesh_indices(f, x...)...)
end

function CartesianIndex_bc(f :: MeshFunction{DD, Q, MT, AT}, x :: Vararg{Union{MeshPoint, <: AbstractValue}, DD}
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return CartesianIndex(_mesh_indices_bc(f, x...)...)
end

# linear index
#----------------------------------------------------------------------------------------------#

"""
    function LinearIndex(f :: MeshFunction{DD, Q, MT, AT}, x :: Vararg{Union{MeshPoint, <: AbstractValue}, DD}
        ) :: Int where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

Returns linear index for access to `f.data`
"""
function LinearIndex(f :: MeshFunction{DD, Q, MT, AT}, x :: Vararg{Union{MeshPoint, <: AbstractValue}, DD}
    ) :: Int where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return LinearIndices(size(f.data))[_mesh_indices(f, x...)...]
end

"""
    function LinearIndex_bc(f :: MeshFunction{DD, Q, MT, AT}, x :: Vararg{Union{MeshPoint, <: AbstractValue}, DD}
        ) :: Int where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

Returns linear index for access to `f.data` under boundary conditions
"""
function LinearIndex_bc(f :: MeshFunction{DD, Q, MT, AT}, x :: Vararg{Union{MeshPoint, <: AbstractValue}, DD}
    ) :: Int where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return LinearIndices(size(f.data))[_mesh_indices_bc(f, x...)...]
end

"""
    function LinearIndex(f :: MeshFunction{DD, Q, MT, AT}, cidx :: CartesianIndex{DD}
        ) :: Int where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

Returns linear index for access to `f.data`
"""
function LinearIndex(f :: MeshFunction{DD, Q, MT, AT}, cidx :: CartesianIndex{DD}
    ) :: Int where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return LinearIndices(size(f.data))[cidx]
end

"""
    function LinearIndex(f :: MeshFunction{DD, Q, MT, AT}, x :: Vararg{Int, DD}
        ) :: Int where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

Returns linear index for access to `f.data`
"""
function LinearIndex(f :: MeshFunction{DD, Q, MT, AT}, x :: Vararg{Int, DD}
    ) :: Int where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return LinearIndices(size(f.data))[x...]
end

# to avoid ambiguities
function LinearIndex(:: MeshFunction{0, Q, MT, AT}) where {Q <: Number, MT <: NTuple{0, Mesh}, AT <: AbstractArray{Q, 0}}
    error("Data dimension of MeshFunction cannot be zero")
end

# conversion to meshes
#----------------------------------------------------------------------------------------------#

"""
    function to_meshes(f :: MeshFunction{DD, Q, MT, AT}, idx :: Int
        ) :: NTuple{DD, Union{MeshPoint}} where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

Returns mesh points
"""
function to_meshes(f :: MeshFunction{DD, Q, MT, AT}, idx :: Int
    ) :: NTuple{DD, Union{MeshPoint}} where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return to_meshes(f, CartesianIndices(size(f.data))[idx])
end

# getindex
#----------------------------------------------------------------------------------------------#

function Base.:getindex(f :: MeshFunction{DD, Q, MT, AT}, x :: Vararg{Union{MeshPoint, <: AbstractValue, Int, UnitRange, Colon}, DD}
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return f.data[_mesh_indices(f, x...)...]
end

function Base.:getindex(f :: MeshFunction{DD, Q, MT, AT}, x :: Vararg{Union{Int, UnitRange, Colon}, DD}
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return f.data[x...]
end

function Base.:getindex(f :: MeshFunction{DD, Q, MT, AT}, cidx :: CartesianIndex{DD}
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}
    
    return f.data[cidx]
end

function Base.:getindex(f :: MeshFunction, idx :: Int)
    return f.data[idx]
end

# to avoid ambiguities
function Base.:getindex(f :: MeshFunction{1, Q, MT, AT}, x :: Int) where {Q <: Number, MT <: NTuple{1, Mesh}, AT <: AbstractArray{Q, 1}}
    return f.data[x]
end

# views
#----------------------------------------------------------------------------------------------#

function Base.:view(f :: MeshFunction{DD, Q, MT, AT}, x :: Vararg{Union{MeshPoint, <: AbstractValue, Int, UnitRange, Colon}, DD}
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return view(f, _mesh_indices(f, x...)...)
end

function Base.:view(f :: MeshFunction{DD, Q, MT, AT}, x :: Vararg{Union{Int, UnitRange, Colon}, DD}
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return view(f.data, x...)
end

# setindex
#----------------------------------------------------------------------------------------------#

function Base.:setindex!(f :: MeshFunction{DD, Q, MT, AT}, val :: Qp, x :: Vararg{Union{MeshPoint, <: AbstractValue, Int, UnitRange, Colon}, DD}
    ) where {DD, Q <: Number, Qp <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    f.data[_mesh_indices(f, x...)...] = val
    return nothing
end

function Base.:setindex!(f :: MeshFunction{DD, Q, MT, AT}, val :: Qp, x :: Vararg{Union{Int, UnitRange, Colon}, DD}
    ) where {DD, Q <: Number, Qp <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    f.data[x...] = val
    return nothing
end

function Base.:setindex!(f :: MeshFunction{DD, Q, MT, AT}, val :: Qp, cidx :: CartesianIndex{DD}
    ) where {DD, Q <: Number, Qp <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    f.data[cidx] = val
    return nothing
end

function Base.:setindex!(f :: MeshFunction{DD, Q, MT, AT}, val :: Qp, idx :: Int
    ) where {DD, Q <: Number, Qp <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    f.data[idx] = val
    return nothing
end

# to avoid ambiguities
function Base.:setindex!(f :: MatsubaraFunctions.MeshFunction{1, Q, MT, AT}, val :: Qp, idx :: Int
    ) where {Q <: Number, Qp <: Number, MT <: NTuple{1, Mesh}, AT <: AbstractVector{Q}}

    f.data[idx] = val
    return nothing
end

# export
#----------------------------------------------------------------------------------------------#

export 
    LinearIndex,
    LinearIndex_bc,
    to_meshes