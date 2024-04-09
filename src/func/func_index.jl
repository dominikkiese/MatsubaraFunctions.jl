# cartesian index
#----------------------------------------------------------------------------------------------#

function Base.:CartesianIndex(f :: MeshFunction{DD, Q, MT, AT}, x :: Vararg{Union{MeshPoint, <: AbstractValue}, DD}
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return CartesianIndex(ntuple(i -> mesh_index(x[i], meshes(f, i)), DD)...)
end

function CartesianIndex_bc(f :: MeshFunction{DD, Q, MT, AT}, x :: Vararg{Union{MeshPoint, <: AbstractValue}, DD}
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return CartesianIndex(ntuple(i -> mesh_index_bc(x[i], meshes(f, i)), DD)...)
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

    return LinearIndices(size(f.data))[ntuple(i -> mesh_index(x[i], meshes(f, i)), DD)...]
end

"""
    function LinearIndex_bc(f :: MeshFunction{DD, Q, MT, AT}, x :: Vararg{Union{MeshPoint, <: AbstractValue}, DD}
        ) :: Int where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

Returns linear index for access to `f.data` under boundary conditions
"""
function LinearIndex_bc(f :: MeshFunction{DD, Q, MT, AT}, x :: Vararg{Union{MeshPoint, <: AbstractValue}, DD}
    ) :: Int where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return LinearIndices(size(f.data))[ntuple(i -> mesh_index_bc(x[i], meshes(f, i)), DD)...]
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
    function to_meshes(f :: MeshFunction{DD, Q, MT, AT}, cidx :: CartesianIndex{DD}
        ) :: NTuple{DD, Union{MeshPoint}} where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

Returns mesh points
"""
function to_meshes(f :: MeshFunction{DD, Q, MT, AT}, cidx :: CartesianIndex{DD}
    ) :: NTuple{DD, Union{MeshPoint}} where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return ntuple(i -> meshes(f, i)[cidx[i]], DD)
end

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

#function Base.:getindex(f :: MeshFunction{DD, Q, MT, AT}, x :: Vararg{Union{MeshPoint, <: AbstractValue, Int, UnitRange, Colon}, DD}
#    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}
#
#    return f.data[ntuple(i -> mesh_index(x[i], meshes(f, i)), DD)...]
#end

@generated function Base.:getindex(f :: MeshFunction{DD, Q, MT, AT}, x :: Vararg{Union{MeshPoint, <: AbstractValue, Int, UnitRange, Colon}, DD}
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}
    return Meta.parse("f.data["*prod(["mesh_index(x[$i], meshes(f, $i)), " for i in 1:DD])[1:end-2]*"]")
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

    return view(f.data, ntuple(i -> mesh_index(x[i], meshes(f, i)), DD)...)
end

function Base.:view(f :: MeshFunction{DD, Q, MT, AT}, x :: Vararg{Union{Int, UnitRange, Colon}, DD}
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return view(f.data, x...)
end

# setindex
#----------------------------------------------------------------------------------------------#

function Base.:setindex!(f :: MeshFunction{DD, Q, MT, AT}, val :: Qp, x :: Vararg{Union{MeshPoint, <: AbstractValue, Int, UnitRange, Colon}, DD}
    ) where {DD, Q <: Number, Qp <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    f.data[ntuple(i -> mesh_index(x[i], meshes(f, i)), DD)...] = val
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