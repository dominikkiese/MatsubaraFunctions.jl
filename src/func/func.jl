# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct MeshFunction{DD, Q <: Number, AT <: AbstractArray{Q, DD}}

MeshFunction type with fields:
* `meshes :: NTuple{DD, Union{<: AbstractMesh}}` : set of meshes
* `data   :: AT`                                 : multidimensional data array
"""
struct MeshFunction{DD, Q <: Number, AT <: AbstractArray{Q, DD}}
    meshes :: NTuple{DD, Union{<: AbstractMesh}}
    data   :: AT

    function MeshFunction(data :: AT, meshes :: Vararg{Union{<: AbstractMesh}, DD}) where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}
        if Q <: Integer || Q <: Complex{Int} error("Integer data type not supported") end
        @DEBUG DD > 0 "Mesh dimension cannot be zero"
        return new{DD, Q, AT}(tuple(meshes...), data)
    end

    function MeshFunction(meshes :: Vararg{Union{<: AbstractMesh}, DD}; data_t :: DataType = ComplexF64) where {DD}
        return MeshFunction(Array{data_t, DD}(undef, length.(meshes)...), meshes...)
    end

    function MeshFunction(f :: MeshFunction) :: MeshFunction
        return MeshFunction(copy(f.data), meshes(f)...)
    end
end

"""
    function meshes(f :: MeshFunction{DD, Q, AT}) :: NTuple{DD, Union{<: AbstractMesh}} where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

Returns `f.meshes`
"""
function meshes(f :: MeshFunction{DD, Q, AT}) :: NTuple{DD, Union{<: AbstractMesh}} where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}
    return f.meshes
end

"""
    function meshes(f :: MeshFunction, idx :: Int64) :: AbstractMesh

Returns `f.meshes[idx]`
"""
function meshes(f :: MeshFunction, idx :: Int64) :: AbstractMesh
    return f.meshes[idx]
end

function Base.:copy(f :: MeshFunction)
    return MeshFunction(f)
end

# load implementation details and export
#-------------------------------------------------------------------------------#

include("func_operations.jl")
include("func_index.jl")
include("func_itp.jl")
include("func_eval.jl")
include("func_io.jl")
include("func_symmetries.jl")

export 
    MeshFunction, 
    meshes