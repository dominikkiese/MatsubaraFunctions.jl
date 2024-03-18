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

    function MeshFunction( # pass by reference
        meshes :: NTuple{DD, Union{<: AbstractMesh}}, 
        data   :: AT
        ) where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

        if Q <: Integer || Q <: Complex{Int} 
            error("Integer data type not supported") 
        end

        @DEBUG DD > 0 "Mesh dimension cannot be zero"
        return new{DD, Q, AT}(meshes, data)
    end

    function MeshFunction(
        mesh  :: MT, 
        data  :: AT
        ) where {DD, Q <: Number, MT <: AbstractMesh, AT <: AbstractArray{Q, DD}}

        return MeshFunction((mesh,), data)
    end

    function MeshFunction( # pass meshes by copy
        meshes :: NTuple{DD, Union{<: AbstractMesh}}
        ;
        data_t :: DataType = ComplexF64
        ) where {DD}
        
        data = Array{data_t, DD}(undef, length.(meshes)...)
        return MeshFunction(Mesh.(meshes), data)
    end

    # pass meshes by copy
    function MeshFunction(meshes :: Vararg{Union{<: AbstractMesh}, DD}; data_t :: DataType = ComplexF64) where {DD}
        return MeshFunction(tuple(meshes...); data_t)
    end

    function MeshFunction(f :: MeshFunction) :: MeshFunction
        return MeshFunction(meshes(f), copy(f.data))
    end

    function MeshFunction(mesh :: MT; data_t :: DataType = ComplexF64) where {MT <: AbstractMesh}
        return MeshFunction((mesh,); data_t)
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