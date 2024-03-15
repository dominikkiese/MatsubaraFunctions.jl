# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct MeshFunction{MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

MeshFunction type with fields:
* `meshes :: NTuple{MD, Union{<: AbstractMesh}}` : set of meshes
* `shape  :: NTuple{SD, Int64}`                  : index structure
* `data   :: AT`                                 : multidimensional data array
"""
struct MeshFunction{MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}
    meshes :: NTuple{MD, Union{<: AbstractMesh}}
    shape  :: NTuple{SD, Int64}          
    data   :: AT

    function MeshFunction( # pass by reference
        meshes :: NTuple{MD, Union{<: AbstractMesh}}, 
        shape  :: NTuple{SD, Int64}, 
        data   :: AT
        ) where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

        if Q <: Integer || Q <: Complex{Int} 
            error("Integer data type not supported") 
        end

        @DEBUG MD > 0 "Mesh dimension cannot be zero"
        @DEBUG MD + SD == DD "Data dimension incompatible with meshes and shape"
        return new{MD, SD, DD, Q, AT}(meshes, shape, data)
    end

    function MeshFunction( # pass by reference
        mesh  :: MT, 
        shape :: NTuple{SD, Int64}, 
        data  :: AT
        ) where {SD, DD, Q <: Number, MT <: AbstractMesh, AT <: AbstractArray{Q, DD}}

        return MeshFunction((mesh,), shape, data)
    end

    # pass by reference
    function MeshFunction(meshes :: NTuple{MD, Union{<: AbstractMesh}}, data :: AT) where {MD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}
        return MeshFunction(meshes, (), data)
    end

    # pass by reference
    function MeshFunction(mesh :: MT, data :: AT) where {DD, Q <: Number, MT <: AbstractMesh, AT <: AbstractArray{Q, DD}}
        return MeshFunction((mesh,), (), data)
    end

    function MeshFunction( # pass meshes by copy
        meshes :: NTuple{MD, Union{<: AbstractMesh}},
        shape  :: Vararg{Int64, SD}
        ;
        data_t :: DataType = ComplexF64
        ) where {MD, SD}
        
        data = Array{data_t, MD + SD}(undef, length.(meshes)..., shape...)
        return MeshFunction(Mesh.(meshes), tuple(shape...), data)
    end

    function MeshFunction( # pass mesh by copy
        mesh   :: MT,
        shape  :: Vararg{Int64, SD}
        ;
        data_t :: DataType = ComplexF64
        ) where {SD, MT <: AbstractMesh}

        return MeshFunction((mesh,), shape...; data_t)
    end

    # pass meshes by copy
    function MeshFunction(meshes :: Vararg{Union{<: AbstractMesh}, MD}; data_t :: DataType = ComplexF64) where {MD}
        return MeshFunction(tuple(meshes...); data_t)
    end

    # pass meshes and data by copy
    function MeshFunction(f :: MeshFunction) :: MeshFunction
        return MeshFunction(Mesh.(meshes(f)), shape(f), copy(f.data))
    end

    # specialization, pass mesh by copy
    function MeshFunction(mesh :: MT; data_t :: DataType = ComplexF64) where {MT <: AbstractMesh}
        return MeshFunction((mesh,); data_t)
    end
end

"""
    function meshes(f :: MeshFunction{MD, SD, DD, Q, AT}) :: NTuple{MD, Union{<: AbstractMesh}} where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

Returns `f.meshes`
"""
function meshes(f :: MeshFunction{MD, SD, DD, Q, AT}) :: NTuple{MD, Union{<: AbstractMesh}} where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}
    return f.meshes
end

"""
    function meshes(f :: MeshFunction, idx :: Int64) :: AbstractMesh

Returns `f.meshes[idx]`
"""
function meshes(f :: MeshFunction, idx :: Int64) :: AbstractMesh
    return f.meshes[idx]
end

"""
    function shape(f :: MeshFunction{MD, SD, DD, Q, AT}) :: NTuple{SD, Int64} where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

Returns `f.shape`
"""
function shape(f :: MeshFunction{MD, SD, DD, Q, AT}) :: NTuple{SD, Int64} where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}
    return f.shape 
end 

"""
    function shape(f :: MeshFunction, idx :: Int64) :: Int64

Returns `f.shape[idx]`
"""
function shape(f :: MeshFunction, idx :: Int64) :: Int64
    return f.shape[idx]
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

export 
    MeshFunction, 
    meshes,
    shape