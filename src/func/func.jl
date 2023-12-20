# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct MeshFunction{MD, SD, DD, Q <: Number}

MeshFunction type with fields:
* `meshes :: NTuple{MD, AbstractMesh}` : set of meshes
* `shape  :: NTuple{SD, Int64}`        : index structure
* `data   :: AbstractArray{Q, DD}`     : multidimensional data array
"""
struct MeshFunction{MD, SD, DD, Q <: Number}
    meshes :: NTuple{MD, AbstractMesh}
    shape  :: NTuple{SD, Int64}          
    data   :: AbstractArray{Q, DD}

    function MeshFunction( # parse by reference
        meshes :: NTuple{MD, AbstractMesh}, 
        shape  :: NTuple{SD, Int64}, 
        data   :: AbstractArray{Q, DD}
        )      :: MeshFunction{MD, SD, DD, Q} where {MD, SD, DD, Q <: Number}

        if Q <: Integer || Q <: Complex{Int} 
            error("Integer data type not supported") 
        end

        @DEBUG MD > 0 "Mesh dimension cannot be zero"
        @DEBUG MD + SD == DD "Data dimension incompatible with meshes and shape"
        return new{MD, SD, DD, Q}(meshes, shape, data)
    end

    function MeshFunction( # parse by reference
        mesh  :: AbstractMesh, 
        shape :: NTuple{SD, Int64}, 
        data  :: AbstractArray{Q, DD}
        )     :: MeshFunction{1, SD, DD, Q} where {SD, DD, Q <: Number}

        return MeshFunction((mesh,), shape, data)
    end

    function MeshFunction( # parse by reference
        meshes :: NTuple{MD, AbstractMesh}, 
        data   :: AbstractArray{Q, DD}
        )      :: MeshFunction{MD, 0, DD, Q} where {MD, DD, Q <: Number}

        return MeshFunction(meshes, (), data)
    end

    function MeshFunction( # parse by reference
        mesh  :: AbstractMesh, 
        data  :: AbstractArray{Q, DD}
        )     :: MeshFunction{1, 0, DD, Q} where {DD, Q <: Number}

        return MeshFunction((mesh,), (), data)
    end

    function MeshFunction( # parse meshes by copy
        meshes :: NTuple{MD, AbstractMesh},
        shape  :: Vararg{Int64, SD}
        ;
        data_t :: DataType = ComplexF64
        )      :: MeshFunction{MD, SD, MD + SD, data_t} where {MD, SD}
        
        data = Array{data_t, MD + SD}(undef, length.(meshes)..., shape...)
        return MeshFunction(Mesh.(meshes), tuple(shape...), data)
    end

    function MeshFunction( # parse mesh by copy
        mesh   :: AbstractMesh,
        shape  :: Vararg{Int64, SD}
        ;
        data_t :: DataType = ComplexF64
        )      :: MeshFunction{1, SD, 1 + SD, data_t} where {SD}

        return MeshFunction((mesh,), shape...; data_t)
    end

    function MeshFunction( # parse meshes by copy
        meshes :: Vararg{AbstractMesh, MD},
        ;
        data_t :: DataType = ComplexF64
        )      :: MeshFunction{MD, 0, MD, data_t} where {MD}

        return MeshFunction(tuple(meshes...); data_t)
    end

    function MeshFunction( # parse meshes and data by copy
        f :: MeshFunction
        ) :: MeshFunction

        return MeshFunction(Mesh.(meshes(f)), shape(f), copy(f.data))
    end

    # specialization
    function MeshFunction( # parse mesh by copy
        mesh   :: AbstractMesh
        ;
        data_t :: DataType = ComplexF64
        )      :: MeshFunction{1, 0, 1, data_t}

        return MeshFunction((mesh,); data_t)
    end
end

"""
    function meshes(
        f :: MeshFunction{MD, SD, DD, Q}
        ) :: NTuple{MD, AbstractMesh} where {MD, SD, DD, Q <: Number}

Returns `f.meshes`
"""
function meshes(
    f :: MeshFunction{MD, SD, DD, Q}
    ) :: NTuple{MD, AbstractMesh} where {MD, SD, DD, Q <: Number}

    return f.meshes
end

"""
    function meshes(
        f   :: MeshFunction,
        idx :: Int64
        )   :: AbstractMesh

Returns `f.meshes[idx]`
"""
function meshes(
    f   :: MeshFunction,
    idx :: Int64
    )   :: AbstractMesh

    return f.meshes[idx]
end

"""
    function shape(
        f :: MeshFunction{MD, SD, DD, Q}
        ) :: NTuple{SD, Int64} where {MD, SD, DD, Q <: Number}

Returns `f.shape`
"""
function shape(
    f :: MeshFunction{MD, SD, DD, Q}
    ) :: NTuple{SD, Int64} where {MD, SD, DD, Q <: Number}

    return f.shape 
end 

"""
    function shape(
        f   :: MeshFunction,
        idx :: Int64
        )   :: Int64

Returns `f.shape[idx]`
"""
function shape(
    f   :: MeshFunction,
    idx :: Int64
    )   :: Int64

    return f.shape[idx]
end 

function Base.:copy(
    f :: MeshFunction
    ) :: MeshFunction

    return MeshFunction(f)
end

# load implementation details and export
#-------------------------------------------------------------------------------#

include("func_operations.jl")
include("func_index.jl")
include("func_itp.jl")
include("func_eval.jl")

export 
    MeshFunction, 
    meshes,
    shape