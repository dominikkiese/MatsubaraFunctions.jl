# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct MeshFunction{MD, SD, DD, Q <: Number}

MeshFunction type with fields:
* `meshes :: NTuple{MD, AbstractMesh}` : set of meshes
* `shape  :: NTuple{SD, Int64}`        : index structure
* `data   :: Array{Q, DD}`             : multidimensional data array
"""
struct MeshFunction{MD, SD, DD, Q <: Number}
    meshes :: NTuple{MD, AbstractMesh}
    shape  :: NTuple{SD, Int64}          
    offset :: NTuple{DD, Int64}          
    data   :: OffsetArray{Q, DD, Array{Q, DD}}

    function MeshFunction(
        meshes :: NTuple{MD, AbstractMesh}, 
        shape  :: NTuple{SD, Int64}, 
        data   :: OffsetArray{Q, DD, Array{Q, DD}}
        )      :: MeshFunction{MD, SD, DD, Q} where {MD, SD, DD, Q <: Number}

        if Q <: Integer || Q <: Complex{Int} 
            error("Integer data type not supported") 
        end

        @DEBUG MD > 0 "Mesh dimension cannot be zero"
        @DEBUG MD + SD == DD "Data dimension incompatible with meshes and shape"
        return new{MD, SD, DD, Q}(meshes, shape, ntuple(i -> firstindex(data, i) - 1, DD), data)
    end

    function MeshFunction(
        meshes :: NTuple{MD, AbstractMesh},
        shape  :: Vararg{Int64, SD}
        ;
        data_t :: DataType = ComplexF64,
        shape_offset :: NTuple{SD,Int64} = ntuple(i -> 0, SD)
        )      :: MeshFunction{MD, SD, MD + SD, data_t} where {MD, SD}
        
        #data = Array{data_t, MD + SD}(undef, length.(meshes)..., shape...)
        data = OffsetArray(Array{data_t, MD + SD}(undef, length.(meshes)..., shape...), (firstindex.(meshes) .- 1) ..., shape_offset...)
        return MeshFunction(Mesh.(meshes), ntuple(i -> shape[i], SD), data)
    end

    function MeshFunction(
        mesh   :: AbstractMesh,
        shape  :: Vararg{Int64, SD}
        ;
        data_t :: DataType = ComplexF64,
        shape_offset :: NTuple{SD,Int64} = ntuple(i -> 0, SD)
        )      :: MeshFunction{1, SD, 1 + SD, data_t} where {SD}

        return MeshFunction((mesh,), shape...; data_t, shape_offset)
    end

    function MeshFunction(
        meshes :: Vararg{AbstractMesh, MD},
        ;
        data_t :: DataType = ComplexF64
        )      :: MeshFunction{MD, 0, MD, data_t} where {MD}

        return MeshFunction(ntuple(i -> meshes[i], MD); data_t)
    end

    function MeshFunction(
        f :: MeshFunction
        ) :: MeshFunction

        return MeshFunction(meshes(f), shape(f), copy(f.data))
    end

    # specialization
    function MeshFunction(
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

"""
    function axes(f :: MeshFunction)

Returns a tuple of valid index ranges for `f.data`
"""
function Base.:axes(f :: MeshFunction)
    return axes(f.data)
end

"""
    function axes(f :: MeshFunction, idx :: Int64)

Returns the range of valid indices along dimension `idx` of `f.data`
"""
function Base.:axes(f :: MeshFunction, idx :: Int64)
    return axes(f.data, idx)
end

function Base.size(
    f :: MeshFunction{GD, SD, DD, Q}
    ) :: NTuple{DD, Int64} where {GD, SD, DD, Q <: Number}
    
    return size(f.data)
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
include("func_symmetries.jl")
include("func_io.jl")

export 
    MeshFunction, 
    meshes,
    shape