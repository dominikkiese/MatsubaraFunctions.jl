# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct MeshFunction{DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

MeshFunction type with fields:
* `meshes :: MT` : set of meshes
* `data   :: AT` : multidimensional data array
"""
struct MeshFunction{DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}
    meshes :: MT
    data   :: AT

    function MeshFunction(meshes :: MT, data :: AT) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}
        if Q <: Integer || Q <: Complex{Int} error("Integer data type not supported") end
        @DEBUG DD > 0 "Mesh dimension cannot be zero"
        return new{DD, Q, MT, AT}(meshes, data)
    end

    function MeshFunction(meshes :: Vararg{Mesh, DD}; data_t :: DataType = ComplexF64) where {DD}
        return MeshFunction(tuple(meshes...), Array{data_t, DD}(undef, length.(meshes)...))
    end

    function MeshFunction(f :: MeshFunction) :: MeshFunction
        return MeshFunction(meshes(f), copy(f.data))
    end
end

"""
    function meshes(f :: MeshFunction{DD, Q, MT, AT}) :: MT where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

Returns `f.meshes`
"""
function meshes(f :: MeshFunction{DD, Q, MT, AT}) :: MT where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}
    return f.meshes
end

"""
    function meshes(f :: MeshFunction{DD, Q, MT, AT}, :: Val{idx}) :: Mesh where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}, idx}

Returns `f.meshes[idx]`. Since idx is static the type of the returned mesh can be inferred at compile time.
"""
function meshes(f :: MeshFunction{DD, Q, MT, AT}, :: Val{idx}) :: Mesh where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}, idx}
    return f.meshes[idx]
end

"""
    function meshes(f :: MeshFunction{DD, Q, MT, AT}, idx :: Int) :: Mesh where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

Returns `f.meshes[idx]`. Since idx is only known at runtime the return mesh type is the union of mesh types in MT.
"""
function meshes(f :: MeshFunction{DD, Q, MT, AT}, idx :: Int) :: Mesh where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}
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