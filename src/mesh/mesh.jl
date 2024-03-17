include("mesh_point.jl")

# abstract types
#-------------------------------------------------------------------------------#

"""
    abstract type AbstractMesh

AbstractMesh type
"""
abstract type AbstractMesh end

# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct Mesh{T <: AbstractMeshPoint} <: AbstractMesh

Mesh type with fields:
* `hash   :: UInt64`    : mesh identifier
* `points :: Vector{T}` : mesh points
* `domain :: Dict`      : implementation details
"""
struct Mesh{T <: AbstractMeshPoint} <: AbstractMesh
    hash   :: UInt64 # no accessor, only for internal use
    points :: Vector{T}  
    domain :: Dict  

    function Mesh(
        hash   :: UInt64,
        points :: Vector{T},
        domain :: Dict  
        ) where {T <: AbstractMeshPoint}

        return new{T}(hash, points, domain)
    end

    # copy constructor
    function Mesh(m :: Mesh)
        return Mesh(m.hash, copy(points(m)), copy(domain(m)))
    end
end

"""
    function points(m :: Mesh{T}) :: Vector{T} where {T <: AbstractMeshPoint}

Returns `m.points`
"""
function points(m :: Mesh{T}) :: Vector{T} where {T <: AbstractMeshPoint}
    return m.points
end

"""
    function points(m :: Mesh{T}, idx :: Int64) :: T where {T <: AbstractMeshPoint}

Returns `m.points[idx]`
"""
function points(m :: Mesh{T}, idx :: Int64) :: T where {T <: AbstractMeshPoint}
    return m.points[idx]
end

"""
    function domain(m :: Mesh) :: Dict

Returns `m.domain`
"""
function domain(m :: Mesh) :: Dict
    return m.domain 
end

function Base.:length(m :: Mesh)
    return length(points(m))
end

# indexing
#-------------------------------------------------------------------------------#

function Base.:eachindex(m :: Mesh)
    return eachindex(points(m))
end

function Base.:firstindex(m :: Mesh)
    return firstindex(points(m))
end

function Base.:lastindex(m :: Mesh)
    return lastindex(points(m))
end

function Base.:getindex(m :: Mesh, idx :: Int64)
    return points(m, idx)
end

function Base.:getindex(m :: Mesh, idxs :: UnitRange{Int64})
    return @view points(m)[idxs]
end

function Base.:copy(m :: Mesh)
    return Mesh(m)
end

# iterate
#-------------------------------------------------------------------------------#
 
function Base.:iterate(m :: Mesh) 
    return m[1], 1 
end 

function Base.:iterate(m :: Mesh, state :: Int64)
    if state < length(m)
        return m[state + 1], state + 1 
    else 
        return nothing 
    end
end

# mapping to mesh index
#-------------------------------------------------------------------------------#

function mesh_index(x :: T, m :: Mesh{T}) where {T <: AbstractMeshPoint}
    @DEBUG x.hash == m.hash "Mesh point invalid"
    return index(x)
end

# load mesh from file 
#-------------------------------------------------------------------------------#

"""
    function load_mesh(h :: HDF5.File, l :: String) :: AbstractMesh

Load mesh with name `l` from file `h`
"""
function load_mesh(h :: HDF5.File, l :: String) :: AbstractMesh
    # direct dispatch on the respective overload
    return load_mesh(h, l, Val(Symbol(read_attribute(h[l], "tag"))))
end

# load implementations and export
#-------------------------------------------------------------------------------#

# for each value type the respective mesh must implement:
# - outer constructor
# - mapping from value type (and if needed plain_value type) to mesh index 
# - comparison operator
# boundary conditions are optional and only need to be implemented if meaningful

include("matsubara/matsubara_mesh.jl")
include("brillouin/brillouin_mesh.jl")
include("index/index_mesh.jl")

export
    AbstractMesh,
    Mesh,
    points,
    domain,
    load_mesh