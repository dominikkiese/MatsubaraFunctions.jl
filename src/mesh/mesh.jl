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
        )      :: Mesh{T} where {T <: AbstractMeshPoint}

        return new{T}(hash, points, domain)
    end

    # copy constructor
    function Mesh(m :: Mesh) :: Mesh
        return Mesh(m.hash, copy(points(m)), copy(domain(m)))
    end
end

"""
    function points(
        m :: Mesh{T}
        ) :: Vector{T} where {T <: AbstractMeshPoint}

Returns `m.points`
"""
function points(
    m :: Mesh{T}
    ) :: Vector{T} where {T <: AbstractMeshPoint}

    return m.points
end

"""
    function points(
        m   :: Mesh{T},
        idx :: Int64
        )   :: T where {T <: AbstractMeshPoint}

Returns `m.points[idx]`
"""
function points(
    m   :: Mesh{T},
    idx :: Int64
    )   :: T where {T <: AbstractMeshPoint}

    return m.points[idx]
end

"""
    function domain(
        m :: AbstractMesh
        ) :: Dict

Returns `m.domain`
"""
function domain(
    m :: AbstractMesh
    ) :: Dict

    return m.domain 
end

function Base.:length(
    m :: AbstractMesh
    ) :: Int64

    return length(points(m))
end

# indexing
#-------------------------------------------------------------------------------#

function Base.:eachindex(
    m :: AbstractMesh
    ) :: Base.OneTo{Int64}

    return eachindex(points(m))
end

function Base.:firstindex(
    m :: AbstractMesh
    ) :: Int64

    return firstindex(points(m))
end

function Base.:lastindex(
    m :: AbstractMesh
    ) :: Int64

    return lastindex(points(m))
end

function Base.:getindex(
    m   :: Mesh{T},
    idx :: Int64
    )   :: T where {T <: AbstractMeshPoint}

    return points(m, idx)
end

function Base.:getindex(
    m    :: Mesh{T},
    idxs :: UnitRange{Int64}
    )    :: SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int64}}, true} where {T <: AbstractMeshPoint}

    return @view points(m)[idxs]
end

function Base.:copy(m :: Mesh) :: Mesh
    return Mesh(m)
end

# iterate
#-------------------------------------------------------------------------------#
 
function Base.:iterate(
    m :: Mesh{T}
    ) :: Tuple{T, Int64} where {T <: AbstractMeshPoint}

    return m[1], 1 
end 

function Base.:iterate(
    m     :: Mesh{T},
    state :: Int64
    )     :: Union{Nothing, Tuple{T, Int64}} where {T <: AbstractMeshPoint}

    if state < length(m)
        return m[state + 1], state + 1 
    else 
        return nothing 
    end
end

# test if mesh point is valid
#-------------------------------------------------------------------------------#

# check if the mesh point is valid
function is_valid(
    m :: Mesh{T}, 
    x :: T
    ) :: Bool where {T <: AbstractMeshPoint}

    return x.hash == m.hash
end

# get index of mesh point if valid
function index_valid(
    m :: Mesh{T}, 
    x :: T
    ) :: Int64 where {T <: AbstractMeshPoint}

    @DEBUG is_valid(m, x) "Mesh point invalid"
    return index(x)
end

# load implementations and export
#-------------------------------------------------------------------------------#

# for each value type the respective mesh must implement:
# - outer constructor
# - comparison operator
# - method to call mesh with value type for mapping to mesh index

include("matsubara/matsubara_freq.jl")
include("matsubara/matsubara_mesh.jl")
include("brillouin/brillouin_pt.jl")
include("brillouin/brillouin_zone.jl")
include("brillouin/brillouin_mesh.jl")

export
    AbstractMesh,
    Mesh,
    points,
    domain