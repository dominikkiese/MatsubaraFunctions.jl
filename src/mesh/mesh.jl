include("mesh_point.jl")

# abstract types
#-------------------------------------------------------------------------------#

"""
    abstract type AbstractMesh

AbstractMesh type
"""
abstract type AbstractMesh end

"""
    abstract type AbstractDomain

AbstractDomain type
"""
abstract type AbstractDomain end

# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct Mesh{T <: AbstractMeshPoint, D <: AbstractDomain} <: AbstractMesh

Mesh type with fields:
* `hash   :: Symbol`    : mesh identifier
* `points :: Vector{T}` : mesh points
* `domain :: D`         : implementation details
"""
struct Mesh{T <: AbstractMeshPoint, D <: AbstractDomain} <: AbstractMesh
    hash   :: Symbol # no accessor, only for internal use
    points :: Vector{T}  
    domain :: D

    function Mesh(hash, points :: Vector{T}, domain :: D) where {T <: AbstractMeshPoint, D <: AbstractDomain}
        return new{T, D}(Symbol(hash), points, domain)
    end
end

"""
    function points(m :: Mesh{T, D}) :: Vector{T} where {T <: AbstractMeshPoint, D <: AbstractDomain}

Returns `m.points`
"""
function points(m :: Mesh{T, D}) :: Vector{T} where {T <: AbstractMeshPoint, D <: AbstractDomain}
    return m.points
end

"""
    function points(m :: Mesh{T, D}, idx :: Int) :: T where {T <: AbstractMeshPoint, D <: AbstractDomain}

Returns `m.points[idx]`
"""
function points(m :: Mesh{T, D}, idx :: Int) :: T where {T <: AbstractMeshPoint, D <: AbstractDomain}
    return m.points[idx]
end

"""
    function domain(m :: Mesh{T, D}) :: D where {T <: AbstractMeshPoint, D <: AbstractDomain}

Returns `m.domain`
"""
function domain(m :: Mesh{T, D}) :: D where {T <: AbstractMeshPoint, D <: AbstractDomain}
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

function Base.:getindex(m :: Mesh, idx :: Int)
    return points(m, idx)
end

function Base.:getindex(m :: Mesh, idxs :: UnitRange{Int})
    return @view points(m)[idxs]
end

# iterate
#-------------------------------------------------------------------------------#
 
function Base.:iterate(m :: Mesh) 
    return m[1], 1 
end 

function Base.:iterate(m :: Mesh, state :: Int)
    if state < length(m)
        return m[state + 1], state + 1 
    else 
        return nothing 
    end
end

# mapping to mesh index
#-------------------------------------------------------------------------------#

function mesh_index(x :: T, m :: Mesh{T, D}) where {T <: AbstractMeshPoint, D <: AbstractDomain}
    @DEBUG x.hash === m.hash "Mesh point invalid"
    return index(x)
end

# load mesh from file 
#-------------------------------------------------------------------------------#

"""
    function load_mesh(h :: HDF5.File, l :: String) :: Mesh

Load mesh with name `l` from file `h`
"""
function load_mesh(h :: HDF5.File, l :: String) :: Mesh
    # direct dispatch on the respective overload
    return load_mesh(h, l, Val(Symbol(read_attribute(h[l], "tag"))))
end

# print 
#-------------------------------------------------------------------------------#

function Base.:show(io :: IO, m :: Mesh{T, D}) where {T <: AbstractMeshPoint, D <: AbstractDomain}
    print(io, CYAN, BOLD, "Mesh ", RESET, "for ", CYAN, BOLD, "$(D) \n", RESET)
    print(io, "=> hash   : $(m.hash) \n")
    print(io, "=> length : $(length(m))")
    return nothing 
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
