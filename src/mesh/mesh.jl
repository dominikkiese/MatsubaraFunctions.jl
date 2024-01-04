include("mesh_point.jl")


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
    points :: OffsetVector{T,Vector{T}}
    domain :: Dict  

    function Mesh(
        hash   :: UInt64,
        points :: OffsetVector{T,Vector{T}},
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
    ) :: OffsetVector{T,Vector{T}} where {T <: AbstractMeshPoint}

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
    ) #:: Base.OneTo{Int64}

    return eachindex(points(m))
end

"""
    function firstindex(
        grid :: AbstractMesh
        )    :: Int64

Returns the index of the first Matsubara frequency in grid
"""
function Base.:firstindex(
    m :: AbstractMesh
    ) :: Int64

    return firstindex(points(m))
end

"""
    function lastindex(
        grid :: AbstractMesh
        )    :: Int64

Returns the index of the last Matsubara frequency in grid
"""
function Base.:lastindex(
    m :: AbstractMesh
    ) :: Int64

    return lastindex(points(m))
end

"""
    function axes(grid :: AbstractMesh)

Returns range of valid indices for Mesh
"""
function Base.:axes(grid :: AbstractMesh)
    return first(axes(grid.points))
end

"""
    function firstvalue(
        grid :: AbstractMesh
        )    :: Float64

Returns the value of the first Matsubara frequency in grid
"""
function firstvalue(
    grid :: AbstractMesh
    )    :: AbstractValue

    return value(first(grid.points))
end

"""
    function lastvalue(
        grid :: AbstractMesh
        )    :: AbstractValue

Returns the value of the last Matsubara frequency in grid
"""
function lastvalue(
    grid :: AbstractMesh
    )    :: AbstractValue
    
    return value(last(grid.points))
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
    )    :: SubArray{T, 1, OffsetVector{T,Vector{T}}, Tuple{UnitRange{Int64}}, true} where {T <: AbstractMeshPoint}

    return @view points(m)[idxs]
end

function Base.:copy(m :: Mesh) :: Mesh
    return Mesh(m)
end


#----------------------------------------------------------------------------------------------#

"""
    function is_inbounds(
        w :: MatsubaraFrequency{PT},
        m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
        ) :: Bool where {PT <: AbstractParticle}

Checks if Matsubara frequency in mesh
"""
function is_inbounds(
    p :: AbstractMeshPoint,
    m :: Mesh{T}
    ) :: Bool where {T}

    return is_inbounds(value(p), m)
end

# returns index to data array corresponding to this frequency if in grid
function (m :: AbstractMesh)(
    p :: Union{AbstractValue,AbstractMeshPoint}
    ) :: Int64

    if is_inbounds(p, m)
        return mesh_index(p, m)
    else 
        error("Point is not in mesh")
    end 
end

# returns index to data array corresponding to closest frequency if in grid
function (m :: AbstractMesh)(
    p :: Float64
    ) :: Int64

    if is_inbounds(p, m)
        delta    = plain_value(m[2]) - plain_value(m[1])
        position = (w - plain_value(f[1])) / delta
        return round(Int64, position) + 1
    else 
        error("Frequency not in grid")
    end 
end


# iterate
#-------------------------------------------------------------------------------#
 
function Base.:iterate(
    m :: Mesh{T}
    ) :: Tuple{T, Int64} where {T <: AbstractMeshPoint}

    return m[firstindex(m)], 1 
end 

function Base.:iterate(
    m     :: Mesh{T},
    state :: Int64
    )     :: Union{Nothing, Tuple{T, Int64}} where {T <: AbstractMeshPoint}

    if state < length(m)
        return m[state + firstindex(m)], state + 1 
    else 
        return nothing 
    end
end

# load implementations and export
#-------------------------------------------------------------------------------#

# for each value type the respective mesh must implement:
# - outer constructor
# - mappings from mesh point and value type to mesh index 
# - boundary conditions
# - comparison operator

include("matsubara/matsubara_mesh.jl")
include("brillouin/brillouin_mesh.jl")

export
    AbstractMesh,
    Mesh,
    points,
    domain