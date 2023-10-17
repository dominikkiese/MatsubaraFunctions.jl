# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct Mesh{T <: AbstractMeshPoint}

Mesh type with fields:
* `hash   :: UInt64`    : mesh identifier
* `points :: Vector{T}` : mesh points
* `domain :: Dict`      : T-specific details
"""
struct Mesh{T <: AbstractMeshPoint}
    hash   :: UInt64
    points :: Vector{T}  
    domain :: Dict  
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
    function point(
        m   :: Mesh{T},
        idx :: Int64
        )   :: T where {T <: AbstractMeshPoint}

Returns `m.points[idx]`
"""
function point(
    m   :: Mesh{T},
    idx :: Int64
    )   :: T where {T <: AbstractMeshPoint}

    # bounds check performed by Base
    return m.points[idx]
end

"""
    function domain(
        m :: Mesh{T}
        ) :: Dict where {T <: AbstractMeshPoint}

Returns `m.domain`
"""
function domain(
    m :: Mesh{T}
    ) :: Dict where {T <: AbstractMeshPoint} 

    return m.domain 
end

function Base.:length(
    m :: Mesh{T}
    ) :: Int64 where {T <: AbstractMeshPoint} 

    return length(points(m))
end

# indexing
#-------------------------------------------------------------------------------#

function Base.:eachindex(
    m :: Mesh{T},
    ) :: Base.OneTo{Int64} where {T <: AbstractMeshPoint}

    return eachindex(points(m))
end

function Base.:firstindex(
    m :: Mesh{T},
    ) :: Int64 where {T <: AbstractMeshPoint}

    return firstindex(points(m))
end

function Base.:lastindex(
    m :: Mesh{T},
    ) :: Int64 where {T <: AbstractMeshPoint}

    return lastindex(points(m))
end

function Base.:getindex(
    m   :: Mesh{T},
    idx :: Int64
    )   :: T where {T <: AbstractMeshPoint}

    return point(m, idx)
end

function Base.:getindex(
    m    :: Mesh{T},
    idxs :: UnitRange{Int64}
    )    :: SubArray{T, 1, Vector{T}, Tuple{UnitRange{Int64}}, true} where {T <: AbstractMeshPoint}

    return @view points(m)[idxs]
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

# call
#-------------------------------------------------------------------------------#

# make Mesh callable with MeshPoint, checks if the point is valid
function (m :: Mesh{T})(
    x :: T
    ) :: Bool where {T <: AbstractMeshPoint}

    return x.hash == m.hash
end

# load implementation and export
#-------------------------------------------------------------------------------#

# for each value type the respective mesh must implement:
# - outer constructor
# - method to call mesh with value type for mapping to mesh index

include("matsubara/matsubara.jl")

export
    Mesh,
    points,
    point,
    domain