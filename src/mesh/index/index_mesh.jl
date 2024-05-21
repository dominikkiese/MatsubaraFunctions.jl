include("index_pt.jl")

# the domain
#-------------------------------------------------------------------------------#

"""
    struct IndexDomain <: AbstractDomain

IndexDomain type with fields:
* `N :: Int`: index count
"""
struct IndexDomain <: AbstractDomain
    N :: Int
end

# outer constructors and accessors
#-------------------------------------------------------------------------------#

"""
    function IndexMesh(N :: Int) :: Mesh{MeshPoint{Index}, IndexDomain}

Construct mesh with indices 1 : N
"""
function IndexMesh(N :: Int) :: Mesh{MeshPoint{Index}, IndexDomain}
    HASH   = hash(N)
    points = Vector{MeshPoint{Index}}(undef, N)

    for lin_idx in 1 : N
        points[lin_idx] = MeshPoint(HASH, lin_idx, Index(lin_idx))
    end 

    return Mesh(HASH, points, IndexDomain(N))
end

"""
    function N(m :: Mesh{MeshPoint{Index}, IndexDomain}) :: Int

Return number of indices in mesh
"""
function N(m :: Mesh{MeshPoint{Index}, IndexDomain}) :: Int
    return domain(m).N
end

"""
    function values(m :: Mesh{MeshPoint{Index}, IndexDomain}) :: Vector{Int}

Return values of all indices in mesh
"""
function Base.:values(m :: Mesh{MeshPoint{Index}, IndexDomain}) :: Vector{Int}
    return plain_value.(points(m))
end

# bounds checking
#-------------------------------------------------------------------------------#

"""
    function is_inbounds(w :: Index, m :: Mesh{MeshPoint{Index}, IndexDomain}) :: Bool

Checks if input in mesh
"""
function is_inbounds(w :: Index, m :: Mesh{MeshPoint{Index}, IndexDomain}) :: Bool
    return 1 <= value(w) <= N(m)
end

"""
    function is_inbounds(w :: Int, m :: Mesh{MeshPoint{Index}, IndexDomain}) :: Bool

Checks if input in mesh
"""
function is_inbounds(w :: Int, m :: Mesh{MeshPoint{Index}, IndexDomain}) :: Bool
    return 1 <= w <= N(m)
end

# mapping to mesh index
#-------------------------------------------------------------------------------#

# from value type
function mesh_index(w :: Index, m :: Mesh{MeshPoint{Index}, IndexDomain})
    @DEBUG is_inbounds(w, m) "Index not in mesh"
    return value(w)
end

# from Int
function mesh_index(w :: Int, m :: Mesh)
    @DEBUG is_inbounds(w, m) "Index not in mesh"
    return w
end

# comparison operator
#-------------------------------------------------------------------------------#

function Base.:(==)(m1 :: Mesh{MeshPoint{Index}, IndexDomain}, m2 :: Mesh{MeshPoint{Index}, IndexDomain})
    if N(m1) != N(m2)
        return false 
    end 

    for idx in eachindex(m1)
        if points(m1, idx) != points(m2, idx)
            return false 
        end 
    end 

    return m1.hash === m2.hash
end

# io
#-------------------------------------------------------------------------------#

"""
    function save!(
        h :: HDF5.File,
        l :: String,
        m :: Mesh{MeshPoint{Index}, IndexDomain}
        ) :: Nothing

Save Matsubara mesh to HDF5 file `h` with label `l`
"""
function save!(
    h :: HDF5.File,
    l :: String,
    m :: Mesh{MeshPoint{Index}, IndexDomain}
    ) :: Nothing

    grp = create_group(h, l)

    # save metadata
    attributes(grp)["tag"] = "IndexMesh"
    attributes(grp)["N"]   = N(m)

    return nothing 
end

"""
    function load_mesh(h :: HDF5.File, l :: String, ::Val{:IndexMesh}) :: Mesh

Overload of load_mesh for IndexMesh
"""
function load_mesh(h :: HDF5.File, l :: String, ::Val{:IndexMesh}) :: Mesh
    @DEBUG read_attribute(h[l], "tag") == "IndexMesh" "Dataset $(l) not tagged as IndexMesh"
    return IndexMesh(read_attribute(h[l], "N"))
end

# export
#-------------------------------------------------------------------------------#

export 
    IndexDomain,
    IndexMesh,
    N,
    values,
    is_inbounds,
    save!,
    load_mesh