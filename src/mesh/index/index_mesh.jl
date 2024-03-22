include("index_pt.jl")

# outer constructors and accessors
#-------------------------------------------------------------------------------#

"""
    function IndexMesh(N :: Int) :: Mesh{MeshPoint{Index}}

Construct mesh with indices 1 : N
"""
function IndexMesh(N :: Int) :: Mesh{MeshPoint{Index}}
    HASH   = hash(N)
    points = Vector{MeshPoint{Index}}(undef, N)

    for lin_idx in 1 : N
        points[lin_idx] = MeshPoint(HASH, lin_idx, Index(lin_idx))
    end 

    domain = Dict(:N => N)
    return Mesh(HASH, points, domain)
end

"""
    function values(m :: Mesh{MeshPoint{Index}}) :: Vector{Int}

Return values of all Indices in mesh
"""
function Base.:values(m :: Mesh{MeshPoint{Index}}) :: Vector{Int}
    return plain_value.(points(m))
end

# bounds checking
#-------------------------------------------------------------------------------#

"""
    function is_inbounds(w :: Index, m :: Mesh{MeshPoint{Index}}) :: Bool

Checks if input in mesh
"""
function is_inbounds(w :: Index, m :: Mesh{MeshPoint{Index}}) :: Bool
    return 1 <= value(w) <= length(m)
end

"""
    function is_inbounds(w :: Int, m :: Mesh{MeshPoint{Index}}) :: Bool

Checks if input in mesh
"""
function is_inbounds(w :: Int, m :: Mesh{MeshPoint{Index}}) :: Bool
    return 1 <= w <= length(m)
end

# mapping to mesh index
#-------------------------------------------------------------------------------#

# from value type
function mesh_index(w :: Index, m :: Mesh{MeshPoint{Index}})
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

function Base.:(==)(m1 :: Mesh{MeshPoint{Index}}, m2 :: Mesh{MeshPoint{Index}})
    if m1.hash != m2.hash
        return false
    end

    if domain(m1)[:N] != domain(m2)[:N]
        return false 
    end 

    for idx in eachindex(m1)
        if points(m1, idx) != points(m2, idx)
            return false 
        end 
    end 

    return true
end

# io
#-------------------------------------------------------------------------------#

"""
    function save!(
        h :: HDF5.File,
        l :: String,
        m :: Mesh{MeshPoint{Index}}
        ) :: Nothing

Save Matsubara mesh to HDF5 file `h` with label `l`
"""
function save!(
    h :: HDF5.File,
    l :: String,
    m :: Mesh{MeshPoint{Index}}
    ) :: Nothing

    grp = create_group(h, l)

    # save metadata
    attributes(grp)["tag"] = "IndexMesh"
    attributes(grp)["N"]   = domain(m)[:N]

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
    IndexMesh,
    values,
    is_inbounds,
    save!,
    load_mesh