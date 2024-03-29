include("brillouin_pt.jl")
include("brillouin_zone.jl")

# outer constructors and accessors
#-------------------------------------------------------------------------------#

"""
    function BrillouinZoneMesh(bz :: BrillouinZone{N}; reverse = false) :: Mesh{MeshPoint{BrillouinPoint{N}}} where {N}

Construct uniform mesh for Brillouin zone. If reverse is true, the mesh is constructed using row mayor ordering.
"""
function BrillouinZoneMesh(bz :: BrillouinZone{N}; reverse = false) :: Mesh{MeshPoint{BrillouinPoint{N}}} where {N}
    HASH     = hash(reverse, hash(bz, hash(bz.L, hash(N))))
    ranges   = ntuple(x -> 0 : bz.L - 1, N)
    lin_idxs = reverse ? transpose(LinearIndices(ranges)) : LinearIndices(ranges)
    itrs     = Iterators.product(ranges...)
    points   = Vector{MeshPoint{BrillouinPoint{N}}}(undef, bz.L^N)

    for idxs in itrs
        lin_idx         = lin_idxs[(idxs .+ 1)...]
        points[lin_idx] = MeshPoint(HASH, lin_idx, BrillouinPoint(idxs...))
    end 

    domain = Dict(:bz => bz, :lin_idxs => lin_idxs, :reverse => reverse)
    return Mesh(HASH, points, domain)
end

# conversion from reciprocal to euclidean coordinates
#-------------------------------------------------------------------------------#

"""
    function euclidean(k :: MeshPoint{BrillouinPoint{N}}, m :: Mesh{MeshPoint{BrillouinPoint{N}}}) :: SVector{N, Float64} where {N}

Convert mesh point to euclidean coordinates
"""
function euclidean(k :: MeshPoint{BrillouinPoint{N}}, m :: Mesh{MeshPoint{BrillouinPoint{N}}}) :: SVector{N, Float64} where {N}
    @DEBUG k.hash === m.hash "Hash must be equal between mesh point and mesh"
    return euclidean(value(k), domain(m)[:bz])
end

"""
    function euclidean(k :: BrillouinPoint{N}, m :: Mesh{MeshPoint{BrillouinPoint{N}}}) :: SVector{N, Float64} where {N}

Convert reciprocal to euclidean coordinates
"""
function euclidean(k :: BrillouinPoint{N}, m :: Mesh{MeshPoint{BrillouinPoint{N}}}) :: SVector{N, Float64} where {N}
    return euclidean(k, domain(m)[:bz])
end

"""
    function euclideans(m :: Mesh{MeshPoint{BrillouinPoint{N}}}) :: Vector{SVector{N, Float64}} where {N}

Returns euclidean coordinates for all momenta in mesh
"""
function euclideans(m :: Mesh{MeshPoint{BrillouinPoint{N}}}) :: Vector{SVector{N, Float64}} where {N}
    return [euclidean(k, m) for k in points(m)]
end

# conversion from euclidean to reciprocal coordinates
#-------------------------------------------------------------------------------#

"""
    function reciprocal(k :: T, m :: Mesh{MeshPoint{BrillouinPoint{N}}}) :: SVector{N, Float64} where {N, T <: AbstractVector{Float64}}

Convert euclidean to reciprocal coordinates
"""
function reciprocal(k :: T, m :: Mesh{MeshPoint{BrillouinPoint{N}}}) :: SVector{N, Float64} where {N, T <: AbstractVector{Float64}}
    @DEBUG length(k) == N "Length mismatch for input vector"
    return reciprocal(k, domain(m)[:bz])
end

"""
    function reciprocals(m :: Mesh{MeshPoint{BrillouinPoint{N}}}) :: Vector{SVector{N, Int}} where {N}

Returns reciprocal coordinates for all momenta in mesh
"""
function reciprocals(m :: Mesh{MeshPoint{BrillouinPoint{N}}}) :: Vector{SVector{N, Int}} where {N}
    return plain_value.(points(m))
end

# bounds checking
#-------------------------------------------------------------------------------#

"""
    function is_inbounds(k :: BrillouinPoint{N}, m :: Mesh{MeshPoint{BrillouinPoint{N}}}) :: Bool where {N}

Checks if input in mesh
"""
function is_inbounds(k :: BrillouinPoint{N}, m :: Mesh{MeshPoint{BrillouinPoint{N}}}) :: Bool where {N}
    return is_inbounds(k, domain(m)[:bz])
end

"""
    function is_inbounds(k :: T, m :: Mesh{MeshPoint{BrillouinPoint{N}}}) :: Bool where {N, T <: AbstractVector{Float64}}

Checks if input in mesh
"""
function is_inbounds(k :: T, m :: Mesh{MeshPoint{BrillouinPoint{N}}}) :: Bool where {N, T <: AbstractVector{Float64}}
    @DEBUG length(k) == N "Length mismatch for input vector"
    return is_inbounds(k, domain(m)[:bz])
end

# overload dummy function
function is_inbounds_bc(idx :: Int, m :: Mesh{MeshPoint{BrillouinPoint{N}}}) where {N}
    error("No inbounds check available for types `Int` and `Mesh{MeshPoint{BrillouinPoint{N}}}`")
end

# periodic boundary conditions
#-------------------------------------------------------------------------------#

"""
    function fold_back(k :: BrillouinPoint{N}, m :: Mesh{MeshPoint{BrillouinPoint{N}}}) :: BrillouinPoint{N} where {N}

Use periodic boundary conditions to fold `k` back into mesh
"""
function fold_back(k :: BrillouinPoint{N}, m :: Mesh{MeshPoint{BrillouinPoint{N}}}) :: BrillouinPoint{N} where {N}
    return fold_back(k, domain(m)[:bz])
end

"""
    function fold_back(k :: T, m :: Mesh{MeshPoint{BrillouinPoint{N}}}) :: SVector{N, Float64} where {N, T <: AbstractVector{Float64}}

Use periodic boundary conditions to fold `k` back into mesh
"""
function fold_back(k :: T, m :: Mesh{MeshPoint{BrillouinPoint{N}}}) :: SVector{N, Float64} where {N, T <: AbstractVector{Float64}}
    @DEBUG length(k) == N "Length mismatch for input vector"
    return fold_back(k, domain(m)[:bz])
end

# mapping to mesh index
#-------------------------------------------------------------------------------#

# from value type
function mesh_index(k :: BrillouinPoint{N}, m :: Mesh{MeshPoint{BrillouinPoint{N}}}) where {N}
    @DEBUG is_inbounds(k, m) "Momentum not in mesh"
    return domain(m)[:lin_idxs][(value(k) .+ 1)...]
end

# from Vector of Float, returns index of closest mesh point
function mesh_index(k :: T, m :: Mesh{MeshPoint{BrillouinPoint{N}}}) where {N, T <: AbstractVector{Float64}}
    @DEBUG length(k) == N "Length mismatch for input vector"
    
    # find surrounding box
    x      = reciprocal(k, m)
    ranges = ntuple(n -> floor(Int, x[n]) : ceil(Int, x[n]), N)
    iters  = collect(Iterators.product(ranges...))

    # determine closest point in box
    min_idx  = 1
    min_dist = norm(euclidean(BrillouinPoint(iters[1]...), m) .- k)

    for i in 2 : length(iters)
        dist = norm(euclidean(BrillouinPoint(iters[i]...), m) .- k)
        
        if dist < min_dist
            min_idx  = i
            min_dist = dist
        end
    end

    # even if k is in the mesh, we need to fold back here
    return mesh_index(fold_back(BrillouinPoint(iters[min_idx]...), m), m)
end

# from value type with bc
function mesh_index_bc(k :: BrillouinPoint{N}, m :: Mesh{MeshPoint{BrillouinPoint{N}}}) where {N}
    return mesh_index(fold_back(k, m), m)
end

# comparison operator
#-------------------------------------------------------------------------------#

function Base.:(==)(m1 :: Mesh{MeshPoint{BrillouinPoint{N}}}, m2 :: Mesh{MeshPoint{BrillouinPoint{N}}}) where {N}
    if !(basis(domain(m1)[:bz]) ≈ basis(domain(m2)[:bz]))
        return false 
    end 

    if domain(m1)[:lin_idxs] != domain(m2)[:lin_idxs]
        return false 
    end 

    if domain(m1)[:reverse] != domain(m2)[:reverse]
        return false 
    end 

    for idx in eachindex(m1)
        if points(m1, idx) != points(m2, idx)
            return false 
        end 
    end 

    return m1.hash === m2.hash
end

# mapping to Wigner-Seitz cell
#-------------------------------------------------------------------------------#

"""
    function to_Wigner_Seitz(m :: Mesh{MeshPoint{BrillouinPoint{N}}}) :: Vector{SVector{N, Float64}} where {N}

Generate list of euclidean points in Wigner Seitz cell at Γ = 0.
"""
function to_Wigner_Seitz(m :: Mesh{MeshPoint{BrillouinPoint{N}}}) :: Vector{SVector{N, Float64}} where {N}
    shifts = get_shifts(domain(m)[:bz])
    return [to_Wigner_Seitz(euclidean(k, m), shifts) for k in m]
end

# io
#-------------------------------------------------------------------------------#

"""
    function save!(
        h :: HDF5.File,
        l :: String,
        m :: Mesh{MeshPoint{BrillouinPoint{N}}}
        ) :: Nothing where {N}

Save Brillouin zone mesh to HDF5 file `h` with label `l`
"""
function save!(
    h :: HDF5.File,
    l :: String,
    m :: Mesh{MeshPoint{BrillouinPoint{N}}}
    ) :: Nothing where {N}

    grp = create_group(h, l)

    # save metadata
    bz                         = domain(m)[:bz]
    attributes(grp)["reverse"] = domain(m)[:reverse]
    attributes(grp)["tag"]     = "BrillouinZoneMesh"
    attributes(grp)["basis"]   = basis(bz)
    attributes(grp)["L"]       = bz.L

    return nothing 
end

"""
    function load_mesh(h :: HDF5.File, l :: String, ::Val{:BrillouinZoneMesh}) :: Mesh

Overload of load_mesh for BrillouinZoneMesh
"""
function load_mesh(h :: HDF5.File, l :: String, ::Val{:BrillouinZoneMesh}) :: Mesh
    @DEBUG read_attribute(h[l], "tag") == "BrillouinZoneMesh" "Dataset $(l) not tagged as BrillouinZoneMesh"

    # load metadata
    reverse = read_attribute(h[l], "reverse")
    basis   = read_attribute(h[l], "basis")
    L       = read_attribute(h[l], "L")

    return BrillouinZoneMesh(BrillouinZone(L, SMatrix{size(basis)..., Float64}(basis)); reverse)
end

# export
#-------------------------------------------------------------------------------#

export 
    BrillouinZoneMesh,
    euclidean,
    euclideans,
    reciprocal,
    reciprocals, 
    is_inbounds,
    fold_back,
    to_Wigner_Seitz,
    save!,
    load_mesh