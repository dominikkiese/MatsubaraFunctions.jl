# outer constructors and accessors
#-------------------------------------------------------------------------------#

"""
    function MatsubaraMesh(
        temperature :: Float64, 
        N           :: Int64,
                    :: Type{Fermion}
        )           :: Mesh{MeshPoint{MatsubaraFrequency{Fermion}}}

Constructs fermionic Matsubara mesh with 2 * N symmetrically spaced frequencies 
"""
function MatsubaraMesh(
    temperature :: Float64, 
    N           :: Int64,
                :: Type{Fermion}
    )           :: Mesh{MeshPoint{MatsubaraFrequency{Fermion}}}

    HASH      = hash(temperature, hash(N, hash(Fermion)))
    idx_range = -N : N - 1
    points    = Vector{MeshPoint{MatsubaraFrequency{Fermion}}}(undef, length(idx_range))

    for lin_idx in eachindex(idx_range)
        points[lin_idx] = MeshPoint(HASH, lin_idx, MatsubaraFrequency(temperature, idx_range[lin_idx], Fermion))
    end 

    domain = Dict(:temperature => temperature, :N => N, :type => "Fermion")
    return Mesh(HASH, points, domain)
end

"""
    function MatsubaraMesh(
        temperature :: Float64, 
        N           :: Int64,
                    :: Type{Boson}
        )           :: Mesh{MeshPoint{MatsubaraFrequency{Boson}}}

Constructs bosonic Matsubara mesh with 2 * N - 1 symmetrically spaced frequencies 
"""
function MatsubaraMesh(
    temperature :: Float64, 
    N           :: Int64,
                :: Type{Boson}
    )           :: Mesh{MeshPoint{MatsubaraFrequency{Boson}}}

    HASH      = hash(temperature, hash(N, hash(Boson)))
    idx_range = -N + 1 : N - 1
    points    = Vector{MeshPoint{MatsubaraFrequency{Boson}}}(undef, length(idx_range))

    for lin_idx in eachindex(idx_range)
        points[lin_idx] = MeshPoint(HASH, lin_idx, MatsubaraFrequency(temperature, idx_range[lin_idx], Boson))
    end 

    domain = Dict(:temperature => temperature, :N => N, :type => "Boson")
    return Mesh(HASH, points, domain)
end

"""
    function temperature(
        m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
        ) :: Float64 where {PT <: AbstractParticle}

Returns temperature for which the mesh is defined
"""
function temperature(
    m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
    ) :: Float64 where {PT <: AbstractParticle}

    return domain(m)[:temperature]
end

"""
    function first_index(
        m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
        ) :: Int64 where {PT <: AbstractParticle}

Returns the index of the first Matsubara frequency in mesh
"""
function first_index(
    m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
    ) :: Int64 where {PT <: AbstractParticle}

    return index(value(m[1]))
end

"""
    function last_index(
        m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
        ) :: Int64 where {PT <: AbstractParticle}

Returns the index of the last Matsubara frequency in mesh
"""
function last_index(
    m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
    ) :: Int64 where {PT <: AbstractParticle}

    return index(value(m[end]))
end

"""
    function first_value(
        m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
        ) :: Float64 where {PT <: AbstractParticle}

Returns the value of the first Matsubara frequency in mesh
"""
function first_value(
    m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
    ) :: Float64 where {PT <: AbstractParticle}

    return value(value(m[1]))
end

"""
    function last_value(
        m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
        ) :: Float64 where {PT <: AbstractParticle}

Returns the value of the last Matsubara frequency in mesh
"""
function last_value(
    m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
    ) :: Float64 where {PT <: AbstractParticle}

    return value(value(m[end]))
end

"""
    function indices(
        m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
        ) :: Vector{Int64} where {PT <: AbstractParticle}

Return indices of all Matsubara frequencies in mesh
"""
function indices(
    m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
    ) :: Vector{Int64} where {PT <: AbstractParticle}

    return index.(value.(points(m)))
end

"""
    function values(
        m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
        ) :: Vector{Float64} where {PT <: AbstractParticle}

Return values of all Matsubara frequencies in mesh
"""
function Base.:values(
    m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
    ) :: Vector{Float64} where {PT <: AbstractParticle}

    return value.(value.(points(m)))
end

# mapping to mesh index
#-------------------------------------------------------------------------------#

"""
    function is_inbounds(
        w :: MatsubaraFrequency{PT},
        m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
        ) :: Bool where {PT <: AbstractParticle}

Checks if Matsubara frequency in mesh
"""
function is_inbounds(
    w :: MatsubaraFrequency{PT},
    m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
    ) :: Bool where {PT <: AbstractParticle}

    @DEBUG temperature(w) ≈ temperature(m) "Temperature must be equal between Matsubara frequency and grid"
    return first_index(m) <= index(w) <= last_index(m)
end

"""
    function is_inbounds(
        w :: Float64, 
        m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
        ) :: Bool where {PT <: AbstractParticle}

Checks if float in mesh
"""
function is_inbounds(
    w :: Float64, 
    m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
    ) :: Bool where {PT <: AbstractParticle}

    return first_value(m) <= w <= last_value(m)
end

# methods for mapping Matsubara frequency to mesh index 
function mesh_index(
    w :: MatsubaraFrequency{PT},
    m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
    ) :: Int64 where {PT <: AbstractParticle}

    @DEBUG temperature(w) ≈ temperature(m) "Temperature must be equal between Matsubara frequency and mesh"
    return index(w) - first_index(m) + 1
end

function mesh_index_extrp(
    w :: MatsubaraFrequency{PT},
    m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
    ) :: Int64 where {PT <: AbstractParticle}

    return max(1, min(mesh_index(w, m), length(m)))
end

# make mesh callable with MatsubaraFrequency
function (m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}})(
    w :: MatsubaraFrequency{PT}
    ) :: Int64 where {PT <: AbstractParticle}

    @DEBUG is_inbounds(w, m) "Matsubara frequency not in mesh"
    return mesh_index(w, m)
end

# make mesh callable with Float64, returns index of closest frequency 
function (m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}})(
    w :: Float64
    ) :: Int64 where {PT <: AbstractParticle}

    @DEBUG is_inbounds(w, m) "Value not in mesh"
    delta    = value(value(m[2])) - value(value(m[1]))
    position = (w - value(value(m[1]))) / delta
    return round(Int64, position) + 1
end

# comparison operator
#-------------------------------------------------------------------------------#

function Base.:(==)(
    m1 :: Mesh{MeshPoint{MatsubaraFrequency{PT}}},
    m2 :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
    )  :: Bool where {PT <: AbstractParticle}

    if m1.hash != m2.hash
        return false
    end

    if !(domain(m1)[:temperature] ≈ domain(m2)[:temperature])
        return false 
    end 

    if domain(m1)[:N] != domain(m2)[:N]
        return false 
    end 

    if domain(m1)[:type] != domain(m2)[:type]
        return false 
    end 

    for idx in eachindex(m1)
        if point(m1, idx) != point(m2, idx)
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
        m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
        ) :: Nothing where {PT <: AbstractParticle}

Save Matsubara mesh to HDF5 file `h` with label `l`
"""
function save!(
    h :: HDF5.File,
    l :: String,
    m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
    ) :: Nothing where {PT <: AbstractParticle}

    grp = create_group(h, l)

    # save metadata
    attributes(grp)["tag"]         = "MatsubaraMesh"
    attributes(grp)["temperature"] = domain(m)[:temperature]
    attributes(grp)["N"]           = domain(m)[:N]
    attributes(grp)["type"]        = domain(m)[:type]

    return nothing 
end

"""
    function load_matsubara_mesh(
        h :: HDF5.File, 
        l :: String
        ) :: AbstractMesh

Load Matsubara mesh with label `l` from HDF5 file `h`
"""
function load_matsubara_mesh(
    h :: HDF5.File, 
    l :: String
    ) :: AbstractMesh

    @DEBUG read_attribute(h[l], "tag") == "MatsubaraMesh" "Dataset $(l) not tagged as MatsubaraMesh"

    temperature = read_attribute(h[l], "temperature")
    N           = read_attribute(h[l], "N")
    type        = read_attribute(h[l], "type")

    if type == "Fermion"
        return MatsubaraMesh(temperature, N, Fermion)
    elseif type == "Boson"
        return MatsubaraMesh(temperature, N, Boson)
    end

    return nothing 
end

# export
#-------------------------------------------------------------------------------#

export 
    MatsubaraMesh,
    temperature,
    first_index,
    last_index,
    first_value,
    last_value,
    indices,
    values,
    is_inbounds,
    save!,
    load_matsubara_mesh