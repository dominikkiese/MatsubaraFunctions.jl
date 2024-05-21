include("matsubara_freq.jl")

# the domain
#-------------------------------------------------------------------------------#

"""
    struct MatsubaraDomain <: AbstractDomain

MatsubaraDomain type with fields:
* `temperature :: Float64` : temperature of the Matsubara mesh
* `N           :: Int`     : frequency count
"""
struct MatsubaraDomain <: AbstractDomain
    temperature :: Float64
    N           :: Int
end

# outer constructors and accessors
#-------------------------------------------------------------------------------#

"""
    function MatsubaraMesh(
        temperature :: Float64, 
        N           :: Int,
                    :: Type{Fermion}
        )           :: Mesh{MeshPoint{MatsubaraFrequency{Fermion}}, MatsubaraDomain}

Constructs fermionic Matsubara mesh with 2 * N symmetrically spaced frequencies 
"""
function MatsubaraMesh(
    temperature :: Float64, 
    N           :: Int,
                :: Type{Fermion}
    )           :: Mesh{MeshPoint{MatsubaraFrequency{Fermion}}, MatsubaraDomain}

    HASH      = hash(temperature, hash(N, hash(Fermion)))
    idx_range = -N : N - 1
    points    = Vector{MeshPoint{MatsubaraFrequency{Fermion}}}(undef, length(idx_range))

    for lin_idx in eachindex(idx_range)
        points[lin_idx] = MeshPoint(HASH, lin_idx, MatsubaraFrequency(temperature, idx_range[lin_idx], Fermion))
    end 

    return Mesh(HASH, points, MatsubaraDomain(temperature, N))
end

"""
    function MatsubaraMesh(
        temperature :: Float64, 
        N           :: Int,
                    :: Type{Boson}
        )           :: Mesh{MeshPoint{MatsubaraFrequency{Boson}}, MatsubaraDomain}

Constructs bosonic Matsubara mesh with 2 * N - 1 symmetrically spaced frequencies 
"""
function MatsubaraMesh(
    temperature :: Float64, 
    N           :: Int,
                :: Type{Boson}
    )           :: Mesh{MeshPoint{MatsubaraFrequency{Boson}}, MatsubaraDomain}

    HASH      = hash(temperature, hash(N, hash(Boson)))
    idx_range = -N + 1 : N - 1
    points    = Vector{MeshPoint{MatsubaraFrequency{Boson}}}(undef, length(idx_range))

    for lin_idx in eachindex(idx_range)
        points[lin_idx] = MeshPoint(HASH, lin_idx, MatsubaraFrequency(temperature, idx_range[lin_idx], Boson))
    end 

    return Mesh(HASH, points, MatsubaraDomain(temperature, N))
end

"""
    function temperature(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) :: Float64 where {PT <: AbstractParticle}

Returns the temperature of the Matsubara mesh
"""
function temperature(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) :: Float64 where {PT <: AbstractParticle}
    return domain(m).temperature
end

"""
    function N(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) :: Int where {PT <: AbstractParticle}

Returns the frequency count used for construction of the Matsubara mesh
"""
function N(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) :: Int where {PT <: AbstractParticle}
    return domain(m).N
end

"""
    function first_index(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) :: Int where {PT <: AbstractParticle}

Returns the index of the first Matsubara frequency in mesh
"""
function first_index(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) :: Int where {PT <: AbstractParticle}
    return index(value(m[1]))
end

"""
    function last_index(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) :: Int where {PT <: AbstractParticle}

Returns the index of the last Matsubara frequency in mesh
"""
function last_index(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) :: Int where {PT <: AbstractParticle}
    return index(value(m[end]))
end

"""
    function first_value(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) :: Float64 where {PT <: AbstractParticle}

Returns the value of the first Matsubara frequency in mesh
"""
function first_value(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) :: Float64 where {PT <: AbstractParticle}
    return plain_value(m[1])
end

"""
    function last_value(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) :: Float64 where {PT <: AbstractParticle}

Returns the value of the last Matsubara frequency in mesh
"""
function last_value(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) :: Float64 where {PT <: AbstractParticle}
    return plain_value(m[end])
end

"""
    function indices(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) :: Vector{Int} where {PT <: AbstractParticle}

Return indices of all Matsubara frequencies in mesh
"""
function indices(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) :: Vector{Int} where {PT <: AbstractParticle}
    return index.(value.(points(m)))
end

"""
    function values(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) :: Vector{Float64} where {PT <: AbstractParticle}

Return values of all Matsubara frequencies in mesh
"""
function Base.:values(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) :: Vector{Float64} where {PT <: AbstractParticle}
    return plain_value.(points(m))
end

# bounds checking
#-------------------------------------------------------------------------------#

"""
    function is_inbounds(w :: MatsubaraFrequency{PT}, m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) :: Bool where {PT <: AbstractParticle}

Checks if Matsubara frequency in mesh
"""
function is_inbounds(w :: MatsubaraFrequency{PT}, m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) :: Bool where {PT <: AbstractParticle}
    @DEBUG temperature(w) ≈ temperature(m) "Temperature must be equal between Matsubara frequency and grid"
    return first_index(m) <= index(w) <= last_index(m)
end

"""
    function is_inbounds(w :: Float64, m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) :: Bool where {PT <: AbstractParticle}

Checks if float in mesh
"""
function is_inbounds(w :: Float64, m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) :: Bool where {PT <: AbstractParticle}
    return first_value(m) <= w <= last_value(m)
end

# overload dummy function
function is_inbounds_bc(w :: Union{MatsubaraFrequency{PT}, Float64}, m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) where {PT <: AbstractParticle}
    return is_inbounds(w, m)
end

function is_inbounds_bc(idx :: Int, m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) where {PT <: AbstractParticle}
    error("No inbounds check available for types `Int` and `Mesh{MeshPoint{MatsubaraFrequency{PT}}`")
end

# mapping to mesh index
#-------------------------------------------------------------------------------#

# from value type
function mesh_index(w :: MatsubaraFrequency{PT}, m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) where {PT <: AbstractParticle}
    @DEBUG is_inbounds(w, m) "Matsubara frequency not in mesh"
    return index(w) - first_index(m) + 1
end

# from Float, returns index of closest frequency
function mesh_index(w :: Float64, m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) where {PT <: AbstractParticle}
    @DEBUG is_inbounds(w, m) "Value not in mesh"
    delta    = plain_value(m[2]) - plain_value(m[1])
    position = (w - plain_value(m[1])) / delta
    return round(Int, position) + 1
end

# comparison operator
#-------------------------------------------------------------------------------#

function Base.:(==)(m1 :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}, m2 :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) where {PT <: AbstractParticle}
    if !(temperature(m1) ≈ temperature(m2))
        return false 
    end 

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
        m :: Mesh{MeshPoint{MatsubaraFrequency{Fermion}}, MatsubaraDomain}
        ) :: Nothing

Save Matsubara mesh to HDF5 file `h` with label `l`
"""
function save!(
    h :: HDF5.File,
    l :: String,
    m :: Mesh{MeshPoint{MatsubaraFrequency{Fermion}}, MatsubaraDomain}
    ) :: Nothing

    grp = create_group(h, l)

    # save metadata
    attributes(grp)["tag"]         = "MatsubaraMesh"
    attributes(grp)["type"]        = "Fermion"
    attributes(grp)["temperature"] = temperature(m)
    attributes(grp)["N"]           = N(m)

    return nothing 
end

"""
    function save!(
        h :: HDF5.File,
        l :: String,
        m :: Mesh{MeshPoint{MatsubaraFrequency{Boson}}, MatsubaraDomain}
        ) :: Nothing

Save Matsubara mesh to HDF5 file `h` with label `l`
"""
function save!(
    h :: HDF5.File,
    l :: String,
    m :: Mesh{MeshPoint{MatsubaraFrequency{Boson}}, MatsubaraDomain}
    ) :: Nothing

    grp = create_group(h, l)

    # save metadata
    attributes(grp)["tag"]         = "MatsubaraMesh"
    attributes(grp)["type"]        = "Boson"
    attributes(grp)["temperature"] = temperature(m)
    attributes(grp)["N"]           = N(m)

    return nothing 
end

"""
    function load_mesh(h :: HDF5.File, l :: String, ::Val{:MatsubaraMesh}) :: Mesh

Overload of load_mesh for MatsubaraMesh
"""
function load_mesh(h :: HDF5.File, l :: String, ::Val{:MatsubaraMesh}) :: Mesh
    @DEBUG read_attribute(h[l], "tag") == "MatsubaraMesh" "Dataset $(l) not tagged as MatsubaraMesh"

    # load metadata
    temperature = read_attribute(h[l], "temperature")
    N           = read_attribute(h[l], "N")
    type        = read_attribute(h[l], "type")

    if type == "Fermion"
        return MatsubaraMesh(temperature, N, Fermion)
    elseif type == "Boson"
        return MatsubaraMesh(temperature, N, Boson)
    else 
        error("Particle type $(type) unknown")
    end
end

# export
#-------------------------------------------------------------------------------#

export 
    MatsubaraDomain,
    MatsubaraMesh,
    temperature,
    N,
    first_index,
    last_index,
    first_value,
    last_value,
    indices,
    values,
    is_inbounds,
    save!,
    load_mesh