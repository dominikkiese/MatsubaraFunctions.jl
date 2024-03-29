include("matsubara_freq.jl")

# outer constructors and accessors
#-------------------------------------------------------------------------------#

"""
    function MatsubaraMesh(
        temperature :: Float64, 
        N           :: Int,
                    :: Type{Fermion}
        )           :: Mesh{MeshPoint{MatsubaraFrequency{Fermion}}}

Constructs fermionic Matsubara mesh with 2 * N symmetrically spaced frequencies 
"""
function MatsubaraMesh(
    temperature :: Float64, 
    N           :: Int,
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
        N           :: Int,
                    :: Type{Boson}
        )           :: Mesh{MeshPoint{MatsubaraFrequency{Boson}}}

Constructs bosonic Matsubara mesh with 2 * N - 1 symmetrically spaced frequencies 
"""
function MatsubaraMesh(
    temperature :: Float64, 
    N           :: Int,
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
    function first_index(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}) :: Int where {PT <: AbstractParticle}

Returns the index of the first Matsubara frequency in mesh
"""
function first_index(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}) :: Int where {PT <: AbstractParticle}
    return index(value(m[1]))
end

"""
    function last_index(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}) :: Int where {PT <: AbstractParticle}

Returns the index of the last Matsubara frequency in mesh
"""
function last_index(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}) :: Int where {PT <: AbstractParticle}
    return index(value(m[end]))
end

"""
    function first_value(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}) :: Float64 where {PT <: AbstractParticle}

Returns the value of the first Matsubara frequency in mesh
"""
function first_value(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}) :: Float64 where {PT <: AbstractParticle}
    return plain_value(m[1])
end

"""
    function last_value(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}) :: Float64 where {PT <: AbstractParticle}

Returns the value of the last Matsubara frequency in mesh
"""
function last_value(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}) :: Float64 where {PT <: AbstractParticle}
    return plain_value(m[end])
end

"""
    function indices(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}) :: Vector{Int} where {PT <: AbstractParticle}

Return indices of all Matsubara frequencies in mesh
"""
function indices(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}) :: Vector{Int} where {PT <: AbstractParticle}
    return index.(value.(points(m)))
end

"""
    function values(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}) :: Vector{Float64} where {PT <: AbstractParticle}

Return values of all Matsubara frequencies in mesh
"""
function Base.:values(m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}) :: Vector{Float64} where {PT <: AbstractParticle}
    return plain_value.(points(m))
end

# bounds checking
#-------------------------------------------------------------------------------#

"""
    function is_inbounds(w :: MatsubaraFrequency{PT}, m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}) :: Bool where {PT <: AbstractParticle}

Checks if Matsubara frequency in mesh
"""
function is_inbounds(w :: MatsubaraFrequency{PT}, m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}) :: Bool where {PT <: AbstractParticle}
    @DEBUG temperature(w) ≈ domain(m)[:temperature] "Temperature must be equal between Matsubara frequency and grid"
    return first_index(m) <= index(w) <= last_index(m)
end

"""
    function is_inbounds(w :: Float64, m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}) :: Bool where {PT <: AbstractParticle}

Checks if float in mesh
"""
function is_inbounds(w :: Float64, m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}) :: Bool where {PT <: AbstractParticle}
    return first_value(m) <= w <= last_value(m)
end

# overload dummy function
function is_inbounds_bc(w :: Union{MatsubaraFrequency{PT}, Float64}, m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}) where {PT <: AbstractParticle}
    return is_inbounds(w, m)
end

function is_inbounds_bc(idx :: Int, m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}) where {PT <: AbstractParticle}
    error("No inbounds check available for types `Int` and `Mesh{MeshPoint{MatsubaraFrequency{PT}}`")
end

# mapping to mesh index
#-------------------------------------------------------------------------------#

# from value type
function mesh_index(w :: MatsubaraFrequency{PT}, m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}) where {PT <: AbstractParticle}
    @DEBUG is_inbounds(w, m) "Matsubara frequency not in mesh"
    return index(w) - first_index(m) + 1
end

# from Float, returns index of closest frequency
function mesh_index(w :: Float64, m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}) where {PT <: AbstractParticle}
    @DEBUG is_inbounds(w, m) "Value not in mesh"
    delta    = plain_value(m[2]) - plain_value(m[1])
    position = (w - plain_value(m[1])) / delta
    return round(Int, position) + 1
end

# comparison operator
#-------------------------------------------------------------------------------#

function Base.:(==)(m1 :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}, m2 :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}) where {PT <: AbstractParticle}
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
    MatsubaraMesh,
    first_index,
    last_index,
    first_value,
    last_value,
    indices,
    values,
    is_inbounds,
    save!,
    load_mesh