# quick and dirty interface to the TRIQS library (https://triqs.github.io/triqs/latest/)

function load_mesh(h :: HDF5.File, l :: String, ::Val{:MeshImFreq}) :: Mesh
    @DEBUG read_attribute(h[l], "Format") == "MeshImFreq" "Dataset $(l) not tagged as MeshImFreq"

    if Bool(read(h[l], "positive_freq_only"))
        error("Positive frequency only meshes are currently not supported")
    end

    β         = read(h[l], "beta")
    size      = read(h[l], "size")   
    statistic = read(h[l], "statistic")

    if statistic == "F" 
        return MatsubaraMesh(1.0 / β, Int(size / 2), Fermion)
    elseif statistic == "B" 
        return MatsubaraMesh(1.0 / β, Int((size + 1) / 2), Boson)
    else 
        error("Statistic $(statistic) unknown")
    end 
end

function load_mesh(h :: HDF5.File, l :: String, ::Val{:MeshBrillouinZone}) :: Mesh
    @DEBUG read_attribute(h[l], "Format") == "MeshBrillouinZone" "Dataset $(l) not tagged as MeshBrillouinZone"

    if length(read(h[l], "brillouin_zone/bravais_lattice/atom_orb_pos")) > 1
        warn("Multiple orbitals in TRIQS mesh are currently ignored")
    end

    dims  = read(h[l], "dims")
    fdims = dims[findall(x -> x > 1, dims)]

    if length(fdims) > 1
        if any(x != fdims[1] for x in fdims)
            error("Non uniform momentum mesh dimensions are currently not supported")
        end
    end

    # transpose units matrix to switch from row to column major ordering
    units = read(h[l], "brillouin_zone/bravais_lattice/units")
    units = transpose(units)
    return BrillouinZoneMesh(BrillouinZone(length(fdims) > 1 ? fdims[1] : 1, SMatrix{size(units)..., Float64}(units)); reverse = true)
end 

load_mesh_triqs(h :: HDF5.File, l :: String) :: Mesh = load_mesh(h, l, Val(Symbol(read_attribute(h[l], "Format"))))

function load_triqs_gf(h :: HDF5.File, l :: String) :: MeshFunction 
    @DEBUG read_attribute(h[l], "Format") == "Gf" "Dataset $(l) not tagged as Gf"
    
    mesh_format = read_attribute(h[l * "/mesh"], "Format")
    meshes      = mesh_format == "MeshProduct" ? (load_mesh_triqs(h, l * "/mesh/$(key)") for key in keys(h[l * "/mesh"])) : [load_mesh_triqs(h, l * "/mesh")]

    # reverse ordering of data dimensions to switch from row to column major ordering
    data = read(h[l], "data")
    data = permutedims(data, collect(ndims(data) : -1 : 1)) 
    idxs = (IndexMesh(size(data, i)) for i in length(meshes) + 1 : ndims(data) - 1)
    tpl  = ntuple(x -> :, ndims(data) - 1)

    return MeshFunction(tuple(meshes..., idxs...), view(data, tpl..., 1) + im .* view(data, tpl..., 2))
end

export 
    load_mesh,
    load_triqs_gf