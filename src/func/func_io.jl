"""
    function save_mesh_function!(
        h :: HDF5.File,
        l :: String,
        f :: MeshFunction
        ) :: Nothing

Save MeshFunction `f` with label `l` to file `h`  
"""
function save_mesh_function!(
    h :: HDF5.File,
    l :: String,
    f :: MeshFunction{MD, SD, DD, Q, Array{Q, DD}}
    ) :: Nothing where{MD, SD, DD, Q <: Number}

    grp = create_group(h, l)

    # save metadata
    attributes(grp)["type"]   = "MeshFunction"
    attributes(grp)["shape"]  = Int64[f.shape...]

    for i in eachindex(meshes(f))
        save!(h, l * "/meshes/mesh_$i", meshes(f, i))
    end

    # save data
    grp["data"] = f.data

    return nothing 
end

"""
    function load_mesh_function(
        h :: HDF5.File,
        l :: String
        ) :: MeshFunction{MD, SD, DD, Q, Array{Q, DD}}

Load MeshFunction with label `l` from file `h`
"""
function load_mesh_function(
    h :: HDF5.File,
    l :: String
    ) :: MeshFunction

    # load the metadata 
    type   = read_attribute(h[l], "type")
    shape  = read_attribute(h[l], "shape")

    @DEBUG type == "MeshFunction" "Type $(l) unknown"

    # load the data
    idxs  = eachindex(keys(h[l * "/meshes"]))
    grid_gpnames = [l * "/meshes/mesh_$i" for i in idxs]
    grids = [load_mesh(h, gname, Val(hash(read_attribute(h[gname], "tag")))) for gname in grid_gpnames]

    return MeshFunction((grids...,), (shape...,), read(h, l * "/data"))
end


export 
    save_mesh_function!,
    load_mesh_function