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
    f :: MeshFunction
    ) :: Nothing

    grp = create_group(h, l)

    # save metadata
    attributes(grp)["type"] = "MeshFunction"

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
        ) :: MeshFunction

Load MeshFunction with label `l` from file `h`
"""
function load_mesh_function(
    h :: HDF5.File,
    l :: String
    ) :: MeshFunction

    # load the metadata 
    type = read_attribute(h[l], "type")
    @DEBUG type == "MeshFunction" "Type $(l) unknown"

    # load the data
    grids = (load_mesh(h, l * "/meshes/mesh_$i") for i in eachindex(keys(h[l * "/meshes"])))
    return MeshFunction(tuple(grids...), read(h, l * "/data"))
end

export 
    save_mesh_function!,
    load_mesh_function