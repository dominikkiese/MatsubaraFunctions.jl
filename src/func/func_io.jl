"""
    function save_matsubara_function!(
        h :: HDF5.File,
        l :: String,
        f :: MeshFunction
        ) :: Nothing

Save MatsubaraFunction `f` with label `l` to file `h`  
"""
function save_mesh_function!(
    h :: HDF5.File,
    l :: String,
    f :: MeshFunction
    ) :: Nothing

    grp = create_group(h, l)

    # save metadata
    attributes(grp)["type"]   = "MeshFunction"
    attributes(grp)["shape"]  = Int64[f.shape...]
    attributes(grp)["offset"] = Int64[f.offset...]

    for i in eachindex(meshes(f))
        save!(h, l * "/meshes/mesh_$i", meshes(f, i))
    end

    # save data
    grp["data"] = f.data

    return nothing 
end

"""
    function load_matsubara_function(
        h :: HDF5.File,
        l :: String
        ) :: MeshFunction

Load MatsubaraFunction with label `l` from file `h`
"""
function load_mesh_function(
    h :: HDF5.File,
    l :: String
    ) :: MeshFunction
    
    function load_mesh(_h, _l)
        if read_attribute(h[_l], "tag") == "MatsubaraMesh"
            return load_matsubara_mesh(_h, _l)
        elseif read_attribute(h[_l], "tag") == "BrillouinZoneMesh"
            return load_brillouin_zone_mesh(_h, _l)
        else
            throw(RuntimeError("unknown tag for mesh at ", _l, " in file ", _h, "."))
        end
    end

    # load the metadata 
    type   = read_attribute(h[l], "type")
    shape  = read_attribute(h[l], "shape")
    offset = read_attribute(h[l], "offset")

    @DEBUG type == "MeshFunction" "Type $(l) unknown"

    # load the data
    #grp = open_group(h, l)
    idxs  = eachindex(keys(h[l * "/meshes"]))
    meshs = [load_mesh(h, l * "/meshes/mesh_$i") for i in idxs]

    return MeshFunction((meshs...,), (shape...,), OffsetArray(read(h, l * "/data"), offset...))
end

#----------------------------------------------------------------------------------------------#

"""
    function save_matsubara_symmetry_group!(
        h  :: HDF5.File,
        l  :: String,
        SG :: MeshSymmetryGroup
        )  :: Nothing

    Save MatsubaraSymmetryGroup `SG` with label `l` to file `h`
"""
function save_mesh_symmetry_group!(
    h  :: HDF5.File,
    l  :: String,
    SG :: MeshSymmetryGroup
    )  :: Nothing

    grp = create_group(h, l)

    # save metadata
    attributes(grp)["type"]        = string(typeof(SG))
    attributes(grp)["speedup"]     = SG.speedup
    attributes(grp)["num_classes"] = length(SG.classes)

    # save data 
    num_in_classes        = [length(class) for class in SG.classes]
    grp["num_in_classes"] = num_in_classes

    # convert classes to matrix for fast write to disk
    mat    = Matrix{Int64}(undef, 3, sum(num_in_classes))
    offset = 0

    for cl_idx in eachindex(SG.classes)
        for e in eachindex(SG.classes[cl_idx])
            mat[1, offset + e] = SG.classes[cl_idx][e][1]
            mat[2, offset + e] = SG.classes[cl_idx][e][2].sgn 
            mat[3, offset + e] = SG.classes[cl_idx][e][2].con 
        end

        offset += num_in_classes[cl_idx]
    end

    grp["classes"] = mat 

    return nothing 
end

"""
    function MeshSymmetryGroup{GD, SD, DD, Q}(
        h :: HDF5.File,
        l :: String,
        ) :: MeshSymmetryGroup{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

Load MeshSymmetryGroup with label `l` from file `h`. If `l` was generated using 
`save_mesh_symmetry_group!`, the type parameters can be obtained from `h` as 
`read_attribute(h[l], "type")`.
"""
function MeshSymmetryGroup{GD, SD, DD, Q}(
    h :: HDF5.File,
    l :: String,
    ) :: MeshSymmetryGroup{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

    # load the metadata 
    type        = read_attribute(h[l], "type")
    speedup     = read_attribute(h[l], "speedup")
    num_classes = read_attribute(h[l], "num_classes")
    @DEBUG startswith(type, "MeshSymmetryGroup") "Type $(l) unknown"

    # load the data 
    num_in_classes = read(h, l * "/num_in_classes")
    mat            = read(h, l * "/classes")
    classes        = Vector{Vector{Tuple{Int64, MeshOperation}}}(undef, num_classes)
    offset         = 0 

    for cl_idx in eachindex(classes)
        class = Vector{Tuple{Int64, MeshOperation}}(undef, num_in_classes[cl_idx])
        
        for e in eachindex(class)
            class[e] = mat[1, offset + e], MeshOperation(Bool(mat[2, offset + e]), Bool(mat[3, offset + e]))
        end

        classes[cl_idx] = class; offset += num_in_classes[cl_idx]
    end 

    return MeshSymmetryGroup{GD, SD, DD, Q}(classes, speedup)
end

#----------------------------------------------------------------------------------------------#

export 
    save_mesh_function!,
    load_mesh_function,
    save_mesh_symmetry_group!