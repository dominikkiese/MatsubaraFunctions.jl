"""
    function save_matsubara_function!(
        h :: HDF5.File,
        l :: String,
        f :: MatsubaraFunction
        ) :: Nothing

Save MatsubaraFunction `f` with label `l` to file `h`  
"""
function save_matsubara_function!(
    h :: HDF5.File,
    l :: String,
    f :: MatsubaraFunction
    ) :: Nothing

    grp = create_group(h, l)

    # save metadata
    attributes(grp)["type"]   = "MatsubaraFunction"
    attributes(grp)["shape"]  = Int64[f.shape...]
    attributes(grp)["offset"] = Int64[f.offset...]

    for i in eachindex(grids(f))
        save_matsubara_grid!(h, l * "/grids/grid_$i", grids(f, i))
    end

    # save data
    grp["data"] = f.data

    return nothing 
end

"""
    function load_matsubara_function(
        h :: HDF5.File,
        l :: String
        ) :: MatsubaraFunction

Load MatsubaraFunction with label `l` from file `h`
"""
function load_matsubara_function(
    h :: HDF5.File,
    l :: String
    ) :: MatsubaraFunction

    # load the metadata 
    type   = read_attribute(h[l], "type")
    shape  = read_attribute(h[l], "shape")
    offset = read_attribute(h[l], "offset")

    @DEBUG type == "MatsubaraFunction" "Type $(l) unknown"

    # load the data
    idxs  = eachindex(keys(h[l * "/grids"]))
    grids = [load_matsubara_grid(h, l * "/grids/grid_$i") for i in idxs]

    return MatsubaraFunction((grids...,), (shape...,), OffsetArray(read(h, l * "/data"), offset...))
end

#----------------------------------------------------------------------------------------------#

"""
    function save_matsubara_symmetry_group!(
        h  :: HDF5.File,
        l  :: String,
        SG :: MatsubaraSymmetryGroup
        )  :: Nothing

    Save MatsubaraSymmetryGroup `SG` with label `l` to file `h`
"""
function save_matsubara_symmetry_group!(
    h  :: HDF5.File,
    l  :: String,
    SG :: MatsubaraSymmetryGroup
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
    function MatsubaraSymmetryGroup{GD, SD, DD, Q}(
        h :: HDF5.File,
        l :: String,
        ) :: MatsubaraSymmetryGroup{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

Load MatsubaraSymmetryGroup with label `l` from file `h`. If `l` was generated using 
`save_matsubara_symmetry_group!`, the type parameters can be obtained from `h` as 
`read_attribute(h[l], "type")`.
"""
function MatsubaraSymmetryGroup{GD, SD, DD, Q}(
    h :: HDF5.File,
    l :: String,
    ) :: MatsubaraSymmetryGroup{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

    # load the metadata 
    type        = read_attribute(h[l], "type")
    speedup     = read_attribute(h[l], "speedup")
    num_classes = read_attribute(h[l], "num_classes")
    @DEBUG startswith(type, "MatsubaraSymmetryGroup") "Type $(l) unknown"

    # load the data 
    num_in_classes = read(h, l * "/num_in_classes")
    mat            = read(h, l * "/classes")
    classes        = Vector{Vector{Tuple{Int64, MatsubaraOperation}}}(undef, num_classes)
    offset         = 0 

    for cl_idx in eachindex(classes)
        class = Vector{Tuple{Int64, MatsubaraOperation}}(undef, num_in_classes[cl_idx])
        
        for e in eachindex(class)
            class[e] = mat[1, offset + e], MatsubaraOperation(Bool(mat[2, offset + e]), Bool(mat[3, offset + e]))
        end

        classes[cl_idx] = class; offset += num_in_classes[cl_idx]
    end 

    return MatsubaraSymmetryGroup{GD, SD, DD, Q}(classes, speedup)
end

#----------------------------------------------------------------------------------------------#

export 
    save_matsubara_function!,
    load_matsubara_function,
    save_matsubara_symmetry_group!