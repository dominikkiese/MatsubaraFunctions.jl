"""
    function save_matsubara_function!(
        h :: HDF5.File,
        l :: String,
        f :: MatsubaraFunction{Dg, Ds, Dt, Q}
        ) :: Nothing where {Dg, Ds, Dt, Q <: Number} 

Save MatsubaraFunction `f` with label `l` to file `h`  
"""
function save_matsubara_function!(
    h :: HDF5.File,
    l :: String,
    f :: MatsubaraFunction{Dg, Ds, Dt, Q}
    ) :: Nothing where {Dg, Ds, Dt, Q <: Number}

    # create new group
    grp = create_group(h, l)

    # add metadata 
    attributes(grp)["shape"] = Int64[f.shape...]

    # add data
    grp["data"] = f.data

    for i in eachindex(f.grids)
        save_matsubara_grid!(h, l * "/grids/grid_$i", f.grids[i])
    end

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

    # read the metadata 
    shape = read_attribute(h[l], "shape")

    # read the data
    idxs  = eachindex(keys(h[l * "/grids"]))
    grids = MatsubaraGrid[load_matsubara_grid(h, l * "/grids/grid_$i") for i in idxs]
    data  = read(h, l * "/data")
    
    return MatsubaraFunction((grids...,), (shape...,), data)
end

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

    # create new group 
    grp = create_group(h, l)

    # add metadata 
    attributes(grp)["speedup"]     = SG.speedup
    attributes(grp)["num_classes"] = length(SG.classes)

    # add data 
    num_in_classes        = Int64[length(class) for class in SG.classes]
    grp["num_in_classes"] = num_in_classes

    # convert classes to big matrix for fast write to disk
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
    function load_matsubara_symmetry_group(
        h :: HDF5.File,
        l :: String
        ) :: MatsubaraSymmetryGroup

Load MatsubaraSymmetryGroup with label `l` from file `h`
"""
function load_matsubara_symmetry_group(
    h :: HDF5.File,
    l :: String
    ) :: MatsubaraSymmetryGroup

    # read the metadata 
    speedup     = read_attribute(h[l], "speedup")
    num_classes = read_attribute(h[l], "num_classes")

    # read the data 
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

    return MatsubaraSymmetryGroup(classes, speedup)
end