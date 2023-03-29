"""
    function save_matsubara_grid!(
        h :: HDF5.File,
        l :: String,
        g :: MatsubaraGrid
        ) :: Nothing

Save MatsubaraGrid g with label l to file h
"""
function save_matsubara_grid!(
    h :: HDF5.File,
    l :: String,
    g :: MatsubaraGrid
    ) :: Nothing

    # create new group
    grp = create_group(h, l)

    # add metadata
    attributes(grp)["T"]    = temperature(g)
    attributes(grp)["type"] = "$(type(g))"

    # add data
    for i in 1 : length(g)
        grp["frequency/$i/val"] = value(g[i])
        grp["frequency/$i/idx"] = index(g[i])
    end

    return nothing 
end 

"""
    function load_matsubara_grid(
        h :: HDF5.File,
        l :: String
        ) :: MatsubaraGrid

Load MatsubaraGrid with label l from file h
"""
function load_matsubara_grid(
    h :: HDF5.File,
    l :: String
    ) :: MatsubaraGrid 

    # read the metadata
    T    = read_attribute(h[l], "T")
    type = read_attribute(h[l], "type")
    data = Vector{MatsubaraFrequency}(undef, length(keys(h[l * "/frequency"])))

    # read the data
    for i in eachindex(data)
        val     = read(h, l * "/frequency/$i/val")
        idx     = read(h, l * "/frequency/$i/idx")
        data[i] = MatsubaraFrequency(T, val, idx, Symbol(type))
    end

    return MatsubaraGrid(T, data, Symbol(type))
end

"""
    function save_matsubara_function!(
        h :: HDF5.File,
        l :: String,
        f :: MatsubaraFunction{Dg, Ds, Dt, Q}
        ) :: Nothing where {Dg, Ds, Dt, Q <: Number} 

Save MatsubaraFunction f with label l to file h   
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

Load MatsubaraFunction with label l from file h
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

    Save MatsubaraSymmetryGroup SG with label l to file h 
"""
function save_matsubara_symmetry_group!(
    h  :: HDF5.File,
    l  :: String,
    SG :: MatsubaraSymmetryGroup
    )  :: Nothing

    # create new group 
    grp = create_group(h, l)

    # add metadata 
    attributes(grp)["num_classes"] = length(SG.classes)

    # add data 
    for cl_idx in eachindex(SG.classes)
        # convert class to matrix for fast write to disk
        class = SG.classes[cl_idx]
        mat   = Matrix{Int64}(undef, 3, length(class))

        for e in eachindex(class)
            mat[1, e] = class[e][1]
            mat[2, e] = class[e][2].sgn 
            mat[3, e] = class[e][2].con 
        end

        grp["class_$(cl_idx)"] = mat 
    end 

    return nothing 
end

"""
function load_matsubara_symmetry_group(
    h :: HDF5.File,
    l :: String
    ) :: MatsubaraSymmetryGroup

Load MatsubaraSymmetryGroup with label l from file h
"""
function load_matsubara_symmetry_group(
    h :: HDF5.File,
    l :: String
    ) :: MatsubaraSymmetryGroup

    # read the metadata 
    num_classes = read_attribute(h[l], "num_classes")

    # read the data 
    classes = Vector{Vector{Tuple{Int64, MatsubaraOperation}}}(undef, num_classes)

    for cl_idx in eachindex(classes)
        mat             = read(h, l * "/class_$(cl_idx)")
        classes[cl_idx] = Tuple{Int64, MatsubaraOperation}[(mat[1, i], MatsubaraOperation(Bool(mat[2, i]), Bool(mat[3, i]))) for i in 1 : size(mat, 2)]
    end 

    return MatsubaraSymmetryGroup(classes)
end