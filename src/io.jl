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

    for cl_idx in eachindex(SG.classes)
        attributes(grp)["num_elements_class_$(cl_idx)"] = length(SG.classes[cl_idx])
    end

    # add data 
    for cl_idx in eachindex(SG.classes)
        for e_idx in eachindex(SG.classes[cl_idx])
            grp["class_$(cl_idx)/element_$(e_idx)/idx"]    = first(SG.classes[cl_idx][e_idx])
            grp["class_$(cl_idx)/element_$(e_idx)/op/sgn"] = last(SG.classes[cl_idx][e_idx]).sgn
            grp["class_$(cl_idx)/element_$(e_idx)/op/con"] = last(SG.classes[cl_idx][e_idx]).con
        end 
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
    num_classes  = read_attribute(h[l], "num_classes")
    num_elements = Int64[read_attribute(h[l], "num_elements_class_$(cl_idx)") for cl_idx in 1 : num_classes]

    # read the data 
    classes = Vector{Vector{Tuple{Int64, MatsubaraOperation}}}(undef, num_classes)

    for cl_idx in eachindex(classes)
        class = Vector{Tuple{Int64, MatsubaraOperation}}(undef, num_elements[cl_idx])

        for e_idx in eachindex(class)
            idx          = read(h, l * "/class_$(cl_idx)/element_$(e_idx)/idx")
            sgn          = read(h, l * "/class_$(cl_idx)/element_$(e_idx)/op/sgn")
            con          = read(h, l * "/class_$(cl_idx)/element_$(e_idx)/op/con")
            class[e_idx] = idx, MatsubaraOperation(sgn, con)
        end 

        classes[cl_idx] = class
    end 

    return MatsubaraSymmetryGroup(classes)
end