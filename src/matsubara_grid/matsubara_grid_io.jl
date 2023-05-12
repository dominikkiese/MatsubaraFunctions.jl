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