"""
    function save_matsubara_grid!(
        h :: HDF5.File,
        l :: String,
        g :: MatsubaraGrid
        ) :: Nothing

Save MatsubaraGrid `g` with label `l` to file `h`
"""
function save_matsubara_grid!(
    h :: HDF5.File,
    l :: String,
    g :: MatsubaraGrid
    ) :: Nothing

    # create new group and tag it
    grp                    = create_group(h, l)
    attributes(grp)["tag"] = "MatsubaraGrid"

    # add metadata
    attributes(grp)["T"]    = temperature(g)
    attributes(grp)["type"] = "$(type(g))"

    # add data
    grp["frequency/vals"] = values(g)
    grp["frequency/idxs"] = indices(g)

    return nothing 
end 

"""
    function load_matsubara_grid(
        h :: HDF5.File,
        l :: String
        ) :: MatsubaraGrid

Load MatsubaraGrid with label `l` from file `h`
"""
function load_matsubara_grid(
    h :: HDF5.File,
    l :: String
    ) :: MatsubaraGrid 

    # check if group has correct tag 
    @check haskey(attributes(h[l]), "tag") "Group $(l) not compatible with MatsubaraFunctions"
    @check read_attribute(h[l], "tag") == "MatsubaraGrid" "Group $(l) not tagged as MatsubaraGrid"

    # read the metadata
    T    = read_attribute(h[l], "T")
    type = read_attribute(h[l], "type")

    # read the data 
    vals = read(h, l * "/frequency/vals")
    idxs = read(h, l * "/frequency/idxs")
    data = MatsubaraFrequency[MatsubaraFrequency(T, val, idx, Symbol(type)) for (val, idx) in zip(vals, idxs)]

    return MatsubaraGrid(T, data, Symbol(type))
end