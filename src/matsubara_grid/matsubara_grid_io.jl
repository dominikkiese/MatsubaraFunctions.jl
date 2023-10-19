"""
    function save_matsubara_grid!(
        h :: HDF5.File,
        l :: String,
        g :: AbstractMatsubaraGrid
        ) :: Nothing

Save MatsubaraGrid `g` with label `l` to file `h`
"""
function save_matsubara_grid!(
    h :: HDF5.File,
    l :: String,
    g :: AbstractMatsubaraGrid
    ) :: Nothing

    grp = create_group(h, l)

    # save only metadata to minimize memory footprint
    attributes(grp)["type"] = string(typeof(g))
    attributes(grp)["T"]    = temperature(g)
    attributes(grp)["N"]    = N(g)

    return nothing 
end 

"""
    function load_matsubara_grid(
        h :: HDF5.File,
        l :: String
        ) :: AbstractMatsubaraGrid

Load MatsubaraGrid with label `l` from file `h`
"""
function load_matsubara_grid(
    h :: HDF5.File,
    l :: String
    ) :: MatsubaraGrid 

    # load the metadata
    type = read_attribute(h[l], "type")
    T    = read_attribute(h[l], "T")
    N    = read_attribute(h[l], "N")

    # generate the grid
    if type == "MatsubaraGrid{Fermion}"
        return MatsubaraGrid(T, N, Fermion)
    elseif type == "MatsubaraGrid{Boson}"
        return MatsubaraGrid(T, N, Boson)
    else 
        error("Type $(type) unknown")
    end
end

#----------------------------------------------------------------------------------------------#

export 
    save_matsubara_grid!,
    load_matsubara_grid