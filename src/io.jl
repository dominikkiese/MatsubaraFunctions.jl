function save_matsubara_grid!(
    g :: MatsubaraGrid, 
    l :: String,
    h :: HDF5.File 
    ) :: Nothing 

    attributes(h)[l * "/T"]    = g.T
    attributes(h)[l * "/type"] = "$(g.type)"
    h[l * "/data"]             = g.data

    return nothing 
end 

function load_matsubara_grid(
    l :: String, 
    h :: HDF5.File 
    ) :: MatsubaraGrid 

    T    = read_attribute(h, l * "/T")
    type = read_attribute(h, l * "/type")
    data = read(h, l * "/data")

    if type == "Fermion"
        return MatsubaraGrid(T, data, :Fermion)
    elseif type == "Boson"
        return MatsubaraGrid(T, data, :Boson)
    else 
        error("Grid type $(type) unknown")
    end 
end

function save_matsubara_function!(
    f :: MatsubaraFunction{Dg, Ds, Dt},
    l :: String,
    h :: HDF5.File
    ) :: Nothing where {Dg, Ds, Dt}

    for i in eachindex(f.grids)
        save_matsubara_grid!(f.grids[i], l * "/grids/grid_$i", h)
    end

    h[l * "/shape"] = Int64[f.shape...] 
    h[l * "/data"]  = f.data

    return nothing 
end

function load_matsubara_function(
    l :: String,
    h :: HDF5.File
    ) :: MatsubaraFunction

    idxs  = eachindex(keys(h[l * "/grids"]))
    grids = MatsubaraGrid[load_matsubara_grid(l * "/grids/grid_$i", h) for i in idxs]
    shape = read(h, l * "/shape")
    data  = read(h, l * "/data")
    
    return MatsubaraFunction((grids...,), (shape...,), data)
end