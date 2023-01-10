function save_matsubara_grid!(
    h :: HDF5.File,
    l :: String,
    g :: MatsubaraGrid{GT}
    ) :: Nothing where {GT <: AbstractGrid}

    attributes(h)[l * "/T"]    = g.T
    attributes(h)[l * "/type"] = "$(g.type)"
    attributes(h)[l * "/GT"]   = "$(GT)"
    h[l * "/data"]             = g.data

    return nothing 
end 

function load_matsubara_grid(
    h :: HDF5.File,
    l :: String
    ) :: MatsubaraGrid 

    T    = read_attribute(h, l * "/T")
    type = read_attribute(h, l * "/type")
    GT   = read_attribute(h, l * "/GT")
    data = read(h, l * "/data")

    return MatsubaraGrid(T, data, Symbol(type), eval(Meta.parse(GT)))
end

function save_matsubara_function!(
    h :: HDF5.File,
    l :: String,
    f :: MatsubaraFunction{Dg, Ds, Dt, GT, Q}
    ) :: Nothing where {Dg, Ds, Dt, GT <: AbstractGrid, Q <: Number}

    for i in eachindex(f.grids)
        save_matsubara_grid!(h, l * "/grids/grid_$i", f.grids[i])
    end

    attributes(h)[l * "/shape"] = Int64[f.shape...] 
    h[l * "/data"]              = f.data

    return nothing 
end

function load_matsubara_function(
    h :: HDF5.File,
    l :: String
    ) :: MatsubaraFunction

    idxs  = eachindex(keys(h[l * "/grids"]))
    grids = MatsubaraGrid[load_matsubara_grid(h, l * "/grids/grid_$i") for i in idxs]
    shape = read_attribute(h, l * "/shape")
    data  = read(h, l * "/data")
    
    return MatsubaraFunction((grids...,), (shape...,), data)
end