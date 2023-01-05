function save_matsubara_function!(
    f :: MatsubaraFunction{Ds, Dg, Dt},
    l :: String,
    h :: HDF5.File
    ) :: Nothing where {Ds, Dg, Dt}

    h[l * "/shape"] = Int64[f.shape...] 
    h[l * "/data"]  = f.data

    for i in eachindex(f.grids)
        h[l * "/grids/grid_$i"] = f.grids[i]
    end 

    return nothing 
end

function load_matsubara_function(
    l :: String,
    h :: HDF5.File
    ) :: MatsubaraFunction

    shape = read(h, l * "/shape")
    data  = read(h, l * "/data")
    grids = Vector{Float64}[read(h, l * "/grids/grid_$i") for i in eachindex(keys(h[l * "/grids"]))]
    
    return MatsubaraFunction((shape...,), (grids...,), data)
end