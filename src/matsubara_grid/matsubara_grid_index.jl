# make Matsubara grid indexable
function Base.:eachindex(
    grid :: MatsubaraGrid
    )    :: Base.OneTo{Int64}

    return eachindex(grid.data)
end

function Base.:lastindex(
    grid :: MatsubaraGrid
    )    :: Int64

    return lastindex(grid.data)
end

function Base.:getindex(
    grid :: MatsubaraGrid,
    idx  :: Int64 
    )    :: MatsubaraFrequency

    # bounds check performed by Base
    return grid.data[idx]
end

function Base.:getindex(
    grid :: MatsubaraGrid,
    idxs :: UnitRange{Int64}
    )    :: SubArray{MatsubaraFrequency, 1, Vector{MatsubaraFrequency}, Tuple{UnitRange{Int64}}, true}

    # bounds check performed by Base
    return @view grid.data[idxs]
end  



# unsafe method for converting MatsubaraFrequency to grid index (no bounds check)
function grid_index(
    w    :: MatsubaraFrequency,
    grid :: MatsubaraGrid
    )    :: Int64 

    return index(w) - index(grid[1]) + 1
end

# safer method for converting MatsubaraFrequency to grid index
# (i.e. frequencies which are out of bounds will be reset to mesh boundaries)
function grid_index_extrp(
    w    :: MatsubaraFrequency,
    grid :: MatsubaraGrid
    )    :: Int64 

    return max(1, min(grid_index(w, grid), length(grid)))
end

# make MatsubaraGrid callable with MatsubaraFrequency
# returns index to data array corresponding to this frequency if in grid
function (f :: MatsubaraGrid)(
    w :: MatsubaraFrequency
    ) :: Int64

    if is_inbounds(w, f)
        return grid_index(w, f)
    else 
        error("Frequency not in grid")
    end 
end

# make MatsubaraGrid callable with Float64
# returns index to data array corresponding to closest frequency if in grid
function (f :: MatsubaraGrid)(
    w :: Float64
    ) :: Int64

    if is_inbounds(w, f)
        delta    = value(f[2]) - value(f[1])
        position = (w - value(f[1])) / delta
        return round(Int64, position) + 1
    else 
        error("Frequency not in grid")
    end 
end



# make Matsubara grid iterable
function Base.:iterate(
    grid :: MatsubaraGrid
    )    :: Tuple{MatsubaraFrequency, Int64}

    return grid[1], 1 
end

function Base.:iterate(
    grid  :: MatsubaraGrid,
    state :: Int64
    )     :: Union{Nothing, Tuple{MatsubaraFrequency, Int64}}

    if state < length(grid)
        return grid[state + 1], state + 1 
    else 
        return nothing 
    end
end