struct Param
    idxs :: NTuple{2, Int64}
    wgts :: NTuple{2, Float64}

    function Param(
        idxs :: NTuple{2, Int64}, 
        wgts :: NTuple{2, Float64}
        )    :: Param

        return new(idxs, wgts)
    end

    # Note: grid is assumed to be sorted & linear and val is assumed to be in grid
    # convenience constructor
    function Param(
        val  :: Float64, 
        grid :: Vector{Float64}
        )    :: Param

        delta    = grid[2] - grid[1]
        position = (val - grid[1]) / delta
        idxs     = floor(Int64, position) + 1, ceil(Int64, position) + 1

        if idxs[1] < idxs[2]
            return Param(idxs, ((grid[idxs[2]] - val) / delta, (val - grid[idxs[1]]) / delta))
        else
            return Param(idxs, (1.0, 0.0))
        end
    end
end