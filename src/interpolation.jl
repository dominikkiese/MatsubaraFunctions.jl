struct Param
    idxs :: NTuple{2, Int64}
    wgts :: NTuple{2, Float64}

    function Param(
        idxs :: NTuple{2, Int64}, 
        wgts :: NTuple{2, Float64}
        )    :: Param

        return new(idxs, wgts)
    end

    # convenience constructor (val is assumed to be in grid)
    function Param(
        val  :: Float64, 
        grid :: MatsubaraGrid
        )    :: Param

        delta    = grid[2] - grid[1]
        position = (val - grid[1]) / delta

        # compute nearest-neighbor indices
        # add min to upper index to improve robustness with respect to rounding errors
        idxs = floor(Int64, position) + 1, min(ceil(Int64, position) + 1, length(grid))

        if idxs[1] < idxs[2]
            return Param(idxs, ((grid[idxs[2]] - val) / delta, (val - grid[idxs[1]]) / delta))
        else
            return Param(idxs, (1.0, 0.0))
        end
    end
end