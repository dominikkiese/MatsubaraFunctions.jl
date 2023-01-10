struct Param
    idxs :: NTuple{2, Int64}
    wgts :: NTuple{2, Float64}

    # basic constructor
    function Param(
        idxs :: NTuple{2, Int64}, 
        wgts :: NTuple{2, Float64}
        )    :: Param

        return new(idxs, wgts)
    end

    # convenience constructor for linear grid
    # we assume that in(val, grid) has been checked and evaluates to true
    function Param(
        val  :: Float64, 
        grid :: MatsubaraGrid{Linear}
        )    :: Param

        delta    = grid[2] - grid[1]
        position = (val - grid[1]) / delta

        # compute nearest-neighbor indices
        # add min to upper index to improve robustness with respect to rounding errors
        dn_idx, up_idx = floor(Int64, position) + 1, min(ceil(Int64, position) + 1, length(grid))

        # compute interpolation weights
        if dn_idx < up_idx
            return Param((dn_idx, up_idx), ((grid[up_idx] - val) / delta, (val - grid[dn_idx]) / delta))
        else
            return Param((dn_idx, up_idx), (1.0, 0.0))
        end
    end

    # convenience constructor for coarse grid
    # we assume that in(val, grid) has been checked and evaluates to true
    function Param(
        val  :: Float64, 
        grid :: MatsubaraGrid{Coarse}
        )    :: Param

        # compute nearest-neighbor indices
        idx = 1 

        while val > grid[idx]
            idx += 1 
        end

        dn_idx, up_idx = idx, idx

        if val < grid[idx]
            dn_idx -= 1 
        end 

        # compute interpolation weights
        if dn_idx < up_idx
            delta = grid[up_idx] - grid[dn_idx]
            return Param((dn_idx, up_idx), ((grid[up_idx] - val) / delta, (val - grid[dn_idx]) / delta))
        else
            return Param((dn_idx, up_idx), (1.0, 0.0))
        end
    end
end