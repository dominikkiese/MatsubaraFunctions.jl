# Note: we do not use the () operator to find the closest index in a MatsubaraGrid
#       in order to to avoid duplicate boundary checks
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
    # we assume that is_inbounds(val, grid) has been checked and evaluates to true
    function Param(
        val  :: Float64, 
        grid :: MatsubaraGrid{Linear}
        )    :: Param

        delta    = value(grid[2]) - value(grid[1])
        position = (val - value(grid[1])) / delta

        # compute nearest-neighbor indices
        # add min to upper index to improve robustness with respect to rounding errors
        dn_idx, up_idx = floor(Int64, position) + 1, min(ceil(Int64, position) + 1, length(grid))

        # compute interpolation weights
        if dn_idx < up_idx
            return Param((dn_idx, up_idx), ((value(grid[up_idx]) - val) / delta, (val - value(grid[dn_idx])) / delta))
        else
            return Param((dn_idx, up_idx), (1.0, 0.0))
        end
    end

    # convenience constructor for coarse grid
    # we assume that is_inbounds(val, grid) has been checked and evaluates to true
    function Param(
        val  :: Float64, 
        grid :: MatsubaraGrid{Coarse}
        )    :: Param

        # compute nearest-neighbor indices
        idx = 1 

        while val > value(grid[idx])
            idx += 1 
        end

        dn_idx, up_idx = idx, idx

        if val < value(grid[idx])
            dn_idx -= 1 
        end 

        # compute interpolation weights
        if dn_idx < up_idx
            delta = value(grid[up_idx]) - value(grid[dn_idx])
            return Param((dn_idx, up_idx), ((value(grid[up_idx]) - val) / delta, (val - value(grid[dn_idx])) / delta))
        else
            return Param((dn_idx, up_idx), (1.0, 0.0))
        end
    end
end