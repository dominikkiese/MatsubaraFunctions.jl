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

    # convenience constructor
    @inline function Param(
        val  :: Float64, 
        grid :: MatsubaraGrid
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
end