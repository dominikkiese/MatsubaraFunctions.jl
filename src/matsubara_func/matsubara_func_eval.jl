#== 
    Evaluation for MatsubaraFunctions:
        -> tail fitting
        -> extrapolation and boundary conditions
        -> interpolation
        -> summation
==#

# compute tail moments in quadratic approximation from upper bound of 1D MatsubaraFunction
# Note: the distance between interpolation nodes is kept constant below T = 1
"""
    function upper_tail_moments(
        f :: MatsubaraFunction{1, SD, DD, Q},
        x :: Vararg{Int64, SD} 
        ) :: SVector{3, Q} where {SD, DD, Q <: Number}

Returns high frequency moments for quadratic model using upper grid bound. Note, that 
the distance between interpolation nodes is kept constant below T = 1.
"""
function upper_tail_moments(
    f :: MatsubaraFunction{1, SD, DD, Q},
    x :: Vararg{Int64, SD} 
    ) :: SVector{3, Q} where {SD, DD, Q <: Number}

    # compute interpolation nodes
    dist       = ceil(Int64, 1.0 / min(temperature(f.grids[1]), 1.0))
    idx1, idx2 = grids_shape(f, 1) - dist, grids_shape(f, 1) - 2 * dist

    # read data
    ydat = SVector{3, Q}(f.data[end, x...], f.data[idx1, x...], f.data[idx2, x...])
    xdat = SVector{3, Float64}(1.0 / value(f.grids[1][end]), 1.0 / value(f.grids[1][idx1]), 1.0 / value(f.grids[1][idx2]))
    
    # generate Vandermonde matrix 
    mat = @SMatrix Float64[1.0 xdat[1] xdat[1] * xdat[1];
                           1.0 xdat[2] xdat[2] * xdat[2];
                           1.0 xdat[3] xdat[3] * xdat[3]]
    
    return inv(mat) * ydat
end

# compute tail moments in quadratic approximation from lower bound of 1D MatsubaraFunction
# Note: the distance between interpolation nodes is kept constant below T = 1
"""
    function lower_tail_moments(
        f :: MatsubaraFunction{1, SD, DD, Q},
        x :: Vararg{Int64, SD} 
        ) :: SVector{3, Q} where {SD, DD, Q <: Number}

Returns high frequency moments for quadratic model using lower grid bound. Note, that 
the distance between interpolation nodes is kept constant below T = 1.
"""
function lower_tail_moments(
    f :: MatsubaraFunction{1, SD, DD, Q},
    x :: Vararg{Int64, SD} 
    ) :: SVector{3, Q} where {SD, DD, Q <: Number}

    # compute interpolation nodes
    dist       = ceil(Int64, 1.0 / min(temperature(f.grids[1]), 1.0))
    idx1, idx2 = 1 + dist, 1 + 2 * dist

    # read data
    ydat = SVector{3, Q}(f.data[1, x...], f.data[idx1, x...], f.data[idx2, x...])
    xdat = SVector{3, Float64}(1.0 / value(f.grids[1][1]), 1.0 / value(f.grids[1][idx1]), 1.0 / value(f.grids[1][idx2]))
    
    # generate Vandermonde matrix 
    mat = @SMatrix Float64[1.0 xdat[1] xdat[1] * xdat[1];
                           1.0 xdat[2] xdat[2] * xdat[2];
                           1.0 xdat[3] xdat[3] * xdat[3]]
    
    return inv(mat) * ydat
end

# extrapolate 1D Matsubara function using tail moments
function extrapolate(
    f :: MatsubaraFunction{1, SD, DD, Q},
    w :: Float64,
    x :: Vararg{Int64, SD} 
    ) :: Q where {SD, DD, Q <: Number}

    if sign(w) < 0.0 
        moments = lower_tail_moments(f, x...)
        return moments[1] + (moments[2] + moments[3] / w) / w
    else 
        moments = upper_tail_moments(f, x...)
        return moments[1] + (moments[2] + moments[3] / w) / w
    end
end



# wrappers for boundary conditions 
struct BC{IP, OP <: Number}
    f :: FunctionWrappers.FunctionWrapper{OP, Tuple{IP}}
end 

function (bc :: BC{IP, OP})(
    x :: IP
    ) :: OP where {IP, OP <: Number}

    return bc.f(x)
end



# call to MatsubaraFunction with MatsubaraFrequency
# Note: in contrast to the [] operator used for indexing with MatsubaraFrequency, 
#       () has well-defined behavior for out of bounds access (bc or extrapolation)
# Note: if extrp = true, we use polynomial extrapolation for 1D grids and 
#       constant extrapolation for higher dimensional grids
@inline function (f :: MatsubaraFunction{GD, SD, DD, Q})(
    w     :: NTuple{GD, MatsubaraFrequency},
    x     :: Vararg{Int64, SD} 
    ; 
    bc    :: Function = x -> 0.0,
    extrp :: Bool     = false
    )     :: Q where{GD, SD, DD, Q <: Number}

    if any(ntuple(i -> !is_inbounds(w[i], f.grids[i]), GD))
        if extrp 
            if GD == 1
                return extrapolate(f, value(w[1]), x...)
            else 
                idxs = ntuple(i -> max(1, min(grid_index(w[i], f.grids[i]), length(f.grids[i]))), GD)
                return f[idxs..., x ...]
            end 
        else 
            bc_t = BC{NTuple{GD, MatsubaraFrequency}, Q}(bc)
            return bc_t(w) 
        end 
    end

    idxs = ntuple(i -> grid_index(w[i], f.grids[i]), GD)
    return f[idxs..., x ...]
end

function (f :: MatsubaraFunction{1, SD, DD, Q})(
    w     :: MatsubaraFrequency,
    x     :: Vararg{Int64, SD} 
    ; 
    bc    :: Function = x -> 0.0,
    extrp :: Bool     = false
    )     :: Q where{SD, DD, Q <: Number}

    return f((w,), x...; bc = x -> bc(x[1]), extrp)
end



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

# call to MatsubaraFunction with Float64 (multilinear interpolation)
# Note: if extrp = true, we use polynomial extrapolation for 1D grids and 
#       constant extrapolation for higher dimensional grids
@inline function (f :: MatsubaraFunction{GD, SD, DD, Q})(
    w     :: NTuple{GD, Float64},
    x     :: Vararg{Int64, SD} 
    ; 
    bc    :: Function = x -> 0.0,
    extrp :: Bool     = false
    )     :: Q where{GD, SD, DD, Q <: Number}

    if any(ntuple(i -> !is_inbounds(w[i], f.grids[i]), GD))
        if extrp 
            if GD == 1
                return extrapolate(f, w[1], x...)
            else 
                p   = ntuple(i -> Param(max(value(f.grids[i][1]), min(w[i], value(f.grids[i][end]))), f.grids[i]), GD)
                val = 0.0

                for cidx in CartesianIndices(ntuple(i -> 2, GD))
                    wgts  = ntuple(i -> p[i].wgts[cidx[i]], GD)
                    idxs  = ntuple(i -> p[i].idxs[cidx[i]], GD)
                    val  += prod(wgts) * f[idxs..., x...]
                end

                return val
            end
        else 
            bc_t = BC{NTuple{GD, Float64}, Q}(bc)
            return bc_t(w)
        end 
    end 

    p   = ntuple(i -> Param(w[i], f.grids[i]), GD)
    val = 0.0

    for cidx in CartesianIndices(ntuple(i -> 2, GD))
        wgts  = ntuple(i -> p[i].wgts[cidx[i]], GD)
        idxs  = ntuple(i -> p[i].idxs[cidx[i]], GD)
        val  += prod(wgts) * f[idxs..., x...]
    end

    return val
end

function (f :: MatsubaraFunction{1, SD, DD, Q})(
    w     :: Float64,
    x     :: Vararg{Int64, SD} 
    ; 
    bc    :: Function = x -> 0.0,
    extrp :: Bool     = false
    )     :: Q where{SD, DD, Q <: Number}

    return f((w,), x...; bc = x -> bc(x[1]), extrp)
end



# compute Matsubara sum for complex-valued MatsubaraFunction on 1D grid
# Note: only viable if f has Laurent series representation with
#       respect to an annulus about the imaginary axis
"""
    function sum_me(
        f :: MatsubaraFunction{1, SD, DD, Q},
        x :: Vararg{Int64, SD}
        ) :: Q where {SD, DD, Q <: Complex}  

Computes the Matsubara sum (with regulator exp(-iw0+)) for a complex valued MatsubaraFunction on 1D grid. 
This is only viable if f has a Laurent series representation with respect to an annulus about the imaginary axis.
"""
function sum_me(
    f :: MatsubaraFunction{1, SD, DD, Q},
    x :: Vararg{Int64, SD}
    ) :: Q where {SD, DD, Q <: Complex}

    # compute tail moments 
    upper_moments = upper_tail_moments(f, x...)
    lower_moments = lower_tail_moments(f, x...)

    # check self-consistency 
    Δ     = norm(upper_moments .- lower_moments)
    scale = 1e-3 + max(norm(lower_moments), norm(upper_moments)) * 1e-2
    err   = Δ / scale

    @assert err <= 1.0 "Tail fits are inconsistent! Try more frequencies or check prerequisites"

    # compute expansion coefficients
    α0 = +0.5 * (upper_moments[1] + lower_moments[1])
    α1 = +0.5 * (upper_moments[2] + lower_moments[2]) * im
    α2 = -0.5 * (upper_moments[3] + lower_moments[3])

    # compute the Matsubara sum using quadratic asymptotic model
    T   = temperature(f.grids[1])
    num = grids_shape(f, 1)
    val = -T * (num * α0 - sum(@view f.data[:, x...])) - 0.5 * (α1 + 0.5 * α2 / T)

    for w in 1 : num
        val += T * α2 / value(f.grids[1][w]) / value(f.grids[1][w])
    end

    return val
end