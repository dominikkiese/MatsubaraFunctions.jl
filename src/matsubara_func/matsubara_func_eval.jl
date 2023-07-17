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
        f  :: MatsubaraFunction{1, SD, DD, Q},
        α0 :: Q,
        x  :: Vararg{Int64, SD} 
        )  :: Tuple{Q, Q} where {SD, DD, Q <: Number}

Returns high frequency moments for quadratic model using upper grid bound. Here, α0 is the
asymptotic limit for large positive frequencies. Note, that the distance between interpolation nodes 
is kept constant below T = 1.
"""
function upper_tail_moments(
    f  :: MatsubaraFunction{1, SD, DD, Q},
    α0 :: Q,
    x  :: Vararg{Int64, SD} 
    )  :: Tuple{Q, Q} where {SD, DD, Q <: Number}

    # compute interpolation nodes
    dist = ceil(Int64, 1.0 / min(temperature(f.grids[1]), 1.0))
    idx  = grids_shape(f, 1) - dist
    @check idx > ceil(Int64, 0.75 * grids_shape(f, 1)) "Grid is too small for extrapolation"

    # read data
    y1, y2     = f.data[end, x...] - α0, f.data[idx, x...] - α0
    x1, x2     = 1.0 / value(f.grids[1][end]), 1.0 / value(f.grids[1][idx])
    x1sq, x2sq = x1 * x1, x2 * x2
    dtinv      = 1.0 / (x1 * x2sq - x2 * x1sq)

    return dtinv * (x2sq * y1 - x1sq * y2), dtinv * (x1 * y2 - x2 * y1)
end

# compute tail moments in quadratic approximation from lower bound of 1D MatsubaraFunction
# Note: the distance between interpolation nodes is kept constant below T = 1
"""
    function lower_tail_moments(
        f  :: MatsubaraFunction{1, SD, DD, Q},
        α0 :: Q,
        x  :: Vararg{Int64, SD} 
        )  :: Tuple{Q, Q} where {SD, DD, Q <: Number}

Returns high frequency moments for quadratic model using lower grid bound. Here, α0 is the
asymptotic limit for large negative frequencies. Note, that the distance between interpolation nodes 
is kept constant below T = 1.
"""
function lower_tail_moments(
    f  :: MatsubaraFunction{1, SD, DD, Q},
    α0 :: Q,
    x  :: Vararg{Int64, SD} 
    )  :: Tuple{Q, Q} where {SD, DD, Q <: Number}

    # compute interpolation nodes
    dist = ceil(Int64, 1.0 / min(temperature(f.grids[1]), 1.0))
    idx  = 1 + dist
    @check idx < floor(Int64, 0.25 * grids_shape(f, 1)) "Grid is too small for extrapolation"

    # read data
    y1, y2     = f.data[1, x...] - α0, f.data[idx, x...] - α0
    x1, x2     = 1.0 / value(f.grids[1][1]), 1.0 / value(f.grids[1][idx])
    x1sq, x2sq = x1 * x1, x2 * x2
    dtinv      = 1.0 / (x1 * x2sq - x2 * x1sq)

    return dtinv * (x2sq * y1 - x1sq * y2), dtinv * (x1 * y2 - x2 * y1)
end

# extrapolate 1D Matsubara function using tail moments
function extrapolate(
    f  :: MatsubaraFunction{1, SD, DD, Q},
    w  :: Float64,
    α0 :: Q,
    x  :: Vararg{Int64, SD} 
    )  :: Q where {SD, DD, Q <: Number}

    if sign(w) < 0.0 
        moments = lower_tail_moments(f, α0, x...)
        return α0 + (moments[1] + moments[2] / w) / w
    else 
        moments = upper_tail_moments(f, α0, x...)
        return α0 + (moments[1] + moments[2] / w) / w
    end
end



# call to MatsubaraFunction with MatsubaraFrequency
# Note: in contrast to the [] operator used for indexing with MatsubaraFrequency, 
#       () has well-defined behavior for out of bounds access (bc or extrapolation)
# Note: if extrp[1] = true, we use polynomial extrapolation for 1D grids and 
#       constant extrapolation for higher dimensional grids. extrp[2] sets the value
#       for the asymptotic limit in the 1D case, but it has no effect for higher dimensional grids
function (f :: MatsubaraFunction{GD, SD, DD, Q})(
    w     :: NTuple{GD, MatsubaraFrequency},
    x     :: Vararg{Int64, SD} 
    ; 
    bc    :: Function       = y -> 0.0,
    extrp :: Tuple{Bool, Q} = (false, Q(0.0))
    )     :: Q where{GD, SD, DD, Q <: Number}

    if any(ntuple(i -> !is_inbounds(w[i], f.grids[i]), GD))
        if extrp[1] 
            if GD == 1
                return extrapolate(f, value(w[1]), extrp[2], x...)
            else 
                return f[CartesianIndex_extrp(f, w, x...)]
            end 
        else 
            return bc(w) 
        end 
    end

    return f[ntuple(i -> grid_index(w[i], f.grids[i]), GD)..., x ...]
end

function (f :: MatsubaraFunction{1, SD, DD, Q})(
    w     :: MatsubaraFrequency,
    x     :: Vararg{Int64, SD} 
    ; 
    bc    :: Function       = y -> 0.0,
    extrp :: Tuple{Bool, Q} = (false, Q(0.0))
    )     :: Q where{SD, DD, Q <: Number}

    return f((w,), x...; bc = y -> bc(y[1]), extrp)
end

function (f :: MatsubaraFunction{GD, 1, DD, Q})(
    w     :: Vararg{MatsubaraFrequency, GD},
    ; 
    bc    :: Function       = y -> 0.0,
    extrp :: Tuple{Bool, Q} = (false, Q(0.0))
    )     :: Q where{GD, DD, Q <: Number}

    @check shape(f, 1) == 1 "MatsubaraFunction is not scalar but vector valued"
    return f((w...,), 1; bc = y -> bc(y), extrp)
end

# specialization to make boundary condition more convenient
function (f :: MatsubaraFunction{1, 1, 2, Q})(
    w     :: MatsubaraFrequency,
    ; 
    bc    :: Function       = y -> 0.0,
    extrp :: Tuple{Bool, Q} = (false, Q(0.0))
    )     :: Q where{Q <: Number}

    @check shape(f, 1) == 1 "MatsubaraFunction is not scalar but vector valued"
    return f((w,), 1; bc = y -> bc(y[1]), extrp)
end



# Note: we do not use the () operator to find the closest index in a MatsubaraGrid
#       in order to to avoid duplicate boundary checks
struct Param
    idxs :: NTuple{2, Int64}
    wgts :: NTuple{2, Float64}

    # default constructor
    function Param(
        idxs :: NTuple{2, Int64}, 
        wgts :: NTuple{2, Float64}
        )    :: Param

        return new(idxs, wgts)
    end

    # convenience constructor
    function Param(
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
function (f :: MatsubaraFunction{GD, SD, DD, Q})(
    w     :: NTuple{GD, Float64},
    x     :: Vararg{Int64, SD} 
    ; 
    bc    :: Function       = y -> 0.0,
    extrp :: Tuple{Bool, Q} = (false, Q(0.0))
    )     :: Q where{GD, SD, DD, Q <: Number}

    if any(ntuple(i -> !is_inbounds(w[i], f.grids[i]), GD))
        if extrp[1] 
            if GD == 1
                return extrapolate(f, w[1], extrp[2], x...)
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
            return bc(w)
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
    bc    :: Function       = y -> 0.0,
    extrp :: Tuple{Bool, Q} = (false, Q(0.0))
    )     :: Q where{SD, DD, Q <: Number}

    return f((w,), x...; bc = y -> bc(y[1]), extrp)
end

function (f :: MatsubaraFunction{GD, 1, DD, Q})(
    w     :: Vararg{Float64, GD},
    ; 
    bc    :: Function       = y -> 0.0,
    extrp :: Tuple{Bool, Q} = (false, Q(0.0))
    )     :: Q where{GD, DD, Q <: Number}

    @check shape(f, 1) == 1 "MatsubaraFunction is not scalar but vector valued"
    return f((w...,), 1; bc = y -> bc(y), extrp)
end

# specialization to make boundary condition more convenient
function (f :: MatsubaraFunction{1, 1, 2, Q})(
    w     :: Float64,
    ; 
    bc    :: Function       = y -> 0.0,
    extrp :: Tuple{Bool, Q} = (false, Q(0.0))
    )     :: Q where{Q <: Number}

    @check shape(f, 1) == 1 "MatsubaraFunction is not scalar but vector valued"
    return f((w,), 1; bc = y -> bc(y[1]), extrp)
end



# compute Matsubara sum for complex-valued MatsubaraFunction on 1D grid
# Note: only viable if f has Laurent series representation with
#       respect to an annulus about the imaginary axis
"""
    function sum_me(
        f  :: MatsubaraFunction{1, SD, DD, Q},
        α0 :: Q,
        x  :: Vararg{Int64, SD}
        )  :: Q where {SD, DD, Q <: Complex}

Computes the Matsubara sum (with regulator exp(-iw0+)) for a complex valued MatsubaraFunction on a 1D grid. Here, `α0`
is the asymptotic limit for large frequencies. This is only viable if `f1` has a Laurent series representation with respect 
to an annulus about the imaginary axis.
"""
function sum_me(
    f  :: MatsubaraFunction{1, SD, DD, Q},
    α0 :: Q,
    x  :: Vararg{Int64, SD}
    )  :: Q where {SD, DD, Q <: Complex}

    # compute tail moments 
    upper_moments = upper_tail_moments(f, α0, x...); upper_max = max(abs.(upper_moments)...)
    lower_moments = lower_tail_moments(f, α0, x...); lower_max = max(abs.(lower_moments)...)

    # check self-consistency 
    diff  = max(abs.(upper_moments .- lower_moments)...);
    scale = 1e-3 + max(upper_max, lower_max) * 1e-2;
    err   = diff / scale

    if err >= 1.0 
        @warn "Tail fits are inconsistent: upper_moments = $(upper_moments), lower_moments = $(lower_moments)" 
    end

    # compute expansion coefficients
    α1 = +0.5 * (upper_moments[1] + lower_moments[1]) * im
    α2 = -0.5 * (upper_moments[2] + lower_moments[2])

    # compute the Matsubara sum using quadratic asymptotic model
    T   = temperature(f.grids[1])
    num = grids_shape(f, 1)
    val = -T * (num * α0 - sum(@view f.data[:, x...])) - 0.5 * (α1 + 0.5 * α2 / T)

    for w in 1 : num
        val += T * α2 / value(f.grids[1][w]) / value(f.grids[1][w])
    end

    return val
end

"""
    function sum_me(
        f  :: MatsubaraFunction{1, 1, 2, Q},
        α0 :: Q
        )  :: Q where {Q <: Complex}

Computes the Matsubara sum (with regulator exp(-iw0+)) for a complex valued MatsubaraFunction on a 1D grid. Here, `α0`
is the asymptotic limit for large frequencies. This is only viable if `f1` has a Laurent series representation with respect 
to an annulus about the imaginary axis. Requires `shape(f, 1) == 1`.
"""
function sum_me(
    f  :: MatsubaraFunction{1, 1, 2, Q},
    α0 :: Q
    )  :: Q where {Q <: Complex}

    @check shape(f, 1) == 1 "MatsubaraFunction is not scalar but vector valued"
    return sum_me(f, α0, 1)
end

"""
    function density(
        f :: MatsubaraFunction{1, SD, DD, Q},
        x :: Vararg{Int64, SD}
        ) :: Q where {SD, DD, Q <: Complex}

Computes the density for a complex valued MatsubaraFunction on a 1D grid. Assumes that `f` decays to zero for large 
frequencies (as a single-particle Green's function would).
"""
function density(
    f :: MatsubaraFunction{1, SD, DD, Q},
    x :: Vararg{Int64, SD}
    ) :: Q where {SD, DD, Q <: Complex}

    return 1.0 + sum_me(f, ComplexF64(0.0), x...)
end

"""
    function density(
        f :: MatsubaraFunction{1, 1, 2, Q}
        ) :: Q where {Q <: Complex}

Computes the density for a complex valued MatsubaraFunction on a 1D grid. Assumes that `f` decays to zero for large 
frequencies (as a single-particle Green's function would). Requires `shape(f, 1) == 1`.
"""
function density(
    f :: MatsubaraFunction{1, 1, 2, Q}
    ) :: Q where {Q <: Complex}

    @check shape(f, 1) == 1 "MatsubaraFunction is not scalar but vector valued"
    return 1.0 + sum_me(f, ComplexF64(0.0), 1)
end