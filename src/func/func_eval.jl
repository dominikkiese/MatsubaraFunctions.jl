# compute tail moments in quadratic approximation from upper bound of 1D MatsubaraFunction
# Note: the distance between interpolation nodes is kept constant below T = 1
"""
    function upper_tail_moments(
        f  :: MeshFunction{1, SD, DD, Q},
        α0 :: Q,
        x  :: Vararg{Int64, SD} 
        )  :: Tuple{Q, Q} where {SD, DD, Q <: Number}

Returns high frequency moments for quadratic model using upper grid bound. Here, α0 is the
asymptotic limit for large positive frequencies.
"""
function upper_tail_moments(
    f  :: MeshFunction{1, SD, DD, Q},
    α0 :: Q,
    x  :: Vararg{Int64, SD} 
    )  :: Tuple{Q, Q} where {SD, DD, Q <: Number}

    # compute nodes
    dist = ceil(Int64, 1.0 / min(temperature(firstvalue(meshes(f, 1))), 1.0))
    idx  = lastindex(meshes(f, 1)) - dist
    @DEBUG idx > ceil(Int64, 0.75 * lastindex(meshes(f, 1))) "Grid is too small for extrapolation"

    # read data
    y1, y2     = f.data[end, x...] - α0, f.data[idx, x...] - α0
    x1, x2     = 1.0 / plain_value(meshes(f, 1)[end]), 1.0 / plain_value(meshes(f, 1)[idx])
    x1sq, x2sq = x1 * x1, x2 * x2
    dtinv      = 1.0 / (x1 * x2sq - x2 * x1sq)

    return dtinv * (x2sq * y1 - x1sq * y2), dtinv * (x1 * y2 - x2 * y1)
end

# compute tail moments in quadratic approximation from lower bound of 1D MatsubaraFunction
# Note: the distance between interpolation nodes is kept constant below T = 1
"""
    function lower_tail_moments(
        f  :: MeshFunction{1, SD, DD, Q},
        α0 :: Q,
        x  :: Vararg{Int64, SD} 
        )  :: Tuple{Q, Q} where {SD, DD, Q <: Number}

Returns high frequency moments for quadratic model using lower grid bound. Here, α0 is the
asymptotic limit for large negative frequencies.
"""
function lower_tail_moments(
    f  :: MeshFunction{1, SD, DD, Q},
    α0 :: Q,
    x  :: Vararg{Int64, SD} 
    )  :: Tuple{Q, Q} where {SD, DD, Q <: Number}

    # compute nodes
    dist = ceil(Int64, 1.0 / min(temperature(meshes(f, 1)), 1.0))
    idx1 = firstindex(meshes(f, 1))
    idx2 = idx1 + dist
    @DEBUG idx2 < floor(Int64, 0.75 * idx1) "Grid is too small for extrapolation"

    # read data
    y1, y2     = f.data[idx1, x...] - α0, f.data[idx2, x...] - α0
    x1, x2     = 1.0 / value(meshes(f, 1)[idx1]), 1.0 / value(meshes(f, 1)[idx2])
    x1sq, x2sq = x1 * x1, x2 * x2
    dtinv      = 1.0 / (x1 * x2sq - x2 * x1sq)

    return dtinv * (x2sq * y1 - x1sq * y2), dtinv * (x1 * y2 - x2 * y1)
end

function extrapolate(
    f  :: MeshFunction{1, SD, DD, Q},
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

#----------------------------------------------------------------------------------------------#

# call to MeshFunction, only value type and mesh point, no interpolation
# Note: in contrast to the [] operator used for indexing with MatsubaraFrequency, 
#       () has well-defined behavior for out of bounds access (extrapolation)
function (f :: MeshFunction{MD, SD, DD, Q})(
    p :: NTuple{MD, Union{AbstractValue, AbstractMeshPoint}},
    x :: Vararg{Int64, SD} 
    ) :: Q where{MD, SD, DD, Q <: Number}

    return f[CartesianIndex_bc(f, p, x...)]
end

# specialization
function (f :: MeshFunction{1, SD, DD, Q})(
    p :: Tuple{AbstractValue},
    x :: Vararg{Int64, SD} 
    ; 
    extrp :: Q = Q(0.0)
    ) :: Q where{SD, DD, Q <: Number}

    if !is_inbounds(p[1], meshes(f, 1))
        return extrapolate(f, value(p[1]), extrp, x...)
    end

    return f[mesh_index(p[1], meshes(f, 1)), x ...]
end

function (f :: MeshFunction{1, SD, DD, Q})(
    p :: Tuple{AbstractMeshPoint},
    x :: Vararg{Int64, SD} 
    ; 
    extrp :: Q = Q(0.0)
    ) :: Q where{SD, DD, Q <: Number}

    return f((value(p[1]),), x...; extrp)
end


function (f :: MeshFunction{1, SD, DD, Q})(
    p :: Union{AbstractValue, AbstractMeshPoint},
    x :: Vararg{Int64, SD} 
    ; 
    extrp :: Q = Q(0.0)
    ) :: Q where{SD, DD, Q <: Number}

    return f((p,), x...; extrp)
end

function (f :: MeshFunction{MD, 0, DD, Q})(
    p :: Vararg{Union{AbstractValue, AbstractMeshPoint}, MD} 
    ) :: Q where{MD, DD, Q <: Number}

    return f((p...,))
end

# specialize 
function (f :: MeshFunction{1, 0, DD, Q})(
    p :: Union{AbstractValue, AbstractMeshPoint}
    ; 
    extrp :: Q = Q(0.0)
    ) :: Q where{DD, Q <: Number}

    return f((p,); extrp)
end

# call to MeshFunction, all types, interpolation
function (f :: MeshFunction{MD, SD, DD, Q})(
    p :: NTuple{MD, Union{AbstractValue, AbstractMeshPoint, Float64, SVector}},
    x :: Vararg{Int64, SD} 
    ) :: Q where{MD, SD, DD, Q <: Number}

    params    = ntuple(n -> InterpolationParam(p[n], meshes(f, n)), MD)
    idx_iters = Iterators.product(ntuple(n -> indices(params[n]), MD)...)
    wgt_iters = Iterators.product(ntuple(n -> weights(params[n]), MD)...)
    val       = 0.0

    for idx in zip(wgt_iters, idx_iters)
        val += prod(first(idx)) * f[last(idx)..., x...]
    end

    return val
end

# extrapolate 1D mesh functions 
function (f :: MeshFunction{1, SD, DD, Q})(
    p     :: Tuple{Float64},
    x     :: Vararg{Int64, SD} 
    ; 
    extrp :: Q = Q(0.0)
    )     :: Q where{SD, DD, Q <: Number}

    if !is_inbounds(p[1], meshes(f, 1))
        return extrapolate(f, p[1], extrp, x...)
    end 

    params = InterpolationParam(p[1], meshes(f, 1))
    return params.weights[1] * f[params.indices[1], x...] + params.weights[2] * f[params.indices[2], x...]
end

function (f :: MeshFunction{1, SD, DD, Q})(
    w     :: Union{AbstractValue, AbstractMeshPoint, Float64, SVector},
    x     :: Vararg{Int64, SD} 
    ; 
    extrp :: Q = Q(0.0)
    )     :: Q where{SD, DD, Q <: Number}

    return f((w,), x...; extrp)
end



function (f :: MeshFunction{MD, 0, DD, Q})(
    p :: Vararg{Union{AbstractValue, AbstractMeshPoint, Float64, SVector}, MD} 
    ) :: Q where{MD, DD, Q <: Number}

    return f((p...,))
end

function (f :: MeshFunction{1, 0, DD, Q})(
    w     :: Union{AbstractValue, AbstractMeshPoint, Float64, SVector}
    ;
    extrp :: Q = Q(0.0)
    )     :: Q where{DD, Q <: Number}

    return f((w,); extrp)
end