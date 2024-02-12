# call to MeshFunction, only value type and mesh point, no interpolation
function (f :: MeshFunction{MD, SD, DD, Q, AT})(
    p   :: NTuple{MD, Union{<: AbstractValue, <: AbstractMeshPoint}},
    x   :: Vararg{Int64, SD}
    ;
    lim :: Q = Q(0.0)
    ) where{MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return any(i -> !is_inbounds_bc(p[i], meshes(f, i)), MD) ? lim : f[CartesianIndex_bc(f, p, x...)]
end

function (f :: MeshFunction{1, SD, DD, Q, AT})(
    p   :: Union{<: AbstractValue, <: AbstractMeshPoint},
    x   :: Vararg{Int64, SD} 
    ;
    lim :: Q = Q(0.0)
    ) where{SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return f((p,), x...; lim)
end

function (f :: MeshFunction{MD, 0, DD, Q, AT})(p :: Vararg{Union{<: AbstractValue, <: AbstractMeshPoint}, MD}; lim :: Q = Q(0.0)
    ) where{MD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return f((p...,); lim)
end

# specialize 
function (f :: MeshFunction{1, 0, DD, Q, AT})(p :: Union{<: AbstractValue, <: AbstractMeshPoint}; lim :: Q = Q(0.0)
    ) where{DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return f((p,); lim)
end

# call to MeshFunction, all types, interpolation
function (f :: MeshFunction{MD, SD, DD, Q, AT})(
    p   :: NTuple{MD, Union{<: AbstractValue, <: AbstractMeshPoint, Float64, <: AbstractVector{Float64}}},
    x   :: Vararg{Int64, SD} 
    ;
    lim :: Q = Q(0.0)
    ) where{MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    if any(i -> !is_inbounds_bc(p[i], meshes(f, i)), MD)
        return lim 
    end

    params    = map((y, m) -> InterpolationParam(y, m), p, meshes(f))
    idx_iters = Iterators.product(indices.(params)...)
    wgt_iters = Iterators.product(weights.(params)...)
    val       = 0.0

    for y in zip(wgt_iters, idx_iters)
        val += prod(first(y)) * f[last(y)..., x...]
    end

    return val
end

function (f :: MeshFunction{1, SD, DD, Q, AT})(
    p   :: Union{<: AbstractValue, <: AbstractMeshPoint, Float64, <: AbstractVector{Float64}},
    x   :: Vararg{Int64, SD} 
    ;
    lim :: Q = Q(0.0)
    ) where{SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return f((p,), x...; lim)
end

function (f :: MeshFunction{MD, 0, DD, Q, AT})(p :: Vararg{Union{<: AbstractValue, <: AbstractMeshPoint, Float64, <: AbstractVector{Float64}}, MD}; lim :: Q = Q(0.0)
    ) where{MD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return f((p...,); lim)
end

# specialize 
function (f :: MeshFunction{1, 0, DD, Q, AT})(p :: Union{<: AbstractValue, <: AbstractMeshPoint, Float64, <: AbstractVector{Float64}}; lim :: Q = Q(0.0)
    ) where{DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return f((p,); lim)
end