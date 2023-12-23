# call to MeshFunction, only value type and mesh point, no interpolation
function (f :: MeshFunction{MD, SD, DD, Q})(
    p :: NTuple{MD, Union{AbstractValue, AbstractMeshPoint}},
    x :: Vararg{Int64, SD} 
    ) :: Q where{MD, SD, DD, Q <: Number}

    return f[CartesianIndex_bc(f, p, x...)]
end

function (f :: MeshFunction{1, SD, DD, Q})(
    p :: Union{AbstractValue, AbstractMeshPoint},
    x :: Vararg{Int64, SD} 
    ) :: Q where{SD, DD, Q <: Number}

    return f((p,), x...)
end

function (f :: MeshFunction{MD, 0, DD, Q})(
    p :: Vararg{Union{AbstractValue, AbstractMeshPoint}, MD} 
    ) :: Q where{MD, DD, Q <: Number}

    return f((p...,))
end

# specialize 
function (f :: MeshFunction{1, 0, DD, Q})(
    p :: Union{AbstractValue, AbstractMeshPoint},
    ) :: Q where{DD, Q <: Number}

    return f((p,))
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

    for x in zip(wgt_iters, idx_iters)
        val += prod(first(x)) * f[last(x)...]
    end

    return val
end

function (f :: MeshFunction{1, SD, DD, Q})(
    p :: Union{AbstractValue, AbstractMeshPoint, Float64, SVector},
    x :: Vararg{Int64, SD} 
    ) :: Q where{SD, DD, Q <: Number}

    return f((p,), x...)
end

function (f :: MeshFunction{MD, 0, DD, Q})(
    p :: Vararg{Union{AbstractValue, AbstractMeshPoint, Float64, SVector}, MD} 
    ) :: Q where{MD, DD, Q <: Number}

    return f((p...,))
end

# specialize 
function (f :: MeshFunction{1, 0, DD, Q})(
    p :: Union{AbstractValue, AbstractMeshPoint, Float64, SVector},
    ) :: Q where{DD, Q <: Number}

    return f((p,))
end