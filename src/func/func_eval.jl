# call to MeshFunction, only value type and mesh point, no interpolation
function (f :: MeshFunction{DD, Q, MT, AT})(p :: Vararg{Union{MeshPoint, <: AbstractValue, Int, UnitRange, Colon}, DD}; lim :: Q = Q(0.0)
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}
    
    return _all_inbounds_bc(f, p...) ? f[_mesh_indices_bc(f, p...)...] : lim
end

# call to MeshFunction, all types, interpolation
function (f :: MeshFunction{DD, Q, MT, AT})(
    p :: Vararg{Union{MeshPoint, <: AbstractValue, Int, Float64, <: AbstractVector{Float64}}, DD}; lim :: Q = Q(0.0)
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    if _all_inbounds_bc(f, p...)
        params    = _get_params(f, p...)
        idx_iters = Iterators.product(indices.(params)...)
        wgt_iters = Iterators.product(weights.(params)...)
        val :: Q  = 0.0

        for y in zip(wgt_iters, idx_iters)
            val += prod(first(y)) * f[last(y)...]
        end

        return val 
    else 
        return lim 
    end
end