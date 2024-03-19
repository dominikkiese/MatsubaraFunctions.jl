# call to MeshFunction, only value type and mesh point, no interpolation
function (f :: MeshFunction{DD, Q, AT})(p :: Vararg{Union{<: AbstractValue, <: AbstractMeshPoint}, DD}; lim :: Q = Q(0.0)
    ) where{DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return all(i -> is_inbounds_bc(p[i], meshes(f, i)), 1 : DD) ? f[map((y, m) -> mesh_index_bc(y, m), p, meshes(f))...] : lim
end

# call to MeshFunction, all types, interpolation
function (f :: MeshFunction{DD, Q, AT})(
    p :: Vararg{Union{<: AbstractValue, <: AbstractMeshPoint, Int, Float64, <: AbstractVector{Float64}}, DD}; lim :: Q = Q(0.0)
    ) where{DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    if all(i -> is_inbounds_bc(p[i], meshes(f, i)), 1 : DD)
        params    = map((y, m) -> InterpolationParam(y, m), p, meshes(f))
        idx_iters = Iterators.product(indices.(params)...)
        wgt_iters = Iterators.product(weights.(params)...)
        val       = 0.0

        for y in zip(wgt_iters, idx_iters)
            val += prod(first(y)) * f[last(y)...]
        end

        return val 
    else 
        return lim 
    end
end