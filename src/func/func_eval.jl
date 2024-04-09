# call to MeshFunction, only value type and mesh point, no interpolation
function (f :: MeshFunction{DD, Q, MT, AT})(p :: Vararg{Union{MeshPoint, <: AbstractValue, Int}, DD}; lim :: Q = Q(0.0)
    ) where{DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}
    
    return all(ntuple(i -> is_inbounds_bc(p[i], meshes(f, i)), DD)) ? _eval_nointerpol(f, p...) : lim
end

@generated function _eval_nointerpol(f :: MeshFunction{DD, Q, MT, AT}, p :: Vararg{Union{MeshPoint, <: AbstractValue, Int}, DD}
    ) where{DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}
    
    # force compiler to dispatch properly for every input type
    return Meta.parse("f["*prod(["mesh_index_bc(p[$i], meshes(f, $i)), " for i in 1 : DD])[1 : end - 2]*"]")
end

# call to MeshFunction, all types, interpolation
function (f :: MeshFunction{DD, Q, MT, AT})(
    p :: Vararg{Union{MeshPoint, <: AbstractValue, Int, Float64, <: AbstractVector{Float64}}, DD}; lim :: Q = Q(0.0)
    ) where{DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    if all(ntuple(i -> is_inbounds_bc(p[i], meshes(f, i)), DD))
        params    = ntuple(i -> InterpolationParam(p[i], meshes(f, i)), DD)
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