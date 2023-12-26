# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct InterpolationParam{N}

InterpolationParam struct with fields:
* `indices :: NTuple{N, Int64}`   : linear indices for reference data
* `weights :: NTuple{N, Float64}` : interpolation weights
"""
struct InterpolationParam{N}
    indices :: NTuple{N, Int64}
    weights :: NTuple{N, Float64}

    function InterpolationParam(
        indices :: NTuple{N, Int64}, 
        weights :: NTuple{N, Float64}
        )       :: InterpolationParam{N} where {N}

        @DEBUG sum(weights) ≈ 1.0 "Weights must add up to 1"
        return new{N}(indices, weights)
    end 

    function InterpolationParam(
        index  :: Int64, 
        weight :: Float64
        )      :: InterpolationParam{1}

        return InterpolationParam((index,), (weight,))
    end 

    # from mesh point or value type 
    function InterpolationParam(
        p :: Union{MeshPoint{T}, T},
        m :: Mesh{MeshPoint{T}}
        ) :: InterpolationParam{1} where {T <: AbstractValue}

        return InterpolationParam(mesh_index_bc(p, m), 1.0)
    end

    # Matsubara mesh
    function InterpolationParam(
        w :: Float64, 
        m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}}
        ) :: InterpolationParam where {PT <: AbstractParticle}

        # calculate mesh spacing and position in mesh
        val   = x -> value(value(x))
        #w     = max(first_value(m), min(w, last_value(m)))
        #delta = val(m[firstindex(m)+1]) - first_value(m)
        #x     = (w - first_value(m)) / delta
        #println("delta: ", delta)
        #println("x: ", x)
        ## calculate mesh indices and interpolation weights
        #indices = floor(Int64, x)   , min(ceil(Int64, x), lastindex(m))
        #println("indices: ", indices)
        #if first(indices) < last(indices)
        #    weights = (val(m[last(indices)]) - w) / delta, (w - val(m[first(indices)])) / delta
        #    return InterpolationParam(indices, weights)
        #else 
        #    return InterpolationParam(first(indices), 1.0)
        #end

        delta    = val(m[2]) - val(m[1])
        position = (w - val(m[0])) / delta

        # nearest-neighbor indices
        # add min to upper index to improve robustness with respect to rounding errors
        dn_idx, up_idx = floor(Int64, position), min(ceil(Int64, position), lastindex(m))

        # interpolation weights
        if dn_idx < up_idx
            return InterpolationParam((dn_idx, up_idx), ((val(m[up_idx]) - w) / delta, (w - val(m[dn_idx])) / delta))
        else
            return InterpolationParam((up_idx, up_idx), (1.0, 0.0))
        end
    end

    # Brillouin zone mesh
    function InterpolationParam(
        k :: SVector{N, Float64},
        m :: Mesh{MeshPoint{BrillouinPoint{N}}}
        ) :: InterpolationParam where {N}

        # calculate mesh spacing and position in mesh
        x      = reciprocal(k, m)
        ranges = ntuple(n -> floor(Int64, x[n]) : ceil(Int64, x[n]), N)
        iters  = collect(Iterators.product(ranges...))

        # calculate mesh indices and interpolation weights
        wgt     = y -> prod(ntuple(n -> 1.0 - abs(y[n] - x[n]), N))
        idx     = y -> mesh_index_bc(BrillouinPoint(y...), m)
        weights = ntuple(n -> wgt(iters[n]), length(iters))
        indices = ntuple(n -> idx(iters[n]), length(iters))

        return InterpolationParam(indices, weights)
    end
end

"""
    function indices(
        p :: InterpolationParam{N}
        ) :: NTuple{N, Int64} where {N} 

Returns `p.indices`
"""
function indices(
    p :: InterpolationParam{N}
    ) :: NTuple{N, Int64} where {N} 

    return p.indices 
end 

"""
    function weights(
        p :: InterpolationParam{N}
        ) :: NTuple{N, Float64} where {N} 

Returns `p.weights`
"""
function weights(
    p :: InterpolationParam{N}
    ) :: NTuple{N, Float64} where {N} 

    return p.weights
end 

# export
#-------------------------------------------------------------------------------#

export 
    InterpolationParam,
    indices,
    weights