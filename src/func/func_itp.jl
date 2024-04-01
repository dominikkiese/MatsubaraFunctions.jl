# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct InterpolationParam{N}

InterpolationParam struct with fields:
* `indices :: NTuple{N, Int}`     : linear indices for reference data
* `weights :: NTuple{N, Float64}` : interpolation weights
"""
struct InterpolationParam{N}
    indices :: NTuple{N, Int}
    weights :: NTuple{N, Float64}

    function InterpolationParam(indices :: NTuple{N, Int}, weights :: NTuple{N, Float64}) where {N}
        @DEBUG sum(weights) ≈ 1.0 "Weights must add up to 1"
        return new{N}(indices, weights)
    end 

    function InterpolationParam(index :: Int, weight :: Float64) 
        return InterpolationParam((index,), (weight,))
    end 

    # from mesh point or value type 
    function InterpolationParam(x :: Union{MeshPoint, <: AbstractValue}, m :: Mesh)
        return InterpolationParam(mesh_index_bc(x, m), 1.0)
    end

    # Matsubara mesh
    function InterpolationParam(w :: Float64, m :: Mesh{MeshPoint{MatsubaraFrequency{PT}}, MatsubaraDomain}) where {PT <: AbstractParticle}
        # calculate mesh spacing and position in mesh
        w     = max(first_value(m), min(w, last_value(m)))
        delta = plain_value(m[2]) - plain_value(m[1])
        x     = (w - plain_value(m[1])) / delta

        # calculate mesh indices and interpolation weights
        indices = floor(Int, x) + 1, min(ceil(Int, x) + 1, length(m))

        if first(indices) < last(indices)
            weights = (plain_value(m[last(indices)]) - w) / delta, (w - plain_value(m[first(indices)])) / delta
            return InterpolationParam(indices, weights)
        else 
            return InterpolationParam(first(indices), 1.0)
        end
    end

    # Brillouin zone mesh
    function InterpolationParam(k :: T, m :: Mesh{MeshPoint{BrillouinPoint{N}}, BrillouinDomain{N}}) where {N, T <: AbstractVector{Float64}}
        @DEBUG length(k) == N "Length mismatch for input vector"

        # calculate mesh spacing and position in mesh
        x      = reciprocal(k, m)
        ranges = ntuple(n -> floor(Int, x[n]) : ceil(Int, x[n]), N)
        iters  = collect(Iterators.product(ranges...))

        # calculate mesh indices and interpolation weights
        wgt     = y -> prod(ntuple(n -> 1.0 - abs(y[n] - x[n]), N))
        idx     = y -> mesh_index_bc(BrillouinPoint(y...), m)
        weights = ntuple(n -> wgt(iters[n]), length(iters))
        indices = ntuple(n -> idx(iters[n]), length(iters))

        return InterpolationParam(indices, weights)
    end

    # Index mesh 
    function InterpolationParam(x :: Int, m :: Mesh{MeshPoint{Index}, IndexDomain})
        return InterpolationParam(x, 1.0)
    end
end

"""
    function indices(p :: InterpolationParam{N}) :: NTuple{N, Int} where {N} 

Returns `p.indices`
"""
function indices(p :: InterpolationParam{N}) :: NTuple{N, Int} where {N} 
    return p.indices 
end 

"""
    function weights(p :: InterpolationParam{N}) :: NTuple{N, Float64} where {N} 

Returns `p.weights`
"""
function weights(p :: InterpolationParam{N}) :: NTuple{N, Float64} where {N}  
    return p.weights
end 

# export
#-------------------------------------------------------------------------------#

export 
    InterpolationParam,
    indices,
    weights