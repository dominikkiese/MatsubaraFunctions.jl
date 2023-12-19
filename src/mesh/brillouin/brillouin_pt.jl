# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct BrillouinPoint{N} <: AbstractValue
    
BrillouinPoint type with fields:
* `index :: SVector{N, Int64}` : coordinates in reciprocal space
"""
struct BrillouinPoint{N} <: AbstractValue where {N}
    index :: SVector{N, Int64}

    function BrillouinPoint(
        index :: SVector{N, Int64}
        )     :: BrillouinPoint{N} where {N}

        return new{N}(index)
    end 

    function BrillouinPoint(
        index :: Vararg{Int64, N}
        )     :: BrillouinPoint{N} where {N}

        return new{N}(SVector{N, Int64}(index...))
    end 
end

"""
    function index(
        k :: BrillouinPoint{N}
        ) :: SVector{N, Int64}

Returns `k.index`
""" 
function index(
    k :: BrillouinPoint{N}
    ) :: SVector{N, Int64} where {N}

    return k.index 
end

# arithmetic operations
#-------------------------------------------------------------------------------#

# addition 
function Base.:+(
    k1 :: BrillouinPoint{N},
    k2 :: BrillouinPoint{N}
    )  :: BrillouinPoint{N} where {N} 

    return BrillouinPoint(index(k1) .+ index(k2))
end

# subtraction 
function Base.:-(
    k1 :: BrillouinPoint{N},
    k2 :: BrillouinPoint{N}
    )  :: BrillouinPoint{N} where {N} 

    return BrillouinPoint(index(k1) .- index(k2))
end

# sign reversal 
function Base.:-(
    k :: BrillouinPoint{N}
    ) :: BrillouinPoint{N} where {N} 

    return BrillouinPoint(-index(k))
end

# comparison operator
#-------------------------------------------------------------------------------#

function Base.:(==)(
    k1 :: BrillouinPoint{N},
    k2 :: BrillouinPoint{N},
    )  :: Bool where {N} 

    return index(k1) == index(k2)
end 

# export
#-------------------------------------------------------------------------------#

export 
    BrillouinPoint,
    index