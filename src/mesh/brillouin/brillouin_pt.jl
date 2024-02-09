# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct BrillouinPoint{N} <: AbstractValue
    
BrillouinPoint type with fields:
* `value :: SVector{N, Int64}` : coordinates in reciprocal space
"""
struct BrillouinPoint{N} <: AbstractValue where {N}
    value :: SVector{N, Int64}

    function BrillouinPoint(
        value :: SVector{N, Int64}
        )     :: BrillouinPoint{N} where {N}

        return new{N}(value)
    end 

    function BrillouinPoint(
        value :: Vararg{Int64, N}
        )     :: BrillouinPoint{N} where {N}

        return new{N}(SVector{N, Int64}(value...))
    end 
end

"""
    function value(
        k :: BrillouinPoint{N}
        ) :: SVector{N, Int64}

Returns `k.value`
""" 
function value(
    k :: BrillouinPoint{N}
    ) :: SVector{N, Int64} where {N}

    return k.value 
end

# arithmetic operations
#-------------------------------------------------------------------------------#

# addition 
function Base.:+(
    k1 :: BrillouinPoint{N},
    k2 :: BrillouinPoint{N}
    )  :: BrillouinPoint{N} where {N} 

    return BrillouinPoint(value(k1) .+ value(k2))
end

# subtraction 
function Base.:-(
    k1 :: BrillouinPoint{N},
    k2 :: BrillouinPoint{N}
    )  :: BrillouinPoint{N} where {N} 

    return BrillouinPoint(value(k1) .- value(k2))
end

# sign reversal 
function Base.:-(
    k :: BrillouinPoint{N}
    ) :: BrillouinPoint{N} where {N} 

    return BrillouinPoint(-value(k))
end

# comparison operator
#-------------------------------------------------------------------------------#

function Base.:(==)(
    k1 :: BrillouinPoint{N},
    k2 :: BrillouinPoint{N},
    )  :: Bool where {N} 

    return value(k1) == value(k2)
end 

# export
#-------------------------------------------------------------------------------#

export 
    BrillouinPoint,
    value