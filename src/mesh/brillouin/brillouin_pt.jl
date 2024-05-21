# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct BrillouinPoint{N} <: AbstractValue
    
BrillouinPoint type with fields:
* `value :: SVector{N, Int}` : coordinates in reciprocal space
"""
struct BrillouinPoint{N} <: AbstractValue where {N}
    value :: SVector{N, Int}

    function BrillouinPoint(value :: SVector{N, Int}) where {N}
        return new{N}(value)
    end 

    function BrillouinPoint(value :: Vararg{Int, N}) where {N}
        return new{N}(SVector{N, Int}(value...))
    end 
end

"""
    function value(k :: BrillouinPoint{N}) :: SVector{N, Int}

Returns `k.value`
""" 
function value(k :: BrillouinPoint{N}) :: SVector{N, Int} where {N}
    return k.value 
end

# getindex
#-------------------------------------------------------------------------------#

function Base.:getindex(k :: BrillouinPoint{N}, i :: Int) where {N}
    return value(k)[i]
end

# arithmetic operations
#-------------------------------------------------------------------------------#

# addition 
function Base.:+(k1 :: BrillouinPoint{N}, k2 :: BrillouinPoint{N}) where {N} 
    return BrillouinPoint(value(k1) .+ value(k2))
end

# subtraction 
function Base.:-(k1 :: BrillouinPoint{N}, k2 :: BrillouinPoint{N}) where {N} 
    return BrillouinPoint(value(k1) .- value(k2))
end

# sign reversal 
function Base.:-(k :: BrillouinPoint{N}) where {N} 
    return BrillouinPoint(-value(k))
end

# comparison operator
#-------------------------------------------------------------------------------#

function Base.:(==)(k1 :: BrillouinPoint{N}, k2 :: BrillouinPoint{N}) where {N} 
    return value(k1) == value(k2)
end 

# print 
#-------------------------------------------------------------------------------#

function Base.:show(io :: IO, k :: BrillouinPoint{N}) where {N}
    print(io, CYAN, BOLD, "BrillouinPoint ", RESET, "of dimension ", CYAN, BOLD, "$(N) \n", RESET)
    print(io, "=> value : $(value(k))")
    return nothing 
end

# export
#-------------------------------------------------------------------------------#

export 
    BrillouinPoint,
    value