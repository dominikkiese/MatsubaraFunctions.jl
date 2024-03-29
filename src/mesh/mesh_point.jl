# abstract types
#-------------------------------------------------------------------------------#

"""
    abstract type AbstractMeshPoint

AbstractMeshPoint type
"""
abstract type AbstractMeshPoint end

"""
    abstract type AbstractValue

AbstractValue type
"""
abstract type AbstractValue end

# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct MeshPoint{T <: AbstractValue} <: AbstractMeshPoint

MeshPoint type with fields:
* `hash  :: Symbol` : mesh identifier
* `index :: Int`    : mesh index
* `value :: T`      : mesh coordinates
"""
struct MeshPoint{T <: AbstractValue} <: AbstractMeshPoint
    hash  :: Symbol # no accessor, only for internal use
    index :: Int 
    value :: T

    function MeshPoint(hash, index :: Int, value :: T) where {T <: AbstractValue}
        return new{T}(Symbol(hash), index, value)
    end
end

"""
    function index(x :: T) :: Int where {T <: AbstractMeshPoint}

Returns `x.index`
"""
function index(x :: T) :: Int where {T <: AbstractMeshPoint}
    return x.index
end 

"""
    function value(x :: MeshPoint{T}) :: T where {T <: AbstractValue}  

Returns `x.value`
"""
function value(x :: MeshPoint{T}) :: T where {T <: AbstractValue} 
    return x.value
end 

"""
    function plain_value(x :: MeshPoint{T}) where {T <: AbstractValue} 

Returns the plain underlying value of the mesh point `value(value(x))`.
Assumes that a function `value` is defined for the value type.
"""
function plain_value(x :: MeshPoint{T}) where {T <: AbstractValue} 
    return value(value(x))
end 

# arithmetic operations
#-------------------------------------------------------------------------------#

# mesh point operations are operations on their values
# each value type must implement +, - and sign reversal

# addition
function Base.:+(x1 :: MeshPoint{T1}, x2 :: MeshPoint{T2}) where {T1 <: AbstractValue, T2 <: AbstractValue} 
    return value(x1) + value(x2)
end

function Base.:+(x1 :: T1, x2 :: MeshPoint{T2}) where {T1 <: AbstractValue, T2 <: AbstractValue} 
    return x1 + value(x2)
end

function Base.:+(x1 :: MeshPoint{T1}, x2 :: T2) where {T1 <: AbstractValue, T2 <: AbstractValue} 
    return value(x1) + x2
end

# subtraction
function Base.:-(x1 :: MeshPoint{T1}, x2 :: MeshPoint{T2}) where {T1 <: AbstractValue, T2 <: AbstractValue} 
    return value(x1) - value(x2)
end

function Base.:-(x1 :: T1, x2 :: MeshPoint{T2}) where {T1 <: AbstractValue, T2 <: AbstractValue} 
    return x1 - value(x2)
end

function Base.:-(x1 :: MeshPoint{T1}, x2 :: T2) where {T1 <: AbstractValue, T2 <: AbstractValue} 
    return value(x1) - x2
end

# sign reversal
function Base.:-(x :: MeshPoint{T}) where {T <: AbstractValue} 
    return -value(x)
end

# comparison operator
#-------------------------------------------------------------------------------#

function Base.:(==)(x1 :: MeshPoint{T}, x2 :: MeshPoint{T}) where {T <: AbstractValue} 
    return (x1.hash === x2.hash) && (index(x1) == index(x2)) && (value(x1) == value(x2))
end

# export
#-------------------------------------------------------------------------------#

export 
    AbstractMeshPoint,
    AbstractValue,
    MeshPoint,
    index,
    value,
    plain_value