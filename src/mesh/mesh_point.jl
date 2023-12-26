# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct MeshPoint{T <: AbstractValue} <: AbstractMeshPoint

MeshPoint type with fields:
* `hash  :: UInt64` : mesh identifier
* `index :: Int64`  : mesh index
* `value :: T`      : mesh coordinates
"""
struct MeshPoint{T <: AbstractValue} <: AbstractMeshPoint
    hash  :: UInt64 # no accessor, only for internal use
    index :: Int64 
    value :: T 
end

"""
    function index(
        x :: AbstractMeshPoint
        ) :: Int64

Returns `x.index`
"""
function index(
    x :: AbstractMeshPoint
    ) :: Int64

    return x.index
end 

"""
    function value(
        x :: MeshPoint{T}
        ) :: T where {T <: AbstractValue}  

Returns `x.value`
"""
function value(
    x :: MeshPoint{T}
    ) :: T where {T <: AbstractValue} 

    return x.value
end 


"""
    function plain_value(
        x :: MeshPoint{T}
        ) where {T <: AbstractValue}  

Returns the plain underlying value of the mesh point `value(x.value)`.
"""
function plain_value(
    x :: MeshPoint{T}
    ) where {T <: AbstractValue} 

    return value(x.value)
end 

"""
    function plain_value(
        x :: AbstractValue
        )

Returns the plain underlying value of the mesh point `value(x.value)`.
"""
function plain_value(
    x :: AbstractValue
    )

    return value(x)
end 

# arithmetic operations
#-------------------------------------------------------------------------------#

# mesh point operations are operations on their values
# each value type must implement +, - and sign reversal

function Base.:+(
    x1 :: MeshPoint{T1},
    x2 :: MeshPoint{T2}
    ) where {T1 <: AbstractValue, T2 <: AbstractValue} 
    
    return value(x1) + value(x2)
end

function Base.:-(
    x1 :: MeshPoint{T1},
    x2 :: MeshPoint{T2}
    ) where {T1 <: AbstractValue, T2 <: AbstractValue} 
    
    return value(x1) - value(x2)
end

function Base.:-(
    x :: MeshPoint{T1},
    ) :: T1 where {T1 <: AbstractValue} 
    
    return -value(x)
end


function Base.:+(
    x1 :: T1,
    x2 :: MeshPoint{T2}
    ) where {T1 <: AbstractValue, T2 <: AbstractValue} 
    
    return x1 + value(x2)
end

function Base.:-(
    x1 :: T1,
    x2 :: MeshPoint{T2}
    ) where {T1 <: AbstractValue, T2 <: AbstractValue} 
    
    return x1 - value(x2)
end

function Base.:+(
    x1 :: MeshPoint{T1},
    x2 :: T2
    ) where {T1 <: AbstractValue, T2 <: AbstractValue} 
    
    return value(x1) + x2
end

function Base.:-(
    x1 :: MeshPoint{T1},
    x2 :: T2
    ) where {T1 <: AbstractValue, T2 <: AbstractValue} 
    
    return value(x1) - x2
end


# comparison operator
#-------------------------------------------------------------------------------#

function Base.:(==)(
    x1 :: MeshPoint{T},
    x2 :: MeshPoint{T}
    )  :: Bool where {T <: AbstractValue} 

    if (x1.hash != x2.hash) || (index(x1) != index(x2)) || (value(x1) != value(x2))
        return false
    end

    return true
end

# export
#-------------------------------------------------------------------------------#

export 
    MeshPoint,
    index,
    value,
    plain_value