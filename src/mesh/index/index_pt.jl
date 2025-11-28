# ----------------------------------------------------------------------------- #
# Index Point Types
# ----------------------------------------------------------------------------- #

"""
    Index Point

Defines the `Index` type for representing an integer index as a mesh point.

Fields:
* `index` : Integer value

Examples:
```julia
# Construct an index point
idx = Index(5)

# Access value
val = value(idx)
```
"""

# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct Index <: AbstractValue

Index type with fields:
* `index :: Int` : Matsubara index
"""
struct Index <: AbstractValue
    index :: Int 
end 

"""
    function value(w :: Index) :: Int

Returns `w.index`
"""
function value(w :: Index) :: Int
    return w.index
end 

# arithmetic operations
#-------------------------------------------------------------------------------#

# import functions for overloading
import Base: +, -

# implement fallbacks, Index does not support arithmetics
-(w :: Index) = error("Arithmetic operations not supported for Index type")

for op in (:+, :-)
    @eval ($op)(w1 :: Index, w2 :: Index) = error("Arithmetic operations not supported for Index type")
end

# comparison operator
#-------------------------------------------------------------------------------#

function Base.:(==)(w1 :: Index, w2 :: Index) 
    return value(w1) == value(w2)
end 

# info
#-------------------------------------------------------------------------------#

function info(w :: Index)
    println(CYAN, BOLD, "Index", RESET)
    println("=> value : $(value(w))")
    return nothing 
end

# export
#-------------------------------------------------------------------------------#

export 
    Index,
    value,
    info