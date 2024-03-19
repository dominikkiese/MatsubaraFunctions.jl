# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct Index <: AbstractValue

Index type with fields:
* `index :: Int64` : Matsubara index
"""
struct Index <: AbstractValue
    index :: Int64 
end 

"""
    function value(w :: Index) :: Int64

Returns `w.index`
"""
function value(w :: Index) :: Int64
    return w.index
end 

# export
#-------------------------------------------------------------------------------#

export 
    Index,
    value