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

# print 
#-------------------------------------------------------------------------------#

function Base.:show(io :: IO, w :: Index)
    print(io, CYAN, BOLD, "Index \n", RESET)
    print(io, "=> value : $(value(w))")
    return nothing 
end

# export
#-------------------------------------------------------------------------------#

export 
    Index,
    value