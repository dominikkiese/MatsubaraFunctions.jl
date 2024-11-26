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