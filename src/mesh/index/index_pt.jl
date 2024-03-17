# abstract types 
#-------------------------------------------------------------------------------#


# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct Index <: AbstractValue

Index type with fields:
* `index       :: Int64`   : Matsubara index
"""
struct Index <: AbstractValue
    index       :: Int64 

    # constructor for index
    function Index(index       :: Int64)
        return new(index)
    end 
end 

"""
    function value(w :: Index) :: Int

Returns `w.index`
"""
function value(w :: Index) :: Int64
    return index(w)
end 


"""
    function index(w :: Index) :: Int

Returns `w.index`
"""
function index(w :: Index) :: Int64
    return w.index
end 


# export
#-------------------------------------------------------------------------------#

export 
    Index,
    value