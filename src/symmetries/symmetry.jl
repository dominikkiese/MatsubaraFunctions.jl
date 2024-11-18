# type def
#-------------------------------------------------------------------------------#

"""
    struct Symmetry{DD}

Symmetry type with fields:
* `f :: Function`

A symmetry is a function that takes a tuple of `DD` MeshPoint value types and 
returns a tuple of `DD` MeshPoint value types and an `Operation`.
"""
struct Symmetry{DD}
    f :: Function
end

# call to Symmetry
#-------------------------------------------------------------------------------#

# explicit return type to enforce proper implementation
function (S :: Symmetry{DD})(w :: NTuple{DD, Union{<: AbstractValue}}
    ) :: Tuple{NTuple{DD, Union{<: AbstractValue}}, Operation} where {DD}

    return S.f(w)
end

# export
#-------------------------------------------------------------------------------#

export 
    Symmetry