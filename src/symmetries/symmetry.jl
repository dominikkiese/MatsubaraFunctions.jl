# type defs
#-------------------------------------------------------------------------------#

"""
    abstract type AbstractSymmetry

AbstractSymmetry type
"""
abstract type AbstractSymmetry end

"""
    struct Symmetry{DD} <: AbstractSymmetry

Symmetry type with fields:
* `f :: Function`

A Symmetry takes a tuple of `DD` MeshPoint value types and 
returns a tuple of `DD` MeshPoint value types and an `Operation`.
"""
struct Symmetry{DD} <: AbstractSymmetry
    f :: Function
end

"""
    struct InitFunction{DD, Q <: Number}

InitFunction type with fields:
* `f :: Function` 

An InitFunction takes a tuple of `DD` MeshPoint value types and 
returns a number of type `Q`.
"""
struct InitFunction{DD, Q <: Number}
    f :: Function
end

# call to Symmetry and InitFunction, force correct input and return type
#-------------------------------------------------------------------------------#

function (S :: Symmetry{DD})(w :: NTuple{DD, Union{<: AbstractValue}}
    ) :: Tuple{NTuple{DD, Union{<: AbstractValue}}, Operation} where {DD}

    return S.f(w)
end

function (I :: InitFunction{DD, Q})(w :: NTuple{DD, Union{<: AbstractValue}}
    ) :: Q where {DD, Q <: Number}
    
    return I.f(w)
end

# export
#-------------------------------------------------------------------------------#

export 
    AbstractSymmetry,
    Symmetry,
    InitFunction