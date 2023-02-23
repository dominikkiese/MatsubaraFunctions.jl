# define particle types for dispatch
"""
    abstract type AbstractParticle

AbstractParticle type
"""
abstract type AbstractParticle end 

"""
    struct Fermion <: AbstractParticle

Fermionic particle type, used for MatsubaraFrequency and MatsubaraGrid constructors
"""
struct Fermion <: AbstractParticle end 

"""
    struct Boson <: AbstractParticle

Bosonic particle type, used for MatsubaraFrequency and MatsubaraGrid constructors
"""
struct Boson <: AbstractParticle end

# Note: we do not allow MatsubaraFrequency to be dispatched on the particle type 
#       to allow for mixed grids in the construction of MatsubaraFunctions
"""
    struct MatsubaraFrequency

MatsubaraFrequency type with fields:
* `T    :: Float64` : physical temperature
* `val  :: Float64` : position on the imaginary axis
* `idx  :: Int64`   : Matsubara index
* `type :: Symbol`  : particle type

Examples:
```julia
# construction
T   = 1.0
idx = 5
v   = MatsubaraFrequency(T, idx, Fermion)
w   = MatsubaraFrequency(T, idx, Boson) 

# usage
w1 = v + v # type(v1) = :Boson
v2 = w - v # type(v2) = :Fermion
v3 = -v    # type(v3) = :Fermion 
```
"""
struct MatsubaraFrequency 
    T    :: Float64 
    val  :: Float64
    idx  :: Int64 
    type :: Symbol 

    # basic constructor
    function MatsubaraFrequency(
        T    :: Float64, 
        val  :: Float64, 
        idx  :: Int64, 
        type :: Symbol 
        )    :: MatsubaraFrequency
        
        @assert type == :Fermion || type == :Boson "MatsubaraFrequency type must be Fermion or Boson"
        return new(T, val, idx, type)
    end 

    # convenience constructor for fermionic frequencies 
    function MatsubaraFrequency(
        T   :: Float64, 
        idx :: Int64, 
            :: Type{Fermion}
        )   :: MatsubaraFrequency 

        return new(T, pi * T * (2 * idx + 1), idx, :Fermion)
    end 

    # convenience constructor for bosonic frequencies 
    function MatsubaraFrequency(
        T   :: Float64, 
        idx :: Int64, 
            :: Type{Boson}
        )   :: MatsubaraFrequency 

        return new(T, 2.0 * pi * T * idx, idx, :Boson)
    end 
end 



# getter functions
"""
    function temperature(
        w :: MatsubaraFrequency
        ) :: Float64

Returns w.T
"""
function temperature(
    w :: MatsubaraFrequency
    ) :: Float64 

    return w.T 
end 

"""
    function value(
        w :: MatsubaraFrequency
        ) :: Float64

Returns w.val
"""
function value(
    w :: MatsubaraFrequency
    ) :: Float64 

    return w.val 
end 

"""
    function index(
        w :: MatsubaraFrequency
        ) :: Int64

Returns w.idx
"""
function index(
    w :: MatsubaraFrequency
    ) :: Int64

    return w.idx
end 

"""
    function type(
        w :: MatsubaraFrequency
        ) :: Symbol

Returns w.type
"""
function type(
    w :: MatsubaraFrequency
    ) :: Symbol

    return w.type
end 



# load methods
include("matsubara_freq_ops.jl")