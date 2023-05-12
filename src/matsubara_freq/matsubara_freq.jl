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

# we do not allow MatsubaraFrequency to be dispatched on the particle type 
# to allow for mixed grids in the construction of MatsubaraFunctions
"""
    struct MatsubaraFrequency

MatsubaraFrequency type with fields:
* `T    :: Float64` : temperature
* `val  :: Float64` : position on the imaginary axis
* `idx  :: Int64`   : Matsubara index
* `type :: Symbol`  : particle type
"""
struct MatsubaraFrequency 
    T    :: Float64 
    val  :: Float64
    idx  :: Int64 
    type :: Symbol 

    # default constructor
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



"""
    function temperature(
        w :: MatsubaraFrequency
        ) :: Float64

Returns `w.T`
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

Returns `w.val`
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

Returns `w.idx`
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

Returns `w.type`
"""
function type(
    w :: MatsubaraFrequency
    ) :: Symbol

    return w.type
end 



# load methods
include("matsubara_freq_ops.jl")