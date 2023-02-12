# define particle types for dispatch
abstract type AbstractParticle end 
struct Fermion <: AbstractParticle end 
struct Boson   <: AbstractParticle end

# Note: we do not allow MatsubaraFrequency to be dispatched on the particle type 
#       to allow for mixed grids in the construction of MatsubaraFunctions
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
function temperature(
    w :: MatsubaraFrequency
    ) :: Float64 

    return w.T 
end 

function value(
    w :: MatsubaraFrequency
    ) :: Float64 

    return w.val 
end 

function index(
    w :: MatsubaraFrequency
    ) :: Int64

    return w.idx
end 

function type(
    w :: MatsubaraFrequency
    ) :: Symbol

    return w.type
end 



# basic operations (dispatching onto the particle type would be faster, but see line 1)
function Base.:+(
    w1 :: MatsubaraFrequency, 
    w2 :: MatsubaraFrequency
    )  :: MatsubaraFrequency 

    @assert temperature(w1) ≈ temperature(w2) "Temperatures must be equal for addition"
    shift = (type(w1) == :Fermion && type(w2) == :Fermion)
    ptype = (type(w1) == type(w2) ? Boson : Fermion)
    return MatsubaraFrequency(temperature(w1), index(w1) + index(w2) + shift, ptype)
end

function Base.:-(
    w1 :: MatsubaraFrequency, 
    w2 :: MatsubaraFrequency
    )  :: MatsubaraFrequency 

    @assert temperature(w1) ≈ temperature(w2) "Temperatures must be equal for subtraction"
    shift = (type(w1) == :Boson && type(w2) == :Fermion)
    ptype = (type(w1) == type(w2) ? Boson : Fermion)
    return MatsubaraFrequency(temperature(w1), index(w1) - index(w2) - shift, ptype)
end