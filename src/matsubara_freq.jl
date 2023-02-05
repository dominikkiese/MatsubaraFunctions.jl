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

    @assert isapprox(temperature(w1), temperature(w2)) "Temperatures must be equal for addition"

    if type(w1) == :Fermion && type(w2) == :Fermion 
        return MatsubaraFrequency(temperature(w1), index(w1) + index(w2) + 1, Boson)
    elseif type(w1) == :Fermion && type(w2) == :Boson
        return MatsubaraFrequency(temperature(w1), index(w1) + index(w2), Fermion)
    elseif type(w1) == :Boson && type(w2) == :Fermion 
        return MatsubaraFrequency(temperature(w1), index(w1) + index(w2), Fermion)
    elseif type(w1) == :Boson && type(w2) == :Boson
        return MatsubaraFrequency(temperature(w1), index(w1) + index(w2), Boson)
    end 
end

function Base.:-(
    w1 :: MatsubaraFrequency, 
    w2 :: MatsubaraFrequency
    )  :: MatsubaraFrequency 

    @assert isapprox(temperature(w1), temperature(w2)) "Temperatures must be equal for subtraction"

    if type(w1) == :Fermion && type(w2) == :Fermion 
        return MatsubaraFrequency(temperature(w1), index(w1) - index(w2), Boson)
    elseif type(w1) == :Fermion && type(w2) == :Boson
        return MatsubaraFrequency(temperature(w1), index(w1) - index(w2), Fermion)
    elseif type(w1) == :Boson && type(w2) == :Fermion 
        return MatsubaraFrequency(temperature(w1), index(w1) - index(w2) - 1, Fermion)
    elseif type(w1) == :Boson && type(w2) == :Boson
        return MatsubaraFrequency(temperature(w1), index(w1) - index(w2), Boson)
    end 
end