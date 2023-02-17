# basic operations
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

function Base.:-(
    w :: MatsubaraFrequency, 
    ) :: MatsubaraFrequency 

    T = temperature(w); n = index(w)
    return type(w) == :Fermion ? MatsubaraFrequency(T, -n - 1, Fermion) : MatsubaraFrequency(T, -n, Boson)
end