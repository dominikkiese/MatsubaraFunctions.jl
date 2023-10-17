# abstract types 
#-------------------------------------------------------------------------------#

"""
    abstract type AbstractParticle

AbstractParticle type
"""
abstract type AbstractParticle end 

"""
    struct Fermion <: AbstractParticle

Fermionic particle type
"""
struct Fermion <: AbstractParticle end 

"""
    struct Boson <: AbstractParticle

Bosonic particle type
"""
struct Boson <: AbstractParticle end

# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct MatsubaraFrequency{PT <: AbstractParticle} <: AbstractValue

MatsubaraFrequency type with fields:
* `temperature :: Float64` : temperature
* `value       :: Float64` : position on the imaginary axis
* `index       :: Int64`   : Matsubara index
"""
struct MatsubaraFrequency{PT <: AbstractParticle} <: AbstractValue
    temperature :: Float64 
    value       :: Float64
    index       :: Int64 

    # constructor for fermionic frequencies 
    function MatsubaraFrequency(
        temperature :: Float64, 
        index       :: Int64, 
                    :: Type{Fermion}
        )           :: MatsubaraFrequency 

        return new{Fermion}(temperature, pi * temperature * (2 * index + 1), index)
    end 

    # constructor for bosonic frequencies 
    function MatsubaraFrequency(
        temperature :: Float64, 
        index       :: Int64, 
                    :: Type{Boson}
        )           :: MatsubaraFrequency 

        return new{Boson}(temperature, 2 * pi * temperature * index, index)
    end 
end 

"""
    function temperature(
        w :: MatsubaraFrequency{PT}
        ) :: Float64 where {PT <: AbstractParticle}

Returns `w.temperature`
"""
function temperature(
    w :: MatsubaraFrequency{PT}
    ) :: Float64 where {PT <: AbstractParticle}

    return w.temperature
end 

"""
    function value(
        w :: MatsubaraFrequency{PT}
        ) :: Float64 where {PT <: AbstractParticle}

Returns `w.value`
"""
function value(
    w :: MatsubaraFrequency{PT}
    ) :: Float64 where {PT <: AbstractParticle}

    return w.value
end 

"""
    function index(
        w :: MatsubaraFrequency{PT}
        ) :: Int64 where {PT <: AbstractParticle}

Returns `w.index`
"""
function index(
    w :: MatsubaraFrequency{PT}
    ) :: Int64 where {PT <: AbstractParticle}

    return w.index
end 

# arithmetic operations
#-------------------------------------------------------------------------------#

# addition
function Base.:+(
    w1 :: MatsubaraFrequency{Fermion}, 
    w2 :: MatsubaraFrequency{Fermion}
    )  :: MatsubaraFrequency{Boson}

    T = temperature(w1)
    @check temperature(w2) ≈ T "Temperatures must be equal for addition"
    return MatsubaraFrequency(T, index(w1) + index(w2) + 1, Boson)
end

function Base.:+(
    w1 :: MatsubaraFrequency{Boson}, 
    w2 :: MatsubaraFrequency{Boson}
    )  :: MatsubaraFrequency{Boson}

    T = temperature(w1)
    @check temperature(w2) ≈ T "Temperatures must be equal for addition"
    return MatsubaraFrequency(T, index(w1) + index(w2), Boson)
end

function Base.:+(
    w1 :: MatsubaraFrequency{Fermion}, 
    w2 :: MatsubaraFrequency{Boson}
    )  :: MatsubaraFrequency{Fermion}

    T = temperature(w1)
    @check temperature(w2) ≈ T "Temperatures must be equal for addition"
    return MatsubaraFrequency(T, index(w1) + index(w2), Fermion)
end

function Base.:+(
    w1 :: MatsubaraFrequency{Boson}, 
    w2 :: MatsubaraFrequency{Fermion}
    )  :: MatsubaraFrequency{Fermion}

    T = temperature(w1)
    @check temperature(w2) ≈ T "Temperatures must be equal for addition"
    return MatsubaraFrequency(T, index(w1) + index(w2), Fermion)
end

# subtraction
function Base.:-(
    w1 :: MatsubaraFrequency{Fermion}, 
    w2 :: MatsubaraFrequency{Fermion}
    )  :: MatsubaraFrequency{Boson}

    T = temperature(w1)
    @check temperature(w2) ≈ T "Temperatures must be equal for subtraction"
    return MatsubaraFrequency(T, index(w1) - index(w2), Boson)
end

function Base.:-(
    w1 :: MatsubaraFrequency{Boson}, 
    w2 :: MatsubaraFrequency{Boson}
    )  :: MatsubaraFrequency{Boson}

    T = temperature(w1)
    @check temperature(w2) ≈ T "Temperatures must be equal for subtraction"
    return MatsubaraFrequency(T, index(w1) - index(w2), Boson)
end

function Base.:-(
    w1 :: MatsubaraFrequency{Fermion}, 
    w2 :: MatsubaraFrequency{Boson}
    )  :: MatsubaraFrequency{Fermion}

    T = temperature(w1)
    @check temperature(w2) ≈ T "Temperatures must be equal for subtraction"
    return MatsubaraFrequency(T, index(w1) - index(w2), Fermion)
end

function Base.:-(
    w1 :: MatsubaraFrequency{Boson}, 
    w2 :: MatsubaraFrequency{Fermion}
    )  :: MatsubaraFrequency{Fermion}

    T = temperature(w1)
    @check temperature(w2) ≈ T "Temperatures must be equal for subtraction"
    return MatsubaraFrequency(T, index(w1) - index(w2) - 1, Fermion)
end

# sign reversal
function Base.:-(
    w :: MatsubaraFrequency{Fermion}, 
    ) :: MatsubaraFrequency{Fermion}

    return MatsubaraFrequency(temperature(w), -index(w) - 1, Fermion)
end

function Base.:-(
    w :: MatsubaraFrequency{Boson}, 
    ) :: MatsubaraFrequency{Boson}

    return MatsubaraFrequency(temperature(w), -index(w), Boson)
end

# export
#-------------------------------------------------------------------------------#

export 
    AbstractParticle,
    Fermion,
    Boson,
    MatsubaraFrequency,
    temperature,
    value, 
    index