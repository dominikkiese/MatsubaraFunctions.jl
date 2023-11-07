"""
    struct MatsubaraFrequency{PT <: AbstractParticle} <: AbstractMatsubaraFrequency

MatsubaraFrequency type with fields:
* `T   :: Float64` : temperature
* `val :: Float64` : position on the imaginary axis
* `idx :: Int64`   : Matsubara index
"""
struct MatsubaraFrequency{PT <: AbstractParticle} <: AbstractMatsubaraFrequency
    T   :: Float64 
    val :: Float64
    idx :: Int64  

    function MatsubaraFrequency(
        T   :: Float64, 
        val :: Float64, 
        idx :: Int64,
            :: Type{PT}
        )   :: MatsubaraFrequency{PT} where {PT <: AbstractParticle}
        
        return new{PT}(T, val, idx)
    end 

    function MatsubaraFrequency(
        T   :: Float64, 
        idx :: Int64, 
            :: Type{Fermion}
        )   :: MatsubaraFrequency{Fermion}

        return new{Fermion}(T, pi * T * (2 * idx + 1), idx)
    end 
 
    function MatsubaraFrequency(
        T   :: Float64, 
        idx :: Int64, 
            :: Type{Boson}
        )   :: MatsubaraFrequency{Boson}

        return new{Boson}(T, 2 * pi * T * idx, idx)
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

#----------------------------------------------------------------------------------------------#

function Base.:+(
    w1 :: MatsubaraFrequency{Fermion}, 
    w2 :: MatsubaraFrequency{Fermion}
    )  :: MatsubaraFrequency{Boson}

    T = temperature(w1)
    @DEBUG temperature(w2) ≈ T "Temperatures must be equal for addition"
    return MatsubaraFrequency(T, index(w1) + index(w2) + 1, Boson)
end

function Base.:+(
    w1 :: MatsubaraFrequency{Boson}, 
    w2 :: MatsubaraFrequency{Boson}
    )  :: MatsubaraFrequency{Boson}

    T = temperature(w1)
    @DEBUG temperature(w2) ≈ T "Temperatures must be equal for addition"
    return MatsubaraFrequency(T, index(w1) + index(w2), Boson)
end

function Base.:+(
    w1 :: MatsubaraFrequency{Fermion}, 
    w2 :: MatsubaraFrequency{Boson}
    )  :: MatsubaraFrequency{Fermion}

    T = temperature(w1)
    @DEBUG temperature(w2) ≈ T "Temperatures must be equal for addition"
    return MatsubaraFrequency(T, index(w1) + index(w2), Fermion)
end

function Base.:+(
    w1 :: MatsubaraFrequency{Boson}, 
    w2 :: MatsubaraFrequency{Fermion}
    )  :: MatsubaraFrequency{Fermion}

    T = temperature(w1)
    @DEBUG temperature(w2) ≈ T "Temperatures must be equal for addition"
    return MatsubaraFrequency(T, index(w1) + index(w2), Fermion)
end

function Base.:-(
    w1 :: MatsubaraFrequency{Fermion}, 
    w2 :: MatsubaraFrequency{Fermion}
    )  :: MatsubaraFrequency{Boson}

    T = temperature(w1)
    @DEBUG temperature(w2) ≈ T "Temperatures must be equal for subtraction"
    return MatsubaraFrequency(T, index(w1) - index(w2), Boson)
end

function Base.:-(
    w1 :: MatsubaraFrequency{Boson}, 
    w2 :: MatsubaraFrequency{Boson}
    )  :: MatsubaraFrequency{Boson}

    T = temperature(w1)
    @DEBUG temperature(w2) ≈ T "Temperatures must be equal for subtraction"
    return MatsubaraFrequency(T, index(w1) - index(w2), Boson)
end

function Base.:-(
    w1 :: MatsubaraFrequency{Fermion}, 
    w2 :: MatsubaraFrequency{Boson}
    )  :: MatsubaraFrequency{Fermion}

    T = temperature(w1)
    @DEBUG temperature(w2) ≈ T "Temperatures must be equal for subtraction"
    return MatsubaraFrequency(T, index(w1) - index(w2), Fermion)
end

function Base.:-(
    w1 :: MatsubaraFrequency{Boson}, 
    w2 :: MatsubaraFrequency{Fermion}
    )  :: MatsubaraFrequency{Fermion}

    T = temperature(w1)
    @DEBUG temperature(w2) ≈ T "Temperatures must be equal for subtraction"
    return MatsubaraFrequency(T, index(w1) - index(w2) - 1, Fermion)
end

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

#----------------------------------------------------------------------------------------------#

"""
    struct MatsubaraIndex{PT <: AbstractParticle} <: AbstractMatsubaraFrequency

MatsubaraIndex type with fields:
* `idx :: Int64` : Matsubara index
"""
struct MatsubaraIndex{PT <: AbstractParticle} <: AbstractMatsubaraFrequency
    idx :: Int64
    
    function MatsubaraIndex(
        idx :: Int64,
            :: Type{PT}
        )   :: MatsubaraIndex{PT} where {PT <: AbstractParticle}

        return new{PT}(idx)
    end

    function MatsubaraIndex(
        w :: MatsubaraFrequency{PT}
        ) :: MatsubaraIndex{PT} where {PT <: AbstractParticle}

        return new{PT}(index(w))
    end
end

function MatsubaraFrequency(
    T :: Float64, 
    x :: MatsubaraIndex{PT}
    ) :: MatsubaraFrequency{PT} where {PT <: AbstractParticle}

    return MatsubaraFrequency(T, index(x), PT)
end

"""
    function index( 
        x :: MatsubaraIndex
        ) :: Int64 

Returns `x.idx`
"""
function index( 
    x :: MatsubaraIndex
    ) :: Int64 

    return x.idx 
end

#----------------------------------------------------------------------------------------------#

function Base.:+(
    w1 :: MatsubaraIndex{Fermion}, 
    w2 :: MatsubaraIndex{Fermion}
    )  :: MatsubaraIndex{Boson}

    return MatsubaraIndex(index(w1) + index(w2) + 1, Boson)
end

function Base.:+(
    w1 :: MatsubaraIndex{Boson}, 
    w2 :: MatsubaraIndex{Boson}
    )  :: MatsubaraIndex{Boson}

    return MatsubaraIndex(index(w1) + index(w2), Boson)
end

function Base.:+(
    w1 :: MatsubaraIndex{Fermion}, 
    w2 :: MatsubaraIndex{Boson}
    )  :: MatsubaraIndex{Fermion}

    return MatsubaraIndex(index(w1) + index(w2), Fermion)
end

function Base.:+(
    w1 :: MatsubaraIndex{Boson}, 
    w2 :: MatsubaraIndex{Fermion}
    )  :: MatsubaraIndex{Fermion}

    return MatsubaraIndex(index(w1) + index(w2), Fermion)
end

function Base.:-(
    w1 :: MatsubaraIndex{Fermion}, 
    w2 :: MatsubaraIndex{Fermion}
    )  :: MatsubaraIndex{Boson}

    return MatsubaraIndex(index(w1) - index(w2), Boson)
end

function Base.:-(
    w1 :: MatsubaraIndex{Boson}, 
    w2 :: MatsubaraIndex{Boson}
    )  :: MatsubaraIndex{Boson}

    return MatsubaraIndex(index(w1) - index(w2), Boson)
end

function Base.:-(
    w1 :: MatsubaraIndex{Fermion}, 
    w2 :: MatsubaraIndex{Boson}
    )  :: MatsubaraIndex{Fermion}

    return MatsubaraIndex(index(w1) - index(w2), Fermion)
end

function Base.:-(
    w1 :: MatsubaraIndex{Boson}, 
    w2 :: MatsubaraIndex{Fermion}
    )  :: MatsubaraIndex{Fermion}

    return MatsubaraIndex(index(w1) - index(w2) - 1, Fermion)
end

function Base.:-(
    w :: MatsubaraIndex{Fermion}, 
    ) :: MatsubaraIndex{Fermion}

    return MatsubaraIndex(-index(w) - 1, Fermion)
end

function Base.:-(
    w :: MatsubaraIndex{Boson}, 
    ) :: MatsubaraIndex{Boson}

    return MatsubaraIndex(-index(w), Boson)
end

#----------------------------------------------------------------------------------------------#

export 
    MatsubaraFrequency,
    temperature,
    value,
    index,
    MatsubaraIndex