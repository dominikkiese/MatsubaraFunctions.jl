#== 
    Basic operations for MatsubaraFunctions:
        -> addition 
        -> subtraction 
        -> multiplication with scalar 
        -> initialization with scalar or MatsubaraFunction
        -> value comparison
        -> flattening of data array into vector
==#

function check_shape_grid!(
    f1 :: MatsubaraFunction{GD, SD, DD, Q}, 
    f2 :: MatsubaraFunction{GD, SD, DD, Q}
    )  :: Nothing where {GD, SD, DD, Q <: Number}

    @check data_shape(f1) == data_shape(f2) "Data shapes must be equal"

    for i in 1 : GD 
        @check type(f1.grids[i]) == type(f2.grids[i]) "Grids must have same particle type" 
        @check temperature(f1.grids[i]) == temperature(f2.grids[i]) "Grids must have same temperature" 
        @check index_range(f1.grids[i]) == index_range(f2.grids[i]) "Grids must have same index range" 
    end

    return nothing 
end

"""
    function add(
        f1 :: MatsubaraFunction{GD, SD, DD, Q}, 
        f2 :: MatsubaraFunction{GD, SD, DD, Q}
        )  :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

Addition of two MatsubaraFunction, returns new MatsubaraFunction
"""
function add(
    f1 :: MatsubaraFunction{GD, SD, DD, Q}, 
    f2 :: MatsubaraFunction{GD, SD, DD, Q}
    ;
    )  :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

    check_shape_grid!(f1, f2)
    return MatsubaraFunction(f1.grids, f1.shape, f1.data .+ f2.data)
end

function Base.:+(
    f1 :: MatsubaraFunction{GD, SD, DD, Q}, 
    f2 :: MatsubaraFunction{GD, SD, DD, Q}
    )  :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

    return add(f1, f2)
end

"""
    function add!(
        f1 :: MatsubaraFunction{GD, SD, DD, Q}, 
        f2 :: MatsubaraFunction{GD, SD, DD, Q}
        )  :: Nothing where {GD, SD, DD, Q <: Number}

Inplace addition of two MatsubaraFunction (`f1 += f2`)
"""
function add!(
    f1 :: MatsubaraFunction{GD, SD, DD, Q}, 
    f2 :: MatsubaraFunction{GD, SD, DD, Q}
    )  :: Nothing where {GD, SD, DD, Q <: Number}

    check_shape_grid!(f1, f2)
    f1.data .+= f2.data

    return nothing 
end

"""
    function subtract(
        f1 :: MatsubaraFunction{GD, SD, DD, Q}, 
        f2 :: MatsubaraFunction{GD, SD, DD, Q}
        )  :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

Subtraction of two MatsubaraFunction, returns new MatsubaraFunction
"""
function subtract(
    f1 :: MatsubaraFunction{GD, SD, DD, Q}, 
    f2 :: MatsubaraFunction{GD, SD, DD, Q}
    )  :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

    check_shape_grid!(f1, f2)
    return MatsubaraFunction(f1.grids, f1.shape, f1.data .- f2.data)
end

function Base.:-(
    f1 :: MatsubaraFunction{GD, SD, DD, Q}, 
    f2 :: MatsubaraFunction{GD, SD, DD, Q}
    )  :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

    return subtract(f1, f2)
end

"""
    function subtract!(
        f1 :: MatsubaraFunction{GD, SD, DD, Q}, 
        f2 :: MatsubaraFunction{GD, SD, DD, Q}
        )  :: Nothing where {GD, SD, DD, Q <: Number}

Inplace subtraction of two MatsubaraFunction (`f1 -= f2`)
"""
function subtract!(
    f1 :: MatsubaraFunction{GD, SD, DD, Q}, 
    f2 :: MatsubaraFunction{GD, SD, DD, Q}
    )  :: Nothing where {GD, SD, DD, Q <: Number}

    check_shape_grid!(f1, f2)
    f1.data .-= f2.data

    return nothing 
end

"""
    function mult(
        f   :: MatsubaraFunction{GD, SD, DD, Q},
        val :: Qp
        )   :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number, Qp <: Number}

Multiplication of MatsubaraFunction with scalar, returns new MatsubaraFunction
"""
function mult(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    val :: Qp
    )   :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number, Qp <: Number}

    # type promotion checked by Base
    return MatsubaraFunction(f.grids, f.shape, val .* f.data)
end

function Base.:*(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    val :: Qp
    )   :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number, Qp <: Number}

    return mult(f, val)
end

function Base.:*(
    val :: Qp,
    f   :: MatsubaraFunction{GD, SD, DD, Q}
    )   :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number, Qp <: Number}

    return mult(f, val)
end

"""
    function mult!(
        f   :: MatsubaraFunction{GD, SD, DD, Q},
        val :: Qp
        )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

Inplace multiplication of MatsubaraFunction with scalar (`f *= val`)
"""
function mult!(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    val :: Qp
    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    # type promotion checked by Base
    f.data .*= val 

    return nothing
end
 
"""
    function set!(
        f   :: MatsubaraFunction{GD, SD, DD, Q},
        val :: Qp,
        )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

Initialize MatsubaraFunction with `val`
"""
function set!(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    val :: Qp,
    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    # type promotion checked by Base
    f.data .= val

    return nothing
end

"""
    function set!(
        f   :: MatsubaraFunction{GD, SD, DD, Q},
        arr :: Array{Qp, DD},
        )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

Initialize MatsubaraFunction with `arr`
"""
function set!(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    arr :: Array{Qp, DD},
    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    # type promotion checked by Base
    f.data .= arr

    return nothing
end

"""
    function set!(
        f1 :: MatsubaraFunction{GD, SD, DD, Q},
        f2 :: MatsubaraFunction{GD, SD, DD, Q},
        )  :: Nothing where {GD, SD, DD, Q <: Number}

Initialize MatsubaraFunction with another MatsubaraFunction (`f1 = f2`)
"""
function set!(
    f1 :: MatsubaraFunction{GD, SD, DD, Q},
    f2 :: MatsubaraFunction{GD, SD, DD, Q},
    )  :: Nothing where {GD, SD, DD, Q <: Number}

    check_shape_grid!(f1, f2)
    f1.data .= f2.data

    return nothing
end

# value comparison
function Base.:(==)(    
    f1 :: MatsubaraFunction{GD, SD, DD, Q},
    f2 :: MatsubaraFunction{GD, SD, DD, Q}
    )  :: Bool where {GD, SD, DD, Q <: Number}

    check_shape_grid!(f1, f2) 
    return f1.data ≈ f2.data 
end

"""
    function flatten(
        f :: MatsubaraFunction{GD, SD, DD, Q}
        ) :: Vector{Q} where {GD, SD, DD, Q <: Number}

Flatten data array of MatsubaraFunction and return vector of the corresponding data type
"""
function flatten(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: Vector{Q} where {GD, SD, DD, Q <: Number}

    x  = Vector{Q}(undef, length(f))
    x .= @view f.data[:]

    return x
end

"""
    function flatten!(
        f :: MatsubaraFunction{GD, SD, DD, Q},
        x :: AbstractVector
        ) :: Nothing where {GD, SD, DD, Q <: Number}

Flatten data array of MatsubaraFunction into vector
"""
function flatten!(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    x :: AbstractVector
    ) :: Nothing where {GD, SD, DD, Q <: Number}

    # dimension checked by Base
    x .= @view f.data[:]
    return nothing
end

"""
    function unflatten!(
        f :: MatsubaraFunction{GD, SD, DD, Q},
        x :: AbstractVector
        ) :: Nothing where {GD, SD, DD, Q <: Number}

Initialize data array of MatsubaraFunction from vector
"""
function unflatten!(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    x :: AbstractVector
    ) :: Nothing where {GD, SD, DD, Q <: Number}
    
    f.data[:] .= x
    return nothing
end