function check_shape_grid!(
    f1 :: MatsubaraFunction{GD, SD, DD, Q}, 
    f2 :: MatsubaraFunction{GD, SD, DD, Q}
    )  :: Nothing where {GD, SD, DD, Q <: Number}

    @DEBUG data_shape(f1) == data_shape(f2) "Data shapes must be equal"

    for i in 1 : GD 
        @DEBUG typeof(grids(f1, i)) == typeof(grids(f2, i)) "Grids must have same particle type" 
        @DEBUG temperature(grids(f1, i)) == temperature(grids(f2, i)) "Grids must have same temperature" 
        @DEBUG index_range(grids(f1, i)) == index_range(grids(f2, i)) "Grids must have same index range" 
    end

    return nothing 
end

"""
    function add(
        f1 :: MatsubaraFunction, 
        f2 :: MatsubaraFunction
        )  :: MatsubaraFunction

Addition of two MatsubaraFunction, returns new MatsubaraFunction. For brevity, use f1 + f2.
"""
function add(
    f1 :: MatsubaraFunction, 
    f2 :: MatsubaraFunction
    )  :: MatsubaraFunction

    check_shape_grid!(f1, f2)
    return MatsubaraFunction(grids(f1), shape(f1), f1.data .+ f2.data)
end

function Base.:+(
    f1 :: MatsubaraFunction, 
    f2 :: MatsubaraFunction
    )  :: MatsubaraFunction

    return add(f1, f2)
end

"""
    function add!(
        f1 :: MatsubaraFunction, 
        f2 :: MatsubaraFunction
        )  :: Nothing

Inplace addition of two MatsubaraFunction (`f1 += f2`)
"""
function add!(
    f1 :: MatsubaraFunction, 
    f2 :: MatsubaraFunction
    )  :: Nothing

    check_shape_grid!(f1, f2)
    f1.data .+= f2.data

    return nothing 
end

"""
    function subtract(
        f1 :: MatsubaraFunction, 
        f2 :: MatsubaraFunction
        )  :: MatsubaraFunction

Subtraction of two MatsubaraFunction, returns new MatsubaraFunction. For brevity, use f1 - f2.
"""
function subtract(
    f1 :: MatsubaraFunction, 
    f2 :: MatsubaraFunction
    )  :: MatsubaraFunction

    check_shape_grid!(f1, f2)
    return MatsubaraFunction(grids(f1), shape(f1), f1.data .- f2.data)
end

function Base.:-(
    f1 :: MatsubaraFunction, 
    f2 :: MatsubaraFunction
    )  :: MatsubaraFunction

    return subtract(f1, f2)
end

"""
    function subtract!(
        f1 :: MatsubaraFunction, 
        f2 :: MatsubaraFunction
        )  :: Nothing

Inplace subtraction of two MatsubaraFunction (`f1 -= f2`)
"""
function subtract!(
    f1 :: MatsubaraFunction, 
    f2 :: MatsubaraFunction
    )  :: Nothing

    check_shape_grid!(f1, f2)
    f1.data .-= f2.data

    return nothing 
end

"""
    function mult(
        f   :: MatsubaraFunction{GD, SD, DD, Q},
        val :: Qp
        )   :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number, Qp <: Number}

Multiplication of MatsubaraFunction with scalar, returns new MatsubaraFunction. For brevity, use val * f or f * val.
"""
function mult(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    val :: Qp
    )   :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number, Qp <: Number}

    return MatsubaraFunction(grids(f), shape(f), val .* f.data)
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

    OffsetArrays.no_offset_view(f.data) .= arr
    return nothing
end

"""
    function set!(
        f1 :: MatsubaraFunction,
        f2 :: MatsubaraFunction
        )  :: Nothing

Initialize MatsubaraFunction with another MatsubaraFunction (`f1 = f2`)
"""
function set!(
    f1 :: MatsubaraFunction,
    f2 :: MatsubaraFunction
    )  :: Nothing

    check_shape_grid!(f1, f2)
    f1.data .= f2.data

    return nothing
end

function Base.:(==)(    
    f1 :: MatsubaraFunction,
    f2 :: MatsubaraFunction
    )  :: Bool

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
    flatten!(f, x)

    return x
end

"""
    function flatten!(
        f :: MatsubaraFunction,
        x :: AbstractVector
        ) :: Nothing

Flatten data array of MatsubaraFunction into vector
"""
function flatten!(
    f :: MatsubaraFunction,
    x :: AbstractVector
    ) :: Nothing

    f_view = @view f.data[:]
    copyto!(x, firstindex(x), f_view, firstindex(f_view), length(f_view))
    return nothing
end

"""
    function unflatten!(
        f :: MatsubaraFunction,
        x :: AbstractVector
        ) :: Nothing

Initialize data array of MatsubaraFunction from vector
"""
function unflatten!(
    f :: MatsubaraFunction,
    x :: AbstractVector
    ) :: Nothing
    
    f_view = @view f.data[:]
    copyto!(f_view, firstindex(f_view), x, firstindex(x))
    return nothing
end

#----------------------------------------------------------------------------------------------#

export 
    add, 
    add!,
    subtract,
    subtract!,
    mult,
    mult!,
    set!,
    flatten,
    flatten!,
    unflatten!