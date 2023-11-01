# ∞-norm
#-------------------------------------------------------------------------------#

"""
    function absmax(
        f :: MeshFunction
        ) :: Float64

Returns ∞-norm of `f.data`
"""
function absmax(
    f :: MeshFunction
    ) :: Float64

    return norm(f.data, inf)
end

"""
    function arg_absmax(
        f :: MeshFunction{MD, SD, DD, Q}
        ) :: CartesianIndex{DD} where {MD, SD, DD, Q <: Number}

Returns cartesian index of `absmax(f)`
"""
function arg_absmax(
    f :: MeshFunction{MD, SD, DD, Q}
    ) :: CartesianIndex{DD} where {MD, SD, DD, Q <: Number}

    return argmax(abs.(f.data))
end

# debugging
#-------------------------------------------------------------------------------#

function debug_f1_f2(
    f1 :: MeshFunction{GD, SD, DD, Q}, 
    f2 :: MeshFunction{GD, SD, DD, Q}
    )  :: Nothing where {GD, SD, DD, Q <: Number}

    @DEBUG size(f1.data) == size(f2.data) "Size of data arrays not equal"
    @DEBUG shape(f1) == shape(f2) "Index structures are not equal"

    for i in 1 : GD 
        @DEBUG meshes(f1, i) == meshes(f2, i) "Meshes are different"
    end

    return nothing 
end

# comparison
#-------------------------------------------------------------------------------#

function Base.:(==)(    
    f1 :: MeshFunction,
    f2 :: MeshFunction
    )  :: Bool

    debug_f1_f2(f1, f2) 
    return f1.data ≈ f2.data 
end

# addition
#-------------------------------------------------------------------------------#

"""
    function add(
        f1 :: MeshFunction, 
        f2 :: MeshFunction
        )  :: MeshFunction

Addition of two MeshFunction, returns new MeshFunction. For brevity, use f1 + f2.
"""
function add(
    f1 :: MeshFunction, 
    f2 :: MeshFunction
    )  :: MeshFunction

    debug_f1_f2(f1, f2)
    return MeshFunction(Mesh.(meshes(f1)), shape(f1), f1.data .+ f2.data)
end

function Base.:+(
    f1 :: MeshFunction, 
    f2 :: MeshFunction
    )  :: MeshFunction

    return add(f1, f2)
end

"""
    function add!(
        f1 :: MeshFunction, 
        f2 :: MeshFunction
        )  :: Nothing

Inplace addition of two MeshFunction (`f1 += f2`)
"""
function add!(
    f1 :: MeshFunction, 
    f2 :: MeshFunction
    )  :: Nothing

    debug_f1_f2(f1, f2)
    f1.data .+= f2.data
    return nothing 
end

# subtraction
#-------------------------------------------------------------------------------#

"""
    function subtract(
        f1 :: MeshFunction, 
        f2 :: MeshFunction
        )  :: MeshFunction

Subtraction of two MeshFunction, returns new MeshFunction. For brevity, use f1 - f2.
"""
function subtract(
    f1 :: MeshFunction, 
    f2 :: MeshFunction
    )  :: MeshFunction

    debug_f1_f2(f1, f2)
    return MeshFunction(Mesh.(meshes(f1)), shape(f1), f1.data .- f2.data)
end

function Base.:-(
    f1 :: MeshFunction, 
    f2 :: MeshFunction
    )  :: MeshFunction

    return subtract(f1, f2)
end

"""
    function subtract!(
        f1 :: MeshFunction, 
        f2 :: MeshFunction
        )  :: Nothing

Inplace subtraction of two MeshFunction (`f1 -= f2`)
"""
function subtract!(
    f1 :: MeshFunction, 
    f2 :: MeshFunction
    )  :: Nothing

    debug_f1_f2(f1, f2)
    f1.data .-= f2.data

    return nothing 
end

# multiplication
#-------------------------------------------------------------------------------#

"""
    function mult(
        f   :: MeshFunction{GD, SD, DD, Q},
        val :: Qp
        )   :: MeshFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number, Qp <: Number}

Multiplication of MeshFunction with scalar, returns new MeshFunction. For brevity, use val * f or f * val.
"""
function mult(
    f   :: MeshFunction{GD, SD, DD, Q},
    val :: Qp
    )   :: MeshFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number, Qp <: Number}

    return MeshFunction(Mesh.(meshes(f)), shape(f), val .* f.data)
end

function Base.:*(
    f   :: MeshFunction{GD, SD, DD, Q},
    val :: Qp
    )   :: MeshFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number, Qp <: Number}

    return mult(f, val)
end

function Base.:*(
    val :: Qp,
    f   :: MeshFunction{GD, SD, DD, Q}
    )   :: MeshFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number, Qp <: Number}

    return mult(f, val)
end

"""
    function mult!(
        f   :: MeshFunction{GD, SD, DD, Q},
        val :: Qp
        )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

Inplace multiplication of MeshFunction with scalar (`f *= val`)
"""
function mult!(
    f   :: MeshFunction{GD, SD, DD, Q},
    val :: Qp
    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    f.data .*= val 
    return nothing
end

# set data
#-------------------------------------------------------------------------------#

"""
    function set!(
        f   :: MeshFunction{GD, SD, DD, Q},
        val :: Qp,
        )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

Initialize MeshFunction with `val`
"""
function set!(
    f   :: MeshFunction{GD, SD, DD, Q},
    val :: Qp,
    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    f.data .= val
    return nothing
end

"""
    function set!(
        f   :: MeshFunction{GD, SD, DD, Q},
        arr :: Array{Qp, DD},
        )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

Initialize MeshFunction with `arr`
"""
function set!(
    f   :: MeshFunction{GD, SD, DD, Q},
    arr :: Array{Qp, DD},
    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    f.data .= arr
    return nothing
end

"""
    function set!(
        f1 :: MeshFunction,
        f2 :: MeshFunction
        )  :: Nothing

Initialize MeshFunction with another MeshFunction (`f1 = f2`)
"""
function set!(
    f1 :: MeshFunction,
    f2 :: MeshFunction
    )  :: Nothing

    debug_f1_f2(f1, f2)
    f1.data .= f2.data
    return nothing
end

# (un-) flatten
#-------------------------------------------------------------------------------#

"""
    function flatten(
        f :: MeshFunction{GD, SD, DD, Q}
        ) :: Vector{Q} where {GD, SD, DD, Q <: Number}

Flatten data array of MeshFunction and return vector of the corresponding data type
"""
function flatten(
    f :: MeshFunction{GD, SD, DD, Q}
    ) :: Vector{Q} where {GD, SD, DD, Q <: Number}

    x  = Vector{Q}(undef, length(f.data))
    x .= @view f.data[:]
    return x
end

"""
    function flatten!(
        f :: MeshFunction,
        x :: AbstractVector
        ) :: Nothing

Flatten data array of MeshFunction into vector
"""
function flatten!(
    f :: MeshFunction,
    x :: AbstractVector
    ) :: Nothing

    x .= @view f.data[:]
    return nothing
end

"""
    function unflatten!(
        f :: MeshFunction,
        x :: AbstractVector
        ) :: Nothing

Initialize data array of MeshFunction from vector
"""
function unflatten!(
    f :: MeshFunction,
    x :: AbstractVector
    ) :: Nothing
    
    f.data[:] .= x
    return nothing
end

# export
#-------------------------------------------------------------------------------#

export 
    absmax,
    arg_absmax,
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