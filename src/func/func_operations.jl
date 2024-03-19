# ∞-norm
#-------------------------------------------------------------------------------#

"""
    function absmax(f :: MeshFunction) :: Float64

Returns ∞-norm of `f.data`
"""
function absmax(f :: MeshFunction) :: Float64
    return norm(f.data, inf)
end

"""
    function arg_absmax(f :: MeshFunction{DD, Q, AT}) :: CartesianIndex{DD} where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

Returns cartesian index of `absmax(f)`
"""
function arg_absmax(f :: MeshFunction{DD, Q, AT}) :: CartesianIndex{DD} where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}
    return argmax(abs.(f.data))
end

# debugging
#-------------------------------------------------------------------------------#

function debug_f1_f2(f1 :: MeshFunction{DD, Q, AT}, f2 :: MeshFunction{DD, Q, BT}
    ) where {DD, Q <: Number, AT <: AbstractArray{Q, DD}, BT <: AbstractArray{Q, DD}}

    @DEBUG size(f1.data) == size(f2.data) "Size of data arrays not equal"
    
    for i in 1 : DD 
        @DEBUG meshes(f1, i) == meshes(f2, i) "Meshes are different"
    end

    return nothing 
end

# comparison
#-------------------------------------------------------------------------------#

function Base.:(==)(f1 :: MeshFunction, f2 :: MeshFunction)
    debug_f1_f2(f1, f2) 
    return f1.data ≈ f2.data 
end

# addition
#-------------------------------------------------------------------------------#

"""
    function add(f1 :: MeshFunction, f2 :: MeshFunction) :: MeshFunction

Addition of two MeshFunction, returns new MeshFunction. For brevity, use f1 + f2.
"""
function add(f1 :: MeshFunction, f2 :: MeshFunction) :: MeshFunction
    debug_f1_f2(f1, f2)
    return MeshFunction(f1.data .+ f2.data, meshes(f1)...)
end

function Base.:+(f1 :: MeshFunction, f2 :: MeshFunction)
    return add(f1, f2)
end

"""
    function add!(f1 :: MeshFunction, f2 :: MeshFunction) :: Nothing

Inplace addition of two MeshFunction (`f1 += f2`)
"""
function add!(f1 :: MeshFunction, f2 :: MeshFunction) :: Nothing
    debug_f1_f2(f1, f2)
    f1.data .+= f2.data
    return nothing 
end

# subtraction
#-------------------------------------------------------------------------------#

"""
    function subtract(f1 :: MeshFunction, f2 :: MeshFunction) :: MeshFunction

Subtraction of two MeshFunction, returns new MeshFunction. For brevity, use f1 - f2.
"""
function subtract(f1 :: MeshFunction, f2 :: MeshFunction) :: MeshFunction
    debug_f1_f2(f1, f2)
    return MeshFunction(f1.data .- f2.data, meshes(f1)...)
end

function Base.:-(f1 :: MeshFunction, f2 :: MeshFunction)
    return subtract(f1, f2)
end

"""
    function subtract!(f1 :: MeshFunction, f2 :: MeshFunction) :: Nothing

Inplace subtraction of two MeshFunction (`f1 -= f2`)
"""
function subtract!(f1 :: MeshFunction, f2 :: MeshFunction) :: Nothing
    debug_f1_f2(f1, f2)
    f1.data .-= f2.data
    return nothing 
end

# multiplication
#-------------------------------------------------------------------------------#

"""
    function mult(f :: MeshFunction{DD, Q, AT}, val :: Qp
        ) :: MeshFunction{DD, Q, AT} where {DD, Q <: Number, Qp <: Number, AT <: AbstractArray{Q, DD}}

Multiplication of MeshFunction with scalar, returns new MeshFunction. For brevity, use val * f or f * val.
"""
function mult(f :: MeshFunction{DD, Q, AT}, val :: Qp
    ) :: MeshFunction{DD, Q, AT} where {DD, Q <: Number, Qp <: Number, AT <: AbstractArray{Q, DD}}

    return MeshFunction(val .* f.data, meshes(f)...)
end

function Base.:*(f :: MeshFunction{DD, Q, AT}, val :: Qp) where {DD, Q <: Number, Qp <: Number, AT <: AbstractArray{Q, DD}}
    return mult(f, val)
end

function Base.:*(val :: Qp, f :: MeshFunction{DD, Q, AT}) where {DD, Q <: Number, Qp <: Number, AT <: AbstractArray{Q, DD}}
    return mult(f, val)
end

"""
    function mult!(f :: MeshFunction{DD, Q, AT}, val :: Qp
        ) :: Nothing where {DD, Q <: Number, Qp <: Number, AT <: AbstractArray{Q, DD}}

Inplace multiplication of MeshFunction with scalar (`f *= val`)
"""
function mult!(f :: MeshFunction{DD, Q, AT}, val :: Qp
    ) :: Nothing where {DD, Q <: Number, Qp <: Number, AT <: AbstractArray{Q, DD}}

    f.data .*= val 
    return nothing
end

# set data
#-------------------------------------------------------------------------------#

"""
    function set!(f :: MeshFunction{DD, Q, AT}, val :: Qp
        ) :: Nothing where {DD, Q <: Number, Qp <: Number, AT <: AbstractArray{Q, DD}}

Initialize MeshFunction with `val`
"""
function set!(f :: MeshFunction{DD, Q, AT}, val :: Qp
    ) :: Nothing where {DD, Q <: Number, Qp <: Number, AT <: AbstractArray{Q, DD}}

    f.data .= val
    return nothing
end

"""
    function set!(f :: MeshFunction{DD, Q, AT}, arr :: Array{Qp, DD}
        ) :: Nothing where {DD, Q <: Number, Qp <: Number, AT <: AbstractArray{Q, DD}}

Initialize MeshFunction with `arr`
"""
function set!(f :: MeshFunction{DD, Q, AT}, arr :: Array{Qp, DD}
    ) :: Nothing where {DD, Q <: Number, Qp <: Number, AT <: AbstractArray{Q, DD}}

    f.data .= arr
    return nothing
end

"""
    function set!(f1 :: MeshFunction, f2 :: MeshFunction) :: Nothing

Initialize MeshFunction with another MeshFunction (`f1 = f2`)
"""
function set!(f1 :: MeshFunction, f2 :: MeshFunction) :: Nothing
    debug_f1_f2(f1, f2)
    f1.data .= f2.data
    return nothing
end

# (un-) flatten
#-------------------------------------------------------------------------------#

"""
    function flatten(f :: MeshFunction{DD, Q, AT}) :: Vector{Q} where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

Flatten data array of MeshFunction and return vector of the corresponding data type
"""
function flatten(f :: MeshFunction{DD, Q, AT}) :: Vector{Q} where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}
    x  = Vector{Q}(undef, length(f.data))
    x .= @view f.data[:]
    return x
end

"""
    function flatten!(f :: MeshFunction{DD, Q, AT}, x :: T
        ) :: Nothing where {DD, Q <: Number, T <: AbstractVector{Q}, AT <: AbstractArray{Q, DD}}

Flatten data array of MeshFunction into vector
"""
function flatten!(f :: MeshFunction{DD, Q, AT}, x :: T
    ) :: Nothing where {DD, Q <: Number, T <: AbstractVector{Q}, AT <: AbstractArray{Q, DD}}

    x .= @view f.data[:]
    return nothing
end

"""
    function unflatten!(f :: MeshFunction{DD, Q, AT}, x :: T
        ) :: Nothing where {DD, Q <: Number, T <: AbstractVector{Q}, AT <: AbstractArray{Q, DD}}

Initialize data array of MeshFunction from vector
"""
function unflatten!(f :: MeshFunction{DD, Q, AT}, x :: T
    ) :: Nothing where {DD, Q <: Number, T <: AbstractVector{Q}, AT <: AbstractArray{Q, DD}}
    
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