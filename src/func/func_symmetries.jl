# operation type
#-------------------------------------------------------------------------------#

"""
    struct Operation 

Operation type with fields:
* `sgn :: Bool` : change sign?
* `con :: Bool` : complex conjugation?
"""
struct Operation 
    sgn :: Bool 
    con :: Bool

    function Operation(sgn :: Bool, con :: Bool)
        return new(sgn, con)
    end 

    function Operation(; sgn :: Bool = false, con :: Bool = false)
        return Operation(sgn, con)
    end
end

"""
    sgn(op :: Operation) :: Bool

Return `op.sgn`
"""
sgn(op :: Operation) :: Bool = op.sgn

"""
    con(op :: Operation) :: Bool

Return `op.con`
"""
con(op :: Operation) :: Bool = op.con

function Base.:*(op1 :: Operation, op2 :: Operation) 
    return Operation(xor(sgn(op1), sgn(op2)), xor(con(op1), con(op2)))
end

function (op :: Operation)(x :: Q) where {Q <: Number}
    if sgn(op); return con(op) ? -conj(x) : -x; end
    return con(op) ? conj(x) : x
end

# symmetry type
#-------------------------------------------------------------------------------#

# a symmetry takes MeshPoint value types & tensor indices and 
# returns MeshPoint value types & tensor indices & an operation
"""
    struct Symmetry{GD, SD}

Symmetry type with fields:
* `f :: Function`
"""
struct Symmetry{GD, SD}
    f :: Function
end

# explicit return type to enforce proper implementation
function (S :: Symmetry{GD, SD})(
    w :: NTuple{GD, Union{<: AbstractValue}},
    x :: NTuple{SD, Int64}
    ) :: Tuple{NTuple{GD, Union{<: AbstractValue}}, NTuple{SD, Int64}, Operation} where {GD, SD}

    return S.f(w, x)
end

function reduce(
    w           :: NTuple{GD, Union{<: AbstractValue}},
    x           :: NTuple{SD, Int64},
    op          :: Operation,
    f           :: MeshFunction{GD, SD, DD, Q, AT},
    checked     :: Array{Bool, DD},
    symmetries  :: Vector{Symmetry{GD, SD}},
    class       :: Vector{Tuple{Int64, Operation}},
    path_length :: Int64
    ;
    max_length  :: Int64 = 0
    ) where {GD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}
 
    for S in symmetries 
        wp, xp, opp = S(w, x)
        new_op      = opp * op

        if all(ntuple(i -> is_inbounds(wp[i], meshes(f, i)), GD))
            idx = LinearIndex(f, ntuple(i -> mesh_index(wp[i], meshes(f, i)), GD)..., xp...)

            if !checked[idx]
                checked[idx] = true 

                # add to symmetry class, reset path length and keep going 
                push!(class, (idx, new_op))
                reduce(wp, xp, new_op, f, checked, symmetries, class, 0; max_length)
            end 

        # if index not valid, increment path length and keep going 
        elseif path_length < max_length
            reduce(wp, xp, new_op, f, checked, symmetries, class, path_length + 1; max_length)
        end 
    end 
end

# symmetry group type
#-------------------------------------------------------------------------------#

"""
    SymmetryGroup{GD, SD, DD, Q <: Number}

SymmetryGroup type with fields:
* `classes :: Vector{Vector{Tuple{Int64, Operation}}}` : collections of symmetry equivalent elements
* `speedup :: Float64`                                 : expected speedup from the symmetry reduction
"""
struct SymmetryGroup{GD, SD, DD, Q <: Number}
    classes :: Vector{Vector{Tuple{Int64, Operation}}}
    speedup :: Float64

    function SymmetryGroup{GD, SD, DD, Q}(
        classes :: Vector{Vector{Tuple{Int64, Operation}}}, speedup :: Float64
        ) where {GD, SD, DD, Q <: Number}

        return new{GD, SD, DD, Q}(classes, speedup)
    end 

    function SymmetryGroup(f :: MeshFunction{GD, SD, DD, Q, AT}) where {GD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}
        return new{GD, SD, DD, Q}([[(idx, Operation())] for idx in eachindex(f.data)], 1.0)
    end
 
    function SymmetryGroup(
        symmetries :: Vector{Symmetry{GD, SD}},
        f          :: MeshFunction{GD, SD, DD, Q, AT}
        ;
        max_length :: Int64 = 0
        ) where {GD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

        checked  = Array{Bool, DD}(undef, size(f.data))
        classes  = Vector{Tuple{Int64, Operation}}[]
        checked .= false

        for idx in eachindex(f.data)
            if !checked[idx]
                checked[idx] = true
                w, x         = to_meshes(f, idx)
                class        = [(idx, Operation())]

                reduce(value.(w), x, Operation(), f, checked, symmetries, class, 0; max_length)
                push!(classes, class)
            end 
        end

        return new{GD, SD, DD, Q}(classes, length(f.data) / length(classes))
    end 
end

# symmetrize data array of the MeshFunction, return error estimate
function (SG :: SymmetryGroup{GD, SD, DD, Q})(
    f :: MeshFunction{GD, SD, DD, Q, AT}
    ) where {GD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    err = 0.0

    for class in SG.classes
        ref = f[class[1][1]]

        for idx in 2 : length(class)
            idx, op = class[idx]
            new_val = op(ref)
            new_err = abs(f[idx] - new_val)

            if new_err > err err = new_err end
            f[idx] = new_val
        end 
    end 

    return err
end

"""
    function get_reduced(
        SG :: SymmetryGroup{GD, SD, DD, Q},
        f  :: MeshFunction{GD, SD, DD, Q, AT}
        )  :: Vector{Q} where {GD, SD, DD, Q <: Number}

Calculate symmetry reduced representation of MeshFunction
"""
function get_reduced(
    SG :: SymmetryGroup{GD, SD, DD, Q},
    f  :: MeshFunction{GD, SD, DD, Q, AT}
    )  :: Vector{Q} where {GD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    return Q[f[first(first(SG.classes[cl]))] for cl in eachindex(SG.classes)] 
end

"""
    function init_from_reduced!(
        SG   :: SymmetryGroup{GD, SD, DD, Q},
        f    :: MeshFunction{GD, SD, DD, Q, AT},
        fvec :: AbstractVector{Q}
        )    :: Nothing where {GD, SD, DD, Q <: Number}

Initialize MeshFunction from symmetry reduced representation
""" 
function init_from_reduced!(
    SG   :: SymmetryGroup{GD, SD, DD, Q},
    f    :: MeshFunction{GD, SD, DD, Q, AT},
    fvec :: AbstractVector{Q}
    )    :: Nothing where {GD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    for cl in eachindex(fvec)
        f[first(first(SG.classes[cl]))] = fvec[cl]
    end 

    SG(f)
    return nothing 
end

# init function type
#-------------------------------------------------------------------------------#

"""
    struct InitFunction{GD, SD, Q <: Number}

InitFunction type with fields:
* `f :: Function` 
"""
struct InitFunction{GD, SD, Q <: Number}
    f :: Function
end

# explicit return type to enforce proper implementation
function (I :: InitFunction{GD, SD, Q})(
    w :: NTuple{GD, Union{<: AbstractValue}},
    x :: NTuple{SD, Int64}
    ) :: Q where {GD, SD, Q <: Number}

    return I.f(w, x)
end

# symmetrize MeshFunction from evaluation of InitFunction
function (SG :: SymmetryGroup{GD, SD, DD, Q})(
    f        :: MeshFunction{GD, SD, DD, Q, AT},
    I        :: InitFunction{GD, SD, Q}
    ;
    mode     :: Symbol = :serial,
    minbatch :: Int64  = 1
    )        :: Nothing where {GD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    if mode === :serial 
        for class in SG.classes
            w, x = to_meshes(f, class[1][1])
            f[class[1][1]] = I(value.(w), x)
        end 

    elseif mode === :polyester
        @batch minbatch = minbatch for class in SG.classes
            w, x = to_meshes(f, class[1][1])
            f[class[1][1]] = I(value.(w), x)
        end

    elseif mode === :threads
        Threads.@threads for class in SG.classes
            w, x = to_meshes(f, class[1][1])
            f[class[1][1]] = I(value.(w), x)
        end

    elseif mode === :hybrid 
        set!(f, 0.0)

        Threads.@threads for clidx in mpi_split(1 : length(SG.classes))
            w, x = to_meshes(f, SG.classes[clidx][1][1])
            f[SG.classes[clidx][1][1]] = I(value.(w), x)
        end

        mpi_allreduce!(f)

    else 
        error("Mode $(mode) unknown. Possible options are: serial, polyester, threads and hybrid")
    end

    SG(f)
    return nothing 
end

# export
#-------------------------------------------------------------------------------#

export
    Operation,
    sgn,
    con,
    Symmetry,
    SymmetryGroup,
    get_reduced,
    init_from_reduced!,
    InitFunction