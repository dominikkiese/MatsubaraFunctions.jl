"""
    struct MatsubaraOperation 

MatsubaraOperation type with fields:
* `sgn :: Bool` : change sign?
* `con :: Bool` : complex conjugation?
"""
struct MatsubaraOperation 
    sgn :: Bool 
    con :: Bool

    function MatsubaraOperation(
        sgn :: Bool,
        con :: Bool
        )   :: MatsubaraOperation 

        return new(sgn, con)
    end 

    function MatsubaraOperation(; 
        sgn :: Bool = false,
        con :: Bool = false
        )   :: MatsubaraOperation
        
        return MatsubaraOperation(sgn, con)
    end
end
 
"""
    sgn(op :: MatsubaraOperation) :: Bool

Return `op.sgn`
"""
sgn(op :: MatsubaraOperation) :: Bool = op.sgn

"""
    con(op :: MatsubaraOperation) :: Bool

Return `op.con`
"""
con(op :: MatsubaraOperation) :: Bool = op.con

function Base.:*(
    op1 :: MatsubaraOperation,
    op2 :: MatsubaraOperation
    )   :: MatsubaraOperation 

    return MatsubaraOperation(xor(sgn(op1), sgn(op2)), xor(con(op1), con(op2)))
end

function (op :: MatsubaraOperation)(
    x :: Q
    ) :: Q where {Q <: Number}

    if sgn(op); return con(op) ? -conj(x) : -x; end
    return con(op) ? conj(x) : x
end

#----------------------------------------------------------------------------------------------#

# a symmetry takes Matsubara frequencies / indices & tensor indices and 
# returns Matsubara frequencies / indices & tensor indices & operation
"""
    struct MatsubaraSymmetry{GD, SD}

MatsubaraSymmetry type with fields:
* `f :: Function`
"""
struct MatsubaraSymmetry{GD, SD}
    f :: Function
end

function (S :: MatsubaraSymmetry{GD, SD})(
    w :: NTuple{GD, MatsubaraFrequency},
    x :: NTuple{SD, Int64}
    ) :: Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}, MatsubaraOperation} where {GD, SD}

    return S.f(w, x)
end

function reduce(
    w           :: NTuple{GD, MatsubaraFrequency},
    x           :: NTuple{SD, Int64},
    op          :: MatsubaraOperation,
    f           :: MatsubaraFunction{GD, SD, DD, Q},
    checked     :: Array{Bool, DD},
    symmetries  :: Vector{MatsubaraSymmetry{GD, SD}},
    class       :: Vector{Tuple{Int64, MatsubaraOperation}},
    path_length :: Int64
    ;
    max_length  :: Int64 = 0
    )           :: Nothing where {GD, SD, DD, Q <: Number}
 
    for S in symmetries 
        wp, xp, opp = S(w, x)
        new_op      = opp * op

        if !any(ntuple(i -> !is_inbounds(wp[i], grids(f, i)), GD))
            idx = LinearIndex(f, ntuple(i -> grid_index(wp[i], grids(f, i)), GD)..., xp...)

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

#----------------------------------------------------------------------------------------------#

"""
    MatsubaraSymmetryGroup{GD, SD, DD, Q <: Number}

MatsubaraSymmetryGroup type with fields:
* `classes :: Vector{Vector{Tuple{Int64, MatsubaraOperation}}}` : collections of symmetry equivalent elements
* `speedup :: Float64`                                          : expected speedup from the symmetry reduction
"""
struct MatsubaraSymmetryGroup{GD, SD, DD, Q <: Number}
    classes :: Vector{Vector{Tuple{Int64, MatsubaraOperation}}}
    speedup :: Float64

    function MatsubaraSymmetryGroup{GD, SD, DD, Q}(
        classes :: Vector{Vector{Tuple{Int64, MatsubaraOperation}}},
        speedup :: Float64
        )       :: MatsubaraSymmetryGroup{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

        return new{GD, SD, DD, Q}(classes, speedup)
    end 

    function MatsubaraSymmetryGroup(
        f :: MatsubaraFunction{GD, SD, DD, Q}
        ) :: MatsubaraSymmetryGroup{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

        return new{GD, SD, DD, Q}([[(idx, MatsubaraOperation())] for idx in eachindex(f.data)], 1.0)
    end
 
    function MatsubaraSymmetryGroup(
        symmetries :: Vector{MatsubaraSymmetry{GD, SD}},
        f          :: MatsubaraFunction{GD, SD, DD, Q}
        ;
        max_length :: Int64 = 0
        )          :: MatsubaraSymmetryGroup{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

        checked  = Array{Bool, DD}(undef, data_shape(f)) 
        classes  = Vector{Tuple{Int64, MatsubaraOperation}}[]
        checked .= false

        for idx in eachindex(f.data)
            if !checked[idx]
                checked[idx] = true
                w, x         = to_Matsubara(f, idx)
                class        = [(idx, MatsubaraOperation())]

                reduce(w, x, MatsubaraOperation(), f, checked, symmetries, class, 0; max_length)
                push!(classes, class)
            end 
        end

        return new{GD, SD, DD, Q}(classes, length(f.data) / length(classes))
    end 
end

# symmetrize data array of the MatsubaraFunction, return error estimate
function (SG :: MatsubaraSymmetryGroup{GD, SD, DD, Q})(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: Float64 where {GD, SD, DD, Q <: Number}

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

#----------------------------------------------------------------------------------------------#

"""
    struct MatsubaraInitFunction{GD, SD, Q <: Number}

MatsubaraInitFunction type with fields:
* `f :: Function` 
"""
struct MatsubaraInitFunction{GD, SD, Q <: Number}
    f :: Function
end

function (I :: MatsubaraInitFunction{GD, SD, Q})(
    w :: NTuple{GD, MatsubaraFrequency},
    x :: NTuple{SD, Int64}
    ) :: Q where {GD, SD, Q <: Number}

    return I.f(w, x)
end

# symmetrize MatsubaraFunction from evaluation of MatsubaraInitFunction
function (SG :: MatsubaraSymmetryGroup{GD, SD, DD, Q})(
    f        :: MatsubaraFunction{GD, SD, DD, Q},
    I        :: MatsubaraInitFunction{GD, SD, Q}
    ;
    mode     :: Symbol = :serial,
    minbatch :: Int64  = 1
    )        :: Nothing where {GD, SD, DD, Q <: Number}

    if mode === :serial 
        for class in SG.classes
            f[class[1][1]] = I(to_Matsubara(f, class[1][1])...)
        end 

    elseif mode === :polyester
        @batch minbatch = minbatch for class in SG.classes
            f[class[1][1]] = I(to_Matsubara(f, class[1][1])...)
        end

    elseif mode === :threads
        Threads.@threads for class in SG.classes
            f[class[1][1]] = I(to_Matsubara(f, class[1][1])...)
        end

    elseif mode === :hybrid 
        set!(f, 0.0)

        Threads.@threads for clidx in mpi_split(1 : length(SG.classes))
            f[SG.classes[clidx][1][1]] = I(to_Matsubara(f, SG.classes[clidx][1][1])...)
        end

        mpi_allreduce!(f)

    else 
        error("Mode $(mode) unknown. Possible options are: serial, polyester, threads and hybrid")
    end

    SG(f)

    return nothing 
end

#----------------------------------------------------------------------------------------------#

export
    MatsubaraOperation,
    sgn,
    con,
    MatsubaraSymmetry,
    MatsubaraSymmetryGroup,
    MatsubaraInitFunction