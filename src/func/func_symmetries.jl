"""
    struct MeshOperation 

MeshOperation type with fields:
* `sgn :: Bool` : change sign?
* `con :: Bool` : complex conjugation?
"""
struct MeshOperation 
    sgn :: Bool 
    con :: Bool

    function MeshOperation(
        sgn :: Bool,
        con :: Bool
        )   :: MeshOperation 

        return new(sgn, con)
    end 

    function MeshOperation(; 
        sgn :: Bool = false,
        con :: Bool = false
        )   :: MeshOperation
        
        return MeshOperation(sgn, con)
    end
end
 
"""
    sgn(op :: MeshOperation) :: Bool

Return `op.sgn`
"""
sgn(op :: MeshOperation) :: Bool = op.sgn

"""
    con(op :: MeshOperation) :: Bool

Return `op.con`
"""
con(op :: MeshOperation) :: Bool = op.con

function Base.:*(
    op1 :: MeshOperation,
    op2 :: MeshOperation
    )   :: MeshOperation 

    return MeshOperation(xor(sgn(op1), sgn(op2)), xor(con(op1), con(op2)))
end

function (op :: MeshOperation)(
    x :: Q
    ) :: Q where {Q <: Number}

    if sgn(op); return con(op) ? -conj(x) : -x; end
    return con(op) ? conj(x) : x
end

#----------------------------------------------------------------------------------------------#

# a symmetry takes Matsubara frequencies / indices & tensor indices and 
# returns Matsubara frequencies / indices & tensor indices & operation
"""
    struct MeshSymmetry{GD, SD}

MeshSymmetry type with fields:
* `f :: Function`
"""
struct MeshSymmetry{GD, SD}
    f :: Function
end

function (S :: MeshSymmetry{GD, SD})(
    w :: NTuple{GD, AbstractValue},
    x :: NTuple{SD, Int64}
    ) :: Tuple{NTuple{GD, Union{AbstractValue, AbstractMeshPoint}}, NTuple{SD, Int64}, MeshOperation} where {GD, SD}

    return S.f(w, x)
end

function (S :: MeshSymmetry{GD, SD})(
    w :: NTuple{GD, Union{AbstractValue, AbstractMeshPoint}},
    x :: NTuple{SD, Int64}
    ) :: Tuple{NTuple{GD, Union{AbstractValue, AbstractMeshPoint}}, NTuple{SD, Int64}, MeshOperation} where {GD, SD}

    return S.f(w, x)
end

function reduce(
    w           :: NTuple{GD, Union{AbstractValue, AbstractMeshPoint}},
    x           :: NTuple{SD, Int64},
    op          :: MeshOperation,
    f           :: MeshFunction{GD, SD, DD, Q},
    checked     :: OffsetArray{Bool, DD, Array{Bool, DD}},
    symmetries  :: Vector{MeshSymmetry{GD, SD}},
    class       :: Vector{Tuple{Int64, MeshOperation}},
    path_length :: Int64
    ;
    max_length  :: Int64 = 0
    )           :: Nothing where {GD, SD, DD, Q <: Number}
 
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

#----------------------------------------------------------------------------------------------#

"""
    MeshSymmetryGroup{GD, SD, DD, Q <: Number}

MeshSymmetryGroup type with fields:
* `classes :: Vector{Vector{Tuple{Int64, MeshOperation}}}` : collections of symmetry equivalent elements
* `speedup :: Float64`                                          : expected speedup from the symmetry reduction
"""
struct MeshSymmetryGroup{GD, SD, DD, Q <: Number}
    classes :: Vector{Vector{Tuple{Int64, MeshOperation}}}
    speedup :: Float64

    function MeshSymmetryGroup{GD, SD, DD, Q}(
        classes :: Vector{Vector{Tuple{Int64, MeshOperation}}},
        speedup :: Float64
        )       :: MeshSymmetryGroup{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

        return new{GD, SD, DD, Q}(classes, speedup)
    end 

    function MeshSymmetryGroup(
        f :: MeshFunction{GD, SD, DD, Q}
        ) :: MeshSymmetryGroup{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

        return new{GD, SD, DD, Q}([[(idx, MeshOperation())] for idx in eachindex(f.data)], 1.0)
    end
    MeshSymmetryGroup
    function MeshSymmetryGroup(
        symmetries :: Vector{MeshSymmetry{GD, SD}},
        f          :: MeshFunction{GD, SD, DD, Q}
        ;
        max_length :: Int64 = 0
        )          :: MeshSymmetryGroup{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

        checked  = OffsetArray(Array{Bool, DD}(undef, size(f)), axes(f)...)
        classes  = Vector{Tuple{Int64, MeshOperation}}[]
        checked .= false

        for idx in eachindex(f.data)
            if !checked[idx]
                checked[idx] = true
                w, x         = to_meshes(f, idx)
                class        = [(idx, MeshOperation())]

                reduce(w, x, MeshOperation(), f, checked, symmetries, class, 0; max_length)
                push!(classes, class)
            end 
        end

        return new{GD, SD, DD, Q}(classes, length(f.data) / length(classes))
    end 
end

# symmetrize data array of the MeshFunction, return error estimate
function (SG :: MeshSymmetryGroup{GD, SD, DD, Q})(
    f :: MeshFunction{GD, SD, DD, Q}
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

"""
    function get_reduced(
        SG :: MeshSymmetryGroup{GD, SD, DD, Q},
        f  :: MeshFunction{GD, SD, DD, Q}
        )  :: Vector{Q} where {GD, SD, DD, Q <: Number}

Calculate symmetry reduced representation of MeshFunction
"""
function get_reduced(
    SG :: MeshSymmetryGroup{GD, SD, DD, Q},
    f  :: MeshFunction{GD, SD, DD, Q}
    )  :: Vector{Q} where {GD, SD, DD, Q <: Number}

    return Q[f[first(first(SG.classes[cl]))] for cl in eachindex(SG.classes)] 
end

"""
    function init_from_reduced!(
        SG   :: MeshSymmetryGroup{GD, SD, DD, Q},
        f    :: MeshFunction{GD, SD, DD, Q},
        fvec :: AbstractVector{Q}
        )    :: Nothing where {GD, SD, DD, Q <: Number}

Initialize MeshFunction from symmetry reduced representation
""" 
function init_from_reduced!(
    SG   :: MeshSymmetryGroup{GD, SD, DD, Q},
    f    :: MeshFunction{GD, SD, DD, Q},
    fvec :: AbstractVector{Q}
    )    :: Nothing where {GD, SD, DD, Q <: Number}

    for cl in eachindex(fvec)
        f[first(first(SG.classes[cl]))] = fvec[cl]
    end 

    SG(f)
    return nothing 
end

#----------------------------------------------------------------------------------------------#

"""
    struct MeshInitFunction{GD, SD, Q <: Number}

MeshInitFunction type with fields:
* `f :: Function` 
"""
struct MeshInitFunction{GD, SD, Q <: Number}
    f :: Function
end

function (I :: MeshInitFunction{GD, SD, Q})(
    w :: NTuple{GD, Union{AbstractValue, AbstractMeshPoint}},
    x :: NTuple{SD, Int64}
    ) :: Q where {GD, SD, Q <: Number}

    return I.f(w, x)
end

# symmetrize MeshFunction from evaluation of MeshInitFunction
function (SG :: MeshSymmetryGroup{GD, SD, DD, Q})(
    f        :: MeshFunction{GD, SD, DD, Q},
    I        :: MeshInitFunction{GD, SD, Q}
    ;
    mode     :: Symbol = :serial,
    minbatch :: Int64  = 1
    )        :: Nothing where {GD, SD, DD, Q <: Number}

    if mode === :serial 
        for class in SG.classes
            f[class[1][1]] = I(to_meshes(f, class[1][1])...)
        end 

    elseif mode === :polyester
        @batch minbatch = minbatch for class in SG.classes
            f[class[1][1]] = I(to_meshes(f, class[1][1])...)
        end

    elseif mode === :threads
        Threads.@threads for class in SG.classes
            f[class[1][1]] = I(to_meshes(f, class[1][1])...)
        end

    elseif mode === :hybrid 
        set!(f, 0.0)

        Threads.@threads for clidx in mpi_split(1 : length(SG.classes))
            f[SG.classes[clidx][1][1]] = I(to_meshes(f, SG.classes[clidx][1][1])...)
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
    MeshOperation,
    sgn,
    con,
    MeshSymmetry,
    MeshSymmetryGroup,
    get_reduced,
    init_from_reduced!,
    MeshInitFunction