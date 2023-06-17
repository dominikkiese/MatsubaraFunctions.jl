# residual modifier after application of symmetry transformation
# sgn -> change of sign | con -> complex conjugation
"""
    struct MatsubaraOperation 

MatsubaraOperation type with fields:
* `sgn :: Bool` : change sign?
* `con :: Bool` : do complex conjugation?
"""
struct MatsubaraOperation 
    sgn :: Bool 
    con :: Bool

    # default constructor 
    function MatsubaraOperation(
        sgn :: Bool,
        con :: Bool
        )   :: MatsubaraOperation 

        return new(sgn, con)
    end 

    # convenience constructor for identity 
    MatsubaraOperation() :: MatsubaraOperation = MatsubaraOperation(false, false)
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

# basic operations 
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



# a symmetry should be defined such that it takes some 
# MatsubaraFrequency arguments and tensor indices and returns a new set 
# of MatsubaraFrequency and tensor indices together with a MatsubaraOperation 
"""
    struct MatsubaraSymmetry{GD, SD}

MatsubaraSymmetry type with fields:
* `f :: FunctionWrappers.FunctionWrapper{Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}, MatsubaraOperation}, Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}}}`
`MatsubaraSymmetry` takes grid coordinates and tensor indices as input and returns a new set of coordinates and indices together with a MatsubaraOperation
"""
struct MatsubaraSymmetry{GD, SD}
    f :: FunctionWrappers.FunctionWrapper{Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}, MatsubaraOperation}, Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}}}
end

# make MatsubaraSymmetry callable
function (S :: MatsubaraSymmetry{GD, SD})(
    w :: NTuple{GD, MatsubaraFrequency},
    x :: NTuple{SD, Int64}
    ) :: Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}, MatsubaraOperation} where {GD, SD}

    return S.f(w, x)
end



# implementation of the symmetry reduction 
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

    # loop over all symmetries 
    for S in symmetries 
        # apply the symmetry operation
        wp, xp, opp = S(w, x)

        # compute new operation
        new_op = opp * op

        # check if index is valid, reset path_length iff we start from valid & unchecked index
        if !any(ntuple(i -> !is_inbounds(wp[i], f.grids[i]), GD))
            # convert to linear index (grid_index method to avoid duplicate inbounds check)
            idx = LinearIndex(f, ntuple(i -> grid_index(wp[i], f.grids[i]), GD)..., xp...)

            # check if index has been used already 
            if !checked[idx]
                # this index is now checked 
                checked[idx] = true 

                # add to symmetry class and keep going 
                push!(class, (idx, new_op))
                reduce(wp, xp, new_op, f, checked, symmetries, class, 0; max_length)
            end 

        # if index is invalid, increment path length and keep going 
        elseif path_length < max_length
            reduce(wp, xp, new_op, f, checked, symmetries, class, path_length + 1; max_length)
        end 
    end 
end

# here, a symmetry group is a list of collections of symmetry equivalent elements
"""
    MatsubaraSymmetryGroup

MatsubaraSymmetryGroup type with fields:
* `classes :: Vector{Vector{Tuple{Int64, MatsubaraOperation}}}` : list of collections of symmetry equivalent elements
* `speedup :: Float64`                                          : expected speedup from the symmetry reduction
"""
struct MatsubaraSymmetryGroup
    classes :: Vector{Vector{Tuple{Int64, MatsubaraOperation}}}
    speedup :: Float64

    # default constructor 
    function MatsubaraSymmetryGroup(
        classes :: Vector{Vector{Tuple{Int64, MatsubaraOperation}}},
        speedup :: Float64
        )       :: MatsubaraSymmetryGroup

        return new(classes, speedup)
    end 

    # identity constructor 
    function MatsubaraSymmetryGroup(
        f :: MatsubaraFunction{GD, SD, DD, Q}
        ) :: MatsubaraSymmetryGroup where {GD, SD, DD, Q <: Number}

        classes = Vector{Tuple{Int64, MatsubaraOperation}}[[(idx, MatsubaraOperation())] for idx in eachindex(f.data)]
        return new(classes, 1.0)
    end

    # convenience constructor from MatsubaraFunction and list of symmetries 
    function MatsubaraSymmetryGroup(
        symmetries :: Vector{MatsubaraSymmetry{GD, SD}},
        f          :: MatsubaraFunction{GD, SD, DD, Q}
        ;
        max_length :: Int64 = 0
        )          :: MatsubaraSymmetryGroup where {GD, SD, DD, Q <: Number}

        # array to check whether index has been sorted into a symmetry class already 
        checked  = Array{Bool, DD}(undef, data_shape(f))
        checked .= false

        # list to store classes of symmetry equivalent elements
        classes = Vector{Tuple{Int64, MatsubaraOperation}}[]

        # loop over elements of f and sort them into symmetry classes 
        for idx in eachindex(f.data)
            if !checked[idx]
                # this index is now checked and generates a new symmetry class
                checked[idx] = true
                w, x         = to_Matsubara(f, idx)
                class        = Tuple{Int64, MatsubaraOperation}[(idx, MatsubaraOperation())]

                # apply all symmetries to current element
                reduce(w, x, MatsubaraOperation(), f, checked, symmetries, class, 0; max_length)

                # add class to list 
                push!(classes, class)
            end 
        end

        return MatsubaraSymmetryGroup(classes, length(f.data) / length(classes))
    end 
end

# make MatsubaraSymmetryGroup callable with MatsubaraFunction. This will iterate 
# over all symmetry classes and symmetrize the data array of the MatsubaraFunction
function (SG :: MatsubaraSymmetryGroup)(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: Nothing where {GD, SD, DD, Q <: Number}

    for class in SG.classes
        ref = f[class[1][1]]

        for idx in 2 : length(class)
            idx, op = class[idx]
            f[idx]  = op(ref)
        end 
    end 

    return nothing 
end



# an init function should be defined such that it takes some 
# MatsubaraFrequency arguments and shape indices and returns a value of type Q
"""
    struct MatsubaraInitFunction{GD, SD, Q <: Number}

MatsubaraInitFunction type with fields:
* `f :: FunctionWrappers.FunctionWrapper{Q, Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}}}` 
`MatsubaraInitFunction` takes grid coordinates and tensor indices as input and returns value of type Q
"""
struct MatsubaraInitFunction{GD, SD, Q <: Number}
    f :: FunctionWrappers.FunctionWrapper{Q, Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}}}
end

# make MatsubaraInitFunction callable
function (I :: MatsubaraInitFunction{GD, SD, Q})(
    w :: NTuple{GD, MatsubaraFrequency},
    x :: NTuple{SD, Int64}
    ) :: Q where {GD, SD, Q <: Number}

    return I.f(w, x)
end

# make MatsubaraSymmetryGroup callable with MatsubaraFunction and MatsubaraInitFunction. This will iterate 
# over all symmetry classes and symmetrize the data array of the MatsubaraFunction starting from an 
# evaluation of the MatsubaraInitFunction
function (SG :: MatsubaraSymmetryGroup)(
    f            :: MatsubaraFunction{GD, SD, DD, Q},
    I            :: MatsubaraInitFunction{GD, SD, Q}
    ;
    mpi_parallel :: Bool = false
    )            :: Nothing where {GD, SD, DD, Q <: Number}

    if mpi_parallel 
        set!(f, 0.0)

        Threads.@threads for clidx in mpi_split(1 : length(SG.classes))
            lidx    = SG.classes[clidx][1][1]
            f[lidx] = I(to_Matsubara(f, lidx)...)
        end

        mpi_allreduce!(f)
    else
        for class in SG.classes
            lidx    = class[1][1]
            f[lidx] = I(to_Matsubara(f, lidx)...)
        end 
    end

    SG(f)

    return nothing 
end