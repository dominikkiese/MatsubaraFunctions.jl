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

    # basic constructor 
    function MatsubaraOperation(
        sgn :: Bool,
        con :: Bool
        )   :: MatsubaraOperation 

        return new(sgn, con)
    end 

    # convenience constructor for identity 
    MatsubaraOperation() :: MatsubaraOperation = MatsubaraOperation(false, false)
end

# getter functions 
"""
    sgn(op :: MatsubaraOperation) :: Bool

Return op.sgn
"""
sgn(op :: MatsubaraOperation) :: Bool = op.sgn

"""
    con(op :: MatsubaraOperation) :: Bool

Return op.con
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
MatsubaraSymmetry takes grid coordinates and tensor indices as input and returns a new set of coordinates and indices together with a MatsubaraOperation
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
                reduce(wp, xp, new_op, f, checked, symmetries, class, 0)
            end 

        # if index is invalid, increment path length and keep going 
        elseif path_length < 10 
            reduce(wp, xp, new_op, f, checked, symmetries, class, path_length + 1)
        end 
    end 
end

# a symmetry group is composed of a list of defining symmetries as well as 
# a list containing collections of symmetry equivalent elements
"""
    MatsubaraSymmetryGroup{GD, SD}

MatsubaraSymmetryGroup type with fields:
* `symmetries :: Vector{MatsubaraSymmetry{GD, SD}}`                : list of symmetry operations
* `classes    :: Vector{Vector{Tuple{Int64, MatsubaraOperation}}}` : list of symmetry classes

Examples:
```julia
# a simple Green's function
ξ = 0.5
T = 1.0
N = 128
g = MatsubaraGrid(T, N, Fermion)
f = MatsubaraFunction(g, 1)

for v in g
    f[v, 1] = 1.0 / (im * value(v) - ξ)
end 

# complex conjugation acting on Green's function
function conj(
    w :: Tuple{MatsubaraFrequency},
    x :: Tuple{Int64}
    ) :: Tuple{Tuple{MatsubaraFrequency}, Tuple{Int64}, MatsubaraOperation}

    return (-w[1],), (x[1],), MatsubaraOperation(false, true)
end 

# compute the symmetry group 
SG = MatsubaraSymmetryGroup([MatsubaraSymmetry{1, 1}(conj)], f)

# obtain another Green's function by symmetrization
function init(
    w :: Tuple{MatsubaraFrequency},
    x :: Tuple{Int64}
    ) :: ComplexF64

    return f[w, x...]
end 

InitFunc = MatsubaraInitFunction{1, 1, ComplexF64}(init)
h = MatsubaraFunction(g, 1)
SG(h, InitFunc)
@assert h.data ≈ f.data
```
"""
struct MatsubaraSymmetryGroup{GD, SD}
    symmetries :: Vector{MatsubaraSymmetry{GD, SD}}
    classes    :: Vector{Vector{Tuple{Int64, MatsubaraOperation}}}

    # basic constructor 
    function MatsubaraSymmetryGroup(
        symmetries :: Vector{MatsubaraSymmetry{GD, SD}},
        classes    :: Vector{Vector{Tuple{Int64, MatsubaraOperation}}}
        )          :: MatsubaraSymmetryGroup{GD, SD} where {GD, SD}

        return new{GD, SD}(symmetries, classes)
    end 

    # convenience constructor from MatsubaraFunction and list of symmetries 
    function MatsubaraSymmetryGroup(
        symmetries :: Vector{MatsubaraSymmetry{GD, SD}},
        f          :: MatsubaraFunction{GD, SD, DD, Q}
        )          :: MatsubaraSymmetryGroup{GD, SD} where {GD, SD, DD, Q <: Number}

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
                reduce(w, x, MatsubaraOperation(), f, checked, symmetries, class, 0)

                # add class to list 
                push!(classes, class)
            end 
        end

        return MatsubaraSymmetryGroup(symmetries, classes)
    end 
end

# make MatsubaraSymmetryGroup callable with MatsubaraFunction. This will iterate 
# over all symmetry classes and symmetrize the data array of the MatsubaraFunction
function (SG :: MatsubaraSymmetryGroup{GD, SD})(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: Nothing where {GD, SD, DD, Q <: Number}

    Threads.@threads for class in SG.classes
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
MatsubaraInitFunction takes grid coordinates and tensor indices as input and returns value of type Q
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
function (SG :: MatsubaraSymmetryGroup{GD, SD})(
    f            :: MatsubaraFunction{GD, SD, DD, Q},
    I            :: MatsubaraInitFunction{GD, SD, Q}
    ;
    mpi_parallel :: Bool = false
    )            :: Nothing where {GD, SD, DD, Q <: Number}

    if mpi_parallel 
        for clidx in 1 : mpi_split(1 : length(SG.classes))
            w, x                       = to_Matsubara(f, SG.classes[clidx][1][1])
            ref                        = I(w, x)
            f[SG.classes[clidx][1][1]] = ref
    
            Threads.@threads for idx in 2 : length(SG.classes[clidx])
                idx, op = SG.classes[clidx][idx]
                f[idx]  = op(ref)
            end 
        end 
    else
        Threads.@threads for class in SG.classes
            w, x           = to_Matsubara(f, class[1][1])
            ref            = I(w, x)
            f[class[1][1]] = ref

            for idx in 2 : length(class)
                idx, op = class[idx]
                f[idx]  = op(ref)
            end 
        end 
    end

    return nothing 
end