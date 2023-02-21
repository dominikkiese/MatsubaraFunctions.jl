# residual modifier after application of symmetry transformation
# sgn -> change of sign | con -> complex conjugation
struct Operation 
    sgn :: Bool 
    con :: Bool

    # basic constructor 
    function Operation(
        sgn :: Bool,
        con :: Bool
        )   :: Operation 

        return new(sgn, con)
    end 

    # convenience constructor for identity 
    Operation() :: Operation = Operation(false, false)
end

# getter functions 
sgn(op :: Operation) :: Bool = op.sgn
con(op :: Operation) :: Bool = op.con

# basic operations 
function Base.:*(
    op1 :: Operation,
    op2 :: Operation
    )   :: Operation 

    return Operation(xor(sgn(op1), sgn(op2)), xor(con(op1), con(op2)))
end

function (op :: Operation)(
    x :: Q
    ) :: Q where {Q <: Number}

    if sgn(op); return con(op) ? -conj(x) : -x; end
    return con(op) ? conj(x) : x
end



# a symmetry should be defined such that it takes some 
# MatsubaraFrequency arguments and shape indices and returns a new set 
# of MatsubaraFrequency and shape indices together with an operation 
struct Symmetry{GD, SD}
    f :: FunctionWrappers.FunctionWrapper{Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}, Operation}, Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}}}
end

# make Symmetry callable
function (S :: Symmetry{GD, SD})(
    w :: NTuple{GD, MatsubaraFrequency},
    x :: NTuple{SD, Int64}
    ) :: Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}, Operation} where {GD, SD}

    return S.f(w, x)
end



# implementation of the symmetry reduction 
function reduce(
    w           :: NTuple{GD, MatsubaraFrequency},
    x           :: NTuple{SD, Int64},
    op          :: Operation,
    f           :: MatsubaraFunction{GD, SD, DD, Q},
    checked     :: Array{Bool, DD},
    symmetries  :: Vector{Symmetry{GD, SD}},
    class       :: Vector{Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}, Operation}},
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
            # convert to linear index 
            idx = LinearIndex(f, ntuple(i -> grid_index(wp[i], f.grids[i]), GD)..., xp...)

            # check if index has been used already 
            if !checked[idx]
                # this index is now checked 
                checked[idx] = true 

                # add to symmetry class and keep going 
                push!(class, (wp, xp, new_op))
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
struct SymmetryGroup{GD, SD}
    symmetries :: Vector{Symmetry{GD, SD}}
    classes    :: Vector{Vector{Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}, Operation}}}

    # basic constructor 
    function SymmetryGroup(
        symmetries :: Vector{Symmetry{GD, SD}},
        classes    :: Vector{Vector{Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}, Operation}}}
        )          :: SymmetryGroup{GD, SD} where {GD, SD}

        return new{GD, SD}(symmetries, classes)
    end 

    # convenience constructor from MatsubaraFunction and list of symmetries 
    function SymmetryGroup(
        symmetries :: Vector{Symmetry{GD, SD}},
        f          :: MatsubaraFunction{GD, SD, DD, Q}
        )          :: SymmetryGroup{GD, SD} where {GD, SD, DD, Q <: Number}

        # array to check whether index has been sorted into a symmetry class already 
        checked  = Array{Bool, DD}(undef, data_shape(f))
        checked .= false

        # list to store classes of symmetry equivalent elements
        classes = Vector{Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}, Operation}}[]

        # loop over elements of f and sort them into symmetry classes 
        for idx in eachindex(f.data)
            if !checked[idx]
                # this index is now checked and generates a new symmetry class
                checked[idx] = true
                w, x         = to_Matsubara(f, idx)
                class        = Tuple{NTuple{GD, MatsubaraFrequency}, NTuple{SD, Int64}, Operation}[(w, x, Operation())]

                # apply all symmetries to current element
                reduce(w, x, Operation(), f, checked, symmetries, class, 0)

                # add class to list 
                push!(classes, class)
            end 
        end

        return SymmetryGroup(symmetries, classes)
    end 
end

# make SymmetryGroup callable with MatsubaraFunction. This will iterate 
# over all symmetry classes and symmetrize the data array of the MatsubaraFunction
function (SG :: SymmetryGroup{GD, SD})(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: Nothing where {GD, SD, DD, Q <: Number}

    for class in SG.classes
        ref = f[class[1][1], class[1][2]...]

        for idx in 2 : length(class)
            w, x, op   = class[idx]
            f[w, x...] = op(ref)
        end 
    end 

    return nothing 
end