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

# a symmetry takes MeshPoint value types and returns MeshPoint value types and an operation
"""
    struct Symmetry{DD}

Symmetry type with fields:
* `f :: Function`
"""
struct Symmetry{DD}
    f :: Function
end

# explicit return type to enforce proper implementation
function (S :: Symmetry{DD})(w :: NTuple{DD, Union{<: AbstractValue}}
    ) :: Tuple{NTuple{DD, Union{<: AbstractValue}}, Operation} where {DD}

    return S.f(w)
end

function reduce(
    w           :: NTuple{DD, Union{<: AbstractValue}},
    op          :: Operation,
    f           :: MeshFunction{DD, Q, MT, AT},
    checked     :: Array{Bool, DD},
    symmetries  :: Vector{Symmetry{DD}},
    data_id     :: Vector{Int},
    data_op     :: Vector{Operation},
    path_length :: Int
    ;
    max_length  :: Int = 0
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}
 
    for S in symmetries 
        wp, opp = S(w)
        new_op  = opp * op

        if _all_inbounds(f, wp...)
            idx = LinearIndex(f, _mesh_indices(f, wp...)...)

            if !checked[idx]
                checked[idx] = true 

                # add to symmetry class, reset path length and keep going 
                push!(data_id, idx)
                push!(data_op, new_op)
                reduce(wp, new_op, f, checked, symmetries, data_id, data_op, 0; max_length)
            end 

        # if index not valid, increment path length and keep going 
        elseif path_length < max_length
            reduce(wp, new_op, f, checked, symmetries, data_id, data_op, path_length + 1; max_length)
        end 
    end 
end

# symmetry group type
#-------------------------------------------------------------------------------#

"""
    SymmetryGroup{DD, Q <: Number}

SymmetryGroup type with fields:
* `data_id    :: Union{Vector{Int}, Base.OneTo{Int}}` : linear MeshFunction indices
* `data_op    :: Vector{Operation}`                   : operations on MeshFunction elements in `data_idx`
* `classes    :: Vector{Int}`                         : indices of irreducible elements in `data_idx`
* `is_trivial :: Bool`                                : is the symmetry group trivial?

If `is_trivial` is `true`, the symmetry group only contains the identity. This allows some optimizations,
such as setting the type of `data_id` from `Vector{Int}` to `Base.OneTo{Int}`.
"""
struct SymmetryGroup{DD, Q <: Number}
    data_id    :: Union{Vector{Int}, Base.OneTo{Int}}
    data_op    :: Vector{Operation}
    classes    :: Vector{Int}
    is_trivial :: Bool

    function SymmetryGroup{DD, Q}(
        data_id    :: Union{Vector{Int}, Base.OneTo{Int}},
        data_op    :: Vector{Operation},
        classes    :: Vector{Int}
        ;
        is_trivial :: Bool = false
        ) where {DD, Q <: Number}  

        return new{DD, Q}(data_id, data_op, classes, is_trivial)
    end 

    # initialize trivial symmetry group from MeshFunction
    function SymmetryGroup(f :: MeshFunction{DD, Q, MT, AT}) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}} 
        return SymmetryGroup{DD, Q}(eachindex(f), Operation[], Int[]; is_trivial = true)
    end
 
    function SymmetryGroup(
        symmetries :: Vector{Symmetry{DD}},
        f          :: MeshFunction{DD, Q, MT, AT}
        ;
        max_length :: Int = 0
        ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

        # fallback to trivial symmetry group if no symmetries are provided
        if isempty(symmetries) return SymmetryGroup(f) end

        # we require dynamic information about length!
        data_id = Int[] 
        data_op = Operation[]
        classes = Int[]
        checked = fill(false, size(f.data))

        for idx in eachindex(f)
            if !checked[idx]
                # start new class
                push!(data_id, idx)
                push!(data_op, Operation())
                push!(classes, length(data_id))
                checked[idx] = true

                # run symmetry reduction
                reduce(value.(to_meshes(f, idx)), Operation(), f, checked, symmetries, data_id, data_op, 0; max_length)
            end 
        end

        # fallback to trivial symmetry group if no reduction was possible
        if length(classes) == length(f)
            @warn "No reduction was possible. Treating symmetry group as trivial."
            return SymmetryGroup(f)
        end

        # if successful, finalize last class and return
        push!(classes, length(data_id) + 1)
        return SymmetryGroup{DD, Q}(data_id, data_op, classes)
    end 
end

"""
    function class(SG :: SymmetryGroup{DD, Q}, n :: Int) :: Union{Tuple{Int, Operation}, Vector{Tuple{Int, Operation}}} where {DD, Q <: Number}

Return the `n`-th class of symmetry equivalent elements.
"""
function class(SG :: SymmetryGroup{DD, Q}, n :: Int) :: Union{Tuple{Int, Operation}, Vector{Tuple{Int, Operation}}} where {DD, Q <: Number}
    if SG.is_trivial
        return SG.data_id[n], Operation() # still use SG.data_id to check if n is valid
    end

    return Tuple{Int, Operation}[(SG.data_id[i], SG.data_op[i]) for i in SG.classes[n] : SG.classes[n + 1] - 1]
end

""" 
    function irreducible(SG :: SymmetryGroup{DD, Q}, n :: Int) :: Int where {DD, Q <: Number}

Return linear index of the irreducible (~reference) element for the `n`-th class.
"""
function irreducible(SG :: SymmetryGroup{DD, Q}, n :: Int) :: Int where {DD, Q <: Number}
    if SG.is_trivial 
        return SG.data_id[n] # still use SG.data_id to check if n is valid
    end

    return SG.data_id[SG.classes[n]]
end

# iteration over SymmetryGroup is iteration over classes
function Base.:length(SG :: SymmetryGroup{DD, Q}) where {DD, Q <: Number}
    return SG.is_trivial ? length(SG.data_id) : length(SG.classes) - 1
end

function Base.:iterate(SG :: SymmetryGroup{DD, Q}) where {DD, Q <: Number}
    return class(SG, 1), 1
end 

function Base.:iterate(SG :: SymmetryGroup{DD, Q}, state :: Int) where {DD, Q <: Number}
    if state < length(SG)
        return class(SG, state + 1), state + 1
    else 
        return nothing
    end
end

# symmetrize data array of the MeshFunction
function (SG :: SymmetryGroup{DD, Q})(f :: MeshFunction{DD, Q, MT, AT}; mode :: Symbol = :serial
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    # if SymmetryGroup is trivial there is nothing to do
    if SG.is_trivial return nothing end

    @inline function symmetrize!(class :: Vector{Tuple{Int, Operation}}, f)
        for tpl in class
            f[first(tpl)] = last(tpl)(f[first(first(class))])
        end 
    end

    if mode === :serial 
        for class in SG
            symmetrize!(class, f)
        end

    elseif (mode === :threads) || (mode === :hybrid)
        Threads.@threads for clidx in 1 : length(SG)
            symmetrize!(class(SG, clidx), f)
        end

    else 
        error("Mode $(mode) unknown. Possible options are: serial, polyester, threads and hybrid (here: hybrid = threads)")
    end
end

"""
    function get_reduced(SG :: SymmetryGroup{DD, Q}, f :: MeshFunction{DD, Q, MT, AT}
        ) :: Vector{Q} where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

Calculate symmetry reduced representation of MeshFunction
"""
function get_reduced(SG :: SymmetryGroup{DD, Q}, f :: MeshFunction{DD, Q, MT, AT}
    ) :: Vector{Q} where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    return Q[f[irreducible(SG, n)] for n in 1 : length(SG)]
end

"""
    function init_from_reduced!(
        SG   :: SymmetryGroup{DD, Q},
        f    :: MeshFunction{DD, Q, MT, AT},
        fvec :: AbstractVector{Q}
        )    :: Nothing where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

Initialize MeshFunction from symmetry reduced representation
""" 
function init_from_reduced!(
    SG   :: SymmetryGroup{DD, Q},
    f    :: MeshFunction{DD, Q, MT, AT},
    fvec :: AbstractVector{Q}
    )    :: Nothing where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    for cl in eachindex(fvec)
        f[irreducible(SG, cl)] = fvec[cl]
    end 

    SG(f)
    return nothing 
end

# init function type
#-------------------------------------------------------------------------------#

"""
    struct InitFunction{DD, Q <: Number}

InitFunction type with fields:
* `f :: Function` 
"""
struct InitFunction{DD, Q <: Number}
    f :: Function
end

# explicit return type to enforce proper implementation
function (I :: InitFunction{DD, Q})(w :: NTuple{DD, Union{<: AbstractValue}}) :: Q where {DD, Q <: Number}
    return I.f(w)
end

# symmetrize MeshFunction from evaluation of InitFunction
function (SG :: SymmetryGroup{DD, Q})(
    f    :: MeshFunction{DD, Q, MT, AT},
    I    :: InitFunction{DD, Q}
    ;
    mode :: Symbol = :serial,
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    @inline function init!(class, f, I)
        idx    = first(first(class))
        f[idx] = I(value.(to_meshes(f, idx)))
    end

    if mode === :serial 
        for class in SG
            init!(class, f, I)
        end 

    elseif mode === :threads
        Threads.@threads for clidx in 1 : length(SG)
            init!(class(SG, clidx), f, I)
        end

    elseif mode === :hybrid 
        set!(f, Q(0))

        Threads.@threads for clidx in mpi_split(1 : length(SG))
            init!(class(SG, clidx), f, I)
        end

        mpi_allreduce!(f)

    else 
        error("Mode $(mode) unknown. Possible options are: serial, polyester, threads and hybrid")
    end

    SG(f; mode)
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
    class,
    irreducible,
    get_reduced,
    init_from_reduced!,
    InitFunction