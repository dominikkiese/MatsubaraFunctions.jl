# type def and accessors
#-------------------------------------------------------------------------------#

"""
    SymmetryGroup{DD, Q <: Number}

SymmetryGroup type with fields:
* `classes    :: Vector{SymmetryClass}`
* `is_trivial :: Bool`
* `range      :: Base.OneTo{Int}`

A SymmetryGroup defines a partitioning of a MeshFunction into symmetry classes. 
Each class contains a list of indices and operations. If `is_trivial` is `true`, 
the symmetry group only contains the identity. This allows some optimizations,
such as relying on `range` instead of `classes` for many operations.
"""
struct SymmetryGroup{DD, Q <: Number}
    classes    :: Vector{SymmetryClass}
    is_trivial :: Bool
    range      :: Base.OneTo{Int}

    function SymmetryGroup{DD, Q}(classes :: Vector{SymmetryClass}; is_trivial :: Bool = false, range :: Base.OneTo{Int} = Base.OneTo(0)
        ) where {DD, Q <: Number}  

        return new{DD, Q}(classes, is_trivial, range)
    end 

    # initialize trivial symmetry group from MeshFunction
    function SymmetryGroup(f :: MeshFunction{DD, Q, MT, AT}) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}} 
        return SymmetryGroup{DD, Q}(SymmetryClass[]; is_trivial = true, range = eachindex(f))
    end
end

"""
    function class(SG :: SymmetryGroup{DD, Q}, n :: Int) :: Union{Tuple{Int, Operation}, Vector{Tuple{Int, Operation}}} where {DD, Q <: Number}

Return the `n`-th class of symmetry equivalent elements.
"""
function class(SG :: SymmetryGroup{DD, Q}, n :: Int) :: SymmetryClass where {DD, Q <: Number}
    if SG.is_trivial
        @DEBUG 0 < n <= length(SG.range) "Invalid class index!"
        return SymmetryClass(Int[n], Operation[Operation()])
    end

    return SG.classes[n]
end

""" 
    function irreducible(SG :: SymmetryGroup{DD, Q}, n :: Int) :: Int where {DD, Q <: Number}

Return linear index of the irreducible element for the `n`-th class.
"""
function irreducible(SG :: SymmetryGroup{DD, Q}, n :: Int) :: Int where {DD, Q <: Number}
    if SG.is_trivial 
        @DEBUG 0 < n <= length(SG.range) "Invalid class index!"
        return n
    end

    # inner first returns tuple of (index, operation)
    return first(first(SG.classes[n]))
end

# construction from Vector{Symmetry}
#-------------------------------------------------------------------------------#

# FIXME: we need something scalable with good performance here!!!

# make iterable
#-------------------------------------------------------------------------------#

function Base.:length(SG :: SymmetryGroup{DD, Q}) where {DD, Q <: Number}
    return SG.is_trivial ? length(SG.range) : length(SG.classes)
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

# symmetrize data array of MeshFunction
#-------------------------------------------------------------------------------#

function (SG :: SymmetryGroup{DD, Q})(f :: MeshFunction{DD, Q, MT, AT}; mode :: Symbol = :serial
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    # if SymmetryGroup is trivial there is nothing to do
    if SG.is_trivial return nothing end

    @inline function symmetrize!(SC :: SymmetryClass, f)
        idx, _ = first(SC)

        for (i, op) in SC
            f[i] = op(f[idx])
        end 
    end

    if mode === :serial 
        for SC in SG
            symmetrize!(SC, f)
        end

    elseif (mode === :threads) || (mode === :hybrid)
        Threads.@threads for clidx in 1 : length(SG)
            symmetrize!(class(SG, clidx), f)
        end

    else 
        error("Mode $(mode) unknown. Possible options are: serial, polyester, threads and hybrid (here: hybrid = threads)")
    end
end

# flattening / unflattening of MeshFunction including symmetries
#-------------------------------------------------------------------------------#

"""
    function reduce!(SG :: SymmetryGroup{DD, Q}, f :: MeshFunction{DD, Q, MT, AT}, fvec :: AbstractVector{Q}
        ) :: Nothing where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

Calculate and store symmetry reduced representation of MeshFunction in `fvec`
"""
function reduce!(SG :: SymmetryGroup{DD, Q}, f :: MeshFunction{DD, Q, MT, AT}, fvec :: AbstractVector{Q}
    ) :: Nothing where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    @DEBUG length(fvec) == length(SG) "Length of fvec does not match the number of classes in the symmetry group!"
    for n in 1 : length(SG)
        fvec[n] = f[irreducible(SG, n)]
    end

    return nothing 
end

"""
    function reduce(SG :: SymmetryGroup{DD, Q}, f :: MeshFunction{DD, Q, MT, AT}
        ) :: Vector{Q} where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

Calculate and return symmetry reduced representation of MeshFunction
"""
function reduce(SG :: SymmetryGroup{DD, Q}, f :: MeshFunction{DD, Q, MT, AT}
    ) :: Vector{Q} where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    fvec = Vector{Q}(undef, length(SG))
    reduce!(SG, f, fvec)
    return fvec
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
    ;
    mode :: Symbol = :serial
    )    :: Nothing where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    for n in eachindex(fvec)
        f[irreducible(SG, n)] = fvec[n]
    end 

    SG(f; mode)
    return nothing 
end

# export
#-------------------------------------------------------------------------------#

export 
    SymmetryGroup,
    class,
    irreducible,
    reduce!,
    reduce,
    init_from_reduced!