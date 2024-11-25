# type def and accessors
#-------------------------------------------------------------------------------#

"""
    SymmetryGroup{DD, Q <: Number}

SymmetryGroup type with fields:
* `classes    :: Vector{SymmetryClass{Q}}`
* `is_trivial :: Bool`
* `range      :: Base.OneTo{Int}`

A SymmetryGroup defines a partitioning of a MeshFunction into symmetry classes. 
Within each class, all elements are symmetry equivalent. If `is_trivial = true`, 
each class contains only one element. This allows some optimizations, such as 
relying on `range` instead of `classes` for many operations.
"""
struct SymmetryGroup{DD, Q <: Number}
    classes    :: Vector{SymmetryClass{Q}}
    is_trivial :: Bool
    range      :: Base.OneTo{Int}

    function SymmetryGroup{DD, Q}(classes :: Vector{SymmetryClass{Q}}; is_trivial :: Bool = false, range :: Base.OneTo{Int} = Base.OneTo(0)
        ) where {DD, Q <: Number}  

        return new{DD, Q}(classes, is_trivial, range)
    end 

    function SymmetryGroup(f :: MeshFunction{DD, Q, MT, AT}) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}} 
        return SymmetryGroup{DD, Q}(SymmetryClass{Q}[]; is_trivial = true, range = eachindex(f))
    end
end

"""
    function classes(SG :: SymmetryGroup{DD, Q}, n :: Int) :: SymmetryClass{Q} where {DD, Q <: Number}

Return the `n`-th class of symmetry equivalent elements.
"""
function classes(SG :: SymmetryGroup{DD, Q}, n :: Int) :: SymmetryClass{Q} where {DD, Q <: Number}
    if SG.is_trivial
        @DEBUG 0 < n <= length(SG.range) "Invalid class index!"
        return SymmetryClass(Int[n], Operation{Q}[Operation{Q}()])
    end

    return SG.classes[n]
end

""" 
    function irreducible(SG :: SymmetryGroup{DD, Q}, n :: Int) :: Int where {DD, Q <: Number}

Return the linear index of the irreducible element for the `n`-th class.
"""
function irreducible(SG :: SymmetryGroup{DD, Q}, n :: Int) :: Int where {DD, Q <: Number}
    if SG.is_trivial 
        @DEBUG 0 < n <= length(SG.range) "Invalid class index!"
        return n
    end

    return first(first(classes(SG, n)))
end

# construction from list of symmetries
#-------------------------------------------------------------------------------#

function symmetry_reduction!(
    w          :: NTuple{DD, Union{<: AbstractValue}},
    op         :: Operation{Q},
    f          :: MeshFunction{DD, Q, MT, AT},
    visited    :: Vector{Bool},
    symmetries :: Vector{<: AbstractSymmetry},
    SC         :: SymmetryClass{Q}
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}
    
    # loop over all symmetries and apply them recursively
    for S in symmetries 
        w_, op_ = S(w) # S can only be called with w if type(S) == Symmetry{DD}

        # add to class if inbounds and not visited yet
        if _all_inbounds(f, w_...)
            idx = LinearIndex(f, _mesh_indices(f, w_...)...)
            
            if !visited[idx]
                visited[idx] = true
                op__         = op_ * op
                push!(SC, idx, op__)
                symmetry_reduction!(w_, op__, f, visited, symmetries, SC)
            end
        end
    end 

    return nothing
end

function SymmetryGroup(symmetries :: Vector{<: AbstractSymmetry}, f :: MeshFunction{DD, Q, MT, AT}
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    # fallback to trivial if no symmetries provided
    if isempty(symmetries) return SymmetryGroup(f) end

    # storage (maximum memory is preallocated) and bookkeeping
    classes = SymmetryClass{Q}[]
    visited = fill(false, length(f))
    sizehint!(classes, length(f))

    # loop over elements of input MeshFunction and perform symmetry reduction
    for idx in eachindex(f)
        if !visited[idx]
            # start new class
            visited[idx] = true
            SC           = SymmetryClass(Int[idx], Operation{Q}[Operation{Q}()])

            # run symmetry reduction
            symmetry_reduction!(value.(to_meshes(f, idx)), Operation{Q}(), f, visited, symmetries, SC)
            push!(classes, SC)
        end
    end
    
    # trim down to actual size and return
    resize!(classes, length(classes))
    return SymmetryGroup{DD, Q}(classes)
end 

# construction with factorized symmetries
#-------------------------------------------------------------------------------#

# use generators to avoid allocations
function get_iters(groups :: Vararg{SymmetryGroup, NG}) where {NG}
    return Iterators.product((x.classes for x in groups)...)
end

function SymmetryGroup(
    symmetries :: Vector{Vector{<: AbstractSymmetry}}, 
    mesh_idxs  :: Vector{Vector{Int}}, 
    f          :: MeshFunction{DD, Q, MT, AT}
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    # all symmetries in list must be assigned to meshes
    @DEBUG length(symmetries) == length(mesh_idxs) "Length of symmetries and mesh indices must match!"
    @DEBUG all(!isempty(x) for x in symmetries) "Symmetries must not be empty!"
    @DEBUG all(!isempty(x) for x in mesh_idxs) "Mesh indices must not be empty!"

    # assignment of symmetries to meshes must be valid, unique and sorted
    all_idxs = vcat(mesh_idxs...)
    @DEBUG all(0 < x <= DD for x in all_idxs) "Mesh indices out of bounds!"
    @DEBUG all_idxs == unique(all_idxs) "Mesh indices must not contain duplicates!"
    @DEBUG issorted(all_idxs) "Mesh indices must be sorted!"

    # build factorized symmetry groups
    groups    = Vector{SymmetryGroup}(undef, length(symmetries))
    mesh_tpls = Tuple[tuple(map(idx -> meshes(f, idx), idxs)...) for idxs in mesh_idxs]
    
    Threads.@threads for i in eachindex(symmetries)
        groups[i] = SymmetryGroup(symmetries[i], MeshFunction(mesh_tpls[i]...; data_t = Q))
    end

    # build product classes
    iters   = collect(get_iters(groups...))
    shapes  = Tuple[length.(mesh_tpl) for mesh_tpl in mesh_tpls]
    classes = Vector{SymmetryClass{Q}}(undef, length(iters))

    Threads.@threads for i in eachindex(iters)
        ids, ops = cartesian_product(iters[i]...)
        data_id  = Vector{Int}(undef, length(ids))

        # id is tuple of linear indices for factorized meshes
        for (j, id) in enumerate(ids) 
            cidx       = CartesianIndex(map(x -> CartesianIndices(shapes[x])[id[x]], eachindex(id))...)
            data_id[j] = LinearIndex(f, cidx)
        end

        classes[i] = SymmetryClass(data_id, ops)
    end

    return SymmetryGroup{DD, Q}(classes)
end 

# make iterable
#-------------------------------------------------------------------------------#

function Base.:length(SG :: SymmetryGroup{DD, Q}) where {DD, Q <: Number}
    return SG.is_trivial ? length(SG.range) : length(SG.classes)
end

function Base.:iterate(SG :: SymmetryGroup{DD, Q}) where {DD, Q <: Number}
    return classes(SG, 1), 1
end 

function Base.:iterate(SG :: SymmetryGroup{DD, Q}, state :: Int) where {DD, Q <: Number}
    if state < length(SG)
        return classes(SG, state + 1), state + 1
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

    # symmetrize by applying operations to irreducible element of the class
    @inline function symmetrize!(SC :: SymmetryClass{Q}, f)
        idx, _ = first(SC)

        for (i, op) in SC
            f[i] = op(f[idx])
        end 
    end

    # serial mode
    if mode === :serial 
        for SC in SG
            symmetrize!(SC, f)
        end
    
    # parallel mode
    elseif (mode === :threads) || (mode === :hybrid)
        Threads.@threads for idx in 1 : length(SG)
            symmetrize!(classes(SG, idx), f)
        end

    # exception handling
    else 
        @warn "Mode $(mode) unknown. Falling back to serial execution ..."
        SG(f; mode = :serial)
    end
end

# symmetrize data array of MeshFunction given InitFunction
#-------------------------------------------------------------------------------#

function (SG :: SymmetryGroup{DD, Q})(f :: MeshFunction{DD, Q, MT, AT}, I :: InitFunction{DD, Q}; mode :: Symbol = :serial,
    ) where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    # initialize by applying InitFunction to irreducible element of the class
    @inline function init!(SC :: SymmetryClass{Q}, f, I)
        idx    = first(first(SC))
        f[idx] = I(value.(to_meshes(f, idx)))
    end

    # serial mode
    if mode === :serial 
        for SC in SG
            init!(SC, f, I)
        end 

    # threaded mode
    elseif mode === :threads
        Threads.@threads for idx in 1 : length(SG)
            init!(classes(SG, idx), f, I)
        end

    # hybrid mode: MPI + threading
    elseif mode === :hybrid 
        set!(f, Q(0))

        Threads.@threads for idx in mpi_split(1 : length(SG))
            init!(classes(SG, idx), f, I)
        end

        mpi_allreduce!(f)

    # exception handling
    else 
        @warn "Mode $(mode) unknown. Falling back to serial execution ..."
        SG(f, I; mode = :serial)
    end

    SG(f; mode)
    return nothing 
end

# flattening / unflattening of MeshFunction including symmetries
#-------------------------------------------------------------------------------#

"""
    function flatten!(SG :: SymmetryGroup{DD, Q}, f :: MeshFunction{DD, Q, MT, AT}, fvec :: AbstractVector{Q}
        ) :: Nothing where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

Calculate and store symmetry reduced representation of MeshFunction in `fvec`
"""
function flatten!(SG :: SymmetryGroup{DD, Q}, f :: MeshFunction{DD, Q, MT, AT}, fvec :: AbstractVector{Q}
    ) :: Nothing where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    @DEBUG length(fvec) == length(SG) "Length of fvec does not match the number of classes in the symmetry group!"

    for n in eachindex(fvec)
        fvec[n] = f[irreducible(SG, n)]
    end

    return nothing 
end

"""
    function flatten(SG :: SymmetryGroup{DD, Q}, f :: MeshFunction{DD, Q, MT, AT}
        ) :: Vector{Q} where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

Calculate and return symmetry reduced representation of MeshFunction
"""
function flatten(SG :: SymmetryGroup{DD, Q}, f :: MeshFunction{DD, Q, MT, AT}
    ) :: Vector{Q} where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    fvec = Vector{Q}(undef, length(SG))
    flatten!(SG, f, fvec)
    return fvec
end

"""
    function unflatten!(
        SG   :: SymmetryGroup{DD, Q},
        f    :: MeshFunction{DD, Q, MT, AT},
        fvec :: AbstractVector{Q}
        )    :: Nothing where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

Initialize MeshFunction from symmetry reduced representation
""" 
function unflatten!(
    SG   :: SymmetryGroup{DD, Q},
    f    :: MeshFunction{DD, Q, MT, AT},
    fvec :: AbstractVector{Q}
    ;
    mode :: Symbol = :serial
    )    :: Nothing where {DD, Q <: Number, MT <: NTuple{DD, Mesh}, AT <: AbstractArray{Q, DD}}

    @DEBUG length(fvec) == length(SG) "Length of fvec does not match the number of classes in the symmetry group!"

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
    classes,
    irreducible,
    flatten!,
    flatten,
    unflatten!