# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct SymmetryClass 

SymmetryClass type with fields:
* `data_id :: Vector{Int}`
* `data_op :: Vector{Operation}`

A symmetry class contains a list of indices and operations. The first index is associated 
with the irreducible element of a MeshFunction and must have the identity as its operation.
"""
struct SymmetryClass 
    data_id :: Vector{Int}
    data_op :: Vector{Operation}

    function SymmetryClass(data_id :: Vector{Int}, data_op :: Vector{Operation})
        @DEBUG sgn(first(data_op)) == false "First operation must be the identity!"
        @DEBUG con(first(data_op)) == false "First operation must be the identity!"
        return new(data_id, data_op)
    end

    function SymmetryClass()
        return new(Int[], Operation[])
    end
end

function Base.:length(SC :: SymmetryClass)
    return length(SC.data_id)
end

function Base.:getindex(SC :: SymmetryClass, idx :: Int)
    return SC.data_id[idx], SC.data_op[idx]
end

function Base.:setindex!(SC :: SymmetryClass, val :: Tuple{Int, Operation}, idx :: Int)
    SC.data_id[idx], SC.data_op[idx] = val
    return nothing
end

function Base.:eachindex(SC :: SymmetryClass)
    return eachindex(SC.data_id)
end

function Base.:first(SC :: SymmetryClass)
    return SC[1]
end

function Base.:popat!(SC :: SymmetryClass, idx :: Int)
    return popat!(SC.data_id, idx), popat!(SC.data_op, idx)
end

function Base.:push!(SC :: SymmetryClass, idx :: Int, op :: Operation)
    push!(SC.data_id, idx)
    push!(SC.data_op, op)
    return SC
end

function Base.:pushfirst!(SC :: SymmetryClass, idx :: Int, op :: Operation)
    pushfirst!(SC.data_id, idx)
    pushfirst!(SC.data_op, op)
    return SC
end

# make SymmetryClass iterable
#-------------------------------------------------------------------------------#

function Base.:iterate(SC :: SymmetryClass)
    return first(SC), 1
end

function Base.:iterate(SC :: SymmetryClass, state :: Int)
    if state < length(SC)
        return SC[state + 1], state + 1
    else 
        return nothing
    end
end

# change irreducible element
#-------------------------------------------------------------------------------#

"""
    function set_irreducible!(SC :: SymmetryClass, idx :: Int) :: Nothing

Set the `idx`-th element as the irreducible element of the symmetry class.
"""
function set_irreducible!(SC :: SymmetryClass, idx :: Int) :: Nothing
    @DEBUG 0 < idx <= length(SC) "Index out of bounds!"

    # move to first position
    pushfirst!(SC, popat!(SC, idx)...)
    _, inv = first(SC)

    # ensure first operation is identity using that operations are self-inverse
    for i in eachindex(SC)
       j, op = SC[i]
       SC[i] = j, op * inv
    end

    return nothing
end

# intersect and merge
#-------------------------------------------------------------------------------#

function Base.:intersect(SC1 :: SymmetryClass, SC2 :: SymmetryClass)
    return intersect(SC1.data_id, SC2.data_id)
end

"""
    function do_intersect(SC1 :: SymmetryClass, SC2 :: SymmetryClass) :: Bool

Check if two symmetry classes intersect.
"""
function do_intersect(SC1 :: SymmetryClass, SC2 :: SymmetryClass) :: Bool
    return any(x -> x in SC1.data_id, SC2.data_id)
end

"""
    function fuse(SC1 :: SymmetryClass, SC2 :: SymmetryClass) :: SymmetryClass

Fuse two symmetry classes into a new one. The irreducible element is chosen as the first
common element between the two classes.
"""
function fuse(SC1 :: SymmetryClass, SC2 :: SymmetryClass) :: SymmetryClass
    @DEBUG do_intersect(SC1, SC2) "Symmetry classes do not intersect!"

    # first element chosen as common irreducible element
    I   = intersection(SC1, SC2)
    idx = first(I)
    
    # redefine irreducible element to idx (valid since choice of irreducible element is arbitrary)
    set_irreducible!(SC1, idx)
    set_irreducible!(SC2, idx)

    # merge data, ensure no duplicates
    fused_data_id = copy(SC1.data_id)
    fused_data_op = copy(SC1.data_op)

    for (i, op) in SC2
        if i in I continue end
        push!(fused_data_id, i)
        push!(fused_data_op, op)
    end

    @DEBUG merged_data_id == unique(merged_data_id) "Merged data has duplicate entries!"
    return SymmetryClass(merged_data_id, merged_data_op)
end

# export
#-------------------------------------------------------------------------------#

export 
    SymmetryClass,
    set_irreducible!,
    do_intersect,
    fuse