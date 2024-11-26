# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct SymmetryClass{Q <: Number} 

SymmetryClass type with fields:
* `data_id :: Vector{Int}`
* `data_op :: Vector{Operation{Q}}`

A SymmetryClass contains a list of indices and operations. The first index is associated 
with the irreducible (~ reference) element of that class and must have the identity 
as its operation.
"""
struct SymmetryClass{Q <: Number}
    data_id :: Vector{Int}
    data_op :: Vector{Operation{Q}}

    function SymmetryClass(data_id :: Vector{Int}, data_op :: Vector{Operation{Q}}) where {Q <: Number}
        @DEBUG !sgn(first(data_op)) && !con(first(data_op)) "First operation must be the identity!"
        return new{Q}(data_id, data_op)
    end

    function SymmetryClass{Q}() where {Q <: Number}
        return new{Q}(Int[], Operation{Q}[])
    end
end

function Base.:length(SC :: SymmetryClass)
    @DEBUG length(SC.data_id) == length(SC.data_op) "Data length mismatch!"
    return length(SC.data_id)
end

function Base.:getindex(SC :: SymmetryClass, idx :: Int)
    return SC.data_id[idx], SC.data_op[idx]
end

function Base.:first(SC :: SymmetryClass)
    return SC[1]
end

function Base.:push!(SC :: SymmetryClass{Q}, idx :: Int, op :: Operation{Q}) where {Q <: Number}
    push!(SC.data_id, idx)
    push!(SC.data_op, op)
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

# cartesian product of SymmetryClass objects 
#-------------------------------------------------------------------------------#

# to avoid ambiguities
function get_product_iters()
    return Iterators.product(), Iterators.product()
end

# use generators to avoid allocations
function get_product_iters(classes :: Vararg{SymmetryClass{Q}, NC}) where {Q <: Number, NC}
    return Iterators.product((x.data_id for x in classes)...), Iterators.product((x.data_op for x in classes)...)
end

# to avoid unbound type parameters
function cartesian_product()
    ids, ops = vec.(collect.(get_product_iters()))
    return ids, ops
end

function cartesian_product(classes :: Vararg{SymmetryClass{Q}, NC}) where {Q <: Number, NC}
    ids, ops = vec.(collect.(get_product_iters(classes...)))
    return ids, prod.(ops)
end

# export
#-------------------------------------------------------------------------------#

export 
    SymmetryClass