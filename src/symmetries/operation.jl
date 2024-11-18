# type def and accessors
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

# multiplication
#-------------------------------------------------------------------------------#

function Base.:*(op1 :: Operation, op2 :: Operation) 
    return Operation(xor(sgn(op1), sgn(op2)), xor(con(op1), con(op2)))
end

# call to Operation
#-------------------------------------------------------------------------------#

function (op :: Operation)(x :: Q) where {Q <: Number}
    if sgn(op); return con(op) ? -conj(x) : -x; end
    return con(op) ? conj(x) : x
end

# export
#-------------------------------------------------------------------------------#

export 
    Operation, 
    sgn,
    con