# type def and accessors
#-------------------------------------------------------------------------------#

"""
    struct Operation{Q <: Number} 

Operation type with fields:
* `sgn :: Bool` : change sign?
* `con :: Bool` : complex conjugation?

If `sgn` is `true`, the operation changes the sign of the input number. 
If `con` is `true`, the operation applies complex conjugation to the input number.
"""
struct Operation{Q <: Number} 
    sgn :: Bool 
    con :: Bool

    function Operation{Q}(; sgn :: Bool = false, con :: Bool = false) where {Q <: Number}
        return new{Q}(sgn, con)
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

function Base.:*(op1 :: Operation{Q}, op2 :: Operation{Q}) where {Q <: Number} 
    return Operation{Q}(sgn = xor(sgn(op1), sgn(op2)), con = xor(con(op1), con(op2)))
end

# call to Operation
#-------------------------------------------------------------------------------#

function (op :: Operation{Q})(x :: Q) where {Q <: Number}
    x_ = sgn(op) ? -x : x
    return con(op) ? conj(x_) : x_
end

# export
#-------------------------------------------------------------------------------#

export 
    Operation, 
    sgn,
    con