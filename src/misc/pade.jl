# algorithm reproduced from https://link.springer.com/article/10.1007/BF00655090
"""
    struct PadeApprox{Q <: Number}

PadeApprox type with fields:
* `coeffs :: Vector{Q}` : coefficients for continued fraction fit
* `xdat   :: Vector{Q}` : evaluation nodes for data
"""
struct PadeApprox{Q <: Number}
    coeffs :: Vector{Q}
    xdat   :: Vector{Q}

    function PadeApprox(xdat :: Vector{Q}, ydat :: Vector{Q}) where {Q <: Number}  
        @DEBUG length(xdat) > 1 "More than one data point must be provided"
        g        = zeros(Q, length(xdat), length(xdat))
        g[1, :] .= ydat
        
        for i in 2 : size(g, 2)
            # terminate if coefficients become too small
            if abs(g[i - 1, i - 1]) < 1e-10
                return @views new{Q}(diag(g)[1 : i - 1], xdat[1 : i - 1])
            end 
            
            for j in 1 : size(g, 1)
                a       = g[i - 1, i - 1] - g[i - 1, j]
                b       = (xdat[j] - xdat[i - 1]) * g[i - 1, j]
                g[i, j] = a / b
            end 
        end
            
        return new{Q}(diag(g), xdat)
    end   
end 

"""
    function coeffs(PA :: PadeApprox{Q}) :: Vector{Q} where {Q <: Number}

Returns `PA.coeffs`
"""
function coeffs(PA :: PadeApprox{Q}) :: Vector{Q} where {Q <: Number}
    return PA.coeffs 
end

"""
    function xdat(PA :: PadeApprox{Q}) :: Vector{Q} where {Q <: Number}

Returns `PA.xdat`
"""
function xdat(PA :: PadeApprox{Q}) :: Vector{Q} where {Q <: Number}
    return PA.xdat
end

function (PA :: PadeApprox{Q})(z :: Q) where {Q <: Number}
    # init recursion for enumerator
    A1 = PA.coeffs[1]
    A2 = A1
    
    # init recursion for denominator
    B1 = Q(1)
    B2 = B1 + (z - PA.xdat[1]) * PA.coeffs[2]
    
    for i in 3 : length(PA.coeffs)
        p = (z - PA.xdat[i - 1]) * PA.coeffs[i]
          
        # update enumerator
        Ai = A2 + p * A1        
        A1 = A2 
        A2 = Ai 
        
        # update denominator
        Bi = B2 + p * B1        
        B1 = B2 
        B2 = Bi 
    end
        
    return A2 / B2
end

# export
#-------------------------------------------------------------------------------#

export 
    PadeApprox,
    coeffs,
    xdat