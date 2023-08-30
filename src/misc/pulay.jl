# algorithm reproduced from https://www.sciencedirect.com/science/article/abs/pii/S0009261416000464?via%3Dihub
"""
    struct PeriodicPulay{Q <: Number}

PeriodicPulay type with fields:
* `x     :: Vector{Q}`       : solution vector
* `Fs    :: Matrix{Q}`       : history for function evaluations 
* `Xs    :: Matrix{Q}`       : history for intermediate solutions
* `aerrs :: Vector{Float64}` : absolute errors
* `rerrs :: Vector{Float64}` : relative errors
"""
struct PeriodicPulay{Q <: Number}
    x     :: Vector{Q}
    Fs    :: Matrix{Q}
    Xs    :: Matrix{Q}
    aerrs :: Vector{Float64}
    rerrs :: Vector{Float64}
    
    function PeriodicPulay(
        x :: Vector{Q}
        ;
        m :: Int64 = 5
        ) :: PeriodicPulay{Q} where {Q <: Number}
        
        @check m >= 1 "Memory size must be >= 1"
        Fs = Matrix{Q}(undef, length(x), m)
        Xs = Matrix{Q}(undef, length(x), m)
        
        return new{Q}(copy(x), Fs, Xs, Float64[], Float64[])
    end 
end

"""
    function solve!(
        f       :: Function,
        P       :: PeriodicPulay{Q}
        ;
        p       :: Int64   = 2,
        iters   :: Int64   = 100,
        α       :: Float64 = 0.25,
        atol    :: Float64 = 1e-8,
        rtol    :: Float64 = 1e-8,
        verbose :: Bool    = false
        )       :: Nothing where {Q <: Number}

Runs the periodic Pulay solver. Here, f(x) = g(x) - x is the residue for the fixed-point equation g(x) = x.
The following keyword arguments are supported:
* `p`       : Pulay period (every p-th iteration Pulay mixing is used)
* `iters`   : maximum number of iterations
* `α`       : mixing factor
* `atol`    : absolute error tolerance
* `rtol`    : relative error tolerance
* `verbose` : show intermediate results?
"""
function solve!(
    f       :: Function,
    P       :: PeriodicPulay{Q}
    ;
    p       :: Int64   = 2,
    iters   :: Int64   = 100,
    α       :: Float64 = 0.25,
    atol    :: Float64 = 1e-8,
    rtol    :: Float64 = 1e-8,
    verbose :: Bool    = false
    )       :: Nothing where {Q <: Number}
    
    # reference fields of P 
    x     = P.x 
    Fs    = P.Fs 
    Xs    = P.Xs
    aerrs = P.aerrs 
    rerrs = P.rerrs
    m     = size(Fs, 2)
    
    # throw warning if not in comfort zone
    mdiv2 = iseven(m) ? Int64(m / 2) : Int64((m + 1) / 2)
    
    if !(2 <= p <= mdiv2)
        @warn "Pulay period not in (2, $(mdiv2)), convergence might be suboptimal"
    end

    verbose && println("Running periodic Pulay solver ...")
    ti = time()
    
    # initial iteration 
    Fp = f(x); xp = copy(x); x .+= α .* Fp

    # init errors, iteration count and memory index
    aerr = Inf 
    rerr = Inf
    iter = 0
    midx = 1
    
    # can we avoid the copies somehow?
    while (aerr > atol) && (rerr > rtol) && (iter < iters)
        if (iter + 1) % p > 0
            verbose && println()
            verbose && println("iteration $(iter) | linear")
            
            F            = f(x)
            Fs[:, midx] .= F .- Fp
            Xs[:, midx] .= x .- xp
            
            aerr = norm(F)
            rerr = aerr / norm(x)
            
            # linear mixing
            Fp .= F; xp .= x; x .+= α .* F
        else 
            verbose && println()
            verbose && println("iteration $(iter) | Pulay")
            
            F            = f(x)
            Fs[:, midx] .= F .- Fp
            Xs[:, midx] .= x .- xp
            
            aerr = norm(F)
            rerr = aerr / norm(x)
            
            # Pulay mixing (use whole history thus far)
            Fmat = view(Fs, :, 1 : midx)
            Xmat = view(Xs, :, 1 : midx)
            
            # use Moore-Penrose pseudoinverse from Base for stability
            Fp .= F; xp .= x; x .+= α .* F .- (Xmat .+ α .* Fmat) * pinv(Fmat) * F
        end
        
        push!(aerrs, aerr)
        push!(rerrs, rerr)
        verbose && println("=> aerr = $(aerr)")
        verbose && println("=> rerr = $(rerr)")
        iter += 1
        midx += 1
        
        # shift solutions in the history kernel once its full (can we circumvent this somehow?)
        if midx > m
            midx = m
            
            for i in 1 : m - 1
                Fs[:, i] .= view(Fs, :, i + 1)
                Xs[:, i] .= view(Xs, :, i + 1)
            end
        end
    end
    
    dt = time() - ti
    verbose && println()
    verbose && println("Done. Time elapsed $(dt)s.")
            
    return nothing
end