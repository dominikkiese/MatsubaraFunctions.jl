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
        
        @DEBUG m >= 1 "Memory size must be >= 1"
        Fs = Matrix{Q}(undef, length(x), m)
        Xs = Matrix{Q}(undef, length(x), m)
        
        return new{Q}(copy(x), Fs, Xs, Float64[], Float64[])
    end 
end

"""
    function solve!(
        f!      :: Function,
        P       :: PeriodicPulay{Q}
        ;
        p       :: Int64   = 3,
        iters   :: Int64   = 100,
        α       :: Float64 = 0.5,
        atol    :: Float64 = 1e-8,
        rtol    :: Float64 = 1e-8,
        verbose :: Bool    = false
        )       :: Nothing where {Q <: Number}

Runs the periodic Pulay solver. Here, f(x) = g(x) - x computes the residue for the fixed-point equation g(x) = x.
`f!` should have the form `(F, x) -> f!(F, x)`, such that the residue can be written into a pre-allocated array `F`. 
The following keyword arguments are supported:
* `p`       : Pulay period (every p-th iteration Pulay mixing is used)
* `iters`   : maximum number of iterations
* `α`       : mixing factor
* `atol`    : absolute error tolerance
* `rtol`    : relative error tolerance
* `verbose` : show intermediate results?
"""
function solve!(
    f!      :: Function,
    P       :: PeriodicPulay{Q}
    ;
    p       :: Int64   = 3,
    iters   :: Int64   = 100,
    α       :: Float64 = 0.5,
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

    # allocate buffers 
    F  = copy(x)
    Fp = copy(x)
    xp = copy(x)
    
    # initial iteration 
    f!(Fp, x)
    x .+= α .* Fp

    # init errors, iteration count and memory index
    aerr = Inf 
    rerr = Inf
    iter = 0
    midx = 1

    verbose && mpi_println("   iter   |   type   |   abs. error   |   rel. error   |   time elapsed (s)   ")
    verbose && mpi_println("------------------------------------------------------------------------------")
    
    while (aerr > atol) && (rerr > rtol) && (iter < iters)
        ti = time()

        f!(F, x)
        Fs[:, midx] .= F .- Fp
        Xs[:, midx] .= x .- xp
        Fp          .= F
        xp          .= x

        aerr = norm(F, Inf)
        rerr = aerr / norm(x, Inf)

        if (iter + 1) % p > 0
            # linear mixing
            x .+= α .* F
        else 
            # Pulay mixing Moore-Penrose pseudoinverse
            Fmat  = view(Fs, :, 1 : midx)
            Xmat  = view(Xs, :, 1 : midx)
            temp  = pinv(Fmat) * F
            x    .= Xmat * temp
            F    .= Fmat * temp 
            x    .= xp .+ α .* Fp .- (x .+ α .* F)
        end
        
        push!(aerrs, aerr)
        push!(rerrs, rerr)

        if mpi_ismain()
            if (iter + 1) % p > 0
                verbose && @printf "    %5d |    LM    |  %5e  |  %5e  |     %5e   \n" iter aerr rerr time() - ti
                verbose && println("------------------------------------------------------------------------------")
            else
                verbose && @printf "    %5d |    PM    |  %5e  |  %5e  |     %5e   \n" iter aerr rerr time() - ti
                verbose && println("------------------------------------------------------------------------------")
            end
        end

        iter += 1
        midx += 1
        
        # shift results in history kernel to the left once its full
        if midx > m
            midx = m

            for i in 1 : m - 1
                Fs[:, i] .= view(Fs, :, i + 1)
                Xs[:, i] .= view(Xs, :, i + 1)
            end
        end
    end
            
    return nothing
end

#----------------------------------------------------------------------------------------------#

export 
    PeriodicPulay,
    solve!