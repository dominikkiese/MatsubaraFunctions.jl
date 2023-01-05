# function to generate grid of fermionic Matsubara frequencies at temperature T 
function mk_fermion_grid(
    T    :: Float64, 
    N    :: Int64
    ; 
    plus :: Bool = false
    )    :: Vector{Float64}

    if plus
        return Float64[pi * T * (2.0 * n + 1) for n in 0 : N]
    else
        return Float64[pi * T * (2.0 * n + 1) for n in (-N - 1) : N]
    end
end

# function to generate grid of bosonic Matsubara frequencies at temperature T 
function mk_boson_grid(    
    T    :: Float64, 
    N    :: Int64
    ; 
    plus :: Bool = false
    )    :: Vector{Float64}

    if plus
        return Float64[2.0 * pi * T * n for n in 0 : N]
    else
        return Float64[2.0 * pi * T * n for n in (-N) : N]
    end
end

# function to generate grid of Matsubara frequencies at temperature T 
function mk_grid(    
    T    :: Float64, 
    N    :: Int64,
    S    :: Symbol
    ; 
    plus :: Bool = false
    )    :: Vector{Float64}

    if S == :Fermion
        return mk_fermion_grid(T, N; plus)
    elseif S == :Boson
        return mk_boson_grid(T, N; plus)
    else
        error("Species must be Fermion or Boson")
    end
end