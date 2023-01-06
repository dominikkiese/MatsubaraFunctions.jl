abstract type MatsubaraGrid end 



# fermionic Matsubara grid
struct FermionGrid <: MatsubaraGrid
    T    :: Float64 
    data :: Vector{Float64}

    # basic constructor 
    function FermionGrid(
        T    :: Float64,
        data :: Vector{Float64}
        )    :: FermionGrid 

        return new(T, data)
    end

    # convenience constructor
    function FermionGrid(
        T    :: Float64,
        N    :: Int64
        ; 
        plus :: Bool = false
        )    :: FermionGrid 

        if plus
            return new(T, Float64[pi * T * (2.0 * n + 1) for n in 0 : N])
        else
            return new(T, Float64[pi * T * (2.0 * n + 1) for n in -N - 1 : N])
        end
    end
end

# make fermionic Matsubara grid indexable
function Base.:getindex(
    fg  :: FermionGrid,
    idx :: Int64 
    )   :: Float64 

    # bounds check performed by Base.Array
    return fg.data[idx]
end

# make fermionic Matsubara grid iterable
function Base.length(
    fg :: FermionGrid
    )  :: Int64

    return length(fg.data)
end

function Base.iterate(
    fg :: FermionGrid
    )  :: Tuple{Float64, Int64}

    return fg[1], 1 
end

function Base.iterate(
    fg    :: FermionGrid,
    state :: Int64
    )     :: Union{Nothing, Tuple{Float64, Int64}}

    if state <= length(fg)
        return fg[state], state + 1 
    else 
        return nothing 
    end
end



# bosonic Matsubara grid
struct BosonGrid <: MatsubaraGrid
    T    :: Float64 
    data :: Vector{Float64}

    # basic constructor 
    function BosonGrid(
        T    :: Float64,
        data :: Vector{Float64}
        )    :: BosonGrid 

        return new(T, data)
    end

    # convenience constructor
    function BosonGrid(
        T    :: Float64,
        N    :: Int64
        ; 
        plus :: Bool = false
        )    :: BosonGrid 

        if plus
            return new(T, Float64[2.0 * pi * T * n for n in 0 : N])
        else
            return new(T, Float64[2.0 * pi * T * n for n in -N : N])
        end
    end
end

# make bosonic Matsubara grid indexable
function Base.:getindex(
    bg  :: BosonGrid,
    idx :: Int64 
    )   :: Float64 

    # bounds check performed by Base.Array
    return bg.data[idx]
end

# make bosonic Matsubara grid iterable
function Base.length(
    bg :: BosonGrid
    )  :: Int64

    return length(bg.data)
end

function Base.iterate(
    bg :: BosonGrid
    )  :: Tuple{Float64, Int64}

    return bg[1], 1 
end

function Base.iterate(
    bg    :: BosonGrid,
    state :: Int64
    )     :: Union{Nothing, Tuple{Float64, Int64}}

    if state <= length(bg)
        return bg[state], state + 1 
    else 
        return nothing 
    end
end