# define particle types for dispatch
abstract type AbstractParticle end 

struct Fermion <: AbstractParticle end 
struct Boson   <: AbstractParticle end



# define grid type
struct MatsubaraGrid
    T    :: Float64 
    data :: Vector{Float64}
    type :: Symbol

    # basic constructor 
    function MatsubaraGrid(
        T    :: Float64,
        data :: Vector{Float64},
        type :: Symbol
        )    :: MatsubaraGrid
    
        return new(T, data, type)
    end

    # convenience constructor for fermionic grid 
    function MatsubaraGrid(
        T    :: Float64,
        N    :: Int64,
             :: Type{Fermion}
        ; 
        plus :: Bool = false
        )    :: MatsubaraGrid 

        if plus
            return MatsubaraGrid(T, Float64[pi * T * (2.0 * n + 1) for n in 0 : N], :Fermion)
        else
            return MatsubaraGrid(T, Float64[pi * T * (2.0 * n + 1) for n in -N - 1 : N], :Fermion)
        end
    end

    # convenience constructor for bosonic grid 
    function MatsubaraGrid(
        T    :: Float64,
        N    :: Int64,
             :: Type{Boson}
        ; 
        plus :: Bool = false
        )    :: MatsubaraGrid 

        if plus
            return MatsubaraGrid(T, Float64[2.0 * pi * T * n for n in 0 : N], :Boson)
        else
            return MatsubaraGrid(T, Float64[2.0 * pi * T * n for n in -N : N], :Boson)
        end
    end
end



# make Matsubara grid indexable
function Base.:getindex(
    grid :: MatsubaraGrid,
    idx  :: Int64 
    )    :: Float64 

    # bounds check performed by Base.Array
    return grid.data[idx]
end



# make Matsubara grid iterable
function Base.length(
    grid :: MatsubaraGrid
    )    :: Int64

    return length(grid.data)
end

function Base.iterate(
    grid :: MatsubaraGrid
    )    :: Tuple{Float64, Int64}

    return grid[1], 1 
end

function Base.iterate(
    grid  :: MatsubaraGrid,
    state :: Int64
    )     :: Union{Nothing, Tuple{Float64, Int64}}

    if state <= length(grid)
        return grid[state], state + 1 
    else 
        return nothing 
    end
end