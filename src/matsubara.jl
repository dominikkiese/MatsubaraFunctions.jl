# Note: we do not allow MatsubaraGrid to be dispatched on the particle type to allow for mixed grids 
#       in the construction of MatsubaraFunctions
struct MatsubaraGrid{GT <: AbstractGrid}
    T    :: Float64 
    data :: Vector{Float64}
    type :: Symbol          

    # basic constructor 
    function MatsubaraGrid(
        T    :: Float64,
        data :: Vector{Float64},
        type :: Symbol,
             :: Type{GT}
        )    :: MatsubaraGrid{GT} where {GT <: AbstractGrid}
    
        return new{GT}(T, data, type)
    end

    # convenience constructor for linear fermionic grid 
    function MatsubaraGrid(
        T    :: Float64,
        N    :: Int64,
             :: Type{Fermion}
        ; 
        plus :: Bool = false
        )    :: MatsubaraGrid{Linear}

        grid = Float64[pi * T * (2 * n + 1) for n in 0 : N]

        if plus
            return MatsubaraGrid(T, grid, :Fermion, Linear)
        else
            return MatsubaraGrid(T, vcat(-reverse(grid), grid), :Fermion, Linear)
        end
    end

    # convenience constructor for coarse fermionic grid 
    function MatsubaraGrid(
        T    :: Float64,
        N    :: Int64,
        z    :: Float64,
             :: Type{Fermion}
        ; 
        plus :: Bool = false
        )    :: MatsubaraGrid{Coarse}

        grid = Float64[pi * T * (2 * round(Int64, z * sinh(n / z)) + 1) for n in 0 : N]

        if plus
            return MatsubaraGrid(T, grid, :Fermion, Coarse)
        else
            return MatsubaraGrid(T, vcat(-reverse(grid), grid), :Fermion, Coarse)
        end
    end

    # convenience constructor for linear bosonic grid 
    function MatsubaraGrid(
        T    :: Float64,
        N    :: Int64,
             :: Type{Boson}
        ; 
        plus :: Bool = false
        )    :: MatsubaraGrid{Linear}

        if plus
            return MatsubaraGrid(T, Float64[2.0 * pi * T * n for n in 0 : N], :Boson, Linear)
        else
            return MatsubaraGrid(T, Float64[2.0 * pi * T * n for n in -N : N], :Boson, Linear)
        end
    end

    # convenience constructor for coarse bosonic grid 
    function MatsubaraGrid(
        T    :: Float64,
        N    :: Int64,
        z    :: Float64,
             :: Type{Boson}
        ; 
        plus :: Bool = false
        )    :: MatsubaraGrid{Coarse}

        if plus
            return MatsubaraGrid(T, Float64[2.0 * pi * T * round(Int64, z * sinh(n / z)) for n in 0 : N], :Boson, Coarse)
        else
            return MatsubaraGrid(T, Float64[2.0 * pi * T * round(Int64, z * sinh(n / z)) for n in -N : N], :Boson, Coarse)
        end
    end
end



# make Matsubara grid indexable
function Base.:lastindex(
    grid :: MatsubaraGrid{GT}
    )    :: Int64 where {GT <: AbstractGrid}

    return lastindex(grid.data)
end

function Base.:getindex(
    grid :: MatsubaraGrid{GT},
    idx  :: Int64 
    )    :: Float64 where {GT <: AbstractGrid}

    # bounds check performed by Base.Array
    return grid.data[idx]
end

function Base.:getindex(
    grid :: MatsubaraGrid{GT},
    idxs :: UnitRange{Int64},
    )    :: Vector{Float64} where {GT <: AbstractGrid}

    # bounds check performed by Base.Array 
    return grid.data[idxs]
end  



# make Matsubara grid iterable
function Base.length(
    grid :: MatsubaraGrid{GT}
    )    :: Int64 where {GT <: AbstractGrid}

    return length(grid.data)
end

function Base.iterate(
    grid :: MatsubaraGrid{GT}
    )    :: Tuple{Float64, Int64} where {GT <: AbstractGrid}

    return grid[1], 1 
end

function Base.iterate(
    grid  :: MatsubaraGrid{GT},
    state :: Int64
    )     :: Union{Nothing, Tuple{Float64, Int64}} where {GT <: AbstractGrid}

    if state <= length(grid)
        return grid[state], state + 1 
    else 
        return nothing 
    end
end