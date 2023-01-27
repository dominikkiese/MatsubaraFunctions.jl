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

    # convenience constructor for linear fermionic grid (total no. frequencies generated is 2 * N)
    function MatsubaraGrid(
        T :: Float64,
        N :: Int64,
          :: Type{Fermion}
        ) :: MatsubaraGrid{Linear}

        grid = Float64[pi * T * (2 * n + 1) for n in 0 : N - 1]
        return MatsubaraGrid(T, vcat(-reverse(grid), grid), :Fermion, Linear)
    end

    # convenience constructor for coarse fermionic grid (total no. frequencies generated is 2 * N)
    # lb is the integer after which only every δth frequency is considered
    function MatsubaraGrid(
        T  :: Float64,
        N  :: Int64,
        lb :: Int64,
        δ  :: Int64,
           :: Type{Fermion}
        )  :: MatsubaraGrid{Coarse}

        grid   = Vector{Float64}(undef, N)
        offset = δ - 1

        for n in 1 : N
            if n <= lb
                grid[n] = pi * T * (2 * (n - 1) + 1)
            else 
                grid[n] = pi * T * (2 * (n - 1 + offset) + 1)
                offset  = offset + δ - 1
            end 
        end

        return MatsubaraGrid(T, vcat(-reverse(grid), grid), :Fermion, Coarse)
    end

    # convenience constructor for linear bosonic grid (total no. frequencies generated is 2 * N - 1)
    function MatsubaraGrid(
        T :: Float64,
        N :: Int64,
          :: Type{Boson}
        ) :: MatsubaraGrid{Linear}

        grid = Float64[2.0 * pi * T * n for n in 0 : N - 1]
        return MatsubaraGrid(T, vcat(-reverse(grid[2 : end]), grid), :Boson, Linear)
    end

    # convenience constructor for coarse bosonic grid (total no. frequencies generated is 2 * N - 1)
    # lb is the integer after which only every δth frequency is considered
    function MatsubaraGrid(
        T  :: Float64,
        N  :: Int64,
        lb :: Int64,
        δ  :: Int64,
           :: Type{Boson}
        )  :: MatsubaraGrid{Coarse}

        grid   = Vector{Float64}(undef, N)
        offset = δ - 1

        for n in 1 : N
            if n <= lb
                grid[n] = 2.0 * pi * T * (n - 1)
            else 
                grid[n] = 2.0 * pi * T * (n - 1 + offset)
                offset  = offset + δ - 1
            end 
        end

        return MatsubaraGrid(T, vcat(-reverse(grid[2 : end]), grid), :Boson, Coarse)
    end
end



# getter functions 
function temperature(
    grid :: MatsubaraGrid{GT}
    )    :: Float64 where {GT <: AbstractGrid}

    return grid.T 
end 

function type(
    grid :: MatsubaraGrid{GT}
    )    :: Symbol where {GT <: AbstractGrid}

    return grid.type 
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

    if state < length(grid)
        return grid[state + 1], state + 1 
    else 
        return nothing 
    end
end