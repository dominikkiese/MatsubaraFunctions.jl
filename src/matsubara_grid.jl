# Note: we do not allow MatsubaraGrid to be dispatched on the particle type 
#       to allow for mixed grids in the construction of MatsubaraFunctions
struct MatsubaraGrid
    T     :: Float64 
    data  :: Vector{MatsubaraFrequency}
    type  :: Symbol          

    # basic constructor 
    function MatsubaraGrid(
        T      :: Float64,
        data   :: Vector{MatsubaraFrequency},
        type   :: Symbol
        ;
        checks :: Bool = true
        )      :: MatsubaraGrid

        if checks 
            for w in data 
                @assert type == w.type "Particle type must be equal between frequencies and grid"
                @assert T ≈ temperature(w) "Temperature must be equal between frequencies and grid"
            end  
        end  
    
        return new(T, data, type)
    end

    # convenience constructor for fermionic grid (total no. frequencies generated is 2 * N)
    function MatsubaraGrid(
        T :: Float64,
        N :: Int64,
          :: Type{Fermion}
        ) :: MatsubaraGrid

        grid = MatsubaraFrequency[MatsubaraFrequency(T, n, Fermion) for n in -N : N - 1]
        return MatsubaraGrid(T, grid, :Fermion; checks = false)
    end

    # convenience constructor for bosonic grid (total no. frequencies generated is 2 * N - 1)
    function MatsubaraGrid(
        T :: Float64,
        N :: Int64,
          :: Type{Boson}
        ) :: MatsubaraGrid

        grid = MatsubaraFrequency[MatsubaraFrequency(T, n, Boson) for n in -N + 1 : N - 1]
        return MatsubaraGrid(T, grid, :Boson; checks = false)
    end
end



# getter functions
function temperature(
    grid :: MatsubaraGrid
    )    :: Float64

    return grid.T 
end 

function type(
    grid :: MatsubaraGrid
    )    :: Symbol

    return grid.type 
end

function Base.:length(
    grid :: MatsubaraGrid
    )    :: Int64

    return length(grid.data)
end

function index_range(
    grid :: MatsubaraGrid
    )    :: NTuple{2, Int64}

    return index(grid.data[1]), index(grid.data[end])
end

function is_inbounds(
    w    :: MatsubaraFrequency,
    grid :: MatsubaraGrid
    )    :: Bool

    # perform checks, otherwise question is ill-defined
    @assert type(w) == type(grid) "Particle type must be equal between frequency and grid"
    @assert temperature(w) ≈ temperature(grid) "Temperature must be equal between frequency and grid"
    idx_range = index_range(grid)
    return idx_range[1] <= index(w) <= idx_range[2]
end

function is_inbounds(
    w    :: Float64,
    grid :: MatsubaraGrid
    )    :: Bool

    return value(grid[1]) <= w <= value(grid[end])
end

function info(
    grid :: MatsubaraGrid
    )    :: Nothing

    println("MatsubaraGrid properties")
    println("------------------------")
    println("Particle type : $(type(grid))")
    println("Temperature   : $(temperature(grid))")
    println("Length        : $(length(grid))")
    println("Index range   : $(index_range(grid))")

    return nothing
end



# make Matsubara grid indexable
function Base.:eachindex(
    grid :: MatsubaraGrid
    )    :: Base.OneTo{Int64}

    return eachindex(grid.data)
end

function Base.:lastindex(
    grid :: MatsubaraGrid
    )    :: Int64

    return lastindex(grid.data)
end

function Base.:getindex(
    grid :: MatsubaraGrid,
    idx  :: Int64 
    )    :: MatsubaraFrequency

    # bounds check performed by Base.Array
    return grid.data[idx]
end

function Base.:getindex(
    grid :: MatsubaraGrid,
    idxs :: UnitRange{Int64}
    )    :: SubArray{MatsubaraFrequency, 1, Vector{MatsubaraFrequency}, Tuple{UnitRange{Int64}}, true}

    # bounds check performed by Base.Array 
    return @view grid.data[idxs]
end  



# unsafe method for converting MatsubaraFrequency to grid index (no bounds check)
function grid_index(
    w    :: MatsubaraFrequency,
    grid :: MatsubaraGrid
    )    :: Int64 

    return index(w) - index(grid[1]) + 1
end

# make MatsubaraGrid callable with MatsubaraFrequency
# returns index to data array corresponding to this frequency if in grid
@inline function (f :: MatsubaraGrid)(
    w :: MatsubaraFrequency
    ) :: Int64

    if is_inbounds(w, f)
        return grid_index(w, f)
    else 
        error("Frequency not in grid")
    end 
end

# make linear MatsubaraGrid callable with Float64
# returns index to data array corresponding to closest frequency if in grid
@inline function (f :: MatsubaraGrid)(
    w :: Float64
    ) :: Int64

    if is_inbounds(w, f)
        delta    = value(f[2]) - value(f[1])
        position = (w - value(f[1])) / delta
        return round(Int64, position) + 1
    else 
        error("Frequency not in grid")
    end 
end



# make Matsubara grid iterable
function Base.:iterate(
    grid :: MatsubaraGrid
    )    :: Tuple{MatsubaraFrequency, Int64}

    return grid[1], 1 
end

function Base.:iterate(
    grid  :: MatsubaraGrid,
    state :: Int64
    )     :: Union{Nothing, Tuple{MatsubaraFrequency, Int64}}

    if state < length(grid)
        return grid[state + 1], state + 1 
    else 
        return nothing 
    end
end