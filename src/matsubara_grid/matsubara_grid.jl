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



# load methods
include("matsubara_grid_index.jl")