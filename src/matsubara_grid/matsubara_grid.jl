# we do not allow MatsubaraGrid to be dispatched on the particle type 
# to allow for mixed grids in the construction of MatsubaraFunctions
"""
    struct MatsubaraGrid

MatsubaraGrid type with fields:
* `T    :: Float64`                    : temperature
* `data :: Vector{MatsubaraFrequency}` : list of MatsubaraFrequency objects
* `type :: Symbol`                     : particle type
"""
struct MatsubaraGrid
    T     :: Float64 
    data  :: Vector{MatsubaraFrequency}
    type  :: Symbol          

    # default constructor 
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



"""
    function temperature(
        grid :: MatsubaraGrid
        )    :: Float64

Returns `grid.T`
"""
function temperature(
    grid :: MatsubaraGrid
    )    :: Float64

    return grid.T 
end 

"""
    function type(
        grid :: MatsubaraGrid
        )    :: Symbol

Returns `grid.type`
"""
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

"""
    function first_index(
        grid :: MatsubaraGrid
        )    :: Int64 

Returns the index of the first Matsubara frequency in grid
"""
function first_index(
    grid :: MatsubaraGrid
    )    :: Int64 

    return index(grid.data[1])
end

"""
    function last_index(
        grid :: MatsubaraGrid
        )    :: Int64 

Returns the index of the last Matsubara frequency in grid
"""
function last_index(
    grid :: MatsubaraGrid
    )    :: Int64 

    return index(grid.data[end])
end

"""
    function index_range(
        grid :: MatsubaraGrid
        )    :: NTuple{2, Int64}

Returns indices of the first and last Matsubara frequency in grid
"""
function index_range(
    grid :: MatsubaraGrid
    )    :: NTuple{2, Int64}

    return first_index(grid), last_index(grid)
end

"""
    function is_inbounds(
        w    :: MatsubaraFrequency,
        grid :: MatsubaraGrid
        )    :: Bool

Checks if `w` is contained in grid
"""
function is_inbounds(
    w    :: MatsubaraFrequency,
    grid :: MatsubaraGrid
    )    :: Bool

    # perform checks, otherwise question is ill-defined
    @assert type(w) == type(grid) "Particle type must be equal between frequency and grid"
    @assert temperature(w) ≈ temperature(grid) "Temperature must be equal between frequency and grid"
    return first_index(grid) <= index(w) <= last_index(grid)
end

"""
    function first_value(
        grid :: MatsubaraGrid
        )    :: Float64

Returns the value of the first Matsubara frequency in grid
"""
function first_value(
    grid :: MatsubaraGrid
    )    :: Float64

    return value(grid.data[1])
end

"""
    function last_value(
        grid :: MatsubaraGrid
        )    :: Float64

Returns the value of the last Matsubara frequency in grid
"""
function last_value(
    grid :: MatsubaraGrid
    )    :: Float64

    return value(grid.data[end])
end

"""
    function value_range(
        grid :: MatsubaraGrid
        )    :: NTuple{2, Float64}

Returns values of the first and last Matsubara frequency in grid
"""
function value_range(
    grid :: MatsubaraGrid
    )    :: NTuple{2, Float64}

    return first_value(grid), last_value(grid)
end

"""
    function is_inbounds(
        w    :: Float64,
        grid :: MatsubaraGrid
        )    :: Bool

Checks if `w` lies within grid bounds
"""
function is_inbounds(
    w    :: Float64,
    grid :: MatsubaraGrid
    )    :: Bool

    return first_value(grid) <= w <= last_value(grid)
end

"""
    function info(
        grid :: MatsubaraGrid
        )    :: Nothing

Prints summary of grid properties
"""
function info(
    grid :: MatsubaraGrid
    )    :: Nothing

    println("MatsubaraGrid properties")
    println("------------------------")
    println("Particle type : $(type(grid))")
    println("Temperature   : $(temperature(grid))")
    println("Length        : $(length(grid))")
    println("Index range   : $(index_range(grid))")
    println("Value range   : $(value_range(grid))")

    return nothing
end



# load methods
include("matsubara_grid_index.jl")
include("matsubara_grid_io.jl")