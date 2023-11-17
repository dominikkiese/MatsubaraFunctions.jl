"""
    struct MatsubaraGrid{PT <: AbstractParticle} <: AbstractMatsubaraGrid

MatsubaraGrid type with fields:
* `T    :: Float64` : temperature
* `data :: OffsetVector{MatsubaraFrequency{PT}, Vector{MatsubaraFrequency{PT}}}` : list of MatsubaraFrequency objects
"""
struct MatsubaraGrid{PT <: AbstractParticle} <: AbstractMatsubaraGrid
    T    :: Float64 
    data :: OffsetVector{MatsubaraFrequency{PT}, Vector{MatsubaraFrequency{PT}}}

    function MatsubaraGrid(
        T    :: Float64,
        data :: OffsetVector{MatsubaraFrequency{PT}, Vector{MatsubaraFrequency{PT}}},
        )    :: MatsubaraGrid{PT} where {PT <: AbstractParticle}

        for w in data 
            @DEBUG T ≈ temperature(w) "Temperature must be equal between frequencies and grid"
        end  
    
        return new{PT}(T, data)
    end

    function MatsubaraGrid(
        T :: Float64,
        N :: Int64,
          :: Type{Fermion}
        ) :: MatsubaraGrid{Fermion}

        return MatsubaraGrid(T, OffsetVector([MatsubaraFrequency(T, n, Fermion) for n in -N : N - 1], -N - 1))
    end

    function MatsubaraGrid(
        T :: Float64,
        N :: Int64,
          :: Type{Boson}
        ) :: MatsubaraGrid{Boson}

        return MatsubaraGrid(T, OffsetVector([MatsubaraFrequency(T, n, Boson) for n in -N + 1 : N - 1], -N))
    end

    function MatsubaraGrid(
        grid :: AbstractMatsubaraGrid
        )    :: AbstractMatsubaraGrid

        return MatsubaraGrid(temperature(grid), copy(grid.data))
    end
end

"""
    function temperature(
        grid :: AbstractMatsubaraGrid
        )    :: Float64

Returns `grid.T`
"""
function temperature(
    grid :: AbstractMatsubaraGrid
    )    :: Float64

    return grid.T 
end 

function Base.:length(
    grid :: AbstractMatsubaraGrid
    )    :: Int64

    return length(grid.data)
end

"""
    function first_index(
        grid :: AbstractMatsubaraGrid
        )    :: Int64

Returns the index of the first Matsubara frequency in grid
"""
function first_index(
    grid :: AbstractMatsubaraGrid
    )    :: Int64

    return firstindex(grid.data)
end

"""
    function last_index(
        grid :: AbstractMatsubaraGrid
        )    :: Int64

Returns the index of the last Matsubara frequency in grid
"""
function last_index(
    grid :: AbstractMatsubaraGrid
    )    :: Int64

    return lastindex(grid.data)
end

"""
    function index_range(
        grid :: AbstractMatsubaraGrid
        )    :: NTuple{2, Int64}

Returns indices of the first and last Matsubara frequency in grid
"""
function index_range(
    grid :: AbstractMatsubaraGrid
    )    :: NTuple{2, Int64}

    return first_index(grid), last_index(grid)
end

function Base.:copy(
    grid :: AbstractMatsubaraGrid
    )    :: AbstractMatsubaraGrid

    return MatsubaraGrid(grid)
end

"""
    function is_inbounds(
        w    :: MatsubaraFrequency{PT},
        grid :: MatsubaraGrid{PT}
        )    :: Bool where {PT <: AbstractParticle}

Checks if `w` is contained in grid
"""
function is_inbounds(
    w    :: MatsubaraFrequency{PT},
    grid :: MatsubaraGrid{PT}
    )    :: Bool where {PT <: AbstractParticle}

    @DEBUG temperature(w) ≈ temperature(grid) "Temperature must be equal between frequency and grid"
    return first_index(grid) <= index(w) <= last_index(grid)
end

"""
    function is_inbounds(
        x    :: MatsubaraIndex{PT},
        grid :: MatsubaraGrid{PT}
        )    :: Bool where {PT <: AbstractParticle}

Checks if `x` is contained in grid
"""
function is_inbounds(
    x    :: MatsubaraIndex{PT},
    grid :: MatsubaraGrid{PT}
    )    :: Bool where {PT <: AbstractParticle}

    return first_index(grid) <= index(x) <= last_index(grid)
end

"""
    function N(
        grid :: AbstractMatsubaraGrid
        )    :: Int64

Returns that value `N` used to construct the grid
"""
function N(
    grid :: AbstractMatsubaraGrid
    )    :: Int64

    return last_index(grid) + 1
end

"""
    function first_value(
        grid :: AbstractMatsubaraGrid
        )    :: Float64

Returns the value of the first Matsubara frequency in grid
"""
function first_value(
    grid :: AbstractMatsubaraGrid
    )    :: Float64

    return value(grid.data[first_index(grid)])
end

"""
    function last_value(
        grid :: AbstractMatsubaraGrid
        )    :: Float64

Returns the value of the last Matsubara frequency in grid
"""
function last_value(
    grid :: AbstractMatsubaraGrid
    )    :: Float64
    
    return value(grid.data[end])
end

"""
    function value_range(
        grid :: AbstractMatsubaraGrid
        )    :: NTuple{2, Float64}

Returns values of the first and last Matsubara frequency in grid
"""
function value_range(
    grid :: AbstractMatsubaraGrid
    )    :: NTuple{2, Float64}

    return first_value(grid), last_value(grid)
end

"""
    function is_inbounds(
        w    :: Float64,
        grid :: AbstractMatsubaraGrid
        )    :: Bool

Checks if `w` lies within grid bounds
"""
function is_inbounds(
    w    :: Float64,
    grid :: AbstractMatsubaraGrid
    )    :: Bool

    return first_value(grid) <= w <= last_value(grid)
end

"""
    function indices(grid :: AbstractMatsubaraGrid)   

Returns list of indices for Matsubara frequencies in grid
"""
function indices(grid :: AbstractMatsubaraGrid)    
    return index.(grid.data)
end 

"""
    Base.:values(grid :: AbstractMatsubaraGrid)

Returns list of values for Matsubara frequencies in grid
"""
function Base.:values(grid :: AbstractMatsubaraGrid)
    return value.(grid.data)
end 

"""
    function info(
        grid :: MatsubaraGrid{PT}
        )    :: Nothing where {PT <: AbstractParticle}

Prints summary of grid properties
"""
function info(
    grid :: MatsubaraGrid{PT}
    )    :: Nothing where {PT <: AbstractParticle}

    println("MatsubaraGrid properties")
    println("------------------------")
    println("Particle type : $(string(PT))")
    println("Temperature   : $(temperature(grid))")
    println("Length        : $(length(grid))")
    println("Index range   : $(index_range(grid))")
    println("Value range   : $(value_range(grid))")

    return nothing
end

#----------------------------------------------------------------------------------------------#

include("matsubara_grid_index.jl")
include("matsubara_grid_io.jl")

#----------------------------------------------------------------------------------------------#

export 
    MatsubaraGrid,
    temperature,
    first_index,
    last_index,
    index_range, 
    is_inbounds,
    N,
    first_value,
    last_value,
    value_range,
    indices, 
    values,
    info