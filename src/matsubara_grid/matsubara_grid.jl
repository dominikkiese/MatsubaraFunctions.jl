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
    function firstindex(
        grid :: AbstractMatsubaraGrid
        )    :: Int64

Returns the index of the first Matsubara frequency in grid
"""
function Base.:firstindex(
    grid :: AbstractMatsubaraGrid
    )    :: Int64

    return firstindex(grid.data)
end

"""
    function lastindex(
        grid :: AbstractMatsubaraGrid
        )    :: Int64

Returns the index of the last Matsubara frequency in grid
"""
function Base.:lastindex(
    grid :: AbstractMatsubaraGrid
    )    :: Int64

    return lastindex(grid.data)
end

"""
    function axes(grid :: AbstractMatsubaraGrid)

Returns range of valid indices for Matsubara grid
"""
function Base.:axes(grid :: AbstractMatsubaraGrid)
    return first(axes(grid.data))
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
    return firstindex(grid) <= index(w) <= lastindex(grid)
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

    return firstindex(grid) <= index(x) <= lastindex(grid)
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

    return lastindex(grid) + 1
end

"""
    function firstvalue(
        grid :: AbstractMatsubaraGrid
        )    :: Float64

Returns the value of the first Matsubara frequency in grid
"""
function firstvalue(
    grid :: AbstractMatsubaraGrid
    )    :: Float64

    return value(grid.data[firstindex(grid)])
end

"""
    function lastvalue(
        grid :: AbstractMatsubaraGrid
        )    :: Float64

Returns the value of the last Matsubara frequency in grid
"""
function lastvalue(
    grid :: AbstractMatsubaraGrid
    )    :: Float64
    
    return value(grid.data[end])
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

    return firstvalue(grid) <= w <= lastvalue(grid)
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
    println("Index range   : $(axes(grid))")
    println("Value range   : ($(firstvalue(grid)), $(lastvalue(grid)))")

    return nothing
end

#----------------------------------------------------------------------------------------------#

include("matsubara_grid_index.jl")
include("matsubara_grid_io.jl")

#----------------------------------------------------------------------------------------------#

export 
    MatsubaraGrid,
    temperature, 
    is_inbounds,
    N,
    firstvalue,
    lastvalue,
    indices, 
    values,
    info