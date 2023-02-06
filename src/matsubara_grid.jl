# Note: we do not allow MatsubaraGrid to be dispatched on the particle type 
#       to allow for mixed grids in the construction of MatsubaraFunctions
struct MatsubaraGrid{GT <: AbstractGrid}
    T     :: Float64 
    data  :: Vector{MatsubaraFrequency}
    type  :: Symbol          

    # basic / safe constructor 
    function MatsubaraGrid(
        T      :: Float64,
        data   :: Vector{MatsubaraFrequency},
        type   :: Symbol,
               :: Type{GT}
        ;
        checks :: Bool = true
        )      :: MatsubaraGrid{GT} where {GT <: AbstractGrid}

        if checks 
            for w in data 
                @assert type == w.type "Particle type must be equal between frequencies and grid"
                @assert isapprox(T, temperature(w)) "Temperature must be equal between frequencies and grid"
            end  
        end  
    
        return new{GT}(T, data, type)
    end

    # convenience constructor for linear fermionic grid (total no. frequencies generated is 2 * N)
    function MatsubaraGrid(
        T :: Float64,
        N :: Int64,
          :: Type{Fermion}
        ) :: MatsubaraGrid{Linear}

        grid = MatsubaraFrequency[MatsubaraFrequency(T, n, Fermion) for n in -N : N - 1]
        return MatsubaraGrid(T, grid, :Fermion, Linear; checks = false)
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

        grid_idxs_p = Vector{Int64}(undef, N)
        grid_idxs_m = Vector{Int64}(undef, N)
        offset      = δ - 1

        for n in 1 : N
            if n <= lb
                grid_idxs_p[n] = n - 1
                grid_idxs_m[n] = -grid_idxs_p[n] - 1
            else 
                grid_idxs_p[n] = n - 1 + offset
                grid_idxs_m[n] = -grid_idxs_p[n] - 1
                offset         = offset + δ - 1
            end 
        end

        grid = MatsubaraFrequency[MatsubaraFrequency(T, n, Fermion) for n in vcat(reverse(grid_idxs_m), grid_idxs_p)]
        return MatsubaraGrid(T, grid, :Fermion, Coarse; checks = false)
    end

    # convenience constructor for linear bosonic grid (total no. frequencies generated is 2 * N - 1)
    function MatsubaraGrid(
        T :: Float64,
        N :: Int64,
          :: Type{Boson}
        ) :: MatsubaraGrid{Linear}

        grid = MatsubaraFrequency[MatsubaraFrequency(T, n, Boson) for n in -N + 1 : N - 1]
        return MatsubaraGrid(T, grid, :Boson, Linear; checks = false)
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

        grid_idxs_p = Vector{Int64}(undef, N)
        grid_idxs_m = Vector{Int64}(undef, N)
        offset      = δ - 1

        for n in 1 : N
            if n <= lb
                grid_idxs_p[n] = n - 1
                grid_idxs_m[n] = -grid_idxs_p[n]
            else 
                grid_idxs_p[n] = n - 1 + offset
                grid_idxs_m[n] = -grid_idxs_p[n]
                offset         = offset + δ - 1
            end 
        end

        grid = MatsubaraFrequency[MatsubaraFrequency(T, n, Boson) for n in vcat(reverse(grid_idxs_m[2 : end]), grid_idxs_p)]
        return MatsubaraGrid(T, grid, :Boson, Coarse; checks = false)
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

function Base.:length(
    grid :: MatsubaraGrid{GT}
    )    :: Int64 where {GT <: AbstractGrid}

    return length(grid.data)
end

function index_range(
    grid :: MatsubaraGrid{GT}
    )    :: NTuple{2, Int64} where {GT <: AbstractGrid}

    return index(grid.data[1]), index(grid.data[end])
end

function is_inbounds(
    w    :: MatsubaraFrequency,
    grid :: MatsubaraGrid{GT}
    )    :: Bool where {GT <: AbstractGrid}

    return value(grid[1]) <= value(w) <= value(grid[end])
end

function is_inbounds(
    w    :: Float64,
    grid :: MatsubaraGrid{GT}
    )    :: Bool where {GT <: AbstractGrid}

    return value(grid[1]) <= w <= value(grid[end])
end

function info(
    grid :: MatsubaraGrid{GT}
    )    :: Nothing where {GT <: AbstractGrid}

    println("MatsubaraGrid properties")
    println("------------------------")
    println("Grid type     : $(GT)")
    println("Particle type : $(type(grid))")
    println("Temperature   : $(temperature(grid))")
    println("Length        : $(length(grid))")
    println("Index range   : $(index_range(grid))")

    return nothing
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
    )    :: MatsubaraFrequency where {GT <: AbstractGrid}

    # bounds check performed by Base.Array
    return grid.data[idx]
end

function Base.:getindex(
    grid :: MatsubaraGrid{GT},
    idxs :: UnitRange{Int64}
    )    :: SubArray{MatsubaraFrequency, 1, Vector{MatsubaraFrequency}, Tuple{UnitRange{Int64}}, true} where {GT <: AbstractGrid}

    # bounds check performed by Base.Array 
    return @view grid.data[idxs]
end  



# make linear Matsubara grid callable with MatsubaraFrequency
# returns index to data array corresponding to this frequency if in grid
function (f :: MatsubaraGrid{Linear})(
    w :: MatsubaraFrequency
    ) :: Int64

    @assert type(w) == type(f) "Particle type must be equal between frequency and grid"
    @assert isapprox(temperature(w), temperature(f)) "Temperature must be equal between frequency and grid"

    if is_inbounds(w, f)
        return index(w) - index(f[1]) + 1
    else 
        error("Frequency not in grid")
    end 
end

# make coarse Matsubara grid callable with MatsubaraFrequency
# returns index to data array corresponding to this frequency if in grid
function (f :: MatsubaraGrid{Coarse})(
    w :: MatsubaraFrequency
    ) :: Int64

    @assert type(w) == type(f) "Particle type must be equal between frequency and grid"
    @assert isapprox(temperature(w), temperature(f)) "Temperature must be equal between frequency and grid"

    if is_inbounds(w, f)
        for i in 1 : length(f)
            if index(f[i]) == index(w)
                return i 
            end 
        end

        error("Frequency not in grid")
    else 
        error("Frequency not in grid")
    end 
end

# make linear Matsubara grid callable with Float64
# returns index to data array corresponding to closest frequency if in grid
function (f :: MatsubaraGrid{Linear})(
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

# make coarse Matsubara grid callable with Float64
# returns index to data array corresponding to closest frequency if in grid
function (f :: MatsubaraGrid{Coarse})(
    w :: Float64
    ) :: Int64

    if is_inbounds(w, f)
        dst = Inf
    
        for i in 1 : length(f) 
            dstp = abs(value(f[i]) - w)
    
            if dstp > dst 
                return i - 1
            else 
                dst = dstp 
            end 
        end

        return length(f)
    else 
        error("Frequency not in grid")
    end 
end



# make Matsubara grid iterable
function Base.:iterate(
    grid :: MatsubaraGrid{GT}
    )    :: Tuple{MatsubaraFrequency, Int64} where {GT <: AbstractGrid}

    return grid[1], 1 
end

function Base.:iterate(
    grid  :: MatsubaraGrid{GT},
    state :: Int64
    )     :: Union{Nothing, Tuple{MatsubaraFrequency, Int64}} where {GT <: AbstractGrid}

    if state < length(grid)
        return grid[state + 1], state + 1 
    else 
        return nothing 
    end
end