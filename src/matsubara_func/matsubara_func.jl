struct MatsubaraFunction{GD, SD, DD, Q <: Number}
    grids :: NTuple{GD, MatsubaraGrid}
    shape :: NTuple{SD, Int64}          
    data  :: Array{Q, DD}

    # basic constructor
    function MatsubaraFunction(
        grids  :: NTuple{GD, MatsubaraGrid}, 
        shape  :: NTuple{SD, Int64}, 
        data   :: Array{Q, DD}
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

        # error for integer data type
        if Q <: Integer || Q <: Complex{Int} error("Integer data type not supported") end
        
        if checks
            # check dimensions
            @assert GD + SD == DD "Dimensions do not match"
            
            # check grids
            for g in grids
                @assert temperature(g) ≈ temperature(grids[1]) "Grids must have same temperature"
                @assert issorted(Float64[value(w) for w in g]) "Grids must be sorted"
            end
        end

        return new{GD, SD, DD, Q}(grids, shape, data)
    end

    # convenience constructors
    function MatsubaraFunction(
        grids  :: NTuple{GD, MatsubaraGrid},
        shape  :: NTuple{SD, Int64},
               :: Type{Q}
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{GD, SD, GD + SD, Q} where {GD, SD, Q <: Number}

        dims = length.(grids)..., shape...
        data = Array{Q, GD + SD}(undef, dims)

        return MatsubaraFunction(grids, shape, data; checks)
    end

    function MatsubaraFunction(
        grid   :: MatsubaraGrid,
        shape  :: NTuple{SD, Int64},
               :: Type{Q}
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{1, SD, 1 + SD, Q} where {SD, Q <: Number}

        return MatsubaraFunction((grid,), shape, Q; checks)
    end

    function MatsubaraFunction(
        grids  :: NTuple{GD, MatsubaraGrid},
        shape  :: Int64,
               :: Type{Q}
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{GD, 1, GD + 1, Q} where {GD, Q <: Number}

        return MatsubaraFunction(grids, (shape,), Q; checks)
    end

    function MatsubaraFunction(
        grid   :: MatsubaraGrid,
        shape  :: Int64,
               :: Type{Q}
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{1, 1, 2, Q} where {Q <: Number}

        return MatsubaraFunction((grid,), (shape,), Q; checks)
    end

    # fallback methods if data type Q is not provided
    function MatsubaraFunction(
        grids  :: NTuple{GD, MatsubaraGrid},
        shape  :: NTuple{SD, Int64}
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{GD, SD, GD + SD, ComplexF64} where {GD, SD}

        return MatsubaraFunction(grids, shape, ComplexF64; checks)
    end

    function MatsubaraFunction(
        grid   :: MatsubaraGrid,
        shape  :: NTuple{SD, Int64},
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{1, SD, 1 + SD, ComplexF64} where {SD}

        return MatsubaraFunction(grid, shape, ComplexF64; checks)
    end

    function MatsubaraFunction(
        grids  :: NTuple{GD, MatsubaraGrid},
        shape  :: Int64,
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{GD, 1, GD + 1, ComplexF64} where {GD}

        return MatsubaraFunction(grids, shape, ComplexF64; checks)
    end

    function MatsubaraFunction(
        grid   :: MatsubaraGrid,
        shape  :: Int64,
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{1, 1, 2, ComplexF64}

        return MatsubaraFunction(grid, shape, ComplexF64; checks)
    end
end



# getter functions 
function grids_shape(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: NTuple{GD, Int64} where {GD, SD, DD, Q <: Number}

    return length.(f.grids)
end

function grids_shape(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    idx :: Int64
    )   :: Int64 where {GD, SD, DD, Q <: Number}

    return length(f.grids[idx])
end

function shape(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: NTuple{SD, Int64} where {GD, SD, DD, Q <: Number}

    return f.shape 
end 

function shape(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    idx :: Int64
    )   :: Int64 where {GD, SD, DD, Q <: Number}

    return f.shape[idx]
end 

function data_shape(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: NTuple{DD, Int64} where {GD, SD, DD, Q <: Number}

    return size(f.data)
end

function data_shape(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    idx :: Int64
    )   :: Int64 where {GD, SD, DD, Q <: Number}

    return size(f.data, idx)
end

function absmax(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: Float64 where {GD, SD, DD, Q <: Number}

    return maximum(abs.(f.data))
end

function argmax(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: Float64 where {GD, SD, DD, Q <: Number}

    return argmax(abs.(f.data))
end

function info(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: Nothing where {GD, SD, DD, Q <: Number}

    println("MatsubaraFunction properties")
    println("----------------------------")
    println("Temperature     : $(temperature(f.grids[1]))")
    println("Grid dimension  : $(GD)")
    println("Shape dimension : $(SD)")
    println("Data dimension  : $(DD)")
    println("Data type       : $(Q)")

    return nothing
end



# load methods 
include("matsubara_func_ops.jl")
include("matsubara_func_index.jl")
include("matsubara_func_eval.jl")
include("matsubara_func_symmetries.jl")