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
                @assert temperature(g) ≈ temperature(grids[1]) "Grids must be defined for the same temperature"
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



# basic operations
function add(
    f1     :: MatsubaraFunction{GD, SD, DD, Q}, 
    f2     :: MatsubaraFunction{GD, SD, DD, Q}
    ;
    checks :: Bool = true
    )      :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

    if checks
        for i in 1 : GD 
            @assert Float64[value(v) for v in f1.grids[i]] ≈ Float64[value(v) for v in f2.grids[i]] "Grids must be equal for addition" 
        end
    end

    return MatsubaraFunction(f1.grids, f1.shape, f1.data .+ f2.data; checks)
end

function add!(
    f1     :: MatsubaraFunction{GD, SD, DD, Q}, 
    f2     :: MatsubaraFunction{GD, SD, DD, Q}
    ;
    checks :: Bool = true
    )      :: Nothing where {GD, SD, DD, Q <: Number}

    if checks
        for i in 1 : GD 
            @assert Float64[value(v) for v in f1.grids[i]] ≈ Float64[value(v) for v in f2.grids[i]] "Grids must be equal for addition" 
        end 
    end

    f1.data .+= f2.data

    return nothing 
end

function subtract(
    f1     :: MatsubaraFunction{GD, SD, DD, Q}, 
    f2     :: MatsubaraFunction{GD, SD, DD, Q}
    ;
    checks :: Bool = true
    )      :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number}

    if checks
        for i in 1 : GD 
            @assert Float64[value(v) for v in f1.grids[i]] ≈ Float64[value(v) for v in f2.grids[i]] "Grids must be equal for subtraction" 
        end 
    end

    return MatsubaraFunction(f1.grids, f1.shape, f1.data .- f2.data; checks)
end

function subtract!(
    f1     :: MatsubaraFunction{GD, SD, DD, Q}, 
    f2     :: MatsubaraFunction{GD, SD, DD, Q}
    ;
    checks :: Bool = true
    )      :: Nothing where {GD, SD, DD, Q <: Number}

    if checks
        for i in 1 : GD 
            @assert Float64[value(v) for v in f1.grids[i]] ≈ Float64[value(v) for v in f2.grids[i]] "Grids must be equal for subtraction" 
        end
    end 

    f1.data .-= f2.data

    return nothing 
end

function mult(
    f      :: MatsubaraFunction{GD, SD, DD, Q},
    val    :: Qp
    ;
    checks :: Bool = true
    )      :: MatsubaraFunction{GD, SD, DD, Q} where {GD, SD, DD, Q <: Number, Qp <: Number}

    # type promotion checked by Base.Array 
    return MatsubaraFunction(f.grids, f.shape, val .* f.data; checks)
end

function mult!(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    val :: Qp
    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    # type promotion checked by Base.Array 
    f.data .*= val 

    return nothing
end
 
function set!(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    val :: Qp,
    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    # type promotion checked by Base.Array 
    f.data .= val

    return nothing
end

function set!(
    f1     :: MatsubaraFunction{GD, SD, DD, Q},
    f2     :: MatsubaraFunction{GD, SD, DD, Q},
    ; 
    checks :: Bool = true
    )      :: Nothing where {GD, SD, DD, Q <: Number}

    if checks
        for i in 1 : GD 
            @assert Float64[value(v) for v in f1.grids[i]] ≈ Float64[value(v) for v in f2.grids[i]] "Grids must be equal for overwrite" 
        end
    end 

    f1.data .= f2.data

    return nothing
end



# getindex methods
function Base.:getindex(
    f :: MatsubaraFunction{GD, SD, DD, Q},
    x :: Vararg{Int64, DD}
    ) :: Q where {GD, SD, DD, Q <: Number}

    # bounds check performed by Base.Array
    return f.data[x...]
end

# setindex! methods
function Base.:setindex!(
    f   :: MatsubaraFunction{GD, SD, DD, Q},
    val :: Qp,
    x   :: Vararg{Int64, DD}
    )   :: Nothing where {GD, SD, DD, Q <: Number, Qp <: Number}

    # bounds check performed by Base.Array
    f.data[x...] = val

    return nothing
end



# compute tail moments in quadratic approximation from upper bound of 1D MatsubaraFunction
function upper_tail_moments(
    f :: MatsubaraFunction{1, SD, DD, Q},
    x :: Vararg{Int64, SD} 
    ) :: SVector{3, Q} where {SD, DD, Q <: Number}

    # read data
    ydat = SVector{3, Q}(f.data[end, x...], f.data[end - 1, x...], f.data[end - 2, x...])
    xdat = SVector{3, Float64}(1.0 / value(f.grids[1][end]), 1.0 / value(f.grids[1][end - 1]), 1.0 / value(f.grids[1][end - 2]))
    
    # generate Vandermonde matrix 
    mat = @SMatrix Float64[1.0 xdat[1] xdat[1] * xdat[1];
                           1.0 xdat[2] xdat[2] * xdat[2];
                           1.0 xdat[3] xdat[3] * xdat[3]]
    
    return inv(mat) * ydat
end

# compute tail moments in quadratic approximation from lower bound of 1D MatsubaraFunction
function lower_tail_moments(
    f :: MatsubaraFunction{1, SD, DD, Q},
    x :: Vararg{Int64, SD} 
    ) :: SVector{3, Q} where {SD, DD, Q <: Number}

    # read data
    ydat = SVector{3, Q}(f.data[1, x...], f.data[2, x...], f.data[3, x...])
    xdat = SVector{3, Float64}(1.0 / value(f.grids[1][1]), 1.0 / value(f.grids[1][2]), 1.0 / value(f.grids[1][3]))
    
    # generate Vandermonde matrix 
    mat = @SMatrix Float64[1.0 xdat[1] xdat[1] * xdat[1];
                           1.0 xdat[2] xdat[2] * xdat[2];
                           1.0 xdat[3] xdat[3] * xdat[3]]
    
    return inv(mat) * ydat
end

# extrapolate 1D Matsubara function using tail moments
function extrapolate(
    f :: MatsubaraFunction{1, SD, DD, Q},
    w :: Float64,
    x :: Vararg{Int64, SD} 
    ) :: Q where {SD, DD, Q <: Number}

    if sign(w) < 0.0 
        moments = lower_tail_moments(f, x...)
        return moments[1] + (moments[2] + moments[3] / w) / w
    else 
        moments = upper_tail_moments(f, x...)
        return moments[1] + (moments[2] + moments[3] / w) / w
    end
end



# call to MatsubaraFunction with MatsubaraFrequencies
@inline function (f :: MatsubaraFunction{GD, SD, DD, Q})(
    w     :: NTuple{GD, MatsubaraFrequency},
    x     :: Vararg{Int64, SD} 
    ; 
    bc    :: Function = x -> 0.0,
    extrp :: Bool     = false
    )     :: Q where{GD, SD, DD, Q <: Number}

    if any(ntuple(i -> !is_inbounds(w[i], f.grids[i]), GD))
        if GD == 1 && extrp 
            return extrapolate(f, value(w[1]), x...)
        else 
            return bc(w) 
        end 
    end

    idxs = ntuple(i -> grid_index(w[i], f.grids[i]), GD)
    return f[idxs..., x ...]
end

function (f :: MatsubaraFunction{1, SD, DD, Q})(
    w     :: MatsubaraFrequency,
    x     :: Vararg{Int64, SD} 
    ; 
    bc    :: Function = x -> 0.0,
    extrp :: Bool     = false
    )     :: Q where{SD, DD, Q <: Number}

    return f((w,), x...; bc = x -> bc(x[1]), extrp)
end



# call to MatsubaraFunction with Float64
@inline function (f :: MatsubaraFunction{GD, SD, DD, Q})(
    w     :: NTuple{GD, Float64},
    x     :: Vararg{Int64, SD} 
    ; 
    bc    :: Function = x -> 0.0,
    extrp :: Bool     = false
    )     :: Q where{GD, SD, DD, Q <: Number}

    if any(ntuple(i -> !is_inbounds(w[i], f.grids[i]), GD))
        if GD == 1 && extrp 
            return extrapolate(f, w[1], x...)
        else 
            return bc(w)
        end 
    end 

    p   = ntuple(i -> Param(w[i], f.grids[i]), GD)
    val = 0.0

    for cidx in CartesianIndices(ntuple(i -> 2, GD))
        wgts  = ntuple(i -> p[i].wgts[cidx[i]], GD)
        idxs  = ntuple(i -> p[i].idxs[cidx[i]], GD)
        val  += prod(wgts) * f[idxs..., x...]
    end

    return val
end

function (f :: MatsubaraFunction{1, SD, DD, Q})(
    w     :: Float64,
    x     :: Vararg{Int64, SD} 
    ; 
    bc    :: Function = x -> 0.0,
    extrp :: Bool     = false
    )     :: Q where{SD, DD, Q <: Number}

    return f((w,), x...; bc = x -> bc(x[1]), extrp)
end



# compute Matsubara sum for MatsubaraFunction on 1D grid
# Note: for real valued data the return value is complex (im * sum), 
#       i.e. we cannot generally infer it from Q -> dynamic dispatch
function sum_me(
    f :: MatsubaraFunction{1, SD, DD, Q},
    x :: Vararg{Int64, SD}
    ) where {SD, DD, Q <: Number}

    # compute tail moments 
    upper_moments = upper_tail_moments(f, x...)
    lower_moments = lower_tail_moments(f, x...)

    # check self-consistency 
    Δ     = norm(upper_moments .- lower_moments)
    scale = 1e-3 + max(norm(lower_moments), norm(upper_moments)) * 1e-2
    err   = Δ / scale

    @assert err <= 1.0 "Tail fits are inconsistent! Try more frequencies or check prerequisites"

    # compute expansion coefficients
    α0 = +0.5 * (upper_moments[1] + lower_moments[1])
    α1 = +0.5 * (upper_moments[2] + lower_moments[2]) * im
    α2 = -0.5 * (upper_moments[3] + lower_moments[3])

    # compute the Matsubara sum using quadratic asymptotic model
    T   = temperature(f.grids[1])
    num = grids_shape(f, 1)
    val = -T * (num * α0 - sum(@view f.data[:, x...])) - 0.5 * (α1 + 0.5 * α2 / T)

    for w in 1 : num
        val += T * α2 / value(f.grids[1][w]) / value(f.grids[1][w])
    end

    return val
end