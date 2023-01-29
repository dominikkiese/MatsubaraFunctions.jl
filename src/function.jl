struct MatsubaraFunction{GD, SD, DD, GT <: AbstractGrid, Q <: Number}
    grids :: NTuple{GD, MatsubaraGrid{GT}}
    shape :: NTuple{SD, Int64}          
    data  :: Array{Q, DD}

    # safe constructor
    function MatsubaraFunction(
        grids  :: NTuple{GD, MatsubaraGrid{GT}}, 
        shape  :: NTuple{SD, Int64}, 
        data   :: Array{Q, DD}
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{GD, SD, DD, GT, Q} where {GD, SD, DD, GT <: AbstractGrid, Q <: Number}

        # error for integer data type
        if Q <: Integer || Q <: Complex{Int} error("Integer data type not supported") end

        # throw warning for GD > 3
        if GD > 3 @warn "Matsubara function not callable on hybercubic grids" end
        
        if checks
            # check dimensions
            @assert GD + SD == DD "Dimensions do not match"
            
            # check grids
            for g in grids
                @assert isapprox(temperature(g), temperature(grids[1])) "Grids must be defined for the same temperature"
                @assert issorted(g) "Grids must be sorted"
            end
        end

        return new{GD, SD, DD, GT, Q}(grids, shape, data)
    end

    # convenience constructors
    function MatsubaraFunction(
        grids  :: NTuple{GD, MatsubaraGrid{GT}},
        shape  :: NTuple{SD, Int64},
               :: Type{Q}
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{GD, SD, GD + SD, GT, Q} where {GD, SD, GT <: AbstractGrid, Q <: Number}

        dims = length.(grids)..., shape...
        data = Array{Q, GD + SD}(undef, dims)

        return MatsubaraFunction(grids, shape, data; checks)
    end

    function MatsubaraFunction(
        grid   :: MatsubaraGrid{GT},
        shape  :: NTuple{SD, Int64},
               :: Type{Q}
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{1, SD, 1 + SD, GT, Q} where {SD, GT <: AbstractGrid, Q <: Number}

        return MatsubaraFunction((grid,), shape, Q; checks)
    end

    function MatsubaraFunction(
        grids  :: NTuple{GD, MatsubaraGrid{GT}},
        shape  :: Int64,
               :: Type{Q}
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{GD, 1, GD + 1, GT, Q} where {GD, GT <: AbstractGrid, Q <: Number}

        return MatsubaraFunction(grids, (shape,), Q; checks)
    end

    function MatsubaraFunction(
        grid   :: MatsubaraGrid{GT},
        shape  :: Int64,
               :: Type{Q}
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{1, 1, 2, GT, Q} where {GT <: AbstractGrid, Q <: Number}

        return MatsubaraFunction((grid,), (shape,), Q; checks)
    end

    # fallback methods if Q is not provided
    function MatsubaraFunction(
        grids  :: NTuple{GD, MatsubaraGrid{GT}},
        shape  :: NTuple{SD, Int64}
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{GD, SD, GD + SD, GT, ComplexF64} where {GD, SD, GT <: AbstractGrid}

        return MatsubaraFunction(grids, shape, ComplexF64; checks)
    end

    function MatsubaraFunction(
        grid   :: MatsubaraGrid{GT},
        shape  :: NTuple{SD, Int64},
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{1, SD, 1 + SD, GT, ComplexF64} where {SD, GT <: AbstractGrid}

        return MatsubaraFunction(grid, shape, ComplexF64; checks)
    end

    function MatsubaraFunction(
        grids  :: NTuple{GD, MatsubaraGrid{GT}},
        shape  :: Int64,
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{GD, 1, GD + 1, GT, ComplexF64} where {GD, GT <: AbstractGrid}

        return MatsubaraFunction(grids, shape, ComplexF64; checks)
    end

    function MatsubaraFunction(
        grid   :: MatsubaraGrid{GT},
        shape  :: Int64,
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{1, 1, 2, GT, ComplexF64} where {GT <: AbstractGrid}

        return MatsubaraFunction(grid, shape, ComplexF64; checks)
    end
end



# getter functions 
function grids_shape(
    f :: MatsubaraFunction{GD, SD, DD, GT, Q}
    ) :: NTuple{GD, Int64} where {GD, SD, DD, GT <: AbstractGrid, Q <: Number}

    return length.(f.grids)
end

function grids_shape(
    f   :: MatsubaraFunction{GD, SD, DD, GT, Q},
    idx :: Int64
    )   :: Int64 where {GD, SD, DD, GT <: AbstractGrid, Q <: Number}

    return length(f.grids[idx])
end

function shape(
    f :: MatsubaraFunction{GD, SD, DD, GT, Q}
    ) :: NTuple{SD, Int64} where {GD, SD, DD, GT <: AbstractGrid, Q <: Number}

    return f.shape 
end 

function data_shape(
    f :: MatsubaraFunction{GD, SD, DD, GT, Q}
    ) :: NTuple{DD, Int64} where {GD, SD, DD, GT <: AbstractGrid, Q <: Number}

    return size(f.data)
end

function absmax(
    f :: MatsubaraFunction{GD, SD, DD, GT, Q}
    ) :: Float64 where {GD, SD, DD, GT <: AbstractGrid, Q <: Number}

    return maximum(abs.(f.data))
end

function argmax(
    f :: MatsubaraFunction{GD, SD, DD, GT, Q}
    ) :: Float64 where {GD, SD, DD, GT <: AbstractGrid, Q <: Number}

    return argmax(abs.(f.data))
end



# basic operations
function add(
    f1     :: MatsubaraFunction{GD, SD, DD, GT, Q}, 
    f2     :: MatsubaraFunction{GD, SD, DD, GT, Q}
    ;
    checks :: Bool = true
    )      :: MatsubaraFunction{GD, SD, DD, GT, Q} where {GD, SD, DD, GT <: AbstractGrid, Q <: Number}

    if checks
        for i in 1 : GD 
            @assert isapprox(f1.grids[i].data, f2.grids[i].data) "Grids must be equal for addition" 
        end
    end

    return MatsubaraFunction(f1.grids, f1.shape, f1.data .+ f2.data; checks)
end

function add!(
    f1     :: MatsubaraFunction{GD, SD, DD, GT, Q}, 
    f2     :: MatsubaraFunction{GD, SD, DD, GT, Q}
    ;
    checks :: Bool = true
    )      :: Nothing where {GD, SD, DD, GT <: AbstractGrid, Q <: Number}

    if checks
        for i in 1 : GD 
            @assert isapprox(f1.grids[i].data, f2.grids[i].data) "Grids must be equal for addition" 
        end 
    end

    f1.data .+= f2.data

    return nothing 
end

function subtract(
    f1     :: MatsubaraFunction{GD, SD, DD, GT, Q}, 
    f2     :: MatsubaraFunction{GD, SD, DD, GT, Q}
    ;
    checks :: Bool = true
    )      :: MatsubaraFunction{GD, SD, DD, GT, Q} where {GD, SD, DD, GT <: AbstractGrid, Q <: Number}

    if checks
        for i in 1 : GD 
            @assert isapprox(f1.grids[i].data, f2.grids[i].data) "Grids must be equal for subtraction" 
        end 
    end

    return MatsubaraFunction(f1.grids, f1.shape, f1.data .- f2.data; checks)
end

function subtract!(
    f1     :: MatsubaraFunction{GD, SD, DD, GT, Q}, 
    f2     :: MatsubaraFunction{GD, SD, DD, GT, Q}
    ;
    checks :: Bool = true
    )      :: Nothing where {GD, SD, DD, GT <: AbstractGrid, Q <: Number}

    if checks
        for i in 1 : GD 
            @assert isapprox(f1.grids[i].data, f2.grids[i].data) "Grids must be equal for subtraction" 
        end
    end 

    f1.data .-= f2.data

    return nothing 
end

function mult(
    f      :: MatsubaraFunction{GD, SD, DD, GT, Q},
    val    :: Qp
    ;
    checks :: Bool = true
    )      :: MatsubaraFunction{GD, SD, DD, GT, Q} where {GD, SD, DD, GT <: AbstractGrid, Q <: Number, Qp <: Number}

    # type promotion checked by Base.Array 
    return MatsubaraFunction(f.grids, f.shape, val .* f.data; checks)
end

function mult!(
    f   :: MatsubaraFunction{GD, SD, DD, GT, Q},
    val :: Qp
    )   :: Nothing where {GD, SD, DD, GT <: AbstractGrid, Q <: Number, Qp <: Number}

    # type promotion checked by Base.Array 
    f.data .*= val 

    return nothing
end
 
function set!(
    f   :: MatsubaraFunction{GD, SD, DD, GT, Q},
    val :: Qp,
    )   :: Nothing where {GD, SD, DD, GT <: AbstractGrid, Q <: Number, Qp <: Number}

    # type promotion checked by Base.Array 
    f.data .= val

    return nothing
end

function set!(
    f1     :: MatsubaraFunction{GD, SD, DD, GT, Q},
    f2     :: MatsubaraFunction{GD, SD, DD, GT, Q},
    ; 
    checks :: Bool = true
    )      :: Nothing where {GD, SD, DD, GT <: AbstractGrid, Q <: Number}

    if checks
        for i in 1 : GD 
            @assert isapprox(f1.grids[i].data, f2.grids[i].data) "Grids must be equal for overwrite" 
        end
    end 

    f1.data .= f2.data

    return nothing
end



# getindex methods
function Base.:getindex(
    f :: MatsubaraFunction{GD, SD, DD, GT, Q},
    x :: Vararg{Int64, DD}
    ) :: Q where {GD, SD, DD, GT <: AbstractGrid, Q <: Number}

    # bounds check performed by Base.Array
    return f.data[x...]
end

# setindex! methods
function Base.:setindex!(
    f   :: MatsubaraFunction{GD, SD, DD, GT, Q},
    val :: Qp,
    x   :: Vararg{Int64, DD}
    )   :: Nothing where {GD, SD, DD, GT <: AbstractGrid, Q <: Number, Qp <: Number}

    # bounds check performed by Base.Array
    f.data[x...] = val

    return nothing
end



# compute tail moments in quadratic approximation from upper bound of 1D MatsubaraFunction
@inline @inbounds function upper_tail_moments(
    f :: MatsubaraFunction{1, SD, DD, GT, Q},
    x :: Vararg{Int64, SD} 
    ) :: SVector{3, Q} where {SD, DD, GT <: AbstractGrid, Q <: Number}

    # read data
    ydat = SVector{3, Q}(f.data[end, x...], f.data[end - 1, x...], f.data[end - 2, x...])
    xdat = SVector{3, Float64}(1.0 / f.grids[1][end], 1.0 / f.grids[1][end - 1], 1.0 / f.grids[1][end - 2])
    
    # generate Vandermonde matrix 
    mat = @SMatrix Float64[1.0 xdat[1] xdat[1] * xdat[1];
                           1.0 xdat[2] xdat[2] * xdat[2];
                           1.0 xdat[3] xdat[3] * xdat[3]]
    
    return inv(mat) * ydat
end

# compute tail moments in quadratic approximation from lower bound of 1D MatsubaraFunction
@inline @inbounds function lower_tail_moments(
    f :: MatsubaraFunction{1, SD, DD, GT, Q},
    x :: Vararg{Int64, SD} 
    ) :: SVector{3, Q} where {SD, DD, GT <: AbstractGrid, Q <: Number}

    # read data
    ydat = SVector{3, Q}(f.data[1, x...], f.data[2, x...], f.data[3, x...])
    xdat = SVector{3, Float64}(1.0 / f.grids[1][1], 1.0 / f.grids[1][2], 1.0 / f.grids[1][3])
    
    # generate Vandermonde matrix 
    mat = @SMatrix Float64[1.0 xdat[1] xdat[1] * xdat[1];
                           1.0 xdat[2] xdat[2] * xdat[2];
                           1.0 xdat[3] xdat[3] * xdat[3]]
    
    return inv(mat) * ydat
end



# call to Matsubara function on 1D grid 
# performs extrapolaton using tail moments if extrp = true (default), else fallback to specified bc
@inbounds function (f :: MatsubaraFunction{1, SD, DD, GT, Q})(
    w     :: Float64,
    x     :: Vararg{Int64, SD} 
    ; 
    bc    :: Float64 = 0.0,
    extrp :: Bool    = true
    )     :: Q where{SD, DD, GT <: AbstractGrid, Q <: Number}

    ax = f.grids[1][1] <= w <= f.grids[1][end]

    if ax
        p = Param(w, f.grids[1]) 
        return p.wgts[1] * f.data[p.idxs[1], x...] + p.wgts[2] * f.data[p.idxs[2], x...]
    else 
        if extrp
            if sign(w) < 0.0 
                moments = lower_tail_moments(f, x...)
                return moments[1] + (moments[2] + moments[3] / w) / w
            else 
                moments = upper_tail_moments(f, x...)
                return moments[1] + (moments[2] + moments[3] / w) / w
            end
        else 
            return bc 
        end 
    end
end

# call to Matsubara function on 2D grid
@inbounds function (f :: MatsubaraFunction{2, SD, DD, GT, Q})(
    w  :: NTuple{2, Float64},
    x  :: Vararg{Int64, SD} 
    ; 
    bc :: Float64 = 0.0
    )  :: Q where{SD, DD, GT <: AbstractGrid, Q <: Number}

    ax1 = f.grids[1][1] <= w[1] <= f.grids[1][end]
    ax2 = f.grids[2][1] <= w[2] <= f.grids[2][end]

    if ax1 && ax2
        p1 = Param(w[1], f.grids[1])
        p2 = Param(w[2], f.grids[2])

        return p1.wgts[1] * p2.wgts[1] * f.data[p1.idxs[1], p2.idxs[1], x...] + 
               p1.wgts[1] * p2.wgts[2] * f.data[p1.idxs[1], p2.idxs[2], x...] + 
               p1.wgts[2] * p2.wgts[1] * f.data[p1.idxs[2], p2.idxs[1], x...] + 
               p1.wgts[2] * p2.wgts[2] * f.data[p1.idxs[2], p2.idxs[2], x...]         
    else 
        return bc 
    end 
end

# call to Matsubara function on 3D grid
@inbounds function (f :: MatsubaraFunction{3, SD, DD, GT, Q})(
    w  :: NTuple{3, Float64},
    x  :: Vararg{Int64, SD} 
    ; 
    bc :: Float64 = 0.0
    )  :: Q where{SD, DD, GT <: AbstractGrid, Q <: Number}

    ax1 = f.grids[1][1] <= w[1] <= f.grids[1][end]
    ax2 = f.grids[2][1] <= w[2] <= f.grids[2][end]
    ax3 = f.grids[3][1] <= w[3] <= f.grids[3][end]

    if ax1 && ax2 && ax3
        p1 = Param(w[1], f.grids[1])
        p2 = Param(w[2], f.grids[2])
        p3 = Param(w[3], f.grids[3])

        return p1.wgts[1] * p2.wgts[1] * p3.wgts[1] * f.data[p1.idxs[1], p2.idxs[1], p3.idxs[1], x...] + 
               p1.wgts[1] * p2.wgts[1] * p3.wgts[2] * f.data[p1.idxs[1], p2.idxs[1], p3.idxs[2], x...] + 
               p1.wgts[1] * p2.wgts[2] * p3.wgts[1] * f.data[p1.idxs[1], p2.idxs[2], p3.idxs[1], x...] + 
               p1.wgts[1] * p2.wgts[2] * p3.wgts[2] * f.data[p1.idxs[1], p2.idxs[2], p3.idxs[2], x...] + 
               p1.wgts[2] * p2.wgts[1] * p3.wgts[1] * f.data[p1.idxs[2], p2.idxs[1], p3.idxs[1], x...] + 
               p1.wgts[2] * p2.wgts[1] * p3.wgts[2] * f.data[p1.idxs[2], p2.idxs[1], p3.idxs[2], x...] + 
               p1.wgts[2] * p2.wgts[2] * p3.wgts[1] * f.data[p1.idxs[2], p2.idxs[2], p3.idxs[1], x...] + 
               p1.wgts[2] * p2.wgts[2] * p3.wgts[2] * f.data[p1.idxs[2], p2.idxs[2], p3.idxs[2], x...]  
    else 
        return bc 
    end 
end



# sum data of Matsubara function on 1D grid for real-valued data
@inline @inbounds function sum_me(
    f :: MatsubaraFunction{1, SD, DD, GT, Q},
    x :: Vararg{Int64, SD} 
    ) :: Q where{SD, DD, GT <: AbstractGrid, Q <: Real}

    slice = @view f.data[:, x...]
    sum   = 0.0

    @turbo for idx in eachindex(slice)
        sum += slice[idx]
    end 

    return sum
end

# sum data of Matsubara function on 1D grid for complex-valued data
@inline @inbounds function sum_me(
    f :: MatsubaraFunction{1, SD, DD, GT, Complex{Q}},
    x :: Vararg{Int64, SD} 
    ) :: Complex{Q} where{SD, DD, GT <: AbstractGrid, Q <: Real}

    slice  = reinterpret(Q, @view f.data[:, x...])
    vw_re  = @view slice[1 : 2 : end - 1]
    vw_im  = @view slice[2 : 2 : end]
    sum_re = 0.0
    sum_im = 0.0

    @turbo for idx in eachindex(vw_re)
        sum_re += vw_re[idx]
    end 

    @turbo for idx in eachindex(vw_im)
        sum_im += vw_im[idx]
    end 

    return sum_re + im * sum_im
end

# compute i * Matsubara sum over Matsubara function on 1D grid for real-valued data  
@inbounds function Base.:sum(
    f :: MatsubaraFunction{1, SD, DD, Linear, Q},
    x :: Vararg{Int64, SD}
    ) :: Complex{Q} where {SD, DD, Q <: Real}

    # compute tail moments 
    upper_moments = upper_tail_moments(f, x...)
    lower_moments = lower_tail_moments(f, x...)

    # check self-consistency 
    @assert norm(upper_moments .- lower_moments) / norm(upper_moments) < 1e-2 "Tail fits are inconsistent! Try more frequencies or check prerequisites."

    # compute expansion coefficients
    α0 = +0.5 * (upper_moments[1] + lower_moments[1])
    α1 = +0.5 * (upper_moments[2] + lower_moments[2]) * im
    α2 = -0.5 * (upper_moments[3] + lower_moments[3])

    # compute the Matsubara sum using quadratic asymptotic model
    T   = temperature(f.grids[1])
    num = grids_shape(f, 1)
    val = -T * (num * α0 - sum_me(f, x...)) - 0.5 * (α1 + 0.5 * α2 / T)
    sum = 0.0

    @turbo for w in 1 : num
        sum += T * α2 / f.grids[1].data[w] / f.grids[1].data[w]
    end

    return val + sum
end

# compute Matsubara sum over Matsubara function on 1D grid for complex-valued data
@inbounds function Base.:sum(
    f :: MatsubaraFunction{1, SD, DD, Linear, Complex{Q}},
    x :: Vararg{Int64, SD}
    ) :: Complex{Q} where {SD, DD, Q <: Real}

    # compute tail moments 
    upper_moments = upper_tail_moments(f, x...)
    lower_moments = lower_tail_moments(f, x...)

    # check self-consistency 
    @assert norm(upper_moments .- lower_moments) / norm(upper_moments) < 1e-2 "Tail fits are inconsistent! Try more frequencies or check prerequisites."

    # compute expansion coefficients
    α0 = +0.5 * (upper_moments[1] + lower_moments[1])
    α1 = +0.5 * (upper_moments[2] + lower_moments[2]) * im
    α2 = -0.5 * (upper_moments[3] + lower_moments[3])

    # compute the Matsubara sum using quadratic asymptotic model
    T   = temperature(f.grids[1])
    num = grids_shape(f, 1)
    val = -T * (num * α0 - sum_me(f, x...)) - 0.5 * (α1 + 0.5 * α2 / T)
    
    sum_re = 0.0; α2_re = real(α2)
    sum_im = 0.0; α2_im = imag(α2)

    @turbo for w in 1 : num
        val     = T / f.grids[1].data[w] / f.grids[1].data[w]
        sum_re += α2_re * val
        sum_im += α2_im * val
    end

    return val + sum_re + im * sum_im
end