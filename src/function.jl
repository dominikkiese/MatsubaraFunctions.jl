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
                @assert isapprox(g.T, grids[1].T) "Grids must be defined for the same temperature"
                @assert issorted(g) "Grids must be sorted"
            end
        end

        return new{GD, SD, DD, GT, Q}(grids, shape, data)
    end

    # convenience constructor
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

    # fallback method if Type{Q} is not provided
    function MatsubaraFunction(
        grids  :: NTuple{GD, MatsubaraGrid{GT}},
        shape  :: NTuple{SD, Int64}
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{GD, SD, GD + SD, GT, ComplexF64} where {GD, SD, GT <: AbstractGrid}

        return MatsubaraFunction(grids, shape, ComplexF64; checks)
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



# basic addition
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

# basic subtraction
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

# basic multiplication with scalar 
function mult(
    f      :: MatsubaraFunction{GD, SD, DD, GT, Q},
    val    :: Qp
    ;
    checks :: Bool = true
    )      :: MatsubaraFunction{GD, SD, DD, GT, Q} where {GD, SD, DD, GT <: AbstractGrid, Q <: Number, Qp <: Number}

    return MatsubaraFunction(f.grids, f.shape, val .* f.data; checks)
end

function mult!(
    f   :: MatsubaraFunction{GD, SD, DD, GT, Q},
    val :: Qp
    )   :: Nothing where {GD, SD, DD, GT <: AbstractGrid, Q <: Number, Qp <: Number}

    f.data .*= val 

    return nothing
end



# getindex methods (bounds check performed by Base.Array)
function Base.:getindex(
    f :: MatsubaraFunction{GD, SD, DD, GT, Q},
    x :: Vararg{Int64, DD}
    ) :: Q where {GD, SD, DD, GT <: AbstractGrid, Q <: Number}

    return f.data[x...]
end

# setindex! method (bounds check performed by Base.Array)
function Base.:setindex!(
    f   :: MatsubaraFunction{GD, SD, DD, GT, Q},
    val :: Qp,
    x   :: Vararg{Int64, DD}
    )   :: Nothing where {GD, SD, DD, GT <: AbstractGrid, Q <: Number, Qp <: Number}

    f.data[x...] = val

    return nothing
end



# call to Matsubara function on 1D grid (allows to pass tail moments for extrapolation)
@inbounds function (f :: MatsubaraFunction{1, SD, DD, GT, Q})(
    w    :: Float64,
    x    :: Vararg{Int64, SD} 
    ; 
    bc   :: Float64 = 0.0,
    tail :: NTuple{P, Float64} = ()
    )    :: Q where{SD, DD, P, GT <: AbstractGrid, Q <: Number}

    bv_dn = f.grids[1][1]
    bv_up = f.grids[1][grids_shape(f, 1)]

    if w < bv_dn 
        val = bc

        for n in eachindex(tail)
            val += tail[n] * (bv_dn / w)^n * (f.data[1, x...] - bc)
        end
            
        return val

    elseif w > bv_up 
        val = bc

        for n in eachindex(tail)
            val += tail[n] * (bv_up / w)^n * (f.data[end, x...] - bc)
        end
            
        return val
         
    else
        p = Param(w, f.grids[1]) 
        return p.wgts[1] * f.data[p.idxs[1], x...] + p.wgts[2] * f.data[p.idxs[2], x...]
    end 
end

# call to Matsubara function on 2D grid
@inbounds function (f :: MatsubaraFunction{2, SD, DD, GT, Q})(
    w  :: NTuple{2, Float64},
    x  :: Vararg{Int64, SD} 
    ; 
    bc :: Float64 = 0.0
    )  :: Q where{SD, DD, GT <: AbstractGrid, Q <: Number}

    ax1 = f.grids[1][1] <= w[1] <= f.grids[1][grids_shape(f, 1)]
    ax2 = f.grids[2][1] <= w[2] <= f.grids[2][grids_shape(f, 2)]

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

    ax1 = f.grids[1][1] <= w[1] <= f.grids[1][grids_shape(f, 1)]
    ax2 = f.grids[2][1] <= w[2] <= f.grids[2][grids_shape(f, 2)]
    ax3 = f.grids[3][1] <= w[3] <= f.grids[3][grids_shape(f, 3)]

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
@inbounds function sum_me(
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
@inbounds function sum_me(
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

# compute Matsubara sum over Matsubara function on 1D grid for real-valued data
@inbounds function Base.:sum(
    f :: MatsubaraFunction{1, SD, DD, Linear, Q},
    x :: Vararg{Int64, SD}
    ) :: Complex{Q} where {SD, DD, Q <: Real}

    # read data
    ydat = SVector{3, Q}(f.data[end, x...], f.data[end - 1, x...], f.data[end - 2, x...])
    xdat = SVector{3, Float64}(1.0 / f.grids[1][end], 1.0 / f.grids[1][end - 1], 1.0 / f.grids[1][end - 2])

    # generate Vandermonde matrix 
    mat = @SMatrix Float64[1.0 xdat[1] xdat[1] * xdat[1];
                           1.0 xdat[2] xdat[2] * xdat[2];
                           1.0 xdat[3] xdat[3] * xdat[3]]

    # compute fitting coefficients 
    coeff = inv(mat) * ydat
    α0    = coeff[1]
    α1    = im * coeff[2]
    α2    = -coeff[3]

    # compute the Matsubara sum using quadratic asymptotic model
    T   = f.grids[1].T
    num = grids_shape(f, 1)
    val = -T * (num * α0 - sum_me(f, x...)) + 0.5 * (α1 - 0.5 * α2 / T)
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

    # read data
    ydat = SVector{3, Complex{Q}}(f.data[end, x...], f.data[end - 1, x...], f.data[end - 2, x...])
    xdat = SVector{3, Float64}(1.0 / f.grids[1][end], 1.0 / f.grids[1][end - 1], 1.0 / f.grids[1][end - 2])

    # generate Vandermonde matrix 
    mat = @SMatrix Float64[1.0 xdat[1] xdat[1] * xdat[1];
                           1.0 xdat[2] xdat[2] * xdat[2];
                           1.0 xdat[3] xdat[3] * xdat[3]]

    # compute fitting coefficients 
    coeff = inv(mat) * ydat
    α0    = coeff[1]
    α1    = im * coeff[2]
    α2    = -coeff[3]

    # compute the Matsubara sum using quadratic asymptotic model
    T   = f.grids[1].T
    num = grids_shape(f, 1)
    val = -T * (num * α0 - sum_me(f, x...)) + 0.5 * (α1 - 0.5 * α2 / T)
    
    sum_re = 0.0; α2_re = real(α2)
    sum_im = 0.0; α2_im = imag(α2)

    @turbo for w in 1 : num
        val     = T / f.grids[1].data[w] / f.grids[1].data[w]
        sum_re += α2_re * val
        sum_im += α2_im * val
    end

    return val + sum_re + im * sum_im
end