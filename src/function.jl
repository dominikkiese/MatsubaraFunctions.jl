struct MatsubaraFunction{Dg, Ds, Dt, Q <: Number}
    grids :: NTuple{Dg, MatsubaraGrid}
    shape :: NTuple{Ds, Int64}          
    data  :: Array{Q, Dt}

    # safe constructor
    function MatsubaraFunction(
        grids  :: NTuple{Dg, MatsubaraGrid}, 
        shape  :: NTuple{Ds, Int64}, 
        data   :: Array{Q, Dt}
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{Dg, Ds, Dt, Q} where {Dg, Ds, Dt, Q <: Number}

        # error for integer data type
        if Q <: Integer || Q <: Complex{Int} error("Integer data type not supported") end

        # throw warning for Dg > 3
        if Dg > 3 @warn "Matsubara function not callable on hybercubic grids" end
        
        if checks
            # check dimensions
            @assert Dg + Ds == Dt "Dimensions do not match"
            
            # check grids
            for g in grids
                @assert isapprox(g.T, grids[1].T) "Grids must be defined for the same temperature"
                @assert issorted(g) "Grids must be sorted"
                delta = g[2] - g[1]

                for i in 2 : (length(g) - 1)
                    @assert isapprox(g[i + 1] - g[i], delta) "Grids must be linear"
                end
            end
        end

        return new{Dg, Ds, Dt, Q}(grids, shape, data)
    end

    # convenience constructor
    function MatsubaraFunction(
        grids  :: NTuple{Dg, MatsubaraGrid},
        shape  :: NTuple{Ds, Int64},
               :: Type{Q}
        ;
        checks :: Bool = true
        )      :: MatsubaraFunction{Dg, Ds, Dg + Ds, Q} where {Dg, Ds, Q <: Number}

        dims = length.(grids)..., shape...
        data = Array{Q, Dg + Ds}(undef, dims)

        return MatsubaraFunction(grids, shape, data; checks)
    end
end



# getter functions 
function grids_shape(
    f :: MatsubaraFunction{Dg, Ds, Dt, Q}
    ) :: NTuple{Dg, Int64} where {Dg, Ds, Dt, Q <: Number}

    return length.(f.grids)
end

function grids_shape(
    f   :: MatsubaraFunction{Dg, Ds, Dt, Q},
    idx :: Int64
    )   :: Int64 where {Dg, Ds, Dt, Q <: Number}

    return length(f.grids[idx])
end

function shape(
    f :: MatsubaraFunction{Dg, Ds, Dt, Q}
    ) :: NTuple{Ds, Int64} where {Dg, Ds, Dt, Q <: Number}

    return f.shape 
end 

function data_shape(
    f :: MatsubaraFunction{Dg, Ds, Dt, Q}
    ) :: NTuple{Dt, Int64} where {Dg, Ds, Dt, Q <: Number}

    return size(f.data)
end



# basic addition
function add(
    f1     :: MatsubaraFunction{Dg, Ds, Dt, Q}, 
    f2     :: MatsubaraFunction{Dg, Ds, Dt, Q}
    ;
    checks :: Bool = true
    )      :: MatsubaraFunction{Dg, Ds, Dt, Q} where {Dg, Ds, Dt, Q <: Number}

    if checks
        for i in 1 : Dg 
            @assert isapprox(f1.grids[i].data, f2.grids[i].data) "Grids must be equal for addition" 
        end
    end

    return MatsubaraFunction(f1.grids, f1.shape, f1.data .+ f2.data; checks)
end

function add!(
    f1     :: MatsubaraFunction{Dg, Ds, Dt, Q}, 
    f2     :: MatsubaraFunction{Dg, Ds, Dt, Q}
    ;
    checks :: Bool = true
    )      :: Nothing where {Dg, Ds, Dt, Q <: Number}

    if checks
        for i in 1 : Dg 
            @assert isapprox(f1.grids[i].data, f2.grids[i].data) "Grids must be equal for addition" 
        end 
    end

    f1.data .+= f2.data

    return nothing 
end

# basic subtraction
function subtract(
    f1     :: MatsubaraFunction{Dg, Ds, Dt, Q}, 
    f2     :: MatsubaraFunction{Dg, Ds, Dt, Q}
    ;
    checks :: Bool = true
    )      :: MatsubaraFunction{Dg, Ds, Dt, Q} where {Dg, Ds, Dt, Q <: Number}

    if checks
        for i in 1 : Dg 
            @assert isapprox(f1.grids[i].data, f2.grids[i].data) "Grids must be equal for subtraction" 
        end 
    end

    return MatsubaraFunction(f1.grids, f1.shape, f1.data .- f2.data; checks)
end

function subtract!(
    f1     :: MatsubaraFunction{Dg, Ds, Dt, Q}, 
    f2     :: MatsubaraFunction{Dg, Ds, Dt, Q}
    ;
    checks :: Bool = true
    )      :: Nothing where {Dg, Ds, Dt, Q <: Number}

    if checks
        for i in 1 : Dg 
            @assert isapprox(f1.grids[i].data, f2.grids[i].data) "Grids must be equal for subtraction" 
        end
    end 

    f1.data .-= f2.data

    return nothing 
end

# basic multiplication with scalar 
function mult(
    f      :: MatsubaraFunction{Dg, Ds, Dt, Q},
    val    :: Qp
    ;
    checks :: Bool = true
    )      :: MatsubaraFunction{Dg, Ds, Dt, Q} where {Dg, Ds, Dt, Q <: Number, Qp <: Number}

    return MatsubaraFunction(f.grids, f.shape, val .* f.data; checks)
end

function mult!(
    f   :: MatsubaraFunction{Dg, Ds, Dt, Q},
    val :: Qp
    )   :: Nothing where {Dg, Ds, Dt, Q <: Number, Qp <: Number}

    f.data .*= val 

    return nothing
end



# getindex method (bounds check performed by Base.Array)
function Base.:getindex(
    f :: MatsubaraFunction{Dg, Ds, Dt, Q},
    x :: Vararg{Int64, Dt}
    ) :: Q where {Dg, Ds, Dt, Q <: Number}

    return f.data[x...]
end

# setindex! method (bounds check performed by Base.Array)
function Base.:setindex!(
    f   :: MatsubaraFunction{Dg, Ds, Dt, Q},
    val :: Qp,
    x   :: Vararg{Int64, Dt}
    )   :: Nothing where {Dg, Ds, Dt, Q <: Number, Qp <: Number}

    f.data[x...] = val

    return nothing
end



# call to Matsubara function on 1D grid (allows to pass tail moments for extrapolation)
@inbounds function (f :: MatsubaraFunction{1, Ds, Dt, Q})(
    w    :: Float64,
    x    :: Vararg{Int64, Ds} 
    ; 
    bc   :: Float64 = 0.0,
    tail :: NTuple{P, Float64} = ()
    )    :: Q where{Ds, Dt, Q <: Number, P}

    bv_dn = f.grids[1][1]
    bv_up = f.grids[1][grids_shape(f, 1)]

    if w < bv_dn 
        val = bc

        for n in eachindex(tail)
            val += tail[n] * (bv_dn / w)^n * f.data[1, x...]
        end
            
        return val 
        
    elseif w > bv_up 
        val = bc

        for n in eachindex(tail)
            val += tail[n] * (bv_up / w)^n * f.data[end, x...]
        end
            
        return val 

    else
        p = Param(w, f.grids[1]) 

        return p.wgts[1] * f.data[p.idxs[1], x...] + p.wgts[2] * f.data[p.idxs[2], x...]
    end 
end

# call to Matsubara function on 2D grid
@inbounds function (f :: MatsubaraFunction{2, Ds, Dt, Q})(
    w  :: NTuple{2, Float64},
    x  :: Vararg{Int64, Ds} 
    ; 
    bc :: Float64 = 0.0
    )  :: Q where{Ds, Dt, Q <: Number}

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
@inbounds function (f :: MatsubaraFunction{3, Ds, Dt, Q})(
    w  :: NTuple{3, Float64},
    x  :: Vararg{Int64, Ds} 
    ; 
    bc :: Float64 = 0.0
    )  :: Q where{Ds, Dt, Q <: Number}

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