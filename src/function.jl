struct MatsubaraFunction{Dg, Ds, Dt}
    grids :: NTuple{Dg, MatsubaraGrid}
    shape :: NTuple{Ds, Int64}          
    data  :: Array{Float64, Dt}

    # safe constructor
    function MatsubaraFunction(
        grids :: NTuple{Dg, MatsubaraGrid}, 
        shape :: NTuple{Ds, Int64}, 
        data  :: Array{Float64, Dt}
        )     :: MatsubaraFunction{Dg, Ds, Dt} where {Dg, Ds, Dt}
        
        # throw warning for Dg > 3
        if Dg > 3 @warn "Matsubara function not callable on hybercubic grids" end

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

        return new{Dg, Ds, Dt}(grids, shape, data)
    end

    # convenience constructor
    function MatsubaraFunction(
        grids :: NTuple{Dg, MatsubaraGrid},
        shape :: NTuple{Ds, Int64}
        )     :: MatsubaraFunction{Dg, Ds, Dg + Ds} where {Dg, Ds}

        dims = length.(grids)..., shape...
        data = Array{Float64, Dg + Ds}(undef, dims)

        return MatsubaraFunction(grids, shape, data)
    end
end



# getter functions 
function grids_shape(
    f :: MatsubaraFunction{Dg, Ds, Dt}
    ) :: NTuple{Dg, Int64} where {Dg, Ds, Dt}

    return length.(f.grids)
end

function grids_shape(
    f   :: MatsubaraFunction{Dg, Ds, Dt},
    idx :: Int64
    )   :: Int64 where {Dg, Ds, Dt}

    return length(f.grids[idx])
end

function shape(
    f :: MatsubaraFunction{Dg, Ds, Dt}
    ) :: NTuple{Ds, Int64} where {Dg, Ds, Dt}

    return f.shape 
end 

function data_shape(
    f :: MatsubaraFunction{Dg, Ds, Dt}
    ) :: NTuple{Dt, Int64} where {Dg, Ds, Dt}

    return size(f.data)
end



# basic addition
function Base.:+(
    f1 :: MatsubaraFunction{Dg, Ds, Dt}, 
    f2 :: MatsubaraFunction{Dg, Ds, Dt}
    )  :: MatsubaraFunction{Dg, Ds, Dt} where {Dg, Ds, Dt}

    for i in 1 : Dg 
        @assert isapprox(f1.grids[i].data, f2.grids[i].data) "Grids must be equal for addition" 
    end

    return MatsubaraFunction(f1.grids, f1.shape, @turbo f1.data .+ f2.data)
end

function add!(
    f1 :: MatsubaraFunction{Dg, Ds, Dt}, 
    f2 :: MatsubaraFunction{Dg, Ds, Dt}
    )  :: Nothing where {Dg, Ds, Dt}

    for i in 1 : Dg 
        @assert isapprox(f1.grids[i].data, f2.grids[i].data) "Grids must be equal for addition" 
    end 

    @turbo f1.data .+= f2.data

    return nothing 
end

# basic subtraction
function Base.:-(
    f1 :: MatsubaraFunction{Dg, Ds, Dt}, 
    f2 :: MatsubaraFunction{Dg, Ds, Dt}
    )  :: MatsubaraFunction{Dg, Ds, Dt} where {Dg, Ds, Dt}

    for i in 1 : Dg 
        @assert isapprox(f1.grids[i].data, f2.grids[i].data) "Grids must be equal for subtraction" 
    end 

    return MatsubaraFunction(f1.grids, f1.shape, @turbo f1.data .- f2.data)
end

function subtract!(
    f1 :: MatsubaraFunction{Dg, Ds, Dt}, 
    f2 :: MatsubaraFunction{Dg, Ds, Dt}
    )  :: Nothing where {Dg, Ds, Dt}

    for i in 1 : Dg 
        @assert isapprox(f1.grids[i].data, f2.grids[i].data) "Grids must be equal for subtraction" 
    end 

    @turbo f1.data .-= f2.data

    return nothing 
end

# basic multiplication with scalar 
function Base.:*(
    val :: Float64, 
    f   :: MatsubaraFunction{Dg, Ds, Dt}
    )   :: MatsubaraFunction{Dg, Ds, Dt} where {Dg, Ds, Dt}

    return MatsubaraFunction(f.grids, f.shape, @turbo val .* f.data)
end

function Base.:*(
    f   :: MatsubaraFunction{Dg, Ds, Dt},
    val :: Float64
    )   :: MatsubaraFunction{Dg, Ds, Dt} where {Dg, Ds, Dt}

    return MatsubaraFunction(f.grids, f.shape, @turbo val .* f.data)
end

function mult!(
    f   :: MatsubaraFunction{Dg, Ds, Dt},
    val :: Float64
    )   :: Nothing where {Dg, Ds, Dt}

    @turbo f.data .*= val 

    return nothing
end



# getindex method (bounds check performed by Base.Array)
function Base.:getindex(
    f :: MatsubaraFunction{Dg, Ds, Dt},
    x :: Vararg{Int64, Dt}
    ) :: Float64 where {Dg, Ds, Dt}

    return f.data[x...]
end

# setindex! method (bounds check performed by Base.Array)
function Base.:setindex!(
    f   :: MatsubaraFunction{Dg, Ds, Dt},
    val :: Float64,
    x   :: Vararg{Int64, Dt}
    )   :: Nothing where {Dg, Ds, Dt}

    f.data[x...] = val

    return nothing
end



# call to Matsubara function on 1D grid
@inbounds @fastmath function (f :: MatsubaraFunction{1, Ds, Dt})(
    w  :: Float64,
    x  :: Vararg{Int64, Ds} 
    ; 
    bc :: Float64 = 0.0
    )  :: Float64 where{Ds, Dt}

    ax = f.grids[1][1] <= w <= f.grids[1][grids_shape(f, 1)]

    if ax
        p = Param(w, f.grids[1]) 
        return p.wgts[1] * f.data[p.idxs[1], x...] + p.wgts[2] * f.data[p.idxs[2], x...]
    else 
        return bc 
    end 
end

# call to Matsubara function on 2D grid
@inbounds @fastmath function (f :: MatsubaraFunction{2, Ds, Dt})(
    w  :: NTuple{2, Float64},
    x  :: Vararg{Int64, Ds} 
    ; 
    bc :: Float64 = 0.0
    )  :: Float64 where{Ds, Dt}

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
@inbounds @fastmath function (f :: MatsubaraFunction{3, Ds, Dt})(
    w  :: NTuple{3, Float64},
    x  :: Vararg{Int64, Ds} 
    ; 
    bc :: Float64 = 0.0
    )  :: Float64 where{Ds, Dt}

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