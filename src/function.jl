struct MatsubaraFunction{Ds, Dg, Dt}
    shape :: NTuple{Ds, Int64}          
    grids :: NTuple{Dg, Vector{Float64}}
    data  :: Array{Float64, Dt}

    # safe constructor
    function MatsubaraFunction(
        shape :: NTuple{Ds, Int64}, 
        grids :: NTuple{Dg, Vector{Float64}}, 
        data  :: Array{Float64, Dt}
        )     :: IndexFunction{Ds, Dg, Dt} where {Ds, Dg, Dt}

        @assert Dg <= 3 "Hybercubic grids currently not supported"
        @assert Ds + Dg == Dt "Dimensions do not match"

        for g in grids
            @assert issorted(g) "Grids must be sorted"
            delta = g[2] - g[1]

            for i in 2 : (length(g) - 1)
                @assert isapprox(g[i + 1] - g[i], delta) "Grids must be linear"
            end
        end

        return new{Ds, Dg, Dt}(shape, grids, data)
    end

    # convenience constructor
    function MatsubaraFunction(
        shape :: NTuple{Ds, Int64}, 
        grids :: NTuple{Dg, Vector{Float64}}
        )     :: MatsubaraFunction{Ds, Dg, Ds + Dg} where {Ds, Dg}

        dims = shape..., length.(grids)...
        data = Array{Float64, Ds + Dg}(undef, dims)

        return MatsubaraFunction(shape, grids, data)
    end
end

# basic addition
function Base.:+(
    f1 :: MatsubaraFunction{Ds, Dg, Dt}, 
    f2 :: MatsubaraFunction{Ds, Dg, Dt}
    )  :: MatsubaraFunction{Ds, Dg, Dt} where {Ds, Dg, Dt}

    for i in 1 : Dg 
        @assert isapprox(f1.grids[i], f2.grids[i]) "Grids must be equal for addition" 
    end

    return MatsubaraFunction(f1.shape, f1.grids, @turbo f1.data .+ f2.data)
end

function add!(
    f1 :: MatsubaraFunction{Ds, Dg, Dt}, 
    f2 :: MatsubaraFunction{Ds, Dg, Dt}
    )  :: Nothing where {Ds, Dg, Dt}

    for i in 1 : Dg 
        @assert isapprox(f1.grids[i], f2.grids[i]) "Grids must be equal for addition" 
    end 

    @turbo f1.data .+= f2.data

    return nothing 
end

# basic subtraction
function Base.:-(
    f1 :: MatsubaraFunction{Ds, Dg, Dt}, 
    f2 :: MatsubaraFunction{Ds, Dg, Dt}
    )  :: MatsubaraFunction{Ds, Dg, Dt} where {Ds, Dg, Dt}

    for i in 1 : Dg 
        @assert isapprox(f1.grids[i], f2.grids[i]) "Grids must be equal for subtraction" 
    end 

    return MatsubaraFunction(f1.shape, f1.grids, @turbo f1.data .- f2.data)
end

function subtract!(
    f1 :: MatsubaraFunction{Ds, Dg, Dt}, 
    f2 :: MatsubaraFunction{Ds, Dg, Dt}
    )  :: Nothing where {Ds, Dg, Dt}

    for i in 1 : Dg 
        @assert isapprox(f1.grids[i], f2.grids[i]) "Grids must be equal for subtraction" 
    end 

    @turbo f1.data .-= f2.data

    return nothing 
end

# basic multiplication with scalar 
function Base.:*(
    val :: Float64, 
    f   :: MatsubaraFunction{Ds, Dg, Dt}
    )   :: MatsubaraFunction{Ds, Dg, Dt} where {Ds, Dg, Dt}

    return MatsubaraFunction(f.shape, f.grids, @turbo val .* f.data)
end

function Base.:*(
    f   :: MatsubaraFunction{Ds, Dg, Dt},
    val :: Float64
    )   :: MatsubaraFunction{Ds, Dg, Dt} where {Ds, Dg, Dt}

    return MatsubaraFunction(f.shape, f.grids, @turbo val .* f.data)
end

function mult!(
    f   :: MatsubaraFunction{Ds, Dg, Dt},
    val :: Float64
    )   :: Nothing where {Ds, Dg, Dt}

    @turbo f.data .*= val 

    return nothing
end

# call to Matsubara function on 1D grid
@inbounds @fastmath function (f :: MatsubaraFunction{Ds, 1, Dt})(
    x  :: NTuple{Ds, Int64}, 
    w  :: Float64
    ; 
    bc :: Float64 = 0.0
    )  :: Float64 where{Ds, Dt}

    ax = first(f.grids[1]) <= w <= last(f.grids[1])

    if ax
        p = Param(w, f.grids[1]) 
        return p.wgts[1] * f.data[x..., p.idxs[1]] + p.wgts[2] * f.data[x..., p.idxs[2]]
    else 
        return bc 
    end 
end

# call to Matsubara function on 2D grid
@inbounds @fastmath function (f :: MatsubaraFunction{Ds, 2, Dt})(
    x  :: NTuple{Ds, Int64}, 
    w  :: NTuple{2, Float64}
    ; 
    bc :: Float64 = 0.0
    )  :: Float64 where{Ds, Dt}

    ax1 = first(f.grids[1]) <= w[1] <= last(f.grids[1])
    ax2 = first(f.grids[2]) <= w[2] <= last(f.grids[2])

    if ax1 && ax2
        p1 = Param(w[1], f.grids[1])
        p2 = Param(w[2], f.grids[2])

        return p1.wgts[1] * p2.wgts[1] * f.data[x..., p1.idxs[1], p2.idxs[1]] + 
               p1.wgts[1] * p2.wgts[2] * f.data[x..., p1.idxs[1], p2.idxs[2]] + 
               p1.wgts[2] * p2.wgts[1] * f.data[x..., p1.idxs[2], p2.idxs[1]] + 
               p1.wgts[2] * p2.wgts[2] * f.data[x..., p1.idxs[2], p2.idxs[2]]         
    else 
        return bc 
    end 
end

# call to Matsubara function on 3D grid
@inbounds @fastmath function (f :: MatsubaraFunction{Ds, 3, Dt})(
    x  :: NTuple{Ds, Int64}, 
    w  :: NTuple{3, Float64}
    ; 
    bc :: Float64 = 0.0
    )  :: Float64 where{Ds, Dt}

    ax1 = first(f.grids[1]) <= w[1] <= last(f.grids[1])
    ax2 = first(f.grids[2]) <= w[2] <= last(f.grids[2])
    ax3 = first(f.grids[3]) <= w[3] <= last(f.grids[3])

    if ax1 && ax2 && ax3
        p1 = Param(w[1], f.grids[1])
        p2 = Param(w[2], f.grids[2])
        p3 = Param(w[3], f.grids[3])

        return p1.wgts[1] * p2.wgts[1] * p3.wgts[1] * f.data[x..., p1.idxs[1], p2.idxs[1], p3.idxs[1]] + 
               p1.wgts[1] * p2.wgts[1] * p3.wgts[2] * f.data[x..., p1.idxs[1], p2.idxs[1], p3.idxs[2]] + 
               p1.wgts[1] * p2.wgts[2] * p3.wgts[1] * f.data[x..., p1.idxs[1], p2.idxs[2], p3.idxs[1]] + 
               p1.wgts[1] * p2.wgts[2] * p3.wgts[2] * f.data[x..., p1.idxs[1], p2.idxs[2], p3.idxs[2]] + 
               p1.wgts[2] * p2.wgts[1] * p3.wgts[1] * f.data[x..., p1.idxs[2], p2.idxs[1], p3.idxs[1]] + 
               p1.wgts[2] * p2.wgts[1] * p3.wgts[2] * f.data[x..., p1.idxs[2], p2.idxs[1], p3.idxs[2]] + 
               p1.wgts[2] * p2.wgts[2] * p3.wgts[1] * f.data[x..., p1.idxs[2], p2.idxs[2], p3.idxs[1]] + 
               p1.wgts[2] * p2.wgts[2] * p3.wgts[2] * f.data[x..., p1.idxs[2], p2.idxs[2], p3.idxs[2]]  
    else 
        return bc 
    end 
end