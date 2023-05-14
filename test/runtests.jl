using Test
using MatsubaraFunctions 

@testset "Frequencies" begin 
    for n in 1 : 10 
        T        = rand()
        wFermion = MatsubaraFrequency(T, rand(-100 : 100), Fermion)
        wBoson   = MatsubaraFrequency(T, rand(-100 : 100), Boson)

        # addition
        @test value(wFermion + wFermion) ≈ value(wFermion) + value(wFermion)
        @test value(wFermion + wBoson) ≈ value(wFermion) + value(wBoson)
        @test value(wBoson + wFermion) ≈ value(wBoson) + value(wFermion)
        @test value(wBoson + wBoson) ≈ value(wBoson) + value(wBoson)

        # subtraction
        @test value(wFermion - wFermion) ≈ value(wFermion) - value(wFermion)
        @test value(wFermion - wBoson) ≈ value(wFermion) - value(wBoson)
        @test value(wBoson - wFermion) ≈ value(wBoson) - value(wFermion)
        @test value(wBoson - wBoson) ≈ value(wBoson) - value(wBoson)

        # reflection 
        @test value(-wFermion) ≈ -value(wFermion)
        @test value(-wBoson) ≈ -value(wBoson)
    end 
end

@testset "Grids" begin 
    T        = rand() 
    wFermion = MatsubaraGrid(T, 100, Fermion)
    wBoson   = MatsubaraGrid(T, 100, Boson)

    # length
    @test length(wFermion) == 200
    @test length(wBoson) == 199

    # iterator
    @test Float64[value(w) for w in wFermion.data] ≈ Float64[value(w) for w in wFermion]
    @test Float64[value(w) for w in wBoson.data] ≈ Float64[value(w) for w in wBoson]

    for n in 1 : 10
        wFermion_idx = rand(1 : length(wFermion))
        wBoson_idx   = rand(1 : length(wBoson))

        # call to fermionic grid
        @test wFermion(wFermion[wFermion_idx]) == wFermion_idx
        @test wFermion(value(wFermion[wFermion_idx])) == wFermion_idx

        # call to bosonic grid
        @test wBoson(wBoson[wBoson_idx]) == wBoson_idx
        @test wBoson(value(wBoson[wBoson_idx])) == wBoson_idx
    end 

    # first element
    @test wFermion(wFermion[1]) == 1
    @test wFermion(value(wFermion[1])) == 1
    @test wBoson(wBoson[1]) == 1
    @test wBoson(value(wBoson[1])) == 1

    # last element
    @test wFermion(wFermion[end]) == length(wFermion)
    @test wFermion(value(wFermion[end])) == length(wFermion)
    @test wBoson(wBoson[end]) == length(wBoson)
    @test wBoson(value(wBoson[end])) == length(wBoson)
end

@testset "Constructors" begin 
    fg = MatsubaraGrid(1.0, 10, Fermion)

    @test typeof(MatsubaraFunction((fg, fg), (2, 2), Float64)) == MatsubaraFunction{2, 2, 4, Float64}
    @test typeof(MatsubaraFunction((fg, fg), (2, 2))) == MatsubaraFunction{2, 2, 4, ComplexF64}
    @test typeof(MatsubaraFunction(fg, (2, 2), Float64)) == MatsubaraFunction{1, 2, 3, Float64}
    @test typeof(MatsubaraFunction(fg, (2, 2))) == MatsubaraFunction{1, 2, 3, ComplexF64}
    @test typeof(MatsubaraFunction(fg, 2, Float64)) == MatsubaraFunction{1, 1, 2, Float64}
    @test typeof(MatsubaraFunction(fg, 2)) == MatsubaraFunction{1, 1, 2, ComplexF64}
end

@testset "Indices" begin 
    g    = MatsubaraGrid(1.0, 10, Fermion)
    f    = MatsubaraFunction((g, g), (10, 10))
    idxs = rand(1 : length(g)), rand(1 : length(g)), rand(1 : 10), rand(1 : 10)
    lidx = LinearIndex(f, (g[idxs[1]], g[idxs[2]]), idxs[3], idxs[4])
    cidx = CartesianIndex(f, (g[idxs[1]], g[idxs[2]]), idxs[3], idxs[4])
    x    = to_Matsubara(f, lidx) 
    y    = to_Matsubara(f, cidx)

    @test CartesianIndex(f, lidx) == cidx 
    @test LinearIndex(f, cidx) == lidx

    @test value(first(x)[1]) ≈ value(g[idxs[1]])
    @test value(first(x)[2]) ≈ value(g[idxs[2]])
    @test last(x)[1] == idxs[3] 
    @test last(x)[2] == idxs[4]

    @test value(first(y)[1]) ≈ value(g[idxs[1]])
    @test value(first(y)[2]) ≈ value(g[idxs[2]])
    @test last(y)[1] == idxs[3] 
    @test last(y)[2] == idxs[4]
end

@testset "Evaluate" begin 
    fg = MatsubaraGrid(1.0, 10, Fermion); nf = length(fg)
    bg = MatsubaraGrid(1.0, 10, Boson); nb = length(bg)

    F1D = MatsubaraFunction((fg,), (1,), rand(nf, 1))
    F2D = MatsubaraFunction((bg, fg), (1,), rand(nf, nf, 1))
    F3D = MatsubaraFunction((bg, fg, fg), (1,), rand(nf, nf, nf, 1))

    f1D = MatsubaraFunction((fg,), (10, 10), rand(nf, 10, 10))
    f2D = MatsubaraFunction((bg, fg), (10, 10), rand(nb, nf, 10, 10))
    f3D = MatsubaraFunction((bg, fg, fg), (10, 10), rand(nb, nf, nf, 10, 10))

    # type checks
    @test typeof(f1D) == MatsubaraFunction{1, 2, 3, Float64}
    @test typeof(f2D) == MatsubaraFunction{2, 2, 4, Float64}
    @test typeof(f3D) == MatsubaraFunction{3, 2, 5, Float64}

    for n in 1 : 10 
        x     = rand(1 : 10), rand(1 : 10)
        idx1D = rand(1 : nf), x...
        idx2D = rand(1 : nb), rand(1 : nf), x...
        idx3D = rand(1 : nb), rand(1 : nf), rand(1 : nf), x...
        w1D   = fg[idx1D[1]]
        w2D   = bg[idx2D[1]], fg[idx2D[2]]
        w3D   = bg[idx3D[1]], fg[idx3D[2]], fg[idx3D[3]]

        # test evaluators for MatsubaraFrequency and Float64
        @test F1D(w1D) ≈ F1D[idx1D[1], 1]; @test F1D(value(w1D)) ≈ F1D[idx1D[1], 1]
        @test F2D(w2D...) ≈ F2D[idx2D[1 : 2]..., 1]; @test F2D(value.(w2D)...) ≈ F2D[idx2D[1 : 2]..., 1]
        @test F3D(w3D...) ≈ F3D[idx3D[1 : 3]..., 1]; @test F3D(value.(w3D)...) ≈ F3D[idx3D[1 : 3]..., 1]

        @test f1D(w1D, x...) ≈ f1D[idx1D...]; @test f1D(value(w1D), x...) ≈ f1D[idx1D...]
        @test f2D(w2D, x...) ≈ f2D[idx2D...]; @test f2D(value.(w2D), x...) ≈ f2D[idx2D...]
        @test f3D(w3D, x...) ≈ f3D[idx3D...]; @test f3D(value.(w3D), x...) ≈ f3D[idx3D...]
    end 
end

@testset "BC" begin 
    T  = 0.5
    ξ  = 0.5
    v  = MatsubaraFrequency(T, 1000, Fermion)
    fg = MatsubaraGrid(T, 512, Fermion)
    f1 = MatsubaraFunction(fg, 1)
    f2 = MatsubaraFunction((fg, fg), 1)

    # default bc
    @test f1(v) ≈ 0.0
    @test f2(v, fg[1]) ≈ 0.0
    @test f1(value(v)) ≈ 0.0
    @test f2(value(v), value(fg[1])) ≈ 0.0

    # constant bc
    @test f1(v; bc = x -> 1.0) ≈ 1.0
    @test f2(v, fg[1]; bc = x -> 1.0) ≈ 1.0
    @test f1(value(v); bc = x -> 1.0) ≈ 1.0
    @test f2(value(v), value(fg[1]); bc = x -> 1.0) ≈ 1.0

    # frequency dependent bc
    @test f1(v; bc = x -> 1.0 / value(x)) ≈ 1.0 / value(v)
    @test f2(v, fg[1]; bc = x -> 1.0 / value(x[1]) / value(x[2])) ≈ 1.0 / value(v) / value(fg[1])
    @test f1(value(v); bc = x -> 1.0 / x) ≈ 1.0 / value(v)
    @test f2(value(v), value(fg[1]); bc = x -> 1.0 / x[1] / x[2]) ≈ 1.0 / value(v) / value(fg[1])
end

@testset "Extrapolation" begin 
    T  = 0.1
    ξ  = 0.5
    fg = MatsubaraGrid(T, 5000, Fermion)
    f1 = MatsubaraFunction(fg, 1)
    f2 = MatsubaraFunction(fg, 1)
    f3 = MatsubaraFunction(fg, 1)
    f4 = MatsubaraFunction(fg, 1, Float64)
    f5 = MatsubaraFunction((fg, fg), 1)

    for v in fg
        f1[v] = 1.0 / (im * value(v))
        f2[v] = 1.0 / (im * value(v) - ξ)
        f3[v] = 1.0 / (im * value(v) - ξ) / (im * value(v) - ξ)
        f4[v] = 1.0 / value(v)

        for vp in fg 
            f5[v, vp] = 1.0 / (im * value(v) - ξ) / (im * value(vp) - ξ)
        end
    end 
    
    w  = fg[end] + fg[end] + fg[end]
    c0 = ComplexF64(0.0)

    # polynomial extrapolation for 1D grids with MatsubaraFrequency argument
    @test isapprox(f1(w; extrp = (true, c0)), 1.0 / (im * value(w)); atol = 1e-6, rtol = 1e-6)
    @test isapprox(f2(w; extrp = (true, c0)), 1.0 / (im * value(w)) - ξ / value(w) / value(w); atol = 1e-6, rtol = 1e-6)
    @test isapprox(f3(w; extrp = (true, c0)), -1.0 / value(w) / value(w); atol = 1e-6, rtol = 1e-6)
    @test isapprox(f4(w; extrp = (true, 0.)), 1.0 / value(w); atol = 1e-6, rtol = 1e-6)

    # polynomial extrapolation for 1D grids with Float64 argument
    @test isapprox(f1(value(w); extrp = (true, c0)), 1.0 / (im * value(w)); atol = 1e-6, rtol = 1e-6)
    @test isapprox(f2(value(w); extrp = (true, c0)), 1.0 / (im * value(w)) - ξ / value(w) / value(w); atol = 1e-6, rtol = 1e-6)
    @test isapprox(f3(value(w); extrp = (true, c0)), -1.0 / value(w) / value(w); atol = 1e-6, rtol = 1e-6)
    @test isapprox(f4(value(w); extrp = (true, 0.)), 1.0 / value(w); atol = 1e-6, rtol = 1e-6)

    # constant extrapolation for higher-dimensional grids with MatsubaraFrequency argument
    @test f5(w, w; extrp = (true, c0)) ≈ 1.0 / (im * value(fg[end]) - ξ) / (im * value(fg[end]) - ξ)
    @test f5(w, fg[1]; extrp = (true, c0)) ≈ 1.0 / (im * value(fg[end]) - ξ) / (im * value(fg[1]) - ξ)
    @test f5(fg[1], w; extrp = (true, c0)) ≈ 1.0 / (im * value(fg[1]) - ξ) / (im * value(fg[end]) - ξ)

    # constant extrapolation for higher-dimensional grids with Float64 argument
    @test f5(value(w), value(w); extrp = (true, c0)) ≈ 1.0 / (im * value(fg[end]) - ξ) / (im * value(fg[end]) - ξ)
    @test f5(value(w), value(fg[1]); extrp = (true, c0)) ≈ 1.0 / (im * value(fg[end]) - ξ) / (im * value(fg[1]) - ξ)
    @test f5(value(fg[1]), value(w); extrp = (true, c0)) ≈ 1.0 / (im * value(fg[1]) - ξ) / (im * value(fg[end]) - ξ)
end

@testset "Summation" begin 
    T  = 0.1
    ξ  = 0.5
    fg = MatsubaraGrid(T, 5000, Fermion)
    f1 = MatsubaraFunction(fg, 1)
    f2 = MatsubaraFunction(fg, 1)
    f3 = MatsubaraFunction(fg, 1)
    f4 = MatsubaraFunction(fg, 1)

    for v in fg
        f1[v] = 1.0 / (im * value(v))
        f2[v] = 1.0 / (im * value(v) - ξ)
        f3[v] = 1.0 / (im * value(v) + ξ)
        f4[v] = 1.0 / (im * value(v) - ξ) / (im * value(v) - ξ)
    end 

    # compute analytic results
    ρ(x, T) = 1.0 / (exp(x / T) + 1.0)
    ρ10     = ρ(+0, T) - 1.0
    ρ1p     = ρ(+ξ, T) - 1.0
    ρ1m     = ρ(-ξ, T) - 1.0
    ρ2      = ρ(ξ, T) * (ρ(ξ, T) - 1.0) / T
    c0      = ComplexF64(0.0)

    # benchmark vs analytic results
    @test isapprox(sum_me(f1, c0, 1), ρ10; atol = 1e-6, rtol = 1e-6)
    @test isapprox(sum_me(f2, c0, 1), ρ1p; atol = 1e-6, rtol = 1e-6)
    @test isapprox(sum_me(f3, c0, 1), ρ1m; atol = 1e-6, rtol = 1e-6)
    @test isapprox(sum_me(f4, c0, 1),  ρ2; atol = 1e-6, rtol = 1e-6)
end

@testset "Symmetries" begin 
    T  = 0.1
    ξ  = 0.5
    g  = MatsubaraGrid(T, 128, Fermion)
    f1 = MatsubaraFunction(g, 1)
    f2 = MatsubaraFunction(g, 1)
    f3 = MatsubaraFunction(g, 1)

    for v in g
        f1[v] = 1.0 / (im * value(v) - ξ)
    end

    # complex conjugation for Green's function
    function conj(
        w :: Tuple{MatsubaraFrequency},
        x :: Tuple{Int64}
        ) :: Tuple{Tuple{MatsubaraFrequency}, Tuple{Int64}, MatsubaraOperation}

        return (-w[1],), (x[1],), MatsubaraOperation(false, true)
    end 

    # compute the symmetry group 
    SG = MatsubaraSymmetryGroup([MatsubaraSymmetry{1, 1}(conj)], f1)

    # symmetrize f2 and compare to f1 
    for class in SG.classes 
        f2[class[1][1]] = f1[class[1][1]]
    end 

    SG(f2)
    @test f2 == f1

    # symmetrize f3 and compare to f1 using MatsubaraInitFunction
    function init(
        w :: Tuple{MatsubaraFrequency},
        x :: Tuple{Int64}
        ) :: ComplexF64

        return f1[w, x...]
    end 

    InitFunc = MatsubaraInitFunction{1, 1, ComplexF64}(init)
    SG(f3, InitFunc)
    @test f3 == f1
end