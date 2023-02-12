using Test
using MatsubaraFunctions 

# check basic operations for MatsubaraFrequency
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
    end 
end

# check basic operations for MatsubaraGrid
@testset "Grids" begin 

    T          = rand() 
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

# check constructors for MatsubaraFunction
@testset "Constructors" begin 
    fg = MatsubaraGrid(1.0, 10, Fermion)

    @test typeof(MatsubaraFunction((fg, fg), (2, 2), Float64)) == MatsubaraFunction{2, 2, 4, Float64}
    @test typeof(MatsubaraFunction((fg, fg), (2, 2))) == MatsubaraFunction{2, 2, 4, ComplexF64}
    @test typeof(MatsubaraFunction(fg, (2, 2), Float64)) == MatsubaraFunction{1, 2, 3, Float64}
    @test typeof(MatsubaraFunction(fg, (2, 2))) == MatsubaraFunction{1, 2, 3, ComplexF64}
    @test typeof(MatsubaraFunction(fg, 2, Float64)) == MatsubaraFunction{1, 1, 2, Float64}
    @test typeof(MatsubaraFunction(fg, 2)) == MatsubaraFunction{1, 1, 2, ComplexF64}
end

# check whether MatsubaraFunctions correctly evaluate on their associated linear grids
@testset "Evaluate" begin 
    fg = MatsubaraGrid(1.0, 10, Fermion); nf = length(fg)
    bg = MatsubaraGrid(1.0, 10, Boson);   nb = length(bg)

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
        @test f1D(w1D, x...) ≈ f1D[idx1D...]; @test f1D(value(w1D), x...)  ≈ f1D[idx1D...]
        @test f2D(w2D, x...) ≈ f2D[idx2D...]; @test f2D(value.(w2D), x...) ≈ f2D[idx2D...]
        @test f3D(w3D, x...) ≈ f3D[idx3D...]; @test f3D(value.(w3D), x...) ≈ f3D[idx3D...]
    end 
end

# check whether boundary conditions work as expected
@testset "BC" begin 
    T  = 0.5
    ξ  = 0.5
    v  = MatsubaraFrequency(T, 1000, Fermion)
    fg = MatsubaraGrid(T, 512, Fermion)
    f1 = MatsubaraFunction((fg,), (1,))
    f2 = MatsubaraFunction((fg, fg), (1,))

    # default bc
    @test f1(v, 1) ≈ 0.0
    @test f2((v, fg[1]), 1) ≈ 0.0
    @test f1(value(v), 1) ≈ 0.0
    @test f2((value(v), value(fg[1])), 1) ≈ 0.0

    # constant bc
    @test f1(v, 1; bc = x -> 1.0) ≈ 1.0
    @test f2((v, fg[1]), 1; bc = x -> 1.0) ≈ 1.0
    @test f1(value(v), 1; bc = x -> 1.0) ≈ 1.0
    @test f2((value(v), value(fg[1])), 1; bc = x -> 1.0) ≈ 1.0

    # frequency dependent bc
    @test f1(v, 1; bc = x -> 1.0 / value(x)) ≈ 1.0 / value(v)
    @test f2((v, fg[1]), 1; bc = x -> 1.0 / value(x[1]) / value(x[2])) ≈ 1.0 / value(v) / value(fg[1])
    @test f1(value(v), 1; bc = x -> 1.0 / x) ≈ 1.0 / value(v)
    @test f2((value(v), value(fg[1])), 1; bc = x -> 1.0 / x[1] / x[2]) ≈ 1.0 / value(v) / value(fg[1])
end

# check whether extrapolation routines work as expected
@testset "Extrapolation" begin 
    T  = 0.5
    ξ  = 0.5
    fg = MatsubaraGrid(T, 512, Fermion)
    f1 = MatsubaraFunction((fg,), (1,))
    f2 = MatsubaraFunction((fg,), (1,))
    f3 = MatsubaraFunction((fg,), (1,))
    f4 = MatsubaraFunction((fg,), (1,), Float64)

    for v in 1 : length(fg)
        f1[v, 1] = 1.0 / (im * value(fg[v]))
        f2[v, 1] = 1.0 / (im * value(fg[v]) - ξ)
        f3[v, 1] = 1.0 / (im * value(fg[v]) - ξ) / (im * value(fg[v]) - ξ)
        f4[v, 1] = 1.0 / value(fg[v])
    end 

    w = 2.0 * value(fg[end])
    @test isapprox(f1(w, 1; extrp = true), 1.0 / (im * w); atol = 1e-10, rtol = 0.0)
    @test isapprox(f2(w, 1; extrp = true), 1.0 / (im * w) - ξ / w / w; atol = 1e-10, rtol = 0.0)
    @test isapprox(f3(w, 1; extrp = true), -1.0 / w / w; atol = 1e-10, rtol = 0.0)
    @test isapprox(f4(w, 1; extrp = true), 1.0 / w; atol = 1e-10, rtol = 0.0)
end

# check whether summation routines work as expected
@testset "Summation" begin 
    T  = 0.5
    ξ  = 0.5
    fg = MatsubaraGrid(T, 512, Fermion)
    f1 = MatsubaraFunction((fg,), (1,))
    f2 = MatsubaraFunction((fg,), (1,))
    f3 = MatsubaraFunction((fg,), (1,))
    f4 = MatsubaraFunction((fg,), (1,))
    f5 = MatsubaraFunction((fg,), (1,), Float64)

    for v in 1 : length(fg)
        f1[v, 1] = 1.0 / (im * value(fg[v]))
        f2[v, 1] = 1.0 / (im * value(fg[v]) - ξ)
        f3[v, 1] = 1.0 / (im * value(fg[v]) + ξ)
        f4[v, 1] = 1.0 / (im * value(fg[v]) - ξ) / (im * value(fg[v]) - ξ)
        f5[v, 1] = 1.0 / value(fg[v])
    end 

    # compute analytic results
    ρ(x, T) = 1.0 / (exp(x / T) + 1.0)
    ρ10     = ρ(+0, T) - 1.0
    ρ1p     = ρ(+ξ, T) - 1.0
    ρ1m     = ρ(-ξ, T) - 1.0
    ρ2      = ρ(ξ, T) * (ρ(ξ, T) - 1.0) / T
    ρ1r     = im * ρ10

    # benchmark vs analytic results
    @test isapprox(sum_me(f1, 1), ρ10; atol = 1e-6, rtol = 0.0)
    @test isapprox(sum_me(f2, 1), ρ1p; atol = 1e-6, rtol = 0.0)
    @test isapprox(sum_me(f3, 1), ρ1m; atol = 1e-6, rtol = 0.0)
    @test isapprox(sum_me(f4, 1),  ρ2; atol = 1e-6, rtol = 0.0)
    @test isapprox(sum_me(f5, 1), ρ1r; atol = 1e-6, rtol = 0.0)
end