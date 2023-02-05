using Test
using MatsubaraFunctions 

# check basic operations for MatsubaraFrequency type
@testset "Frequencies" begin 

    for n in 1 : 10 
        T        = rand()
        wFermion = MatsubaraFrequency(T, rand(-100 : 100), Fermion)
        wBoson   = MatsubaraFrequency(T, rand(-100 : 100), Boson)

        @test isapprox(value(wFermion + wFermion), value(wFermion) + value(wFermion))
        @test isapprox(value(wFermion + wBoson), value(wFermion) + value(wBoson))
        @test isapprox(value(wBoson + wFermion), value(wBoson) + value(wFermion))
        @test isapprox(value(wBoson + wBoson), value(wBoson) + value(wBoson))
        @test isapprox(value(wFermion - wFermion), value(wFermion) - value(wFermion))
        @test isapprox(value(wFermion - wBoson), value(wFermion) - value(wBoson))
        @test isapprox(value(wBoson - wFermion), value(wBoson) - value(wFermion))
        @test isapprox(value(wBoson - wBoson), value(wBoson) - value(wBoson))
    end 
end

# check basic operations for MatsubaraGrid
@testset "Grids" begin 

    T          = rand() 
    wFermion_l = MatsubaraGrid(T, 100, Fermion)
    wBoson_l   = MatsubaraGrid(T, 100, Boson)
    wFermion_c = MatsubaraGrid(T, 100, 25, 4, Fermion)
    wBoson_c   = MatsubaraGrid(T, 100, 25, 4, Boson)

    @test isapprox(length(wFermion_l), 200)
    @test isapprox(length(wBoson_l), 199)
    @test isapprox(length(wFermion_c), 200)
    @test isapprox(length(wBoson_c), 199)

    @test isapprox(Float64[value(w) for w in wFermion_l.data], Float64[value(w) for w in wFermion_l])
    @test isapprox(Float64[value(w) for w in wBoson_l.data], Float64[value(w) for w in wBoson_l])
    @test isapprox(Float64[value(w) for w in wFermion_c.data], Float64[value(w) for w in wFermion_c])
    @test isapprox(Float64[value(w) for w in wBoson_c.data], Float64[value(w) for w in wBoson_c])

    for n in 1 : 10
        wFermion_l_idx = rand(1 : length(wFermion_l))
        wBoson_l_idx   = rand(1 : length(wBoson_l))
        wFermion_c_idx = rand(1 : length(wFermion_c))
        wBoson_c_idx   = rand(1 : length(wBoson_c))

        @test isapprox(wFermion_l(wFermion_l[wFermion_l_idx]), wFermion_l_idx)
        @test isapprox(wFermion_l(value(wFermion_l[wFermion_l_idx])), wFermion_l_idx)
        @test isapprox(wBoson_l(wBoson_l[wBoson_l_idx]), wBoson_l_idx)
        @test isapprox(wBoson_l(value(wBoson_l[wBoson_l_idx])), wBoson_l_idx)
        @test isapprox(wFermion_c(wFermion_c[wFermion_c_idx]), wFermion_c_idx)
        @test isapprox(wFermion_c(value(wFermion_c[wFermion_c_idx])), wFermion_c_idx)
        @test isapprox(wBoson_c(wBoson_c[wBoson_c_idx]), wBoson_c_idx)
        @test isapprox(wBoson_c(value(wBoson_c[wBoson_c_idx])), wBoson_c_idx)
    end 

    @test isapprox(wFermion_l(wFermion_l[1]), 1)
    @test isapprox(wFermion_l(value(wFermion_l[1])), 1)
    @test isapprox(wBoson_l(wBoson_l[1]), 1)
    @test isapprox(wBoson_l(value(wBoson_l[1])), 1)
    @test isapprox(wFermion_c(wFermion_c[1]), 1)
    @test isapprox(wFermion_c(value(wFermion_c[1])), 1)
    @test isapprox(wBoson_c(wBoson_c[1]), 1)
    @test isapprox(wBoson_c(value(wBoson_c[1])), 1)

    @test isapprox(wFermion_l(wFermion_l[end]), length(wFermion_l))
    @test isapprox(wFermion_l(value(wFermion_l[end])), length(wFermion_l))
    @test isapprox(wBoson_l(wBoson_l[end]), length(wBoson_l))
    @test isapprox(wBoson_l(value(wBoson_l[end])), length(wBoson_l))
    @test isapprox(wFermion_c(wFermion_c[end]), length(wFermion_c))
    @test isapprox(wFermion_c(value(wFermion_c[end])), length(wFermion_c))
    @test isapprox(wBoson_c(wBoson_c[end]), length(wBoson_c))
    @test isapprox(wBoson_c(value(wBoson_c[end])), length(wBoson_c))
end

# check constructors for MatsubaraFunction
@testset "Constructors" begin 
    fg = MatsubaraGrid(1.0, 10, Fermion)

    @test typeof(MatsubaraFunction((fg, fg), (2, 2), Float64)) == MatsubaraFunction{2, 2, 4, Linear, Float64}
    @test typeof(MatsubaraFunction((fg, fg), (2, 2))) == MatsubaraFunction{2, 2, 4, Linear, ComplexF64}
    @test typeof(MatsubaraFunction(fg, (2, 2), Float64)) == MatsubaraFunction{1, 2, 3, Linear, Float64}
    @test typeof(MatsubaraFunction(fg, (2, 2))) == MatsubaraFunction{1, 2, 3, Linear, ComplexF64}
    @test typeof(MatsubaraFunction(fg, 2, Float64)) == MatsubaraFunction{1, 1, 2, Linear, Float64}
    @test typeof(MatsubaraFunction(fg, 2)) == MatsubaraFunction{1, 1, 2, Linear, ComplexF64}
end

# check whether MatsubaraFunctions correctly evaluate on their associated linear grids
@testset "Linear" begin 
    fg = MatsubaraGrid(1.0, 10, Fermion); nf = length(fg)
    bg = MatsubaraGrid(1.0, 10, Boson);   nb = length(bg)

    f1D = MatsubaraFunction((fg,), (10, 10), rand(nf, 10, 10))
    f2D = MatsubaraFunction((bg, fg), (10, 10), rand(nb, nf, 10, 10))
    f3D = MatsubaraFunction((bg, fg, fg), (10, 10), rand(nb, nf, nf, 10, 10))

    @test typeof(f1D) == MatsubaraFunction{1, 2, 3, Linear, Float64}
    @test typeof(f2D) == MatsubaraFunction{2, 2, 4, Linear, Float64}
    @test typeof(f3D) == MatsubaraFunction{3, 2, 5, Linear, Float64}

    for n in 1 : 10 
        x     = rand(1 : 10), rand(1 : 10)
        idx1D = rand(1 : nf), x...
        idx2D = rand(1 : nb), rand(1 : nf), x...
        idx3D = rand(1 : nb), rand(1 : nf), rand(1 : nf), x...
        w1D   = fg[idx1D[1]]
        w2D   = bg[idx2D[1]], fg[idx2D[2]]
        w3D   = bg[idx3D[1]], fg[idx3D[2]], fg[idx3D[3]]

        @test isapprox(f1D(w1D, x...), f1D[idx1D...], atol = 1e-14, rtol = 0.0)
        @test isapprox(f2D(w2D, x...), f2D[idx2D...], atol = 1e-14, rtol = 0.0)
        @test isapprox(f3D(w3D, x...), f3D[idx3D...], atol = 1e-14, rtol = 0.0)
    end 
end

# check whether MatsubaraFunctions correctly evaluate on their associated coarse grids
@testset "Coarse" begin 
    fg = MatsubaraGrid(1.0, 10, 5, 5, Fermion); nf = length(fg)
    bg = MatsubaraGrid(1.0, 10, 5, 5, Boson);   nb = length(bg)

    f1D = MatsubaraFunction((fg,), (10, 10), rand(nf, 10, 10))
    f2D = MatsubaraFunction((bg, fg), (10, 10), rand(nb, nf, 10, 10))
    f3D = MatsubaraFunction((bg, fg, fg), (10, 10), rand(nb, nf, nf, 10, 10))

    @test typeof(f1D) == MatsubaraFunction{1, 2, 3, Coarse, Float64}
    @test typeof(f2D) == MatsubaraFunction{2, 2, 4, Coarse, Float64}
    @test typeof(f3D) == MatsubaraFunction{3, 2, 5, Coarse, Float64}

    for n in 1 : 10 
        x     = rand(1 : 10), rand(1 : 10)
        idx1D = rand(1 : nf), x...
        idx2D = rand(1 : nb), rand(1 : nf), x...
        idx3D = rand(1 : nb), rand(1 : nf), rand(1 : nf), x...
        w1D   = fg[idx1D[1]]
        w2D   = bg[idx2D[1]], fg[idx2D[2]]
        w3D   = bg[idx3D[1]], fg[idx3D[2]], fg[idx3D[3]]

        @test isapprox(f1D(w1D, x...), f1D[idx1D...], atol = 1e-14, rtol = 0.0)
        @test isapprox(f2D(w2D, x...), f2D[idx2D...], atol = 1e-14, rtol = 0.0)
        @test isapprox(f3D(w3D, x...), f3D[idx3D...], atol = 1e-14, rtol = 0.0)
    end 
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

    ρ(x, T) = 1.0 / (exp(x / T) + 1.0)
    ρ10     = ρ(+0, T) - 1.0
    ρ1p     = ρ(+ξ, T) - 1.0
    ρ1m     = ρ(-ξ, T) - 1.0
    ρ2      = ρ(ξ, T) * (ρ(ξ, T) - 1.0) / T
    ρ1r     = im * ρ10

    @test isapprox(sum_me(f1, 1), ρ10; atol = 1e-6, rtol = 0.0)
    @test isapprox(sum_me(f2, 1), ρ1p; atol = 1e-6, rtol = 0.0)
    @test isapprox(sum_me(f3, 1), ρ1m; atol = 1e-6, rtol = 0.0)
    @test isapprox(sum_me(f4, 1),  ρ2; atol = 1e-6, rtol = 0.0)
    @test isapprox(sum_me(f5, 1), ρ1r; atol = 1e-6, rtol = 0.0)
end