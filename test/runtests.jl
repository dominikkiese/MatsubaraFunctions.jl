using Test
using MatsubaraFunctions 

# check whether MatsubaraFunctions correctly evaluate on their associated linear grids
@testset "Linear" begin 
    fg = MatsubaraGrid(1.0, 10, Fermion); nf = length(fg)
    bg = MatsubaraGrid(1.0, 10, Boson);   nb = length(bg)

    @test isapprox(fg.data, Float64[v for v in fg], atol = 1e-14, rtol = 0.0)
    @test isapprox(bg.data, Float64[w for w in bg], atol = 1e-14, rtol = 0.0)

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

    @test isapprox(fg.data, Float64[v for v in fg], atol = 1e-14, rtol = 0.0)
    @test isapprox(bg.data, Float64[w for w in bg], atol = 1e-14, rtol = 0.0)

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
        f1[v, 1] = 1.0 / (im * fg[v])
        f2[v, 1] = 1.0 / (im * fg[v] - ξ)
        f3[v, 1] = 1.0 / (im * fg[v] - ξ) / (im * fg[v] - ξ)
        f4[v, 1] = 1.0 / fg[v]
    end 

    w = 2.0 * fg[end]
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
        f1[v, 1] = 1.0 / (im * fg[v])
        f2[v, 1] = 1.0 / (im * fg[v] - ξ)
        f3[v, 1] = 1.0 / (im * fg[v] + ξ)
        f4[v, 1] = 1.0 / (im * fg[v] - ξ) / (im * fg[v] - ξ)
        f5[v, 1] = 1.0 / fg[v]
    end 

    ρ(x, T) = 1.0 / (exp(x / T) + 1.0)
    ρ10     = ρ(+0, T) - 1.0
    ρ1p     = ρ(+ξ, T) - 1.0
    ρ1m     = ρ(-ξ, T) - 1.0
    ρ2      = ρ(ξ, T) * (ρ(ξ, T) - 1.0) / T
    ρ1r     = im * (ρ(0.0, T) - 1.0)

    @test isapprox(sum(f1, 1), ρ10; atol = 1e-6, rtol = 0.0)
    @test isapprox(sum(f2, 1), ρ1p; atol = 1e-6, rtol = 0.0)
    @test isapprox(sum(f3, 1), ρ1m; atol = 1e-6, rtol = 0.0)
    @test isapprox(sum(f4, 1),  ρ2; atol = 1e-6, rtol = 0.0)
    @test isapprox(sum(f5, 1), ρ1r; atol = 1e-6, rtol = 0.0)
end