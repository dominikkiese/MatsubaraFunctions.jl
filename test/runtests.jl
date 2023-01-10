using Test
using MatsubaraFunctions 

@testset "Linear" begin 
    fg = MatsubaraGrid(1.0, 10, Fermion)
    bg = MatsubaraGrid(1.0, 10, Boson)

    nf = length(fg)
    nb = length(bg)

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

        @test isapprox(f1D(w1D, x...), f1D[idx1D...])
        @test isapprox(f2D(w2D, x...), f2D[idx2D...])
        @test isapprox(f3D(w3D, x...), f3D[idx3D...])
    end 
end

@testset "Coarse" begin 
    fg = MatsubaraGrid(1.0, 10, 2.0, Fermion)
    bg = MatsubaraGrid(1.0, 10, 2.0, Boson)

    nf = length(fg)
    nb = length(bg)

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

        @test isapprox(f1D(w1D, x...), f1D[idx1D...])
        @test isapprox(f2D(w2D, x...), f2D[idx2D...])
        @test isapprox(f3D(w3D, x...), f3D[idx3D...])
    end 
end