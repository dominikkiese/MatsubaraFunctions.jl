using Test
using MatsubaraFunctions 

@testset "Evaluation" begin 
    gf = mk_grid(1.0, 10, :Fermion)
    gb = mk_grid(1.0, 10, :Boson)

    nf = length(gf)
    nb = length(gb)

    f1D = MatsubaraFunction((10, 10), (gf,), rand(10, 10, nf))
    f2D = MatsubaraFunction((10, 10), (gb, gf), rand(10, 10, nb, nf))
    f3D = MatsubaraFunction((10, 10), (gb, gf, gf), rand(10, 10, nb, nf, nf))

    for n in 1 : 10 
        x     = rand(1 : 10), rand(1 : 10)
        idx1D = rand(1 : nf)
        idx2D = rand(1 : nb), rand(1 : nf) 
        idx3D = rand(1 : nb), rand(1 : nf), rand(1 : nf)

        w1D = gf[idx1D[1]]
        w2D = gb[idx2D[1]], gf[idx2D[2]]
        w3D = gb[idx3D[1]], gf[idx3D[2]], gf[idx3D[3]]

        @test isapprox(f1D(x, w1D), f1D.data[x..., idx1D])
        @test isapprox(f2D(x, w2D), f2D.data[x..., idx2D...])
        @test isapprox(f3D(x, w3D), f3D.data[x..., idx3D...])
    end 
end