@testset "Extrapolation" begin 
    g  = MatsubaraGrid(0.1, 5000, Fermion)
    f1 = MatsubaraFunction(g, 1)
    f2 = MatsubaraFunction(g, 1)
    f3 = MatsubaraFunction((g, g), 1)

    for v in g
        f1[v] = 1.0 / (im * value(v) - 0.5)
        f2[v] = 1.0 / (im * value(v) - 0.5) / (im * value(v) - 0.5)

        for vp in g
            f3[v, vp] = 1.0 / (im * value(v) - 0.5) / (im * value(vp) - 0.5)
        end
    end 
    
    w = g[end] + g[end] + g[end]

    # polynomial extrapolation for 1D grids with MatsubaraFrequency argument
    @test isapprox(f1(w), 1.0 / (im * value(w) - 0.5); atol = 1e-6, rtol = 1e-6)
    @test isapprox(f2(w), 1.0 / (im * value(w) - 0.5) / (im * value(w) - 0.5); atol = 1e-6, rtol = 1e-6)

    # polynomial extrapolation for 1D grids with Float64 argument
    @test isapprox(f1(value(w)), 1.0 / (im * value(w) - 0.5); atol = 1e-6, rtol = 1e-6)
    @test isapprox(f2(value(w)), 1.0 / (im * value(w) - 0.5) / (im * value(w) - 0.5); atol = 1e-6, rtol = 1e-6)

    # constant extrapolation for higher-dimensional grids with MatsubaraFrequency argument
    @test f3(w, w)    ≈ 1.0 / (im * value(g[end]) - 0.5) / (im * value(g[end]) - 0.5)
    @test f3(w, g[1]) ≈ 1.0 / (im * value(g[end]) - 0.5) / (im * value(g[1]) - 0.5)
    @test f3(g[1], w) ≈ 1.0 / (im * value(g[1]) - 0.5) / (im * value(g[end]) - 0.5)

    # constant extrapolation for higher-dimensional grids with Float64 argument
    @test f3(value(w), value(w))    ≈ 1.0 / (im * value(g[end]) - 0.5) / (im * value(g[end]) - 0.5)
    @test f3(value(w), value(g[1])) ≈ 1.0 / (im * value(g[end]) - 0.5) / (im * value(g[1]) - 0.5)
    @test f3(value(g[1]), value(w)) ≈ 1.0 / (im * value(g[1]) - 0.5) / (im * value(g[end]) - 0.5)
end