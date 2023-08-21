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

    # benchmark vs analytic results
    @test isapprox(sum_me(f1), ρ10; atol = 1e-6, rtol = 1e-6)
    @test isapprox(sum_me(f2), ρ1p; atol = 1e-6, rtol = 1e-6)
    @test isapprox(sum_me(f3), ρ1m; atol = 1e-6, rtol = 1e-6)
    @test isapprox(sum_me(f4),  ρ2; atol = 1e-6, rtol = 1e-6)

    @test isapprox(density(f1), ρ(+0, T); atol = 1e-6, rtol = 1e-6)
    @test isapprox(density(f2), ρ(+ξ, T); atol = 1e-6, rtol = 1e-6)
    @test isapprox(density(f3), ρ(-ξ, T); atol = 1e-6, rtol = 1e-6)
end