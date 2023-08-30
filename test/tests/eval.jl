function eval_test(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: Nothing where {GD, SD, DD, Q <: Number}

    g_idxs = ntuple(i -> rand(eachindex(grids(f, i))), GD)
    x_idxs = ntuple(i -> rand(1 : shape(f, i)), SD)
    freqs  = ntuple(i -> grids(f, i)[g_idxs[i]], GD)
    val    = f[g_idxs..., x_idxs...]

    @test f[freqs, x_idxs...]         ≈ val 
    @test f(freqs, x_idxs...)         ≈ val 
    @test f(value.(freqs), x_idxs...) ≈ val 

    return nothing 
end

@testset "Evaluate" begin 
    g = MatsubaraGrid(1.0, 10, Fermion)
    n = length(g)

    eval_test(MatsubaraFunction((g,), (1,), rand(n, 1)))
    eval_test(MatsubaraFunction((g,), (10,), rand(n, 10)))
    eval_test(MatsubaraFunction((g, g), (10,), rand(n, n, 10)))
    eval_test(MatsubaraFunction((g, g), (10, 10), rand(n, n, 10, 10)))
end