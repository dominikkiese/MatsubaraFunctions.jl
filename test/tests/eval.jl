function eval_test(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: Nothing where {GD, SD, DD, Q <: Number}

    g_idxs = ntuple(i -> rand(eachindex(grids(f, i))), GD)
    x_idxs = ntuple(i -> rand(axes(f.data, GD+i)), SD)
    freqs  = ntuple(i -> grids(f, i)[g_idxs[i]], GD)
    val    = f[g_idxs..., x_idxs...]

    @test f[freqs, x_idxs...]         ≈ val 
    @test f(freqs, x_idxs...)         ≈ val 
    @test f(value.(freqs), x_idxs...) ≈ val 

    return nothing 
end

@testset "Evaluate" begin 
    gf = MatsubaraGrid(1.0, 5, Fermion)
    gb = MatsubaraGrid(1.0, 5, Boson)
    nf = length(gf)
    nb = length(gb)

    # scalar-valued
    f1 = MatsubaraFunction(gf)
    set!(f1, rand(nf))
    eval_test(f1)

    f2 = MatsubaraFunction((gf, gb))
    set!(f2, rand(nf, nb))
    eval_test(f2)

    # tensor-valued
    f3 = MatsubaraFunction(gf, 5)
    set!(f3, rand(nf, 5))
    eval_test(f3)

    f4 = MatsubaraFunction(gf, 5, 5)
    set!(f4, rand(nf, 5, 5))
    eval_test(f4)

    f5 = MatsubaraFunction((gf, gb), 5)
    set!(f5, rand(nf, nb, 5))
    eval_test(f5)

    f6 = MatsubaraFunction((gf, gb), 5, 5)
    set!(f6, rand(nf, nb, 5, 5))
    eval_test(f6)
end