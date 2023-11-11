function Base.:isapprox(
    t1 :: Tuple,
    t2 :: Tuple
    )  :: Bool 

    for (a1, a2) in zip(t1, t2)
        if !(a1 ≈ a2)
            return false 
        end 
    end 

    return true 
end

function index_test(
    f :: MatsubaraFunction{GD, SD, DD, Q}
    ) :: Nothing where {GD, SD, DD, Q <: Number}

    # test for MatsubaraFrequency
    g_idxs = ntuple(i -> rand(eachindex(grids(f, i))), GD)
    x_idxs = ntuple(i -> rand(1 : shape(f, i)), SD)
    freqs  = ntuple(i -> grids(f, i)[g_idxs[i]], GD)
    l_idx  = LinearIndex(f, freqs, x_idxs...)
    c_idx  = CartesianIndex(f, freqs, x_idxs...)

    @test CartesianIndex(f, l_idx) == c_idx 
    @test LinearIndex(f, c_idx)    == l_idx
 
    w1 = to_Matsubara(f, l_idx)
    w2 = to_Matsubara(f, c_idx)

    @test value.(first(w1)) ≈ value.(freqs)
    @test value.(first(w2)) ≈ value.(freqs)

    @test last(w1) == x_idxs 
    @test last(w2) == x_idxs   

    # test for MatsubaraIndex
    idxs  = MatsubaraIndex.(freqs)
    l_idx = LinearIndex(f, idxs, x_idxs...)
    c_idx = CartesianIndex(f, idxs, x_idxs...)

    @test CartesianIndex(f, l_idx) == c_idx 
    @test LinearIndex(f, c_idx)    == l_idx

    w1 = to_Matsubara(f, l_idx)
    w2 = to_Matsubara(f, c_idx)

    @test value.(first(w1)) ≈ value.(freqs)
    @test value.(first(w2)) ≈ value.(freqs)

    @test last(w1) == x_idxs 
    @test last(w2) == x_idxs  

    return nothing 
end

@testset "Index" begin 
    gf = MatsubaraGrid(1.0, 5, Fermion)
    gb = MatsubaraGrid(1.0, 5, Boson)

    # scalar-valued
    f1 = MatsubaraFunction(gf)
    index_test(f1)

    f2 = MatsubaraFunction((gf, gb))
    index_test(f2)

    # tensor-valued
    f3 = MatsubaraFunction(gf, 5)
    index_test(f3)

    f4 = MatsubaraFunction(gf, 5, 5)
    index_test(f4)

    f5 = MatsubaraFunction((gf, gb), 5)
    index_test(f5)

    f6 = MatsubaraFunction((gf, gb), 5, 5)
    index_test(f6)
end