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
    g = MatsubaraGrid(1.0, 100, Fermion)
   
    index_test(MatsubaraFunction(g, 1))
    index_test(MatsubaraFunction(g, 10))
    index_test(MatsubaraFunction((g, g), 10))
    index_test(MatsubaraFunction((g, g), (10, 10)))
end