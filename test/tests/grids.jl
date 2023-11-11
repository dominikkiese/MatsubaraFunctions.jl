@testset "Grids" begin 
    wFermion = MatsubaraGrid(0.1, 100, Fermion)
    wBoson   = MatsubaraGrid(0.1, 100, Boson)

    # length
    @test length(wFermion) == 200
    @test length(wBoson)   == 199

    # iterator
    @test [values(wFermion)...] ≈ Float64[value(w) for w in wFermion]
    @test [values(wBoson)...]   ≈ Float64[value(w) for w in wBoson]

    wFermion_idx = rand(eachindex(wFermion))
    wBoson_idx   = rand(eachindex(wBoson))

    # call to fermionic grid
    @test wFermion(wFermion[wFermion_idx])                 == wFermion_idx
    @test wFermion(MatsubaraIndex(wFermion[wFermion_idx])) == wFermion_idx
    @test wFermion(value(wFermion[wFermion_idx]))          == wFermion_idx

    # call to bosonic grid
    @test wBoson(wBoson[wBoson_idx])                 == wBoson_idx
    @test wBoson(MatsubaraIndex(wBoson[wBoson_idx])) == wBoson_idx
    @test wBoson(value(wBoson[wBoson_idx]))          == wBoson_idx
end