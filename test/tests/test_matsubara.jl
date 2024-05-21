@testset "MatsubaraFrequencies" begin 
    for trial in 1 : 10
        wFermion = MatsubaraFrequency(1.0, rand(-100 : 100), Fermion)
        wBoson   = MatsubaraFrequency(1.0, rand(-100 : 100), Boson)

        # addition
        @test value(wFermion + wFermion) ≈ value(wFermion) + value(wFermion)
        @test value(wFermion + wBoson)   ≈ value(wFermion) + value(wBoson)
        @test value(wBoson + wFermion)   ≈ value(wBoson) + value(wFermion)
        @test value(wBoson + wBoson)     ≈ value(wBoson) + value(wBoson)

        # subtraction
        @test value(wFermion - wFermion) ≈ value(wFermion) - value(wFermion)
        @test value(wFermion - wBoson)   ≈ value(wFermion) - value(wBoson)
        @test value(wBoson - wFermion)   ≈ value(wBoson) - value(wFermion)
        @test value(wBoson - wBoson)     ≈ value(wBoson) - value(wBoson)

        # reflection 
        @test value(-wFermion) ≈ -value(wFermion)
        @test value(-wBoson)   ≈ -value(wBoson)
    end
end

@testset "MatsubaraMesh" begin 
    mFermion = MatsubaraMesh(1.0, 10, Fermion)
    mBoson   = MatsubaraMesh(1.0, 10, Boson)

    # length
    @test length(mFermion) == 20
    @test length(mBoson)   == 19

    # iterator
    @test indices(mFermion) ≈ [index(value(w)) for w in mFermion]
    @test indices(mBoson)   ≈ [index(value(w)) for w in mBoson]

    @test values(mFermion) ≈ [plain_value(w) for w in mFermion]
    @test values(mBoson)   ≈ [plain_value(w) for w in mBoson]

    # call to grid
    for trial in 1 : 10
        pFermion = rand(points(mFermion))
        @test MatsubaraFunctions.mesh_index(value(pFermion), mFermion) == index(pFermion)
        @test MatsubaraFunctions.mesh_index(plain_value(pFermion), mFermion) == index(pFermion)

        pBoson = rand(points(mBoson))
        @test MatsubaraFunctions.mesh_index(value(pBoson), mBoson) == index(pBoson)
        @test MatsubaraFunctions.mesh_index(plain_value(pBoson), mBoson) == index(pBoson)
    end 
    
    # io
    file = h5open(dirname(@__FILE__) * "/test.h5", "w")

    save!(file, "testFermion", mFermion)
    mFermion_p = load_mesh(file, "testFermion")
    @test mFermion == mFermion_p

    save!(file, "testBoson", mBoson)
    mBoson_p = load_mesh(file, "testBoson")
    @test mBoson == mBoson_p

    close(file)
    rm(dirname(@__FILE__) * "/test.h5")
end