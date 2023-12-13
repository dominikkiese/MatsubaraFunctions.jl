@testset "MatsubaraFrequencies" begin 
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

@testset "MatsubaraMesh" begin 
    mFermion = MatsubaraMesh(1.0, 10, Fermion)
    mBoson   = MatsubaraMesh(1.0, 10, Boson)

    # iterator
    @test indices(mFermion) ≈ [index(value(w)) for w in mFermion]
    @test indices(mBoson)   ≈ [index(value(w)) for w in mBoson]

    @test values(mFermion) ≈ [value(value(w)) for w in mFermion]
    @test values(mBoson)   ≈ [value(value(w)) for w in mBoson]

    # call to grid
    for trial in 1 : 10
        pFermion = rand(points(mFermion))
        @test mFermion(value(pFermion)) == index(pFermion)
        @test mFermion(value(value(pFermion))) == index(pFermion)

        pBoson = rand(points(mBoson))
        @test mBoson(value(pBoson)) == index(pBoson)
        @test mBoson(value(value(pBoson))) == index(pBoson)
    end 
    
    # io
    file = h5open("test.h5", "w")

    save!(file, "testFermion", mFermion)
    mFermion_p = load_matsubara_mesh(file, "testFermion")
    @test mFermion == mFermion_p

    save!(file, "testBoson", mBoson)
    mBoson_p = load_matsubara_mesh(file, "testBoson")
    @test mBoson == mBoson_p

    close(file)
    rm("test.h5")
end