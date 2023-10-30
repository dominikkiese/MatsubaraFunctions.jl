@testset "Matsubara" begin 
    mFermion = MatsubaraMesh(1.0, 10, Fermion)
    mBoson   = MatsubaraMesh(1.0, 10, Boson)

    for trial in 1 : 10
        pFermion = rand(points(mFermion))
        @test mFermion(value(pFermion)) == index(pFermion)
        @test mFermion(value(value(pFermion))) == index(pFermion)

        pBoson = rand(points(mBoson))
        @test mBoson(value(pBoson)) == index(pBoson)
        @test mBoson(value(value(pBoson))) == index(pBoson)
    end 
    
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