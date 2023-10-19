@testset "Frequencies" begin 
    wFermion = MatsubaraFrequency(0.1, rand(-100 : 100), Fermion)
    wBoson   = MatsubaraFrequency(0.1, rand(-100 : 100), Boson)
    iFermion = MatsubaraIndex(wFermion)
    iBoson   = MatsubaraIndex(wBoson)

    # addition
    @test value(wFermion + wFermion) ≈ value(wFermion) + value(wFermion)
    @test value(wFermion + wBoson)   ≈ value(wFermion) + value(wBoson)
    @test value(wBoson + wFermion)   ≈ value(wBoson) + value(wFermion)
    @test value(wBoson + wBoson)     ≈ value(wBoson) + value(wBoson)

    @test value(MatsubaraFrequency(0.1, iFermion + iFermion)) ≈ value(wFermion) + value(wFermion)
    @test value(MatsubaraFrequency(0.1, iFermion + iBoson))   ≈ value(wFermion) + value(wBoson)
    @test value(MatsubaraFrequency(0.1, iBoson + iFermion))   ≈ value(wBoson) + value(wFermion)
    @test value(MatsubaraFrequency(0.1, iBoson + iBoson))     ≈ value(wBoson) + value(wBoson)

    # subtraction
    @test value(wFermion - wFermion) ≈ value(wFermion) - value(wFermion)
    @test value(wFermion - wBoson)   ≈ value(wFermion) - value(wBoson)
    @test value(wBoson - wFermion)   ≈ value(wBoson) - value(wFermion)
    @test value(wBoson - wBoson)     ≈ value(wBoson) - value(wBoson)

    @test value(MatsubaraFrequency(0.1, iFermion - iFermion)) ≈ value(wFermion) - value(wFermion)
    @test value(MatsubaraFrequency(0.1, iFermion - iBoson))   ≈ value(wFermion) - value(wBoson)
    @test value(MatsubaraFrequency(0.1, iBoson - iFermion))   ≈ value(wBoson) - value(wFermion)
    @test value(MatsubaraFrequency(0.1, iBoson - iBoson))     ≈ value(wBoson) - value(wBoson)

    # reflection 
    @test value(-wFermion) ≈ -value(wFermion)
    @test value(-wBoson)   ≈ -value(wBoson)

    @test value(MatsubaraFrequency(0.1, -iFermion)) ≈ -value(wFermion)
    @test value(MatsubaraFrequency(0.1, -iBoson))   ≈ -value(wBoson)
end