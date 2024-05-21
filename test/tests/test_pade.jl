@testset "Pade" begin
    # generate some dummy data
    as   = ntuple(x -> rand(BigFloat), 4)
    f(x) = as[1] / (1.0 + as[2] * x / (1.0  + as[3] * x / (1.0 + as[4] * x)))

    # build Pade approximant
    xdata = Vector{BigFloat}(0.01 : 0.01 : 1.0)
    ydata = f.(xdata)
    PA    = PadeApprox(xdata, ydata)

    @test length(coeffs(PA)) == 5
    @test PA.(xdata) ≈ ydata
end