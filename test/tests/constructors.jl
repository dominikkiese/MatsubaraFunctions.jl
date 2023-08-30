@testset "Constructors" begin 
    g = MatsubaraGrid(1.0, 10, Fermion)

    @test typeof(MatsubaraFunction((g, g), (2, 2), Float64)) == MatsubaraFunction{2, 2, 4, Float64}
    @test typeof(MatsubaraFunction((g, g), (2, 2))) == MatsubaraFunction{2, 2, 4, ComplexF64}
    @test typeof(MatsubaraFunction(g, (2, 2), Float64)) == MatsubaraFunction{1, 2, 3, Float64}
    @test typeof(MatsubaraFunction(g, (2, 2))) == MatsubaraFunction{1, 2, 3, ComplexF64}
    @test typeof(MatsubaraFunction(g, 2, Float64)) == MatsubaraFunction{1, 1, 2, Float64}
    @test typeof(MatsubaraFunction(g, 2)) == MatsubaraFunction{1, 1, 2, ComplexF64}
    @test typeof(MatsubaraFunction(MatsubaraFunction(g, (2, 2)))) == MatsubaraFunction{1, 2, 3, ComplexF64}
end