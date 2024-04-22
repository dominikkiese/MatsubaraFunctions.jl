@testset "Constructors" begin 
    gf = MatsubaraGrid(1.0, 1, Fermion)
    gb = MatsubaraGrid(1.0, 1, Boson)

    # scalar-valued
    @test typeof(MatsubaraFunction(gf)) == MatsubaraFunction{1, 0, 1, ComplexF64}
    @test typeof(MatsubaraFunction((gf, gb))) == MatsubaraFunction{2, 0, 2, ComplexF64}

    # tensor-valued
    @test typeof(MatsubaraFunction(gf, 1)) == MatsubaraFunction{1, 1, 2, ComplexF64}
    @test typeof(MatsubaraFunction(gf, 1, 1)) == MatsubaraFunction{1, 2, 3, ComplexF64}
    @test typeof(MatsubaraFunction((gf, gb), 1)) == MatsubaraFunction{2, 1, 3, ComplexF64}
    @test typeof(MatsubaraFunction((gf, gb), 1, 1)) == MatsubaraFunction{2, 2, 4, ComplexF64}

    # change data type 
    @test typeof(MatsubaraFunction(gf; data_t = Float64)) == MatsubaraFunction{1, 0, 1, Float64}
end