@testset "FuncConstructors" begin 
    k1 = (2.0 * pi / 3) .* SVector{2, Float64}(1, +sqrt(3.0))
    k2 = (2.0 * pi / 3) .* SVector{2, Float64}(1, -sqrt(3.0))
    m1 = BrillouinZoneMesh(BrillouinZone(6, k1, k2))
    m2 = MatsubaraMesh(1.0, 10, Fermion)

    # scalar-valued
    @test typeof(MeshFunction(m1)) == MeshFunction{1, 0, 1, ComplexF64}
    @test typeof(MeshFunction((m1, m2))) == MeshFunction{2, 0, 2, ComplexF64}

    # tensor-valued
    @test typeof(MeshFunction(m1, 1)) == MeshFunction{1, 1, 2, ComplexF64}
    @test typeof(MeshFunction(m1, 1, 1)) == MeshFunction{1, 2, 3, ComplexF64}
    @test typeof(MeshFunction((m1, m2), 1)) == MeshFunction{2, 1, 3, ComplexF64}
    @test typeof(MeshFunction((m1, m2), 1, 1)) == MeshFunction{2, 2, 4, ComplexF64}

    # change data type 
    @test typeof(MeshFunction(m1; data_t = Float64)) == MeshFunction{1, 0, 1, Float64}
end