@testset "FuncConstructors" begin 
    k1    = (2.0 * pi / 3) .* SVector{2, Float64}(1, +sqrt(3.0))
    k2    = (2.0 * pi / 3) .* SVector{2, Float64}(1, -sqrt(3.0))
    m1    = BrillouinZoneMesh(BrillouinZone(6, k1, k2))
    m2    = MatsubaraMesh(1.0, 10, Fermion)
    data1 = rand(length(m1))
    data2 = rand(length(m1), 5, 5)
    data3 = rand(length(m1), length(m2), 5, 5)

    # from data, check if parsed by reference
    f1     = MeshFunction(m1, data1)
    f2     = MeshFunction(m1, (5, 5), data2)
    f3     = MeshFunction((m1, m2), (5, 5), data3)
    data1 .= rand(length(m1))
    data2 .= rand(length(m1), 5, 5)
    data3 .= rand(length(m1), length(m2), 5, 5)

    @test typeof(f1) == MeshFunction{1, 0, 1, Float64}
    @test typeof(f2) == MeshFunction{1, 2, 3, Float64}
    @test typeof(f3) == MeshFunction{2, 2, 4, Float64}
    @test f1.data ≈ data1
    @test f2.data ≈ data2
    @test f3.data ≈ data3
    
    # from meshes
    @test typeof(MeshFunction(m1)) == MeshFunction{1, 0, 1, ComplexF64}
    @test typeof(MeshFunction((m1, m2))) == MeshFunction{2, 0, 2, ComplexF64}
    @test typeof(MeshFunction(m1, 1)) == MeshFunction{1, 1, 2, ComplexF64}
    @test typeof(MeshFunction(m1, 1, 1)) == MeshFunction{1, 2, 3, ComplexF64}
    @test typeof(MeshFunction((m1, m2), 1)) == MeshFunction{2, 1, 3, ComplexF64}
    @test typeof(MeshFunction((m1, m2), 1, 1)) == MeshFunction{2, 2, 4, ComplexF64}
    @test typeof(MeshFunction(m1; data_t = Float64)) == MeshFunction{1, 0, 1, Float64}

    # copy constructor 
    @test typeof(MeshFunction(MeshFunction((m1, m2), 1, 1))) == MeshFunction{2, 2, 4, ComplexF64}
end