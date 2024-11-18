@testset "FuncConstructors" begin 
    k1    = (2.0 * pi / 3) .* SVector{2, Float64}(1, +sqrt(3.0))
    k2    = (2.0 * pi / 3) .* SVector{2, Float64}(1, -sqrt(3.0))
    m1    = BrillouinZoneMesh(BrillouinZone(6, k1, k2))
    m2    = MatsubaraMesh(1.0, 10, Fermion)
    m3    = IndexMesh(7)
    data1 = rand(length(m1))
    data2 = rand(length(m1), length(m2), length(m3))

    # from data, check if passed by reference
    f1     = MeshFunction((m1,), data1)
    f2     = MeshFunction((m1, m2, m3), data2)
    data1 .= rand(length(m1))
    data2 .= rand(length(m1), length(m2), length(m3))
    @test f1.data ≈ data1
    @test f2.data ≈ data2

    # from data view
    f1     = MeshFunction((m1,), view(data1, :))
    f2     = MeshFunction((m1, m2, m3), view(data2, :, :, :))
    data1 .= rand(length(m1))
    data2 .= rand(length(m1), length(m2), length(m3))
    @test f1.data ≈ data1
    @test f2.data ≈ data2
    
    # from meshes
    @test try
        MeshFunction(m1)
        true 
        catch 
        false 
    end

    @test try
        MeshFunction(m1, m2, m3)
        true 
        catch 
        false 
    end

    @test try
        MeshFunction(m1; data_t = Float64)
        true 
        catch 
        false 
    end

    # copy constructor 
    @test try
        MeshFunction(MeshFunction(m1, m2, m3))
        true 
        catch 
        false 
    end
end