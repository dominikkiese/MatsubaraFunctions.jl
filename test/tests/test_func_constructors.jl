@testset "FuncConstructors" begin 
    k1    = (2.0 * pi / 3) .* SVector{2, Float64}(1, +sqrt(3.0))
    k2    = (2.0 * pi / 3) .* SVector{2, Float64}(1, -sqrt(3.0))
    m1    = BrillouinZoneMesh(BrillouinZone(6, k1, k2))
    m2    = MatsubaraMesh(1.0, 10, Fermion)
    m3    = IndexMesh(7)
    data1 = rand(length(m1))
    data2 = rand(length(m1), 5, 5)

    # from data, check if passed by reference
    f1     = MeshFunction(m1, data1)
    f2     = MeshFunction(m1, (5, 5), data2)
    data1 .= rand(length(m1))
    data2 .= rand(length(m1), 5, 5)

    @test typeof(f1) == MeshFunction{1, 0, 1, Float64, Array{Float64, 1}}
    @test typeof(f2) == MeshFunction{1, 2, 3, Float64, Array{Float64, 3}}
    @test f1.data ≈ data1
    @test f2.data ≈ data2

    # from data view
    f1     = MeshFunction(m1, view(data1, :))
    f2     = MeshFunction(m1, (5, 5), view(data2, :, :, :))
    data1 .= rand(length(m1))
    data2 .= rand(length(m1), 5, 5)

    @test typeof(f1) == MeshFunction{1, 0, 1, Float64, SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}}
    @test typeof(f2) == MeshFunction{1, 2, 3, Float64, SubArray{Float64, 3, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}}, true}}
    @test f1.data ≈ data1
    @test f2.data ≈ data2
    
    # from meshes
    @test typeof(MeshFunction(m1)) == MeshFunction{1, 0, 1, ComplexF64, Array{ComplexF64, 1}}
    @test typeof(MeshFunction((m1, m2))) == MeshFunction{2, 0, 2, ComplexF64, Array{ComplexF64, 2}}
    @test typeof(MeshFunction((m1, m3))) == MeshFunction{2, 0, 2, ComplexF64, Array{ComplexF64, 2}}
    @test typeof(MeshFunction(m1, 1)) == MeshFunction{1, 1, 2, ComplexF64, Array{ComplexF64, 2}}
    @test typeof(MeshFunction(m1, 1, 1)) == MeshFunction{1, 2, 3, ComplexF64, Array{ComplexF64, 3}}
    @test typeof(MeshFunction((m1, m2), 1)) == MeshFunction{2, 1, 3, ComplexF64, Array{ComplexF64, 3}}
    @test typeof(MeshFunction((m1, m2), 1, 1)) == MeshFunction{2, 2, 4, ComplexF64, Array{ComplexF64, 4}}
    @test typeof(MeshFunction((m1, m3), 1)) == MeshFunction{2, 1, 3, ComplexF64, Array{ComplexF64, 3}}
    @test typeof(MeshFunction((m1, m3), 1, 1)) == MeshFunction{2, 2, 4, ComplexF64, Array{ComplexF64, 4}}
    @test typeof(MeshFunction(m1; data_t = Float64)) == MeshFunction{1, 0, 1, Float64, Array{Float64, 1}}

    # copy constructor 
    @test typeof(MeshFunction(MeshFunction((m1, m2, m3), 1, 1))) == MeshFunction{3, 2, 5, ComplexF64, Array{ComplexF64, 5}}
end