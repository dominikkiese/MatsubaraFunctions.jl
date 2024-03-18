function index_test(f :: MeshFunction{DD, Q, AT}) where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}
    
    for trial in 1 : 10
        m_idxs = ntuple(i -> rand(eachindex(meshes(f, i))), DD)
        pts    = ntuple(i -> meshes(f, i)[m_idxs[i]], DD)
        l_idx  = LinearIndex(f, pts...)
        c_idx  = CartesianIndex(f, pts...)

        @test CartesianIndex(f, l_idx) == c_idx 
        @test LinearIndex(f, c_idx)    == l_idx
    
        args1 = to_meshes(f, l_idx)
        args2 = to_meshes(f, c_idx)

        @test args1 == pts
        @test args2 == pts

    end

    return nothing 
end

@testset "FuncIndex" begin 
    k1 = (2.0 * pi / 3) .* SVector{2, Float64}(1, +sqrt(3.0))
    k2 = (2.0 * pi / 3) .* SVector{2, Float64}(1, -sqrt(3.0))
    m1 = BrillouinZoneMesh(BrillouinZone(6, k1, k2))
    m2 = MatsubaraMesh(1.0, 10, Fermion)
    m3 = IndexMesh(7)

    index_test(MeshFunction(m1))
    index_test(MeshFunction(m1, m2))
    index_test(MeshFunction(m1, m3))
    index_test(MeshFunction((m1, m2, m3, m3)))
    index_test(MeshFunction(m3))
end