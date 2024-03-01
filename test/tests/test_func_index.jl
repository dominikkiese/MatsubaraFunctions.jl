function index_test(f :: MeshFunction{MD, SD, DD, Q, AT}) where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}
    
    for trial in 1 : 10
        m_idxs = ntuple(i -> rand(eachindex(meshes(f, i))), MD)
        x_idxs = ntuple(i -> rand(1 : shape(f, i)), SD)
        pts    = ntuple(i -> meshes(f, i)[m_idxs[i]], MD)
        l_idx  = LinearIndex(f, pts, x_idxs...)
        c_idx  = CartesianIndex(f, pts, x_idxs...)

        @test CartesianIndex(f, l_idx) == c_idx 
        @test LinearIndex(f, c_idx)    == l_idx
    
        args1 = to_meshes(f, l_idx)
        args2 = to_meshes(f, c_idx)

        @test first(args1) == pts
        @test first(args2) == pts

        @test last(args1) == x_idxs 
        @test last(args2) == x_idxs   
    end

    return nothing 
end

@testset "FuncIndex" begin 
    k1 = (2.0 * pi / 3) .* SVector{2, Float64}(1, +sqrt(3.0))
    k2 = (2.0 * pi / 3) .* SVector{2, Float64}(1, -sqrt(3.0))
    m1 = BrillouinZoneMesh(BrillouinZone(6, k1, k2))
    m2 = MatsubaraMesh(1.0, 10, Fermion)

    index_test(MeshFunction(m1))
    index_test(MeshFunction(m1, m2))
    index_test(MeshFunction(m1, 5))
    index_test(MeshFunction((m1, m2), 5, 5))
end