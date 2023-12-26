function eval_test(
    f :: MeshFunction{MD, SD, DD, Q}
    ) :: Nothing where {MD, SD, DD, Q <: Number}

    for trial in 1 : 10
        m_idxs = ntuple(i -> rand(eachindex(meshes(f, i))), MD)
        x_idxs = ntuple(i -> rand(axes(f.data, MD + i)), SD)
        pts    = ntuple(i -> meshes(f, i)[m_idxs[i]], MD)
        val    = f[m_idxs..., x_idxs...]
        #println("m_idxs: ", m_idxs)
        @test f[pts, x_idxs...]         ≈ val
        @test f[value.(pts), x_idxs...] ≈ val
        @test f(pts, x_idxs...)         ≈ val
        @test f(value.(pts), x_idxs...) ≈ val
    end

    return nothing 
end

@testset "FuncEval" begin 
begin
    k1 = (2.0 * pi / 3) .* SVector{2, Float64}(1, +sqrt(3.0))
    k2 = (2.0 * pi / 3) .* SVector{2, Float64}(1, -sqrt(3.0))
    m1 = BrillouinZoneMesh(BrillouinZone(6, k1, k2))
    m2 = MatsubaraMesh(1.0, 10, Fermion)
end

    # scalar-valued
    f1 = MeshFunction(m1)
    set!(f1, rand(size(f1.data)...))
    eval_test(f1)

    f2 = MeshFunction(m1, m2)
    set!(f2, rand(size(f2.data)...))
    eval_test(f2)

    # tensor-valued
    f3 = MeshFunction(m1, 5)
    set!(f3, rand(size(f3.data)...))
    eval_test(f3)

    f4 = MeshFunction((m1, m2), 5, 5)
    set!(f4, rand(size(f4.data)...))
    eval_test(f4)

    f5 = MeshFunction((m1, m2), 5)
    set!(f5, rand(size(f5.data)...))
    eval_test(f5)

    f6 = MeshFunction((m1, m2), 5, 5)
    set!(f6, rand(size(f6.data)...))
    eval_test(f6)

    # test interpolation
    for trial in 1 : 10
        idxs = ntuple(n -> rand(eachindex(meshes(f2, n))), 2)
        pts  = ntuple(n -> meshes(f2, n)[idxs[n]], 2)
        val  = f2[idxs...]
        
        @test f2(euclidean(pts[1], m1), plain_value(pts[2])) ≈ val
    end
end