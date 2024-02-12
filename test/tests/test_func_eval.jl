# TO DO: add test for out of bounds access
function eval_test(f :: MeshFunction{MD, SD, DD, Q, AT}) where {MD, SD, DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    for trial in 1 : 10
        m_idxs = ntuple(i -> rand(eachindex(meshes(f, i))), MD)
        x_idxs = ntuple(i -> rand(1 : shape(f, i)), SD)
        pts    = ntuple(i -> meshes(f, i)[m_idxs[i]], MD)
        val    = f[m_idxs..., x_idxs...]

        @test f[pts, x_idxs...]         ≈ val
        @test f[value.(pts), x_idxs...] ≈ val
        @test f(pts, x_idxs...)         ≈ val
        @test f(value.(pts), x_idxs...) ≈ val
    end

    return nothing 
end

@testset "FuncEval" begin 
    k1 = (2.0 * pi / 3) .* SVector{2, Float64}(1, +sqrt(3.0))
    k2 = (2.0 * pi / 3) .* SVector{2, Float64}(1, -sqrt(3.0))
    m1 = BrillouinZoneMesh(BrillouinZone(6, k1, k2))
    m2 = MatsubaraMesh(1.0, 10, Fermion)

    f1 = MeshFunction(m1)
    set!(f1, rand(size(f1.data)...))
    eval_test(f1)

    f2 = MeshFunction(m1, m2)
    set!(f2, rand(size(f2.data)...))
    eval_test(f2)

    f3 = MeshFunction(m1, 5)
    set!(f3, rand(size(f3.data)...))
    eval_test(f3)

    f4 = MeshFunction((m1, m2), 5, 5)
    set!(f4, rand(size(f4.data)...))
    eval_test(f4)

    # test interpolation
    for trial in 1 : 10
        idxs = ntuple(n -> rand(eachindex(meshes(f2, n))), 2)
        pts  = ntuple(n -> meshes(f2, n)[idxs[n]], 2)
        val  = f2[idxs...]

        @test f2(euclidean(pts[1], m1), plain_value(pts[2])) ≈ val
    end

    # test out of bounds access 
    w = MatsubaraFrequency(1.0, 100, Fermion)
    @test f2(first(m1), w) ≈ ComplexF64(0.0) 
    @test f2(first(m1), w; lim = ComplexF64(1.0)) ≈ ComplexF64(1.0)
    @test f4((first(m1), w), 1, 2) ≈ ComplexF64(0.0) 
    @test f4((first(m1), w), 3, 5; lim = ComplexF64(1.0)) ≈ ComplexF64(1.0)
end