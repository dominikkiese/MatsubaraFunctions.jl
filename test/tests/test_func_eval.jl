# TO DO: add test for out of bounds access
function eval_test(f :: MeshFunction{DD, Q, AT}) where {DD, Q <: Number, AT <: AbstractArray{Q, DD}}

    for trial in 1 : 10
        m_idxs = ntuple(i -> rand(eachindex(meshes(f, i))), DD)
        #x_idxs = ntuple(i -> rand(1 : shape(f, i)), SD)
        pts    = ntuple(i -> meshes(f, i)[m_idxs[i]], DD)
        val    = f[m_idxs...]

        @test f[pts...]         ≈ val
        @test f[value.(pts)...] ≈ val
        @test f(pts...)         ≈ val
        @test f(value.(pts)...) ≈ val
    end

    return nothing 
end

@testset "FuncEval" begin 
    k1 = (2.0 * pi / 3) .* SVector{2, Float64}(1, +sqrt(3.0))
    k2 = (2.0 * pi / 3) .* SVector{2, Float64}(1, -sqrt(3.0))
    m1 = BrillouinZoneMesh(BrillouinZone(6, k1, k2))
    m2 = MatsubaraMesh(1.0, 10, Fermion)
    m3 = IndexMesh(7)

    f1 = MeshFunction(m1)
    set!(f1, rand(size(f1.data)...))
    eval_test(f1)

    f2 = MeshFunction(m1, m2)
    set!(f2, rand(size(f2.data)...))
    eval_test(f2)

    f3 = MeshFunction(m1, m3)
    set!(f3, rand(size(f3.data)...))
    eval_test(f3)

    f4 = MeshFunction(m1, m2, m3, m3)
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
    @test f4(first(m1), w, MatsubaraFunctions.Index(1), MatsubaraFunctions.Index(2)) ≈ ComplexF64(0.0) 
    @test f4(first(m1), w, MatsubaraFunctions.Index(3), MatsubaraFunctions.Index(5); lim = ComplexF64(1.0)) ≈ ComplexF64(1.0)
end