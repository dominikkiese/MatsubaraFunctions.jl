function index_test(
    f :: MeshFunction{MD, SD, DD, Q}
    ) :: Nothing where {MD, SD, DD, Q <: Number}
    
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

    return nothing 
end

@testset "FuncIndex" begin 
    mFermion = MatsubaraMesh(1.0, 10, Fermion)
    mBoson   = MatsubaraMesh(1.0, 10, Boson)

    # scalar-valued
    index_test(MeshFunction(mFermion))
    index_test(MeshFunction(mFermion, mBoson))

    # tensor-valued
    index_test(MeshFunction(mFermion, 5))
    index_test(MeshFunction((mFermion, mBoson), 5, 5))
end