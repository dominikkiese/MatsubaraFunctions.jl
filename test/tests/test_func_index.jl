function Base.:isapprox(
    t1 :: Tuple,
    t2 :: Tuple
    )  :: Bool 

    for (a1, a2) in zip(t1, t2)
        if !(a1 ≈ a2)
            return false 
        end 
    end 

    return true 
end

function Base.:isapprox(
    t1 :: BrillouinPoint{N},
    t2 :: BrillouinPoint{N}
    )  :: Bool where{N}

    return all(t1.index ≈ t2.index)
end


function Base.:isapprox(
    t1 :: MatsubaraFrequency{PT},
    t2 :: MatsubaraFrequency{PT}
    )  :: Bool where{PT}

    return t1.temperature ≈ t2.temperature && t1.index == t2.index && t1.value ≈ t2.value
end

function index_test(
    f :: MeshFunction{MD, SD, DD, Q}
    ) :: Nothing where {MD, SD, DD, Q <: Number}
    
    for trial in 1 : 10
        # get random Meshpoints and miscellaneous indices
        m_idxs = ntuple(i -> rand(eachindex(meshes(f, i))), MD)
        x_idxs = ntuple(i -> rand(axes(f, MD + i)), SD)
        
        # test for MeshPoint
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
        
        @test value.(first(args1)) ≈ value.(pts)
        @test value.(first(args2)) ≈ value.(pts)

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
    index_test(MeshFunction(m1, 5, 5))
    index_test(MeshFunction((m1, m2), 5))
    index_test(MeshFunction((m1, m2), 5, 5))
end