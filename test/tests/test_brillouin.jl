@testset "BrillouinPoints" begin 
    k1 = (2.0 * pi / 3) .* SVector{2, Float64}(1, +sqrt(3.0))
    k2 = (2.0 * pi / 3) .* SVector{2, Float64}(1, -sqrt(3.0))
    m  = BrillouinZoneMesh(BrillouinZone(6, k1, k2))

    for trial in 1 : 10 
        x1    = m[rand(eachindex(m))]
        x2    = m[rand(eachindex(m))]
        q1    = value(x1)
        q2    = value(x2)
        q1vec = euclidean(q1, m)
        q2vec = euclidean(q2, m)

        # addition 
        @test euclidean(x1 + x2, m) ≈ q1vec + q2vec
        @test euclidean(x1 + q2, m) ≈ q1vec + q2vec
        @test euclidean(q1 + x2, m) ≈ q1vec + q2vec
        @test euclidean(q1 + q2, m) ≈ q1vec + q2vec

        # subtraction
        @test euclidean(x1 - x2, m) ≈ q1vec - q2vec
        @test euclidean(x1 - q2, m) ≈ q1vec - q2vec
        @test euclidean(q1 - x2, m) ≈ q1vec - q2vec
        @test euclidean(q1 - q2, m) ≈ q1vec - q2vec

        # reflection
        @test euclidean(-x1, m) ≈ -q1vec
        @test euclidean(-q1, m) ≈ -q1vec
        @test euclidean(-x2, m) ≈ -q2vec
        @test euclidean(-q2, m) ≈ -q2vec
    end
end

@testset "BrillouinZoneMesh" begin 
    k1 = (2.0 * pi / 3) .* SVector{2, Float64}(1, +sqrt(3.0))
    k2 = (2.0 * pi / 3) .* SVector{2, Float64}(1, -sqrt(3.0))
    m  = BrillouinZoneMesh(BrillouinZone(6, k1, k2))

    # iterator
    @test reciprocals(m) ≈ [plain_value(k) for k in m]
    @test euclideans(m)  ≈ [euclidean(k, m) for k in m]

    # mapping to mesh index
    for trial in 1 : 10
        p  = rand(points(m))
        n  = SVector{2, Int}(rand(-4 : 4), rand(-4 : 4))
        q1 = value(p) + BrillouinPoint(domain(m)[:bz].L .* n)
        q2 = euclidean(p, m) + basis(domain(m)[:bz]) * n

        @test MatsubaraFunctions.mesh_index(value(p), m) == index(p)
        @test MatsubaraFunctions.mesh_index(euclidean(p, m), m) == index(p)
        @test MatsubaraFunctions.mesh_index_bc(q1, m) == index(p)
        @test MatsubaraFunctions.mesh_index(q2, m) == index(p) # call with Vector of Float does backfolding anyways
    end 

    # io
    file = h5open(dirname(@__FILE__) * "/test.h5", "w")
    save!(file, "testBrillouin", m)

    mp = load_mesh(file, "testBrillouin")
    @test m == mp

    close(file)
    rm(dirname(@__FILE__) * "/test.h5")
end