@testset "BrillouinZoneMesh" begin 
    k1 = (2.0 * pi / 3) .* SVector{2, Float64}(1, +sqrt(3.0))
    k2 = (2.0 * pi / 3) .* SVector{2, Float64}(1, -sqrt(3.0))
    m  = BrillouinZoneMesh(BrillouinZone(6, k1, k2))

    # iterator
    @test reciprocals(m) ≈ [index(value(k)) for k in m]
    @test euclideans(m)  ≈ [euclidean(k, m) for k in m]

    # mapping to mesh index
    for trial in 1 : 10
        p  = rand(points(m))
        n  = SVector{2, Int64}(rand(-4 : 4), rand(-4 : 4))
        q1 = value(p) + BrillouinPoint(domain(m)[:bz].L .* n)
        q2 = euclidean(p, m) + basis(domain(m)[:bz]) * n

        @test MatsubaraFunctions.mesh_index(value(p), m) == index(p)
        @test MatsubaraFunctions.mesh_index(euclidean(p, m), m) == index(p)
        @test MatsubaraFunctions.mesh_index_bc(q1, m) == index(p)
        @test MatsubaraFunctions.mesh_index(q2, m) == index(p) # call with Vector of Float needs to do backfolding anyways
    end 

    # io
    file = h5open("test.h5", "w")

    save!(file, "testBrillouin", m)
    mp = load_brillouin_zone_mesh(file, "testBrillouin")
    @test m == mp

    close(file)
    rm("test.h5")
end