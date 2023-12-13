@testset "BrillouinZoneMesh" begin 
    k1 = (2.0 * pi / 3) .* SVector{2, Float64}(1, +sqrt(3.0))
    k2 = (2.0 * pi / 3) .* SVector{2, Float64}(1, -sqrt(3.0))
    m  = BrillouinZoneMesh(BrillouinZone(6, k1, k2))

    # iterator
    @test reciprocals(m) ≈ [index(value(k)) for k in m]
    @test euclideans(m)  ≈ [euclidean(k, m) for k in m]

    # call to grid
    for trial in 1 : 10
        p = rand(points(m))
        @test m(value(p)) == index(p)
        @test m(euclidean(p, m)) == index(p)
    end 

    # periodic boundary conditions 
    for trial in 1 : 10 
        p  = rand(points(m))
        n  = SVector{2, Int64}(rand(-4 : 4), rand(-4 : 4))
        q1 = value(p) + BrillouinPoint(domain(m)[:bz].L .* n)
        q2 = euclidean(p, m) + basis(domain(m)[:bz]) * n
        @test m(q1) == index(p)
        @test m(q2) == index(p)
    end

    # io
    file = h5open("test.h5", "w")

    save!(file, "testBrillouin", m)
    mp = load_brillouin_zone_mesh(file, "testBrillouin")
    @test m == mp

    close(file)
    rm("test.h5")
end