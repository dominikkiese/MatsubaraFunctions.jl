@testset "IO" begin 
    g = MatsubaraMesh(0.1, 128, Fermion)
    f = MeshFunction(g)

    for v in g
        f[v] = 1.0 / (im * plain_value(v) - 0.5)
    end

    # write f and SG to file 
    file = h5open(dirname(@__FILE__) * "/test.h5", "w")
    save_mesh_function!(file, "f", f)

    # load f from file 
    fp = load_mesh_function(file, "f")
    @test f == fp

    # close and remove test file 
    close(file)
    rm(dirname(@__FILE__) * "/test.h5")
end