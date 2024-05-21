@testset "IndexMesh" begin 
    idxs = IndexMesh(10)

    # length
    @test length(idxs) == 10

    # iterator
    @test values(idxs) ≈ [plain_value(w) for w in idxs]

    # call to grid
    for trial in 1 : 10
        idx = rand(points(idxs))
        @test MatsubaraFunctions.mesh_index(value(idx), idxs) == plain_value(idx)
        @test MatsubaraFunctions.mesh_index(plain_value(idx), idxs) == plain_value(idx)
    end 
    
    # io
    file = h5open(dirname(@__FILE__) * "/test.h5", "w")

    save!(file, "testIndices", idxs)
    idxs_p = load_mesh(file, "testIndices")
    @test idxs == idxs_p

    close(file)
    rm(dirname(@__FILE__) * "/test.h5")
end