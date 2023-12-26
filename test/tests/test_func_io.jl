@testset "IO" begin 
    g = MatsubaraMesh(0.1, 128, Fermion)
    f = MeshFunction(g)

    for v in g
        f[v] = 1.0 / (im * plain_value(v) - 0.5)
    end

    function conj(
        w :: Tuple{Union{MeshPoint{MatsubaraFrequency{PT}},MatsubaraFrequency{PT}}},
        x :: Tuple{}
        ) :: Tuple{Tuple{MatsubaraFrequency{PT}}, Tuple{}, MeshOperation} where{PT}

        return (-w[1],), (), MeshOperation(false, true)
    end 

    SG = MeshSymmetryGroup([MeshSymmetry{1, 0}(conj)], f)

    # write f and SG to file 
    file = h5open(dirname(@__FILE__) * "/test.h5", "w")
    save_mesh_function!(file, "f", f)
    save_mesh_symmetry_group!(file, "SG", SG)
    close(file)
    
    # load f and SG from file 
    file = h5open(dirname(@__FILE__) * "/test.h5", "r")
    fp  = load_mesh_function(file, "f")
    SGp = MeshSymmetryGroup{1, 0, 1, ComplexF64}(file, "SG")
    
    @test f == fp
    t = true 

    for (cl, clp) in zip(SG.classes, SGp.classes)
        if cl != clp 
            t = false 
            break 
        end
    end 

    @test t

    # close and remove test file 
    close(file)
    rm(dirname(@__FILE__) * "/test.h5")
end