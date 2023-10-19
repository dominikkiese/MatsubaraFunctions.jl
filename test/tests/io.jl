@testset "IO" begin 
    g = MatsubaraGrid(0.1, 128, Fermion)
    f = MatsubaraFunction(g)

    for v in g
        f[v] = 1.0 / (im * value(v) - 0.5)
    end

    function conj(
        w :: Tuple{MatsubaraFrequency},
        x :: Tuple{}
        ) :: Tuple{Tuple{MatsubaraFrequency}, Tuple{}, MatsubaraOperation}

        return (-w[1],), (), MatsubaraOperation(false, true)
    end 

    SG = MatsubaraSymmetryGroup([MatsubaraSymmetry{1, 0}(conj)], f)

    # write f and SG to file 
    file = h5open(dirname(@__FILE__) * "/test.h5", "w")
    save_matsubara_function!(file, "f", f)
    save_matsubara_symmetry_group!(file, "SG", SG)

    # load f and SG from file 
    fp  = load_matsubara_function(file, "f")
    SGp = MatsubaraSymmetryGroup{1, 0, 1, ComplexF64}(file, "SG")
    
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