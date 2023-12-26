@testset "Symmetries" begin 
    g  = MatsubaraMesh(0.1, 128, Fermion)
    f1 = MeshFunction(g)
    f2 = MeshFunction(g)
    f3 = MeshFunction(g)

    for v in g
        f1[v] = 1.0 / (im * value(value(v)) - 0.5)
    end

    # complex conjugation for Green's function
    function conj(
        w :: Tuple{MeshPoint{MatsubaraFrequency{PT}}},
        x :: Tuple{}
        ) :: Tuple{Tuple{MeshPoint{MatsubaraFrequency{PT}}}, Tuple{}, MeshOperation} where{PT}

        return (-w[1],), (), MeshOperation(false, true)
    end 


    # compute the symmetry group 
    SG = MeshSymmetryGroup([MeshSymmetry{1, 0}(conj)], f1)

    # symmetrize f2 and compare to f1 
    for class in SG.classes 
        f2[class[1][1]] = f1[class[1][1]]
    end 

    SG(f2)
    @test f2 == f1
    @test SG(f1) < 1e-14

    # reduce f1, then init f2 and compare
    f1vec = get_reduced(SG, f1) 
    init_from_reduced!(SG, f2, f1vec)
    @test f2 == f1

    # symmetrize f3 and compare to f1 using MeshInitFunction
    function init(
        w :: Tuple{MeshPoint{MatsubaraFrequency{PT}}},
        x :: Tuple{}
        ) :: ComplexF64 where{PT}

        return f1[w]
    end 

    InitFunc = MeshInitFunction{1, 0, ComplexF64}(init)
    SG(f3, InitFunc; mode = :serial)
    @test f3 == f1

    SG(f3, InitFunc; mode = :polyester)
    @test f3 == f1

    SG(f3, InitFunc; mode = :threads)
    @test f3 == f1

    SG(f3, InitFunc; mode = :hybrid)
    @test f3 == f1
end