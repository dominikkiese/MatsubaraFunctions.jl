@testset "Symmetries" begin 
    g  = MatsubaraGrid(0.1, 128, Fermion)
    f1 = MatsubaraFunction(g)
    f2 = MatsubaraFunction(g)
    f3 = MatsubaraFunction(g)

    for v in g
        f1[v] = 1.0 / (im * value(v) - 0.5)
    end

    # complex conjugation for Green's function
    function conj(
        w :: Tuple{MatsubaraFrequency},
        x :: Tuple{}
        ) :: Tuple{Tuple{MatsubaraFrequency}, Tuple{}, MatsubaraOperation}

        return (-w[1],), (), MatsubaraOperation(false, true)
    end 

    # compute the symmetry group 
    SG = MatsubaraSymmetryGroup([MatsubaraSymmetry{1, 0}(conj)], f1)

    # symmetrize f2 and compare to f1 
    for class in SG.classes 
        f2[class[1][1]] = f1[class[1][1]]
    end 

    SG(f2)
    @test f2 == f1
    @test SG(f1) < 1e-14

    # symmetrize f3 and compare to f1 using MatsubaraInitFunction
    function init(
        w :: Tuple{MatsubaraFrequency},
        x :: Tuple{}
        ) :: ComplexF64

        return f1[w]
    end 

    InitFunc = MatsubaraInitFunction{1, 0, ComplexF64}(init)
    SG(f3, InitFunc; mode = :serial)
    @test f3 == f1

    SG(f3, InitFunc; mode = :polyester)
    @test f3 == f1

    SG(f3, InitFunc; mode = :threads)
    @test f3 == f1

    SG(f3, InitFunc; mode = :hybrid)
    @test f3 == f1
end