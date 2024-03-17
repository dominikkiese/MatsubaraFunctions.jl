@testset "Symmetries" begin 
    # build BrillouinZone for square lattice
    k1 = (2.0 * pi) .* SVector{2, Float64}(1, 0)
    k2 = (2.0 * pi) .* SVector{2, Float64}(0, 1)
    bz = BrillouinZone(6, k1, k2)

    # build dummy mesh functions
    g1 = BrillouinZoneMesh(bz)
    g2 = MatsubaraMesh(1.0, 128, Fermion)
    f1 = MeshFunction(g1, g2)
    f2 = MeshFunction(g1, g2)
    f3 = MeshFunction(g1, g2)

    # init f1
    for k in meshes(f1, 1), v in meshes(f1, 2)
        p        = euclidean(k, g1)
        f1[k, v] = 1.0 / (im * plain_value(v) - cos(p[1]) - cos(p[2]))
    end

    # complex conjugation for Green's function
    function conj(w :: Tuple{BrillouinPoint{2}, MatsubaraFrequency}, x :: Tuple{})
        return (w[1], -w[2]), (), Operation(sgn = false, con = true)
    end 

    # reflection along zone diagonal
    function ref(w :: Tuple{BrillouinPoint{2}, MatsubaraFrequency}, x :: Tuple{})
        return (BrillouinPoint(value(w[1])[2], value(w[1])[1]), w[2]), (), Operation()
    end 

    # rotation by π/2
    function rot(w :: Tuple{BrillouinPoint{2}, MatsubaraFrequency}, x :: Tuple{}, g :: Mesh{MeshPoint{BrillouinPoint{2}}})
        kp = fold_back(BrillouinPoint(value(w[1])[2], -value(w[1])[1]), g)
        return (kp, w[2]), (), Operation()
    end 

    # compute the symmetry group 
    SG = SymmetryGroup([Symmetry{2, 0}(conj), Symmetry{2, 0}(ref), Symmetry{2, 0}((w, x) -> rot(w, x, g1))], f1)

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

    # symmetrize f3 and compare to f1 using InitFunction
    init(w :: Tuple{BrillouinPoint{2}, MatsubaraFrequency}, x :: Tuple{}) :: ComplexF64 = f1[w]
    InitFunc = InitFunction{2, 0, ComplexF64}(init)

    SG(f3, InitFunc; mode = :serial)
    @test f3 == f1

    SG(f3, InitFunc; mode = :polyester)
    @test f3 == f1

    SG(f3, InitFunc; mode = :threads)
    @test f3 == f1

    SG(f3, InitFunc; mode = :hybrid)
    @test f3 == f1
end