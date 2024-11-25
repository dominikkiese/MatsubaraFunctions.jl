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
    function conj(w :: Tuple{BrillouinPoint{2}, MatsubaraFrequency})
        return (w[1], -w[2]), Operation{ComplexF64}(sgn = false, con = true)
    end 

    function conj(w :: Tuple{MatsubaraFrequency})
        return (-w[1],), Operation{ComplexF64}(sgn = false, con = true)
    end 

    # reflection along zone diagonal
    function ref(w :: Tuple{BrillouinPoint{2}, MatsubaraFrequency})
        return (BrillouinPoint(value(w[1])[2], value(w[1])[1]), w[2]), Operation{ComplexF64}()
    end 

    function ref(w :: Tuple{BrillouinPoint{2}})
        return (BrillouinPoint(value(w[1])[2], value(w[1])[1]),), Operation{ComplexF64}()
    end 

    # rotation by π/2
    function rot(w :: Tuple{BrillouinPoint{2}, MatsubaraFrequency}, g :: Mesh{MeshPoint{BrillouinPoint{2}}})
        kp = fold_back(BrillouinPoint(value(w[1])[2], -value(w[1])[1]), g)
        return (kp, w[2]), Operation{ComplexF64}()
    end 

    function rot(w :: Tuple{BrillouinPoint{2}}, g :: Mesh{MeshPoint{BrillouinPoint{2}}})
        kp = fold_back(BrillouinPoint(value(w[1])[2], -value(w[1])[1]), g)
        return (kp,), Operation{ComplexF64}()
    end 

    # compute the symmetry group with and without factorization
    symmetries = Vector{<:AbstractSymmetry}[[Symmetry{1}(ref), Symmetry{1}(x -> rot(x, g1))], [Symmetry{1}(conj)]]
    mesh_idxs  = Vector{Int}[[1], [2]]
    SG         = SymmetryGroup([Symmetry{2}(conj), Symmetry{2}(ref), Symmetry{2}((w) -> rot(w, g1))], f1)
    SG_        = SymmetryGroup(symmetries, mesh_idxs, f1)

    # test if equivalent (then continue using SG)
    idxs = vcat([x.data_id for x in SG.classes]...)
    @test idxs == unique(idxs)
    @test length(idxs) == length(f1)

    idxs_ = vcat([x.data_id for x in SG_.classes]...)
    @test idxs_ == unique(idxs_)
    @test length(idxs_) == length(f1)

    for SC in SG 
        found_match = false
        p           = sortperm(SC.data_id)

        for SC_ in SG_
            p_ = sortperm(SC_.data_id)

            if (SC.data_id[p] == SC_.data_id[p_]) && (SC.data_op[p] == SC_.data_op[p_])
                found_match = true
                break
            end
        end

        @test found_match
    end

    # symmetrize f2 and compare to f1 
    for class in SG
        idx     = first(first(class))
        f2[idx] = f1[idx]
    end 

    SG(f2)
    @test f2 == f1

    # reduce f1, then init f2 and compare
    f1vec = flatten(SG, f1) 
    unflatten!(SG, f2, f1vec)
    @test f2 == f1

    # symmetrize f3 and compare to f1 using InitFunction
    init(w :: Tuple{BrillouinPoint{2}, MatsubaraFrequency}) :: ComplexF64 = f1[w...]
    InitFunc = InitFunction{2, ComplexF64}(init)

    SG(f3, InitFunc; mode = :serial)
    @test f3 == f1

    SG(f3, InitFunc; mode = :threads)
    @test f3 == f1

    SG(f3, InitFunc; mode = :hybrid)
    @test f3 == f1
end