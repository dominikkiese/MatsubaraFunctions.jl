@testset "FuncConstructors" begin 
    mFermion = MatsubaraMesh(1.0, 10, Fermion)
    mBoson   = MatsubaraMesh(1.0, 10, Boson)

    # scalar-valued
    @test typeof(MeshFunction(mFermion)) == MeshFunction{1, 0, 1, ComplexF64}
    @test typeof(MeshFunction((mFermion, mBoson))) == MeshFunction{2, 0, 2, ComplexF64}

    # tensor-valued
    @test typeof(MeshFunction(mFermion, 1)) == MeshFunction{1, 1, 2, ComplexF64}
    @test typeof(MeshFunction(mFermion, 1, 1)) == MeshFunction{1, 2, 3, ComplexF64}
    @test typeof(MeshFunction((mFermion, mBoson), 1)) == MeshFunction{2, 1, 3, ComplexF64}
    @test typeof(MeshFunction((mFermion, mBoson), 1, 1)) == MeshFunction{2, 2, 4, ComplexF64}

    # change data type 
    @test typeof(MeshFunction(mFermion; data_t = Float64)) == MeshFunction{1, 0, 1, Float64}
end