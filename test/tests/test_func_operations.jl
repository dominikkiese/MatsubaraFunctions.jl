@testset "FuncOperations" begin 
    mFermion = MatsubaraMesh(1.0, 10, Fermion)
    mBoson   = MatsubaraMesh(1.0, 10, Boson)
    f1       = MeshFunction((mFermion, mBoson), 2, 2)
    f2       = MeshFunction((mFermion, mBoson), 2, 2)
    f3       = MeshFunction((mFermion, mBoson), 2, 2)
    data1    = rand(ComplexF64, size(f1.data)...)
    data2    = rand(ComplexF64, size(f2.data)...)
    data3    = data1 .+ data2

    # init MeshFunction
    set!(f1, data1)
    set!(f2, data2)
    set!(f3, data3)

    # addition
    f4 = f1 + f2 
    @test f3 == f4 
    
    add!(f1, f2)
    @test f3 == f1

    set!(f1, data1)

    # subtraction
    f4 = f3 - f2 
    @test f4 == f1

    subtract!(f3, f2)
    @test f3 == f1

    set!(f3, data3)

    # multiplication
    set!(f2, 2.0 .* data1)

    f4 = 2.0 * f1 
    @test f4 == f2 

    set!(f2, data2)

    # flatten / unflatten
    x = flatten(f1)
    y = copy(x)
    flatten!(f2, y)

    unflatten!(f3, x)
    unflatten!(f4, y)

    @test f1 == f3 
    @test f2 == f4
end