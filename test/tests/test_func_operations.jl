@testset "FuncOperations" begin 
    k1    = (2.0 * pi / 3) .* SVector{2, Float64}(1, +sqrt(3.0))
    k2    = (2.0 * pi / 3) .* SVector{2, Float64}(1, -sqrt(3.0))
    m1    = BrillouinZoneMesh(BrillouinZone(6, k1, k2))
    m2    = MatsubaraMesh(1.0, 10, Fermion)
    m3    = IndexMesh(10)
    f1    = MeshFunction(m1, m2, m3, m3)
    f2    = MeshFunction(m1, m2, m3, m3)
    f3    = MeshFunction(m1, m2, m3, m3)
    data1 = rand(ComplexF64, size(f1.data)...)
    data2 = rand(ComplexF64, size(f2.data)...)
    data3 = data1 .+ data2

    # init MeshFunction
    set!(f1, data1)
    set!(f2, data2)
    set!(f3, data3)

    # addition
    f4 = f1 + π
    @test f4.data ≈ data1 .+ π

    add!(f4, π)
    @test f4.data ≈ data1 .+ 2.0 .* π

    f4 = f1 + f2 
    @test f3 == f4 
    
    add!(f1, f2)
    @test f3 == f1

    set!(f1, data1)

    # subtraction
    f4 = f1 - π
    @test f4.data ≈ data1 .- π

    subtract!(f4, π)
    @test f4.data ≈ data1 .- 2.0 .* π

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

    # multiplication + addition 
    set!(f4, data1)
    mult_add!(f4, f2, 2.0)
    @test f4.data == data1 .+ 2.0 .* data2

    # flatten / unflatten
    x = flatten(f1)
    y = copy(x)
    flatten!(f2, y)

    unflatten!(f3, x)
    unflatten!(f4, y)

    @test f1 == f3 
    @test f2 == f4
end