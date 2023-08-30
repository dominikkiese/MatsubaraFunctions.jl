@testset "Pulay" begin 

    MPI.Init()
    
    # example borrowed from NLsolve
    function f!(F, x)
        F[1] = (x[1] + 3.0) * (x[2]^3 - 7.0) + 18.0 - x[1]
        F[2] = sin(x[2] * exp(x[1]) - 1.0) - x[2]
        return nothing
    end

    P = PeriodicPulay([0.1, 1.2])
    solve!(f!, P)

    @test P.x ≈ [-0.6653795670474839, -0.9983433268474549]
end