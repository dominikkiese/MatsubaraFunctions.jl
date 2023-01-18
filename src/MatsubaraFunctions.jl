module MatsubaraFunctions

    using Test 
    using HDF5
    using LoopVectorization
    using StaticArrays

    include("types.jl")
    include("matsubara.jl")
    include("interpolation.jl")
    include("function.jl")
    include("io.jl")

    export
        # types.jl 
        AbstractParticle,
        Fermion,
        Boson,
        AbstractGrid,
        Linear,
        Coarse,

        # matsubara.jl
        MatsubaraGrid, 

        # function.jl
        MatsubaraFunction, 
        grids_shape,
        shape, 
        data_shape,
        add,
        add!, 
        subtract,
        subtract!, 
        mult,
        mult!,
        tail_moments,
        sum,

        # io.jl
        save_matsubara_grid!, 
        load_matsubara_grid, 
        save_matsubara_function!, 
        load_matsubara_function
end
