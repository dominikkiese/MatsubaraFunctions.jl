module MatsubaraFunctions

    using Test 
    using LoopVectorization
    using HDF5

    include("matsubara.jl")
    include("interpolation.jl")
    include("function.jl")
    include("io.jl")

    export
        # matsubara.jl
        MatsubaraGrid, 
        FermionGrid, 
        BosonGrid,

        # function.jl
        MatsubaraFunction, 
        grids_shape,
        shape, 
        data_shape,
        add!, 
        subtract!, 
        mult!,

        # io.jl
        save_matsubara_grid!, 
        load_matsubara_grid, 
        save_matsubara_function!, 
        load_matsubara_function
end
