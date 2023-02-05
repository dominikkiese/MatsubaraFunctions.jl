module MatsubaraFunctions

    using LinearAlgebra 
    using StaticArrays
    using HDF5
    using Test 

    include("types.jl")
    include("matsubara_freq.jl")
    include("matsubara_grid.jl")
    include("interpolation.jl")
    include("matsubara_func.jl")
    include("io.jl")

    export
        # types.jl 
        AbstractParticle,
        Fermion,
        Boson,
        AbstractGrid,
        Linear,
        Coarse,

        # matsubara_freq.jl 
        MatsubaraFrequency,
        temperature,
        value, 
        index,
        type,

        # matsubara_grid.jl
        MatsubaraGrid, 
        index_range,
        is_inbounds,
        info,

        # matsubara_func.jl
        MatsubaraFunction, 
        grids_shape,
        shape, 
        data_shape,
        absmax,
        argmax,
        add,
        add!, 
        subtract,
        subtract!, 
        mult,
        mult!,
        set!,
        upper_tail_moments,
        lower_tail_moments,
        sum_me,

        # io.jl
        save_matsubara_grid!, 
        load_matsubara_grid, 
        save_matsubara_function!, 
        load_matsubara_function
end