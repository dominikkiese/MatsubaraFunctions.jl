module MatsubaraFunctions

    using LinearAlgebra 
    using StaticArrays
    using FunctionWrappers
    using MPI
    using HDF5 
    using Documenter
    
    include("matsubara_freq/matsubara_freq.jl")
    include("matsubara_grid/matsubara_grid.jl")
    include("matsubara_func/matsubara_func.jl")

    export
        # matsubara_freq.jl 
        AbstractParticle,
        Fermion,
        Boson,
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
        save_matsubara_grid!, 
        load_matsubara_grid,

        # matsubara_func.jl
        MatsubaraFunction, 
        grids_shape,
        shape, 
        data_shape,
        absmax,
        argmax,
        mpi_comm,
        mpi_rank,
        mpi_size,
        mpi_split,
        mpi_allreduce!,
        mpi_ismain,
        mpi_println,
        mpi_info,
        mpi_barrier,
        add,
        add!, 
        subtract,
        subtract!, 
        mult,
        mult!,
        set!,
        LinearIndex,
        to_Matsubara,
        upper_tail_moments,
        lower_tail_moments,
        sum_me,
        MatsubaraOperation,
        sgn,
        con,
        MatsubaraSymmetry,
        MatsubaraSymmetryGroup,
        speedup,
        MatsubaraInitFunction,
        save_matsubara_function!, 
        load_matsubara_function,
        save_matsubara_symmetry_group!,
        load_matsubara_symmetry_group
end