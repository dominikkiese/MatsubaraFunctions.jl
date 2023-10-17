module MatsubaraFunctions

    using PrecompileTools

    @recompile_invalidations begin
        using Printf
        using LinearAlgebra
        using FunctionWrappers
        using MPI
        using Polyester
        using HDF5
        using Documenter
        using Aqua
    end

    # macro for sanity checks
    sanity_checks() = true 

    macro check(expr, msgs)
        esc(:(if $(@__MODULE__).sanity_checks() @assert($expr, $msgs...) end))
    end

    include("mesh/mesh_point.jl")
    include("mesh/mesh.jl")

    # include("matsubara_freq/matsubara_freq.jl")
    # include("matsubara_grid/matsubara_grid.jl")
    # include("matsubara_func/matsubara_func.jl")
    # include("misc/mpi_helpers.jl")
    # include("misc/pade.jl")
    # include("misc/pulay.jl")
    
    #==
    export
        # mpi_helpers
        mpi_comm,
        mpi_rank,
        mpi_size,
        mpi_split,
        mpi_allreduce!,
        mpi_ismain,
        mpi_println,
        mpi_info,
        mpi_barrier,
        
        # pade
        PadeApprox,
        coeffs, 
        xdat,

        # pulay 
        PeriodicPulay,
        solve!,

        # matsubara_freq
        AbstractParticle,
        Fermion,
        Boson,
        MatsubaraFrequency,
        temperature,
        value, 
        index,
        type,

        # matsubara_grid
        MatsubaraGrid, 
        first_index,
        last_index,
        index_range,
        is_inbounds,
        N,
        first_value,
        last_value,
        value_range,
        indices,
        values,
        info,
        MatsubaraIndex,
        save_matsubara_grid!, 
        load_matsubara_grid,

        # matsubara_func
        MatsubaraFunction,
        grids, 
        grids_shape,
        shape, 
        data_shape,
        absmax,
        add,
        add!, 
        subtract,
        subtract!, 
        mult,
        mult!,
        set!,
        flatten,
        flatten!,
        unflatten!,
        LinearIndex,
        to_Matsubara,
        upper_tail_moments,
        lower_tail_moments,
        sum_me,
        density,
        MatsubaraOperation,
        sgn,
        con,
        MatsubaraSymmetry,
        MatsubaraSymmetryGroup,
        MatsubaraInitFunction,
        save_matsubara_function!, 
        load_matsubara_function,
        save_matsubara_symmetry_group!,
        load_matsubara_symmetry_group
    ==#
end